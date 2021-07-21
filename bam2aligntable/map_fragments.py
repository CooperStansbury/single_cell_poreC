import numpy as np
import pandas as pd
import dask
from importlib import reload

import pysam
from pysam import FastaFile, FastxFile
from pysam import AlignmentFile
from Bio import Restriction
from Bio.Seq import Seq
import pyranges as pr
import networkx as nx

from enum import Enum
from typing import Dict, List, NewType, Optional
from intake.source.base import DataSource, Schema
from pydantic import BaseModel, confloat, conint, constr

import nanopore_datamodels as ndm


def map_fragments(pore_c_table, fragment_df, min_overlap_length: int, containment_cutoff: float):
    """A function to join alignments and the fragments from a virtually digested reference genome.
    
    Joining is done using CURRENT version of PyRanges `.join()` with `.new_position("intersection")`
    logic. This takes the alignment file and maps it to boundaries of the virtually
    digested reference with the following logic: From pyranges/pyranges/methods/new_position.py:
    
    ```
    if new_pos == "intersection":

    new_starts = pd.Series(
        np.where(start1.values > start2.values, start1, start2),
        index=df.index,
        dtype=in_dtype)

    new_ends = pd.Series(
        np.where(end1.values < end2.values, end1, end2),
        index=df.index,
        dtype=in_dtype)
    ```
    
    Parameters:
    -----------------------------
        : pore_c_table (pd.DataFrame): alignments from a `.bam` conforming to the nanopore PoreCRecord model
        : fragment_df (pd.DataFrame): the digested reference genome conforming to nanopore FragmentRecord
        : min_overlap_length (int): the minimun number of base pairs to consider when filtering overlaps
        : containment_cutoff (float): cut-off for the percentage of the fragment covered by the aligned read
        
    Returns:
    -----------------------------
        : align_df (pd.DataFrame): a dataframe with aligned reads and their mapped fragments on the 
        reference genome
    """
    pore_c_table['align_start'] = pore_c_table['start'].copy()
    pore_c_table['align_end'] = pore_c_table['end'].copy()
    
    align_range = pr.PyRanges(
        pore_c_table[["chrom", "start", "end", "align_idx", "align_start", "align_end"]].rename(
            columns={"chrom": "Chromosome", "start": "Start", "end": "End"}
        )
    )
    
    fragment_range = pr.PyRanges(
        fragment_df[["chrom", "start", "end", "fragment_id"]].rename(
            columns={"chrom": "Chromosome", "start": "Start", "end": "End"}
        )
    )
    # all overlaps, one to many
    overlaps = align_range.join(fragment_range).new_position("intersection")
    

    if len(overlaps) == 0:
        raise ValueError("No overlaps found between alignments and fragments, this shouldn't happen")
    
    overlaps = (
        overlaps.df.rename(
              columns={
                "Start": "start",
                "End": "end",
                "Start_b": "fragment_start",
                "End_b": "fragment_end",
            }
        ) 
        .eval("overlap_length = (end - start)")
        .eval("perc_of_alignment = (100.0 * overlap_length) / (align_end - align_start)")
        .eval("perc_of_fragment = (100.0 * overlap_length) / (fragment_end - fragment_start)")
        .eval(f"is_contained = (perc_of_fragment >= {containment_cutoff})")
    )
    
    # per-alignment statistics
    by_align = overlaps.groupby("align_idx", sort=True)

    rank = by_align["overlap_length"].rank(method="first", ascending=False).astype(int)
    overlaps["overlap_length_rank"] = rank

    best_overlap = overlaps[overlaps.overlap_length_rank == 1].set_index(["align_idx"])

    contained_fragments = (
        by_align["is_contained"]
        .agg(["size", "sum"])
        .astype({"sum": int})
        .rename(columns={"size": "num_overlapping_fragments", "sum": "num_contained_fragments"})
    )
    align_df = contained_fragments.join(
        best_overlap[
            [
                "fragment_id",
                "fragment_start",
                "fragment_end",
                "overlap_length",
                "perc_of_alignment",
                "perc_of_fragment",
                "is_contained",
            ]
        ]
    )
    dtype = {col: dtype for col, dtype in pore_c_table.dtypes.items() if col in align_df.columns}
    align_df = align_df.astype(dtype)
    return align_df



def assign_fragments(align_table: ndm.AlignmentRecordDf,
                     fragment_df: ndm.FragmentRecordDf,
                     mapping_quality_cutoff: int = 1,
                     min_overlap_length: int = 10,
                     containment_cutoff: float = 99.0,
                    ) -> ndm.PoreCRecordDf:
    """A function to wrap the mapping of fragments
    
    Parameters:
    -----------------------------
        : align_table (ndm.AlignmentRecordDf): alignments from a `.bam` conforming to the nanopore PoreCRecord model
        : fragment_df (pd.DataFrame): the digested reference genome conforming to nanopore FragmentRecord
        : mapping_quality_cutoff (int): threshold for mapping quality
        : min_overlap_length (int): the minimun number of base pairs to consider when filtering overlaps
        : containment_cutoff (float): cut-off for the percentage of the fragment covered by the aligned read
        
    Returns:
    -----------------------------
        : pore_c_table (pd.DataFrame): a dataframe with filters and alignment reference fragements
    """
    # initialise everything as having passed filter
    pore_c_table = ndm.PoreCRecord.init_dataframe(align_table)
    dtype = pore_c_table.dtypes

    align_types = pore_c_table.align_type.value_counts()
    some_aligns = align_types["unmapped"] != align_types.sum()

    if not some_aligns:
        pore_c_table["pass_filter"] = False
        pore_c_table["filter_reason"] = "unmapped"
        return pore_c_table

    fragment_assignments = map_fragments(pore_c_table, 
                                         fragment_df, 
                                         min_overlap_length, 
                                         containment_cutoff)

    pore_c_table = pore_c_table.set_index("align_idx", drop=False)
    pore_c_table.update(fragment_assignments)
    pore_c_table = pore_c_table.astype(dtype)
    
    # run the filtering operations
    pore_c_table = run_filters(pore_c_table, dtype,
                               mapping_quality_cutoff, 
                               min_overlap_length)
    
    return pore_c_table.reset_index(drop=True).sort_values(["read_idx", "align_idx"], ascending=True)


def run_filters(pore_c_table, dtype, mapping_quality_cutoff, min_overlap_length):
    """A function to filter the mapped fragmennts based on the mapping quality and
    the numbver of overlapping base pairs
    
    Parameters:
    -----------------------------
        : pore_c_table (pd.DataFrame): a dataframe with filters and alignment reference fragements
        : dtype (dict): data types
        : mapping_quality_cutoff (int): threshold for mapping quality
        : min_overlap_length (int): the minimun number of base pairs to consider when filtering overlaps
        
    Returns:
    -----------------------------
        : pore_c_table (pd.DataFrame): a dataframe with filters and alignment reference fragements
    """
     # apply alignment-level filters
    unmapped_mask = pore_c_table.align_type == "unmapped"
    pore_c_table.loc[unmapped_mask, "pass_filter"] = False
    pore_c_table.loc[unmapped_mask, "filter_reason"] = "unmapped"

    # if not unmapped, but mapping quality below cutoff then fail
    fail_mq_mask = ~unmapped_mask & (pore_c_table.mapping_quality <= mapping_quality_cutoff)
    pore_c_table.loc[fail_mq_mask, "pass_filter"] = False
    pore_c_table.loc[fail_mq_mask, "filter_reason"] = "low_mq"

    # short overlap filter
    short_overlap_mask = ~(unmapped_mask | fail_mq_mask) & (pore_c_table.overlap_length < min_overlap_length)
    pore_c_table.loc[short_overlap_mask, "pass_filter"] = False
    pore_c_table.loc[short_overlap_mask, "filter_reason"] = "short_overlap"

    # no need to do other checks if nothing left
    if pore_c_table["pass_filter"].any():
        # for the remaining alignments filtering happens on a per-read basis
        by_read_res = (
            pore_c_table[pore_c_table.pass_filter]
            .groupby("read_idx", sort=False, as_index=False)
            .apply(apply_per_read_filters)
        )
        pore_c_table.update(by_read_res[["pass_filter", "filter_reason"]])
        pore_c_table = pore_c_table.astype(dtype)

    pore_c_table.loc[pore_c_table.pass_filter, "filter_reason"] = "pass"
    return pore_c_table


def apply_per_read_filters(read_df):
    """A wrapper function to construct a set of filters for each read which passes the 
    minimum overlap and maping quality threhsolds. These are read-level filters.
        
    Parameters:
    -----------------------------
        : read_df (pd.DataFrame): a dataframe of alignments and mapped fragments
        
    Returns:
    -----------------------------
        : read_df (pd.DataFrame): a dataframe with filters and alignment reference fragements
    """
    return read_df.pipe(filter_singleton).pipe(filter_exact_overlap_on_query).pipe(filter_shortest_path)


def filter_singleton(read_df):
    """A filter to drop reads that have only a single alignment
    
    i.e., if you have a single alignment at this point you fail
    
    Parameters:
    -----------------------------
        : read_df (pd.DataFrame): a dataframe of alignments and mapped fragments
        
    Returns:
    -----------------------------
        : read_df (pd.DataFrame): a dataframe with filters and alignment reference fragements
    """
    if len(read_df) == 1:  
        read_df.loc[:, "pass_filter"] = False
        read_df.loc[:, "filter_reason"] = "singleton"
    return read_df


def filter_exact_overlap_on_query(read_df):
    """A filter select best reads by their alignment score, if there are 
    duplicates

    Parameters:
    -----------------------------
        : read_df (pd.DataFrame): a dataframe of alignments and mapped fragments
        
    Returns:
    -----------------------------
        : read_df (pd.DataFrame): a dataframe with filters and alignment reference fragements
    """
    overlap_on_read = read_df.duplicated(subset=["read_start", "read_end"], keep=False)
    if overlap_on_read.any():
        best_align_idx = read_df.loc[overlap_on_read, :].groupby(["read_start", "read_end"])["align_score"].idxmax()
        overlap_on_read[best_align_idx.values] = False
        read_df.loc[overlap_on_read, "pass_filter"] = False
        read_df.loc[overlap_on_read, "filter_reason"] = "overlap_on_read"
    return read_df


def filter_shortest_path(read_df, aligner="minimap2"):
    """A filter function to wrap the creation of an alignment graph and then 
    find read on the shortest path via bellman_ford methof (networkx implementation). 
    
    Parameters:
    -----------------------------
        : read_df (pd.DataFrame): a dataframe of alignments and mapped fragments
        : aligner (str): can be minimap or bwa
        
    Returns:
    -----------------------------
        : read_df (pd.DataFrame): a dataframe with filters and alignment reference fragements
    """
    aligns = read_df[read_df.pass_filter]
    num_aligns = len(aligns)
    if num_aligns < 2:
        # can't build a graph, so nothing is filtered by this method
        return read_df

    if aligner == "minimap2":
        gap_fn = minimap_gapscore
    elif aligner == "bwa":
        gap_fn = bwa_gapscore
    else:
        raise ValueError(f"Unrecognised aligner: {aligner}")
    graph = create_align_graph(aligns, gap_fn)
    distance, shortest_path = nx.single_source_bellman_ford(graph, "ROOT", "SINK")
    for idx in aligns.index:
        if idx not in shortest_path:
            read_df.at[idx, "pass_filter"] = False
            read_df.at[idx, "filter_reason"] = "not_on_shortest_path"
    return read_df


def create_align_graph(aligns, gap_fn):
    """A function to construct a read graph to find reads along the shortest path. 
    
    Parameters:
    -----------------------------
        : aligns (pd.DataFrame): passing reads (whose filter value=='T')
        : gap_fn (callable): a function to compute the gap score
        
    Returns:
    -----------------------------
        : graph (nx.DiGraph): a directed read graph
    """
    # we'll visit the alignments in order of increasing endpoint on the read, need to keep
    # the ids as the index in the original list of aligns for filtering later
    aligns = aligns[["read_start", "read_end", "read_length", "align_score"]].copy().sort_values(["read_end"])
    node_ids = list(aligns.index)
    graph = nx.DiGraph()
    # initialise graph with root and sink node, and one for each alignment
    # edges in the graph will represent transitions from one alignment segment
    # to the next
    graph.add_nodes_from(["ROOT", "SINK"] + node_ids)
    for align in aligns.itertuples():
        align_idx = align.Index
        align_score = align.align_score
        gap_penalty_start = gap_fn(align.read_start)
        graph.add_edge("ROOT", align_idx, weight=gap_penalty_start - align_score)
        gap_penalty_end = gap_fn(int(align.read_length - align.read_end))
        graph.add_edge(align_idx, "SINK", weight=gap_penalty_end)

    # for each pair of aligned segments add an edge
    for idx_a, align_idx_a in enumerate(node_ids[:-1]):
        align_a_end = aligns.at[align_idx_a, "read_end"]
        for align_idx_b in node_ids[idx_a + 1 :]:  # noqa: E203 black does this
            align_b_score = aligns.at[align_idx_b, "align_score"]
            align_b_read_start = aligns.at[align_idx_b, "read_start"]
            gap_penalty = gap_fn(abs(int(align_b_read_start) - int(align_a_end)))
            graph.add_edge(align_idx_a, align_idx_b, weight=gap_penalty - align_b_score)

    return graph


def minimap_gapscore(length, o1=4, o2=24, e1=2, e2=1):
    """Minimap gap scoring function
    From this minimap working paper:
    
    "Minimap2 performs DP-based global alignment between adjacent anchors in a
    chain. It uses a 2-piece affine gap cost \citep{Gotoh:1990aa}"
    
    Parameters:
    -----------------------------
        : length (int): number of read starts?
        : o1 (int): function param
        : o2 (int): function param
        : e1 (int): function param
        : e2 (int): function param
        
    Returns:
    -----------------------------
        : graph (nx.DiGraph): a directed read graph
    """
    return min([o1 + int(length) * e1, o2 + int(length) * e2])