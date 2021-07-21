import numpy as np
import pandas as pd
import dask

import pysam
from pysam import FastaFile, FastxFile
from pysam import AlignmentFile
from Bio import Restriction
from Bio.Seq import Seq

import pyranges as pr
from enum import Enum
from typing import Dict, List, NewType, Optional
from intake.source.base import DataSource, Schema
from pydantic import BaseModel, confloat, conint, constr

import nanopore_datamodels as ndm


def find_fragment_intervals(digest_param: str, seq: str) -> List[int]:
    """Finds the start positions of all matches of the regex in the sequence
    
    Parameters:
    -----------------------------
        : digest_param (str): the digestion enzyme
        : seq (str): the sequence 
        
    Returns:
    -----------------------------
        : intervals (dict): start and end of digested fragements along reference genome
    """
    enz = getattr(Restriction, digest_param, None)
    if enz is None:
        raise ValueError("Enzyme not found: {}".format(digest_param))
    
    s = Seq(seq)
    positions = [_ - 1 for _ in enz.search(s)]
    intervals = to_intervals(positions, len(seq))
    return intervals


def create_fragment_dataframe(seqid: str, seq: str, digest_param: str) -> pd.DataFrame:
    """Iterate over the sequences in a fasta file and find the
    match positions for the restriction fragment
    
    Parameters:
    -----------------------------
        : seqid (str): the sequence id
        : seq (str): the sequence 
        : digest_param (str): the digestion enzyme
        
    Returns:
    -----------------------------
        : intervals (pd.DataFrame): the distances between cut points
    """
    intervals = (
        pd.DataFrame(find_fragment_intervals(digest_param, seq))
        .assign(chrom=seqid)
        .eval("fragment_length = end - start")
    )
    return intervals


def to_intervals(positions: List[int], chrom_length: int):
    """A utility function get start end intervals from reference"""
    prefix, suffix = [], []
    if (len(positions) == 0) or positions[0] != 0:
        prefix = [0]
    if (len(positions) == 0) or positions[-1] != chrom_length:
        suffix = [chrom_length]
    endpoints = np.array(prefix + positions + suffix)
    return {"start": endpoints[:-1], "end": endpoints[1:]}


def create_virtual_digest(reference_fasta, digest_param):
    """A function to iterate over the sequences in a fasta file 
    and find the match positions for the restriction fragment
    
    Parameters:
    -----------------------------
        : reference_fasta (str): the input file name (.fna.gz reference file)
        
    Returns:
    -----------------------------
        : frag_df (pd.DataFrame): dataframe with digested fragments from the reference genome
    """
    # convert the sequences to a dask bag
    ff =  ndm.IndexedFasta(reference_fasta)
    seq_bag = ff.to_dask()
    dtype = ndm.FragmentRecord.pandas_dtype(overrides={"chrom": pd.CategoricalDtype(ff._chroms, ordered=True)})

    frag_df = (
        pd.concat(
            seq_bag.map(lambda x: (x["seqid"], x["seq"], digest_param))
            .starmap(create_fragment_dataframe)
            .compute()
        )
        .astype(dict(chrom=dtype["chrom"]))
        .sort_values(["chrom", "start"])
        .assign(fragment_id=lambda x: np.arange(len(x), dtype=int) + 1)
        .astype(dtype)
    )

    return frag_df