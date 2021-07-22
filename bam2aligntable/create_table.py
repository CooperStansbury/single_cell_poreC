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


def reformat_bam(input_bam, output_bam):
    """Reformat query_name in INPUT_SAM  and write to OUTPUT_SAM
    This tool reformats an alignment file so that it works with downstream
    steps in the Pore-C pipeline. For both files you can supply '-' if you want
    to read/write from/to stdin/stdout. The 'query_name' field of the alignment
    file will be reformatted so that each alignment in the SAM file has a
    unique query name:
    \b
        <read_id> -> <read_id>:<read_idx>:<align_idx>
    Where 'read_idx' is a unique integer id for each read within the file and
    'align_idx' is a unique integer id for each alignment within the file. The
    tool also adds a 'BX' tag consisting of the 'read_id' to each record.
    """
    
    input_bam = pysam.AlignmentFile(input_bam, "rb")
    output_bam = pysam.AlignmentFile(output_bam, "rb")

    read_indices = {}
    for align_idx, align in enumerate(input_bam.fetch(until_eof=True)):
        read_id = align.query_name
        read_idx = read_indices.get(read_id, None)
        if read_idx is None:
            read_idx = len(read_indices)
            read_indices[read_id] = read_idx
        align.query_name = f"{read_id}:{read_idx}:{align_idx}"
        output_bam.write(align)
        
    output_bam.close()
    
    align_idx += 1
    num_reads = len(read_indices)


def create_table(alignments_bam):
    """A function to Convert a BAM file to a tabular format sorted 
    by read for downstream analysis.
    
    Adapted from: https://github.com/nanoporetech/pore-c/
    
    Parameters:
    -----------------------------
        : alignments_bam (str): the input file name (.bam file from alignment tool)
        
    Returns:
    -----------------------------
        : align_df (pd.DataFrame): dataframe with aligned reads
    """
    af = AlignmentFile(alignments_bam)
    chrom_order = list(af.references)
    assert "NULL" not in chrom_order
    chrom_order.append("NULL")

    align_df = ndm.AlignmentRecord.to_dataframe(
        [ndm.AlignmentRecord.from_aligned_segment(a) for a in af], chrom_order=chrom_order
    )
    align_df = align_df.sort_values(["read_name"])
    num_aligns, num_reads = len(align_df), align_df.read_idx.nunique()
    return align_df