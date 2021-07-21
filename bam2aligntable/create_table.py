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