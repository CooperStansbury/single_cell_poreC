import pandas as pd
import numpy as np
import os
import sys

def get_output_name(in_name):
    """A function to return a string based on an input file
    for output file naming
    
    Parameters:
    -----------------------------
        : in_name (str): the input file name
        
    Returns:
    -----------------------------
        : out_name (str): string of for output naming
    """
    out_name = os.path.splitext(in_name)[0]
    return out_name


def read_parquet_dir(dirpath):
    """A function to ingest a directory of parquet files into a pandas 
    dataframe
    
    Parameters:
    -----------------------------
        : dirpath (str): location of the directory
        
    Returns:
    -----------------------------
        : df (pd.DataFrame): the 'raw' alignment table    
    """
    # full path resolution
    dirpath = os.path.abspath(dirpath)
    
    df_list = []
    
    for parq_file in os.listdir(dirpath):
        parq_path = f"{dirpath}/{parq_file}" # full file path
        
        if parq_file.endswith('parquet'):
            df = pd.read_parquet(parq_path, engine='fastparquet')
            df_list.append(df)
        
    df = pd.concat(df_list, ignore_index=True)
    return df


def read_assembly(assmebly_path):
    """A function to read an assembly path for chromosome mapping
    
    NOTE: This function expects the column headers to contain EXACTLY:
    
    >>> ['Chromosome', 'Total length (bp)', 'GenBank accession',
       'RefSeq accession']
    
    
    Parameters:
    -----------------------------
        : assmebly_path (str): location of the mapping file
        
    Returns:
    -----------------------------
        : assembly (pd.DataFrame): the 'raw' chromosome mappings
    """
    assembly = pd.read_csv(assmebly_path)
    
    cols = ['Chromosome', 
            'Total length (bp)', 
            'GenBank accession',
            'RefSeq accession']
    
    if not set(assembly.columns) == set(cols):
        raise ValueError("Error: columns of assembly file need to be: 'Chromosome', \
         'Total length (bp)', 'GenBank accession','RefSeq accession'")
    
    return assembly


def load_alignments_from_csv(filepath):
    """A function to load processed BAM files (from pore-C snakemake)
    into a dataframe for hyperedge construction. 
    
    NOTES:
        (1) the snakemake pipeline is here: 
        https://github.com/nanoporetech/Pore-C-Snakemake. 
        We are working with the alignment table (.parquet) that comes from 
        Pore-C-Snakemake/rules/mapping.smk
        
        (2) this function should not add any new dervied columns, as multiple different file
        formats should be 'ingestable'
    
    Parameters:
    -----------------------------
        : filepath (str): location of the .csv file
        
    Returns:
    -----------------------------
        : df (pd.DataFrame): the 'raw' alignment table
    """
    df = pd.read_csv(filepath)
    return df
