import pandas as pd
import numpy as np
import sys
import os


def map_chromosome_names(df, assembly):
    """A function to map the accession to a chromosome number and other metadata
    based on the input assembly.
    
    
    Parameters:
    -----------------------------
        : df (pd.DataFrame): the 'raw' alignment table    
        : assembly (pd.DataFrame): the assmebly table
        
    Returns:
    -----------------------------
        : df (pd.DataFrame): the alignment table with chromosome information
    """
    df['chrom'] = df['chrom'].astype(str).str.strip()
    assembly['RefSeq accession'] = assembly['RefSeq accession'].astype(str).str.strip()
    
    df = pd.merge(df, assembly, 
                  left_on='chrom', 
                  right_on='RefSeq accession', 
                  how='left')
    
    return df


def add_fragment_midpoints(df):
    """A function to add a column of fragment midpoints. These coordinates are reference genome 
    positions halfway (rounding up) along the mapped fragment. This decision was made by Dr. Can Chen
    based on the Nanopore logic in pore-C snakemake during generation of the pairs table.
    
    Parameters:
    -----------------------------
        :  df (pd.DataFrame): the 'raw' alignment table
        
    Returns:
    -----------------------------
        : df (pd.DataFrame): the alignment table with 1 additional column
    """
    df['fragment_mean'] = np.ceil((df['fragment_start'] + df['fragment_end']) /2 )
    df['fragment_mean'] = df['fragment_mean'].astype(int)
    return df


def subset_chromosome(df, chrom, verbose=True):
    """A function to filter the alignment table to fragments within a single 
    chromosome.
    
    NOTES:
        (1) There may be reads that contain mapped fragments in different chromsomes. THESE ARE DROPPED.
    
    Parameters:
    -----------------------------
        :  df (pd.DataFrame): the alignment table 
        :  chrom (int): the chomosome to select
        :  verbose (bool): if True, print numner of rows
        
    Returns:
    -----------------------------
        : df (pd.DataFrame): the alignment table for fragments in a single chromosome
    """    
    if verbose:
        print(f"\nnumber of total mapped fragments {chrom}: {len(df)}")
    df = df[df['Chromosome'] == str(chrom)]
    if verbose:
        print(f"number of mapped fragments in chromosome {chrom}: {len(df)}")
    return df


def drop_low_fragment_count_reads(df, n=1, verbose=True):
    """A function to drop reads that have n or fewer fragments 
        
    Parameters:
    -----------------------------
        :  df (pd.DataFrame): the alignment table 
        :  n (int): reads with n or fewer fragments will be dropped
        :  verbose (bool): if True, print numner of rows
        
    Returns:
    -----------------------------
        : df (pd.DataFrame): the alignment table without low fragement reads
    """
    grped  = df.groupby('read_idx', as_index=False).agg({
        'fragment_id' : 'count',
    })
    if verbose:
        print(f"\nn reads BEFORE striping {n}-fragment reads: {len(grped)}")
    
    # filter 
    grped = grped[grped['fragment_id'] > n]
    df = df[df['read_idx'].isin(grped['read_idx'])]
    
    if verbose:
        print(f"n reads AFTER striping {n}-fragment reads: {len(grped)}")
        
    return df


def per_read_filter(df, criterion, verbose=0):
    """A function to filter individual reads and remove duplicate with-in read 
    fragments
    
    Parameters:
    -----------------------------
      :  df (pd.DataFrame): the alignment table
      :  criterion (str): column name of the column to use for determining which 
      fragment (if duplciated), should be kept. The logic takes the fragment with 
      the HIGHEST value for the criterion specified. 
      :  verbose (int): one of 0, 1, 2
          0: show nothing
          1: print before and after alignment dataframe sizes (rows)
          2: print individual read filter results (debugging only) 
        
    Returns:
    -----------------------------
        : df (pd.DataFrame): the 'raw' alignment table
    """
    df_list = []
    
    if verbose > 0:
        print(f"\nalignment table total fragments: {len(df)}")
    
    for read_idx in df['read_idx'].unique():
        read_df = df[df['read_idx'] == read_idx].reset_index()
        
        read_df['read_min'] = read_df['fragment_start'].min() # add first fragment
        read_df['read_max'] = read_df['fragment_end'].max() # add last fragment
        read_df = read_df.sort_values(by=['fragment_id', criterion], ascending=False)
        
        before = len(read_df)
        read_df = read_df.drop_duplicates(subset=['fragment_id'], keep='first')
        after = len(read_df)
        read_df['n_fragments'] = len(read_df)
        
        first_fragment = read_df[read_df['fragment_start'] == read_df['read_min']]['fragment_id'].unique()[0]
        last_fragment = read_df[read_df['fragment_end'] == read_df['read_max']]['fragment_id'].unique()[0]
        read_df['first_fragment'] = first_fragment
        read_df['last_fragment'] = last_fragment
        
        if verbose == 2:
            print(read_idx, before, after)
        
        df_list.append(read_df)
        
    res = pd.concat(df_list, ignore_index=True)
    
    if verbose > 0:
        print(f"alignment table fragments after filtering reads: {len(res)}")
    return res


def get_maximal_reads(df, n=2):
    """A function to return the n reads that share bookending 
    fragment IDS, sorted by mean mapping quality (read)
    
    Parameters:
    -----------------------------
        :  df (pd.DataFrame): the alignment table
        :  n (int): the number of top reads to take (biologically infeasible to take > 4)
        
    Returns:
    -----------------------------
        : df (pd.DataFrame): the alignment table with maximal reads
    """
    
    # reduce rows with unique fragments 
    # - since we perform this excersise on reads, not fragments
    grped = df.groupby('read_idx', as_index=False).agg({
        'first_fragment' : 'first', 
        'last_fragment' : 'first', 
        'n_fragments' : 'first',
        'perc_of_alignment' : np.mean
    })
    
    # NOTE: we may need to do some cross-read fragment
    # decision logic here
    
    # sort by number of fragments and mapping quality
    grped = grped.sort_values(by=['first_fragment', 'last_fragment', 'n_fragments', 'perc_of_alignment'], ascending=False)
    
    # add an identifier for reads that have the same start and end fragment 
    N_TOP = grped.groupby(['first_fragment', 'last_fragment']).cumcount()
    grped.loc[:, 'N_TOP'] = N_TOP + 1
    
    # compute the maximal number of fragments for reads sharing bookended 
    # fragment IDS
    grped['matching_group_max'] = df.groupby(['first_fragment', 'last_fragment'])["n_fragments"].transform(np.max)
    
    # set a flag to take the top N reads with the highest number of fragements
    mask = (grped['n_fragments'] == grped['matching_group_max']) & (grped['N_TOP'] <= n)
    grped['SELECT'] = np.where(mask, 1, 0)

    
    grped = grped[grped['SELECT'] == 1]
    read_ids = grped['read_idx']
    
    # filter the original data frame for only those reads
    df = df[df['read_idx'].isin(read_ids)]
    return df


def build_paohvis_output(df):
    """A function to format the output for paohviz
    
    Parameters:
    -----------------------------
        :  df (pd.DataFrame): the alignment table 
        
    Returns:
    -----------------------------
        : pao_df (pd.DataFrame): paohvis structured CSV
    """
    pao_df = df[['read_idx', 'fragment_id', 'chrom']].reset_index(drop=True)
    pao_df['slot_1'] = ""
    pao_df['slot_2'] = ""
    pao_df = pao_df[['read_idx', 'fragment_id', 'slot_1', 'slot_2', 'chrom']]
    return pao_df 


def report_alignments(df):
    """A function to generate a dataframe of diagnostic information
    from the alignments table
    
    Parameters:
    -----------------------------
        :  df (pd.DataFrame): the alignment table 
        
    Returns:
    -----------------------------
        :  results (pd.DataFrame): diagnostic information
    """
    grped  = df.groupby('read_idx', as_index=False).agg({
        'fragment_id' : 'count',
        'mapping_quality' : np.mean,
        'num_contained_fragments' : np.sum,
    })
    
    params = {
        'n_unique_reads' : df['read_idx'].nunique(),
        'n_fragments' : len(df),
        'n_unique_fragments' : df['fragment_id'].nunique(),
        'n_unique_fragments' : df['fragment_id'].nunique(),
        'mean_read_length' : df['read_length'].mean(),
        'std_read_length' : df['read_length'].std(),
        'max_read_length' : df['read_length'].max(),
        'min_read_length' : df['read_length'].min(),
        'mean_read_mapping_qual' : df['mapping_quality'].mean(),
        'std_read_mapping_qual' : df['mapping_quality'].std(),
        'max_read_mapping_qual' : df['mapping_quality'].max(),
        'min_read_mapping_qual' : df['mapping_quality'].min(),
        'mean_read_perc_of_alignment' : df['perc_of_alignment'].mean(),
        'std_read_perc_of_alignment' : df['perc_of_alignment'].std(),
        'max_read_perc_of_alignment' : df['perc_of_alignment'].max(),
        'min_read_perc_of_alignment' : df['perc_of_alignment'].min(),
        'mean_fragments_per_read' : grped['fragment_id'].mean(),
        'std_fragments_per_read' : grped['fragment_id'].std(),
        'max_fragments_per_read' : grped['fragment_id'].max(),
        'min_fragments_per_read' : grped['fragment_id'].min(),
        'mean_contained_fragments_per_read' : grped['num_contained_fragments'].mean(),
        'std_contained_fragments_per_read' : grped['num_contained_fragments'].std(),
        'max_contained_fragments_per_read' : grped['num_contained_fragments'].max(),
        'min_contained_fragments_per_read' : grped['num_contained_fragments'].min(),
    }
    
    new_rows = []
    for k, v in params.items():
        row = {
            'parameter' : k,
            'value' : v
        }
        
        new_rows.append(row)
        
    results = pd.DataFrame(new_rows)
    return results
    
    