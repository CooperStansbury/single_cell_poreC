import sys 
import os
import pandas as pd
import numpy as np
import argparse
from datetime import datetime

# local import
import file_io as _io
import read_filters as _rf
    
    
if __name__ == '__main__':
    
    ############################################################################################
    # INPUT ARGUMENT DEFINITIONS
    ############################################################################################
    
    desc = """A Python3 commandline tool process Pore-C-SnakeMake outputs for PaohVis"""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-parquet_dir", 
                        help="The path to a directory containing `.parquet' files.")
    parser.add_argument("-assembly_path", 
                        help="The path to assembly file for chromosome mapping.")
    parser.add_argument("-output_dir", 
                    help="The path to a directory for all output files.")
    parser.add_argument("-chrom", nargs='?', default=19,
                        help="If not `None`, will subset to a specific chromosome. Defaults to 19.")
    parser.add_argument("-filter_frag", nargs='?', default=1,
                        help="Reads with `filter_frag` or fewer fragments will be dropped.")
    parser.add_argument("-frag_criterion", nargs='?', default='perc_of_alignment',
                        help="A valid column name to use asa tie-breaking metric for sequences\
                        that mapped to the same location on the reference genome.")
    parser.add_argument("-keep_n_reads", nargs='?', default=2,
                    help="The number of reads that share beginning and ending fragments to keep. It is\
                         biologically feasible for up to 4, with low probability (parental and M phase cells).")
    
    ############################################################################################
    # INPUT ARGUMENT PARSING
    ############################################################################################
    args = parser.parse_args()
    
    # argument parsing - may need to handle more robustly
    PARQUET_DIRECTORY = args.parquet_dir
    ASSEMBLY_FILEPATH = args.assembly_path
    OUTPUT_DIR = args.output_dir
    CHROMOSOME = args.chrom
    FILTER_N_FRAGMENT = args.filter_frag
    FRAGMENT_CRITERION = args.frag_criterion
    N_TOP_READS = args.keep_n_reads
    TODAY = datetime.today().strftime('%Y-%m-%d')
    
    if PARQUET_DIRECTORY is None:
        raise ValueError("Error: must specify an input directory via `-parquet_dir'")
        
    if ASSEMBLY_FILEPATH is None:
        raise ValueError("Error: must specify an assembly file via `-assembly_path'")
        
    if OUTPUT_DIR is None:
        raise ValueError("Error: must specify an output directory `-output_dir'")      

    ############################################################################################
    # DATA PROCESSING
    ############################################################################################
    
    diagnstics = []
    assembly = _io.read_assembly(ASSEMBLY_FILEPATH)
    df = _io.read_parquet_dir(PARQUET_DIRECTORY)
    
    # get intial diagnostics
    res = _rf.report_alignments(df)
    res['stage'] = 'initial'
    diagnstics.append(res)

    df = _rf.map_chromosome_names(df, assembly)
    df = _rf.add_fragment_midpoints(df) # compute mid points
    
    if not CHROMOSOME == 'None':
        df = _rf.subset_chromosome(df, CHROMOSOME, verbose=False) # get a single chromosome
    df = _rf.drop_low_fragment_count_reads(df, n=FILTER_N_FRAGMENT, verbose=False) # drop reads with low fragment counts
    df = _rf.per_read_filter(df, FRAGMENT_CRITERION, verbose=0) # filter for duplicate within read fragments 
    df = _rf.drop_low_fragment_count_reads(df, n=FILTER_N_FRAGMENT, verbose=False) # drop reads with low fragment counts
    df = _rf.get_maximal_reads(df, N_TOP_READS) # filter for reads spanning the same fragments
    
    # get final diagnostics
    res = _rf.report_alignments(df)
    res['stage'] = 'final'
    diagnstics.append(res)
    
    
    ############################################################################################
    # OUTPUTS
    ############################################################################################

    # save diagnostics
    diagnostic_df = pd.concat(diagnstics, ignore_index=True)
    filename = f"{OUTPUT_DIR}diagnstics_{TODAY}.csv"
    diagnostic_df.to_csv(filename, index=False)
    print(f"done saving: {filename}")

    # save prepared data
    filename = f"{OUTPUT_DIR}cleaned_reads_{TODAY}.csv"
    df.to_csv(filename, index=False)
    print(f"done saving: {filename}")
    
    # save paohviz data
    pao_df = _rf.build_paohvis_output(df)
    filename = f"{OUTPUT_DIR}paohvis_input_{TODAY}.csv"
    pao_df.to_csv(filename, index=False)
    print(f"done saving: {filename}")
    
    
