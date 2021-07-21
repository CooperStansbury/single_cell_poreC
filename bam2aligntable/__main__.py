import sys 
import os
import pandas as pd
import numpy as np
import argparse
from datetime import datetime

# local imports
import create_table
import virtual_digest
import map_fragments
    
if __name__ == '__main__':
    
    ############################################################################################
    # INPUT ARGUMENT DEFINITIONS
    ############################################################################################
    
    desc = """A Python3 commandline tool process bam files for downstream analysis"""
    parser = argparse.ArgumentParser(description=desc)

    parser.add_argument("-bam", help="The path to a bam file.")
    parser.add_argument("-ref", help="The path to a refenerence genome file.")
    parser.add_argument("-output_dir", 
                    help="The path to a directory for all output files.")
   
    ############################################################################################
    # INPUT ARGUMENT PARSING
    ############################################################################################
    args = parser.parse_args()
    
    bam_path = args.bam
    ref_path = args.ref
    output_dir = args.output_dir

    ############################################################################################
    # DATA PROCESSING
    ############################################################################################
    
    # load the bam file
    align_df = create_table.create_table(bam_path)
    print(f"align_df.shape: {align_df.shape}")
    
    # load the reference file
    frag_df = virtual_digest.create_virtual_digest(ref_path, digest_param='NlaIII')
    print(f"frag_df.shape: {frag_df.shape}")
    
    # assign fragments
    
    pore_c = map_fragments.assign_fragments(align_df, frag_df)
    print(f"pore_c.shape: {pore_c.shape}")
    
    
    ############################################################################################
    # OUTPUTS
    ############################################################################################
    
    TODAY = datetime.today().strftime('%Y-%m-%d')
        
    # save align_df data
    filename = f"{output_dir}raw_align_table_{TODAY}.csv"
    align_df.to_csv(filename, index=False)
    print(f"done saving: {filename}")
    
    # save frag_df data
    filename = f"{output_dir}digested_fragments_table_{TODAY}.csv"
    frag_df.to_csv(filename, index=False)
    print(f"done saving: {filename}")
    
    # save frag_df data
    filename = f"{output_dir}align_table_{TODAY}.csv"
    pore_c.to_csv(filename, index=False)
    print(f"done saving: {filename}")
    
    

    