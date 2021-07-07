# Kung PaohVis

This is a tool to bridge data output from the Nanopore SnakeMake pipeline (https://github.com/nanoporetech/Pore-C-Snakemake) to the hyperedge visualization tool PaohVis (https://aviz.fr/paohvis). The tool takes a directory of `.parquet` files (output from some alignment tooling like bwa), and converts them into a `.csv` file and a cleaned data file. The processing does two important QC tasks:

1. **Unique within-read fragment selection.** The tool will filter out any sequence templates that are aligned to identical locations on the reference genome within a single read. 

1. **Unique 'bookending' read selection.** The tool will filter out any reads that span the same starting and ending fragments. There are parameters to control this selection.

## Dependencies

Currently, the tools relies on `fastparquet` (https://pypi.org/project/fastparquet/) to read `.parquet` files using `pandas`.

## Usage

```
usage: kung_paohvis [-h] [-parquet_dir PARQUET_DIR] [-assembly_path ASSEMBLY_PATH] [-chrom [CHROM]] [-filter_frag [FILTER_FRAG]] [-frag_criterion [FRAG_CRITERION]] [-keep_n_reads [KEEP_N_READS]]

A Python3 commandline tool process Pore-C-SnakeMake outputs for PaohVis

optional arguments:
  -h, --help            show this help message and exit
  -parquet_dir PARQUET_DIR
                        The path to a directory containing `.parquet` files.
  -assembly_path ASSEMBLY_PATH
                        The path to assembly file for chromosome mapping.
  -chrom [CHROM]        If not `None`, will subset to a specific chromosome. Defaults to 19.
  -filter_frag [FILTER_FRAG]
                        Reads with `filter_frag` or fewer fragments will be dropped.
  -frag_criterion [FRAG_CRITERION]
                        A valid column name to use asa tie-breaking metric for sequences that mapped to the same location on the reference genome.
  -keep_n_reads [KEEP_N_READS]
                        The number of reads that share beginning and ending fragments to keep. It is biologically feasible for up to 4, with low probability (parental and M phase cells).
```