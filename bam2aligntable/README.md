# bam2aligntable

This is a tool to read `.bam` files, convert them to tabular structures based on the Nanopore data model, and perform a virtual digest of the referenece genome to map aligned fragments to regions on the reference genome between digestion sites. 

## Dependencies

Currently, the tools relies on the following libraries:

- 

## Usage
Example invocation:

```
python3 bam2aligntable/ -bam development/run02_merged.bam -ref ./development/GCF_000001635.27_GRCm39_genomic.fna -output_dir outputs/
```

Flags:

```
usage: [-h] [-bam BAM] [-ref REF] [-output_dir OUTPUT_DIR]

A Python3 commandline tool process bam files for downstream analysis

optional arguments:
  -h, --help            show this help message and exit
  -bam BAM              The path to a bam file.
  -ref REF              The path to a refenerence genome file.
  -output_dir OUTPUT_DIR
                        The path to a directory for all output files.
```