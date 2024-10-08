## Introduction

methfast is a program for extracting methylation values from a bed-like file for a target bed file. It can accept gzipped or plain text files and by default assumes the format for the bed file with methylation is like this:

chrom start end methylation coverage

Other formats can be used by specifying the column numbers for the methylation and coverage values.

methfast is written in C using Heng Li's cgranges library to perform overlap queries between methylation and bed intervals.

## Dependencies and Installation

methfast can be compiled via:

gcc methfast.c cgranges.c -o methfast -lz -lm

## Usage

methfast <methylation_bed(.gz)> <target_bed> [-f <frac_col>] [-c <cov_col>] [-m <meth_col>] [-u <unmeth_col>]

Options:
  -f <int>   column number for methylation fraction (default is 4)
  -c <int>   column number for total coverage (default is 5)
  -m <int>   column number for methylated coverage (no default)
  -u <int>   column number for unmethylated coverage (no default)
	
### Test data

A small test dataset is included, which can be tested via:
	
   methfast meth.bed.gz regions.bed


https://github.com/samtools/htslib
	
https://github.com/lh3/cgranges/
