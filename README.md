## Introduction

tagbam is a program that modifies a bam file so that it contain a tag for each read with the best-fitting overlapping interval. The intended use is to assign reads to amplicons from amplicon sequencing (eg, Haloplex). The program accepts an input bam file and a bed file with amplicon intervals and creates a new bam file with the name of the best-fitting, strand-specific amplicon in the XN bam tag field. The amplicon bed file must be in the proper bed format (chr, start, end, name, score, strand). The best-fitting amplicon is determined by the minimum sum of the distances between the start and end of a paired-end read fragment and each end of the amplicon. A maximum distance threshold can be specified (default=6 bp). 

tagbam is written in C using the htslib library to parse and manipulate the Bam files and Heng Li's cgranges library to perform overlap query between reads and bed intervals.

## Dependencies and Installation

The htslib must be installed. tagbam can be compiled via:

gcc tagbam.c cgranges.c -o tagbam -lhts -lz

## Usage

tagbam [options] <input.bam> <amplicon.bed> <output.bam>

Options:
  -d <int>   maximum sum of the distance bt amplicon and read ends (default=6)
  -v         verbose
	
### Test data

A small test dataset is included, which can be tested via:
	
   tagbam sub.bam amps.bed out.bam 


https://github.com/samtools/htslib
	
https://github.com/lh3/cgranges/
