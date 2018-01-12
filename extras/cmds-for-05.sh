#!/bin/bash

set -e

# This script downloads the data and runs the analysis for the BEDTools
# component of the workshop

mkdir -p 12_jan_2018
cd 12_jan_2018

# Download data and unpack
wget https://github.com/lcdb/genomics-workshop-2018/raw/master/data/package.tar.gz
tar -xf package.tar.gz

cd data

# Download DESeq2 results from GEO and uncompress
wget -O deseq-results.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE77625&format=file&file=GSE77625%5FmRNA%5FCD%5Fvs%5F16wkHFD%5FDESeq2%5Fresults%2Etxt%2Egz"
gunzip deseq-results.txt.gz

# Inspection of files
head GSE77625/GSE77625_h3k4me3_chow.bed
head deseq-results.txt
head extra/transcripts.bed

# Only on helix:
# module load bedtools


# Show help
bedtools

# Move to extra directory and run some example bedtools commands
cd extra
bedtools intersect -a x.bed -b y.bed
bedtools intersect -a x.bed -b y.bed -u
bedtools intersect -a y.bed -b x.bed -u
bedtools intersect -a x.bed -b y.bed -v
bedtools intersect -a y.bed -b x.bed -v
bedtools flank -r 0 -l 10 -i x.bed -g ex.chromsizes

# Go back up a dir

cd ..

# Find enhancer-like regions in chow condition
bedtools intersect \
  -a GSE77625/GSE77625_h3k4me1_chow.bed \
  -b GSE77625/GSE77625_h3k27ac_chow.bed \
  > enhancer-like_chow.bed

# Only keep enhancer-like regions that are NOT in transcripts.
bedtools intersect \
  -a enhancer-like_chow.bed \
  -b extra/transcripts.bed \
  -v \
  > intergenic_enhancer-like_chow.bed

# Get the TSS of each transcript
bedtools flank \
  -l 1 \
  -r 0 \
  -s \
  -g extra/mm10.chromsizes \
  -i extra/transcripts.bed \
  > tsses.bed

# Get the regions that had H3K27ac in HFD but not chow (so, gained H3K27ac
# regions)
bedtools intersect \
  -a GSE77625/GSE77625_h3k27ac_hfd.bed \
  -b GSE77625/GSE77625_h3k27ac_chow.bed \
  -v \
  > gained_h3k27ac.bed

# Get the TSSes that gained H3K27ac.
bedtools intersect \
  -a tsses.bed \
  -b gained_h3k27ac.bed \
  -u \
  > tsses_with_gained_h3k27ac.bed


cd ../..
tar -cvf 12_jan_2018.tar.gz 12_jan_2018
