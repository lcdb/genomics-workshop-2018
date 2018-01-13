#!/bin/bash

# This file illustrates a minimal set of steps to prepare for an RNA-seq
# analysis, starting from a directory full of FASTQ files as you might get from
# a sequencing facility and ending with files of read counts in each gene that
# can be imported into R for use with DESeq2.
#
# This does not includes the steps of downloading and preparing genome fasta
# files, aligner indexes, and annotation files. For some geneome assemblies,
# these are already available on Helix. See
# https://hpc.nih.gov/apps/db.php#genomes for details. Otherwise you will need
# to download FASTA files from your provider of choice (typicially UCSC,
# GENCODE, or Ensembl) and prepare them appropriately, e.g., using hisat2-build
#
# If running on helix, you will need to load the following modules to get the
# necessary software available:
#
#    module load cutadapt
#    module load fastq_screen
#    module load samtools
#    module load hisat
#    module load fastqc
#    module load subread
#    module load picard
#    module load multiqc

### operates on all .fastq.gz files in the directory in which this is executed.
for file in `ls *fastq.gz | sed 's/.fastq.gz//g'` ;
do
    ### remove adapter sequences and other garbage content from your raw reads
    ###    Requires configuration and reference data
    cutadapt -a file:references/adapters.fa -q 20 --minimum-length=25 $file.fastq.gz \
         -o $file.cutadapt.fastq.gz > $file.cutadapt.fastq.gz.log 2>&1

    ### quality control: test for species contamination in your sequences.
    ###    Requires configuration and reference data
    fastq_screen --outdir tmp_workspace/ --force --aligner bowtie2 --conf fastq_screen.config \
         --subset 100000 $file.fastq.gz > $file.fastq_screen.log 2> $file.fastq_screen.error
    ### MORE QUALITY CONTROL GOES HERE! :)

    ### align your sequences to a reference genome.
    ###    Requires configuration and reference data
    hisat2 -x references/dm6-hisat-index -U $file.cutadapt.fastq.gz --threads 6 -S $file.cutadapt.sam \
       > $file.cutadapt.bam.log 2>&1

    ### remove unmapped reads from the SAM file and convert to BAM
    samtools view -Sb -F 0x04 $file.cutadapt.sam > $file.cutadapt.unsorted.bam

    ### Remove the original SAM file, which is quite large
    rm $file.cutadapt.sam

    ### Sort the BAM file for later duplicates marking and indexing
    samtools sort -o $file.cutadapt.bam -O BAM $file.cutadapt.unsorted.bam

    ### Remove the unsorted BAM to conserve space.
    rm $file.cutadapt.unsorted.bam

    ### perform quality control on your aligned data
    fastqc --noextract --quiet --outdir fastqc/ $file.cutadapt.bam
    ### EVEN MORE QUALITY CONTROL GOES HERE! :):)

    ### mark duplicate sequences in your aligned data
    picard  -Xmx4g MarkDuplicates INPUT=$file.cutadapt.bam OUTPUT=$file.cutadapt.markdups.bam \
        METRICS_FILE=$file.cutadapt.markdups.bam.metrics &> $file.cutadapt.markdups.log

    ### get gene counts based on reference genome annotations for your aligned data
    featurecounts -T 1 -a references/dm6.gtf -o $file.cutadapt.bam.featurecounts.txt $file.cutadapt.bam \
    &> $file.cutadapt.bam.featurecounts.txt.log

    ### summarize all that wonderful quality control you've done!
    LC_ALL=en_US.UTF.8 LC_LANG=en_US.UTF-8 multiqc --outdir data/aggregation/ \
      --force --filename multiqc.html --config multiqc_config.yaml data/rnaseq_samples/ \
      data/aggregration/ &> data/aggregation/multiqc.html.log
done
