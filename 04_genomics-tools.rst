
.. _integration:

Integration
===========
Today we will be working on data from
https://www.ncbi.nlm.nih.gov/pubmed/29241556. These data are from mouse,
looking at the effects of high-fat diet on gene expression (RNA-seq) and
histone modifications (H3K4me1, H3K4me3, H3K27ac).

We will:

- identify candidate intergenic enhancer regions
- identify differentially expressed genes that have enhancer marks nearby
- run a Fisher's exact test to see if differentially upregulated genes are
  enriched for nearby enhancer marks.
- run a Fisher's exact test to see if differentially upregulated genes are
  enriched for increased H3K4me3 at their promoters.


By the end of today, you will:

- Understand what BEDTools is used for and how to find help

- Recognize that BEDTools is just one command-line software tool and that
  it has lots of idiosyncracies, but it does have things in common with
  other tools.

- Be reminded of how to move data back and forth between helix and laptop: Here
  we will do all BEDTools ops on helix, all else on laptop.

- Understand the steps to integrate genomic regions described by *location*
  (here, histone mods identified by ChIP-seq) with experimental results
  described by *gene name* (here, differential expression results)

- Be able to follow a protocol for performing basic bioinformatics analyses
  using BEDTools and R

- Be reminded of how to plot with ggplot2 in R

- Know how to perform a Fisher's exact test in R and interpret the results


Preparation
-----------

Download example files
~~~~~~~~~~~~~~~~~~~~~~
These are data I've prepared ahead of time that we will be using today. You can
see how they were prepared `here
<https://github.com/lcdb/genomics-workshop-2018/blob/master/data/Snakefile>`_.

.. code-block:: bash

    wget https://github.com/lcdb/genomics-workshop-2018/raw/master/data/package.tar.gz
    tar -xf package.tar.gz

.. code-block:: bash

    wget -O - "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE77625&format=file&file=GSE77625%5FmRNA%5FCD%5Fvs%5F16wkHFD%5FDESeq2%5Fresults%2Etxt%2Egz" > GSE77625.txt.gz

Introspection
~~~~~~~~~~~~~

Data are from this paper, published last month:
https://www.ncbi.nlm.nih.gov/pubmed/29241556. They used mm9 coordinates. To
save time, I've already lifted them over to mm10 coordinates.

The directory structure looks like this::

    ├── extra                           # directory of extra files I've created
    │   ├── mm10.chromsizes             # "chromsizes" file for mm10
    │   ├── transcripts.bed             # BED file of transcripts in mm10
    │   ├── x.bed                       # example BED file for teaching
    │   └── y.bed                       # example BED file for teaching
    └── GSE77625                        # directory of files downloaded from GEO
        ├── GSE77625_h3k27ac_chow.bed
        ├── GSE77625_h3k27ac_hfd.bed
        ├── GSE77625_h3k4me1_chow.bed
        ├── GSE77625_h3k4me1_hfd.bed
        └── GSE77625_h3k4me3_chow.bed

Use `head` on each file. You can learn more about the BED format on the `UCSC
page <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_.

If we look closely at the BED files from GEO, they are from the MACS peak caller::

    > head GSE77625/GSE77625_h3k4me3_chow.bed
    chr1    3670401 3672727 MACS_filtered_peak_1    1035.15
    chr1    4491528 4493999 MACS_filtered_peak_2    1440.85
    chr1    4571176 4572360 MACS_filtered_peak_3    1393.38
    chr1    4784173 4786416 MACS_filtered_peak_4    3100.00
    chr1    4807096 4809645 MACS_filtered_peak_5    3100.00
    chr1    4856979 4858869 MACS_filtered_peak_6    3100.00
    chr1    5017846 5021206 MACS_filtered_peak_7    3100.00
    chr1    5082648 5084029 MACS_filtered_peak_8    3100.00
    chr1    6213756 6215799 MACS_filtered_peak_9    3100.00
    chr1    6382408 6383469 MACS_filtered_peak_10   1113.67

Side note on the 5th column
---------------------------

What is that last column? After digging around on the GEO page, I found methods
info in one of the `sample pages for that GEO
entry <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM2055366>`_. In the
"data processing section, they say they used MACS 1.4.0rc2. This is an old
version of MACS, but searching for it I found the `original site has
a README <http://liulab.dfci.harvard.edu/MACS/README.html>`_. At the end of that
README is a description of "Output files". It says::

    Output files

        NAME_peaks.xls is a tabular file which contains information about
        called peaks. You can open it in excel and sort/filter using excel
        functions. Information include: chromosome name, start position of
        peak, end position of peak, length of peak region, peak summit position
        related to the start position of peak region, number of tags in peak
        region, -10*log10(pvalue) for the peak region (e.g. pvalue is 1e-10,
        then this value should be 100), fold enrichment for this region against
        random Poisson distribution with local lambda, FDR in percentage.
        Coordinates in XLS is 1-based which is different with BED format.

        NAME_peaks.bed is BED format file which contains the peak
        locations. You can load it to UCSC genome browser or Affymetrix IGB
        software.

        NAME_summits.bed is in BED format, which contains the peak
        summits locations for every peaks. The 5th column in this file
        is the summit height of fragment pileup. If you want to find
        the motifs at the binding sites, this file is recommended.

I don't think they've converted ``NAME_peaks.xls``, because we don't have that
many columns. I don't think ``NAME_summits.bed`` is what we're looking at,
because I would expect that to be 1-bp peaks. Looking at our BED files, they
are definitely larger. I then downloaded the `tarball package of MACS
<https://github.com/downloads/taoliu/MACS/MACS-1.4.2-1.tar.gz>`_ unpacked it,
and read the README there. It was different! Near the bottom of that page,
I found this::

     2. NAME_peaks.bed is BED format file which contains the peak
     locations. You can load it to UCSC genome browser or Affymetrix IGB
     software. The 5th column in this file is the -10*log10pvalue of peak
     region.

     3. NAME_summits.bed is in BED format, which contains the peak summits
     locations for every peaks. The 5th column in this file is the summit
     height of fragment pileup. If you want to find the motifs at the
     binding sites, this file is recommended.

So I **think** that the 5th column is the -10*log10(pval) of each peak region.


Recap on data
-------------

- Demonstrate that peaks (or domains since this is histone mod data) don't have
  gene IDs

- Demonstrate that we don't have a tsses file

- Demonstrate that the deseq results don't have genomic coords

- Talk about the annoyances in this dataset:

    - peaks are in mm9 coords
    - DESeq2 output is keyed by gene symbol
    - The R data packages that map gene ID to coordinate use Ensembl IDs, not symbol
    - We need to map gene symbol to Ensembl ID, then use that new Ensembl ID to
      lookup the coordinates.
    - Talk about transcripts and genes. What we want is a file of TSSes of
      transcripts for each gene, labeled by that gene.


- Point to the snakefile needed to prep these data. It was more work to prep
  the data than it will be to do this analysis. Also point out that this is
  usually the case.

    - ``transcripts.bed`` has been created for you
    - BED files have been lifted over from mm9 to mm10
    - We don't need to lift over DESeq2 results. Why?


What is bedtools
----------------
- One example of a command-line program. Runs on Mac and Linux, not windows.
  Only way to use it is on the command line. Why?
- Other examples: aligning reads, extracting sequences, counting reads in regions
- Many also have web sites but those sites often are limited. BLAST, multiple
  alignment (clusal, muscle), HMMER


Learning a new tool
-------------------
- Read the docs, try to get it to run. Inspect input and output.
- Getting help (no args; ``-h``, and how this is a convention for
  arbitrary programs)

- Go through the commands, focusing on the "genome arithmetic" section,
  highlighting the big ones (intersect, genomecov, subtract, merge). Don't
  worry about details, the point is to understand what kinds of things BEDTools
  can do, so you know where to look later.

Exercise: which command could we use for getting TSSes?


Example data
------------

.. image:: extras/bedtools/images/bedtools_intersect_-a_x.bed_-b_y.bed.png


Use ``x.bed`` and ``y.bed``. Draw them on the board.

- explain ``intersect``
- explain ``-a`` and ``-b``, and how the order matters. Will likely need smaller files
  for experimenting with
- explain ``-v``

.. code-block:: bash

    cat x.bed
    cat y.bed
    bedtools intersect -a x.bed -b y.bed
    bedtools intersect -a x.bed -b y.bed -u
    bedtools intersect -a x.bed -b y.bed -v
    bedtools intersect -a y.bed -b x.bed -v


Enhancer-like: had H3K4me1 and H3K27ac
--------------------------------------
When we have files with meaningful information in them, we can get interesting regions.

.. code-block:: bash

    bedtools intersect -a GSE77625/GSE77625_h3k4me1_chow.bed -b GSE77625/GSE77625_h3k27ac_chow.bed > enhancer-like_chow.bed
    bedtools intersect -a GSE77625/GSE77625_h3k4me1_hfd.bed -b GSE77625/GSE77625_h3k27ac_hfd.bed > enhancer-like_hfd.bed

    # Intergenic
    bedtools intersect -a enhancer-like_chow.bed -b extra/transcripts.bed -v > intergenic_enhancer-like_chow.bed

    # Closest gene to each enhancer
    bedtools closest -a intergenic_enhancer-like_chow.bed -b extra/transcripts.bed -D a -io -d > closest_genes_to_enhancer_chow.bed

Gotchas
-------

Sorting is important! In fact, we get a hidden error in the "closest" call.

.. code-block:: bash

    bedtools intersect -a GSE77625/GSE77625_h3k4me1_chow.bed -b GSE77625/GSE77625_h3k27ac_chow.bed | bedtools sort -i - > enhancer-like_chow.bed
    bedtools intersect -a GSE77625/GSE77625_h3k4me1_hfd.bed -b GSE77625/GSE77625_h3k27ac_hfd.bed | bedtools sort -i - > enhancer-like_hfd.bed

    # Intergenic
    bedtools intersect -a enhancer-like_chow.bed -b extra/transcripts.bed -v | bedtools sort -i - > intergenic_enhancer-like_chow.bed

    # Closest gene to each enhancer
    bedtools closest -a intergenic_enhancer-like_chow.bed -b extra/transcripts.bed -D a -io -d > closest_genes_to_enhancer_chow.bed

Flank to get tsses
------------------

- go through the flags for flank
- explain chromsizes
- explain bed file
- explain each argument -l, -r, -s, -g, -i
- point out that it's not documented what will happen with a zero-length flank
  -- and highlight that undocumented features are common. Best way to figure it
  out is to experiment.

.. code-block:: bash

    bedtools flank -l 1 -r 0 -s -g extra/mm10.chromsizes -i extra/transcripts.bed > tsses.bed

Gained H3K4me1
--------------

Find gained H3K4me1

- explain that this is NOT the best way to do differential peak calling, but 1)
  people do it anyway, and 2) it will suffice for now. The paper did things
  differently, and found very few differential regions.

.. code-block:: bash

    bedtools intersect -a GSE77625/GSE77625_h3k4me1_hfd.bed -b GSE77625/GSE77625_h3k4me1_chow.bed -v > gained_h3k4me1.bed


We won't do this due to time constraints, but how would you further restrict
gained H3K4me1 sites to only keep those that *also* have gained H3K27ac sites?

How would you get *lost* H3K4me1 sites? And those that also lost H3K27ac?

TSSes with gained H3K4me1
-------------------------

- explain ``-u``
- reminder that this is how we're connecting peaks to gene IDs

.. code-block:: bash

    bedtools intersect -a tsses.bed -b gained_h3k4me1.bed -u > tsses_with_gained_h3k4me1.bed

Move to R
---------

This will likely be on laptops. So we need to set up file transfer from helix.

See https://hpc.nih.gov/docs/transfer.html, we should require Filezilla to be
installed on laptops.

.. code-block:: r

    df <- read.table('GSE77625/GSE77625_chow-vs-HFD-deseq2_results.txt')
    gained <- read.table('tsses_with_gained_h3k4me1.bed')
    closest_to_en <- read.table('closest_genes_to_enhancer_chow.bed')

    head(df)
    head(gained)
    head(closest_to_en)

    df$gained <- FALSE
    df$gained[rownames(df) %in% gained$V4] <- TRUE

    df$closest_to_en <- FALSE
    df$closest_to_en[rownames(df) %in% closest_to_en$V9] <- TRUE


    df$up <- FALSE
    df$dn <- FALSE
    valid <- !is.na(df$padj)
    sig <- valid & df$padj < 0.1
    df[sig & df$log2FoldChange > 0, 'up'] <- TRUE
    df[sig & df$log2FoldChange < 0, 'dn'] <- TRUE

    table(df$up)
    table(df$dn)
    table(df$closest_to_en)

    # which genes went up AND gained h3k4me1?
    idx <- df$up & df$gained
    rownames(df)[idx]

    write.table(rownames(df)[idx], file='upregulated_that_gained_h3k4me1.txt', quote=FALSE, col.names=FALSE, row.names=FALSE)

    # do we want to go here? Maybe just demonstrate; this is a whole 'nother
    # workshop.
    library(ggplot2)
    ggplot(df) + aes(x=log2FoldChange) + geom_histogram(aes(y=..density..)) + facet_grid(gained~.)

    # up- or down-regulated foldchanges are no different in gained or not
    wilcox.test(df$log2FoldChange[df$gained & df$up], df$log2FoldChange[!df$gained & df$up])
    wilcox.test(df$log2FoldChange[df$gained & df$dn], df$log2FoldChange[!df$gained & df$dn])

    # both up- and downregulated genes are enriched for gain in H3K4me1.
    #
    fisher.test(
        matrix(
          c(
             sum(df$up & df$gained),
             sum(df$up & !df$gained),
             sum(!df$up & df$gained),
             sum(!df$up & !df$gained)
          ),
          nrow=2)
    )

    fisher.test(
        matrix(
          c(
             sum(df$dn & df$gained),
             sum(df$dn & !df$gained),
             sum(!df$dn & df$gained),
             sum(!df$dn & !df$gained)
          ),
          nrow=2)
    )
