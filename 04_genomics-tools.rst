
.. _integration:

Integration
===========

.. container:: goals

    - Understand what BEDTools is used for and how to find help

    - Recognize that BEDTools is just one command-line software tool and that
      it has lots of idiosyncracies, but it does have things in common with
      other tools.

    - Refresher on moving data back and forth between helix and laptop: Here we
      will do all BEDTools ops on helix, all else on laptop.

    - Rather than do basic analysis in bash (wc, grep, cut, awk, etc) we bring
      it into R as soon as possible.

    - Identify genes whose TSSes gained H3K4me1 in a high-fat diet and that
      were significantly upregulated.

    - Perform a Fisher's exact test in R to see if upregulated genes are
      significantly associated with gain in H3K4me1.

Preparation
-----------

(draw x.bed and y.bed on the board ahead of time)

Downloading arbitrary files
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Download zip file from github

Explain (prob on whiteboard) the flow of data. I.e, it's not going to your
laptop. It's going straight to helix. We're not using bandwidth to the room.

.. code-block:: bash

    wget <url here>


Unzipping files
~~~~~~~~~~~~~~~

Unzip into a working directory.

.. code-block:: bash

    unzip <filename>


Uncompressing tarballs
~~~~~~~~~~~~~~~~~~~~~~

Same thing; maybe we should supply a zip file and a tarball?


Notes on compressed files
~~~~~~~~~~~~~~~~~~~~~~~~~

Maybe just say that there are different kinds, and the way you unpack them
depends on the kind.

- ``.gz``, ``.zip``, ``.tar.gz``, ``.tar``, ``.bz2``


Introspection
~~~~~~~~~~~~~

Explain the data and where it came from, and our challenges.

Just published last month: https://www.ncbi.nlm.nih.gov/pubmed/29241556

- Head on each file. Explain BED format, and how really cols 1, 2, 3 are the
  important ones. Other cols are wildcards and depend on the tool that created
  them. https://github.com/taoliu/MACS is the docs for MACS2; these files are
  5-columns, which don't match any of the descriptions. Turns out they used the
  old MACS.

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

    - transcripts.bed has been creatd for you
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
