
.. _integration:

Integration
===========
Today we will be working on data from
https://www.ncbi.nlm.nih.gov/pubmed/29241556. These data are from mouse,
looking at the effects of high-fat diet on gene expression (RNA-seq) and
histone modifications (H3K4me1, H3K4me3, H3K27ac). The data are available on
GEO, `GSE77625 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=gse77625>`_.

We will:

- identify candidate intergenic enhancer regions
- identify regions with gained H3K27ac at their promoter
- plot up/down genes with/without gained H3K27ac


By the end of today, you will:

- Understand what BEDTools is used for and how to find help

- Recognize that BEDTools is just one command-line software tool and that
  it has lots of idiosyncracies, but it does have things in common with
  other tools.

- Be reminded of how to move data back and forth between helix and laptop: Here
  we will do all BEDTools ops on helix, all else on laptop.

- Be able to follow a protocol for performing basic bioinformatics analyses
  using BEDTools and R

- Be reminded of how to plot with ggplot2 in R


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


Change to the newly-created ``data`` directory:

.. code-block:: bash

    cd data

While we're in this directory, download the DESeq2 results from GEO into this
directory. Note that special characters in GEO URLs need quote around them.
We use the ``-O`` flag to ``wget``:

.. code-block:: bash

    wget -O deseq-results.txt.gz "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE77625&format=file&file=GSE77625%5FmRNA%5FCD%5Fvs%5F16wkHFD%5FDESeq2%5Fresults%2Etxt%2Egz"

Then gunzip:

.. code-block:: bash

    gunzip deseq-results.txt.gz

Introspection
~~~~~~~~~~~~~

Data are from this paper, published last month:
https://www.ncbi.nlm.nih.gov/pubmed/29241556. They used mm9 coordinates. To
save time, I've already lifted them over to mm10 coordinates.

The directory structure looks like this::

    ├── deseq-results.txt               # file we just downloaded from GEO
    ├── extra                           # directory of extra files I've created
    │   ├── mm10.chromsizes                 # "chromsizes" file for mm10
    │   ├── transcripts.bed                 # BED file of transcripts in mm10
    │   ├── genes.bed                       # BED file of genes in mm10
    │   ├── x.bed                           # example BED file for teaching
    │   └── y.bed                           # example BED file for teaching
    └── GSE77625                        # directory of files downloaded from GEO
        ├── GSE77625_h3k27ac_chow.bed       # H3K27ac domains in chow
        ├── GSE77625_h3k27ac_hfd.bed        # H3K27ac domains in HFD
        ├── GSE77625_h3k4me1_chow.bed       # H3K4me3 domains in chow
        ├── GSE77625_h3k4me1_hfd.bed        # H3K4me1 domains in HFD
        └── GSE77625_h3k4me3_chow.bed       # H3K4me3 domains in chow

Use `head` on each file. You can learn more about the BED format on the `UCSC
page <https://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_.

If we look closely at the BED files from GEO, they are from the MACS peak
caller.  Note that peaks (or domains since this is histone mod data) have
genomic coordinates but don't have gene IDs::

    $ head GSE77625/GSE77625_h3k4me3_chow.bed
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

:Question: How many peaks are there? Which condition and which mark has the
           most peaks?

Note that DESeq2 results have gene IDs, but don't have genomic coordinates::

    $ head deseq-results.txt
              baseMean          log2FoldChange     lfcSE               pvalue                 padj
    Serpina6  5895.82500928936  2.48928902278076   0.0545379886307599  0                      0
    Rhobtb1   3291.54687137     1.95276508740858   0.0611612877537507  1.08731956604379e-223  9.72389887912965e-220
    Saa4      21111.1219005361  2.96047167002528   0.123787400517557   2.09907006812668e-126  1.25146557461712e-122
    Asl       42410.5484534983  -1.72142049473088  0.0773954122626814  1.351328300561e-109    6.04246449595849e-106
    Bhlhe40   2310.29138629314  1.99643457257362   0.0910106893881505  1.17135999139523e-106  4.190188961219e-103
    Aacs      1422.67899510803  3.27241537853794   0.155903781676187   8.10004134319361e-98   2.41462232440602e-94
    Got1      14865.1943802654  -2.53245801431311  0.122703727971087   1.23073925012224e-94   3.14471460395519e-91
    Ccnd1     1305.62849727339  2.48414252966812   0.12291203459522    7.87666962994332e-91   1.76102641251458e-87
    Dact2     579.546268731826  -2.71692983532472  0.136127448792337   1.25892024134677e-88   2.50189415963648e-85

:Question: Is this data organized by transcript or gene?
:Question: How many lines? How many transcripts/genes?
:Question: Why don't we need to lift over DESeq2 results to mm10?

Often we want to know "which genes are bound by a protein", and that's what
we'll be figuring out. To do this, we need gene coordinates, or better,
transcript coordinates. There are many ways of doing this, none of them
straightforward. Most coordinates are provided for Ensembl or RefSeq IDs, but
the authors only provided gene symbol which complicates things.

Common sources for coordinates:

- The `UCSC Table Browser <https://genome.ucsc.edu/goldenPath/help/hgTablesHelp.html>`_
  (requires navigating the interface, and finding by trial-and-error one of the
  table that has gene IDs in the right format)

- `GENCODE <https://www.gencodegenes.org>`_ (data are in GTF format, which can
  be quite difficult to parse)

- `Ensembl BioMart <http://ensembl.org/biomart/martview>`_ (requires navigating
  the interface; download data require reformatting to be useful)

- `BioConductor AnnotationHub <https://bioconductor.org/packages/release/bioc/html/AnnotationHub.html>`_
  (requires quite a bit of R knowledge)

To save time, I've done this in advance (in `this file
<https://github.com/lcdb/genomics-workshop-2018/blob/master/data/Snakefile>`_,
if you're interested). In fact, the preparation may be about as much effort as
the actual analysis! This is not uncommon. The results are in the
``extra/transcripts.bed`` file::

    $ head extra/transcripts.bed
    chr1    3205901 3216344 Xkr4    0       -       ENSMUST00000162897      ENSMUSG00000051951
    chr1    3206523 3215632 Xkr4    0       -       ENSMUST00000159265      ENSMUSG00000051951
    chr1    3214482 3671498 Xkr4    0       -       ENSMUST00000070533      ENSMUSG00000051951
    chr1    4343507 4360314 Rp1     0       -       ENSMUST00000027032      ENSMUSG00000025900
    chr1    4490928 4496413 Sox17   0       -       ENSMUST00000027035      ENSMUSG00000025902
    chr1    4491713 4496363 Sox17   0       -       ENSMUST00000116652      ENSMUSG00000025902
    chr1    4773206 4785710 Mrpl15  0       -       ENSMUST00000130201      ENSMUSG00000033845
    chr1    4773211 4785739 Mrpl15  0       -       ENSMUST00000156816      ENSMUSG00000033845
    chr1    4774436 4785698 Mrpl15  0       -       ENSMUST00000045689      ENSMUSG00000033845
    chr1    4776377 4785739 Mrpl15  0       -       ENSMUST00000115538      ENSMUSG00000033845

:Question: What are the columns? Is this a standard BED file?

What is BEDTools?
-----------------
BEDTools is a "Swiss-army knife of tools for a wide-range of genomics analysis
tasks", especially "genome arithmetic".  Anything that has to do with genomic
coordinates (peaks, gene regions, genomic regions of any kind) can usually be
answered with BEDTools. Using BEDTools is sort of like running a gel. It's a
general tool that's commonly used, and can give you some very interesting
results -- but you have to put the right information into it and make sure
you're getting out what you expect.

- bedtools docs: http://bedtools.readthedocs.io/en/latest/index.html
- extended tutorial: http://quinlanlab.org/tutorials/bedtools/bedtools.html

BEDTools in context
-------------------
BEDTools is one example of a command-line bioinformatics program. It runs on
Mac and Linux, but not Windows. Only way to use it is on the command line,
hence needing to know how to get around in Bash.

:Question: Why do you think the only way to use most bioinformatics programs is
           from the command line?

Other command line tools align reads, extract sequences, count reads in
regions. Still others have companion web servers, though such sites often are
limited. BLAST, multiple alignment (clusal, muscle), HMMER are examples of
this.

Working at the command line puts you in the drivers seat, the same drivers seat
that other bioinformaticians and the tool authors themselves use.


Learning a new tool
-------------------
Learning a new tool is not trivial. You need to read the documentation (which
may be poor or non-existent), try to get it to run. Run it on some small test
data to get a feel for what it wants as input and what it wants as output.

We saw ``man`` as a way of getting help. This is usually for built-in Linux
command line tools. Bioinformatics tools rarely integrate into the ``man``
system. So instead, try getting help by running the program with no args, or
try ``--help`` or ``-h``. This is just a convention; some programs do not
behave nicely!

We will start learning BEDTools by briefly go through the commands. The point
is *not* for you to remember what command does what, but to get a feel for what
*kinds of things* it can do. Then the next time you run across a problem,
you'll think "that seems like something BEDTools could do" and that will give
you a starting point for your searches. It may also give you ideas about what
you can do with your own data.

On Helix, many tools are installed, but we have to enable them first. They are
in "modules", and we need to load the module we want:

.. code-block::

    module load bedtools

This will be enabled as long as we are still connected to Helix during this
session, or we explicitly say ``module unload bedtools``.

See https://hpc.nih.gov/apps for available programs. For example, `here's
the page for bedtools <https://hpc.nih.gov/apps/bedtools.html>`_.

.. code-block:: bash

    bedtools

::

    bedtools: flexible tools for genome arithmetic and DNA sequence analysis.
    usage:    bedtools <subcommand> [options]

    The bedtools sub-commands include:

    [ Genome arithmetic ]
        intersect     Find overlapping intervals in various ways.
        window        Find overlapping intervals within a window around an interval.
        closest       Find the closest, potentially non-overlapping interval.
        coverage      Compute the coverage over defined intervals.
        map           Apply a function to a column for each overlapping interval.
        genomecov     Compute the coverage over an entire genome.
        merge         Combine overlapping/nearby intervals into a single interval.
        cluster       Cluster (but don't merge) overlapping/nearby intervals.
        complement    Extract intervals _not_ represented by an interval file.
        shift         Adjust the position of intervals.
        subtract      Remove intervals based on overlaps b/w two files.
        slop          Adjust the size of intervals.
        flank         Create new intervals from the flanks of existing intervals.
        sort          Order the intervals in a file.
        random        Generate random intervals in a genome.
        shuffle       Randomly redistrubute intervals in a genome.
        sample        Sample random records from file using reservoir sampling.
        spacing       Report the gap lengths between intervals in a file.
        annotate      Annotate coverage of features from multiple files.

    [ Multi-way file comparisons ]
        multiinter    Identifies common intervals among multiple interval files.
        unionbedg     Combines coverage intervals from multiple BEDGRAPH files.

    [ Paired-end manipulation ]
        pairtobed     Find pairs that overlap intervals in various ways.
        pairtopair    Find pairs that overlap other pairs in various ways.

    [ Format conversion ]
        bamtobed      Convert BAM alignments to BED (& other) formats.
        bedtobam      Convert intervals to BAM records.
        bamtofastq    Convert BAM records to FASTQ records.
        bedpetobam    Convert BEDPE intervals to BAM records.
        bed12tobed6   Breaks BED12 intervals into discrete BED6 intervals.

    [ Fasta manipulation ]
        getfasta      Use intervals to extract sequences from a FASTA file.
        maskfasta     Use intervals to mask sequences from a FASTA file.
        nuc           Profile the nucleotide content of intervals in a FASTA file.

    [ BAM focused tools ]
        multicov      Counts coverage from multiple BAMs at specific intervals.
        tag           Tag BAM alignments based on overlaps with interval files.

    [ Statistical relationships ]
        jaccard       Calculate the Jaccard statistic b/w two sets of intervals.
        reldist       Calculate the distribution of relative distances b/w two files.
        fisher        Calculate Fisher statistic b/w two feature files.

    [ Miscellaneous tools ]
        overlap       Computes the amount of overlap from two intervals.
        igv           Create an IGV snapshot batch script.
        links         Create a HTML page of links to UCSC locations.
        makewindows   Make interval "windows" across a genome.
        groupby       Group by common cols. & summarize oth. cols. (~ SQL "groupBy")
        expand        Replicate lines based on lists of values in columns.
        split         Split a file into multiple files with equal records or base pairs.

    [ General help ]
        --help        Print this help menu.
        --version     What version of bedtools are you using?.
        --contact     Feature requests, bugs, mailing lists, etc.


:Exercise: Which command could we use for getting upstream and downstream
           regions of each gene?

:Exercise: Assuming two files `tsses.bed` and `peaks.bed`, how would you
           get promoters with a peak 1kb upstream of TSSes?

Example data
------------

Change to the ``data/extra`` directory.

To get a feel for the BEDTools commands we'll be using, we will be using the
following example files:

.. code-block:: bash

    $ cat x.bed
    chr1    1       100     feature1
    chr1    100     200     feature2
    chr1    150     500     feature3
    chr1    900     950     feature4

.. code-block:: bash

    $ cat y.bed
    chr1    155     200
    chr1    800     901


Intersection is probably the most commonly-used tool. However, note the number
of regions we get back in the result.


.. image:: extras/bedtools/images/bedtools_intersect_-a_x.bed_-b_y.bed.png

:Question: Why do you think there are two regions returned near the 200 bp mark
           in the image above?

Using ``-u`` keeps things in ``a`` that intersect with ``b``. Quoting from the
help::

    -u      Write the original A entry _once_ if _any_ overlaps found in B.
            - In other words, just report the fact >=1 hit was found.
            - Overlaps restricted by -f and -r.

.. image:: extras/bedtools/images/bedtools_intersect_-a_x.bed_-b_y.bed_-u.png

Using ``-u`` is not symmetrical: it matters which file is provided as ``a`` and
which one as ``b``. Here we've switched them, and you can compare with the
previous results:

.. image:: extras/bedtools/images/bedtools_intersect_-a_y.bed_-b_x.bed_-u.png

``-v`` means NOT. Here, "regions in ``a`` that do not intersect ``b``". From the help::

    -v      Only report those entries in A that have _no overlaps_ with B.
            - Similar to "grep -v" (an homage).

.. image:: extras/bedtools/images/bedtools_intersect_-a_x.bed_-b_y.bed_-v.png

``-v`` is asymmetrical as well:

.. image:: extras/bedtools/images/bedtools_intersect_-a_y.bed_-b_x.bed_-v.png

Here is one we can use for getting promoters. Note that a value of zero  (``-r
0``) does not report anything to the right. This is not actually in the
documentation, it is something discovered by experimenting on test files!

.. image:: extras/bedtools/images/bedtools_flank_-r_0_-l_10_-i_x.bed_-g_genome.chromsizes.png


Working with real data
----------------------
When we have files with meaningful information in them, we can get interesting
regions.

:Question: What does the following code do, in biologically-meaningful terms?

.. code-block:: bash

    bedtools intersect -a GSE77625/GSE77625_h3k4me1_chow.bed -b GSE77625/GSE77625_h3k27ac_chow.bed

These commands are about to get long. Here's the same command, but wrapped on
separate lines with a backslash. It's a way of formatting commands: bash will
glue the lines together. It's important to have the spaces right before the
backslashes! If you're typing this in, you can put it all in one line and skip
using the backslashes. This is mostly formatting for display.

.. code-block:: bash

    bedtools intersect \
      -a GSE77625/GSE77625_h3k4me1_chow.bed \
      -b GSE77625/GSE77625_h3k27ac_chow.bed

We need to name the output something useful so we can refer to it later. As we
will see, naming things can get surpisingly annoying.

Let's name the output ``enhancer-like_chow.bed``;

.. code-block:: bash

    bedtools intersect \
      -a GSE77625/GSE77625_h3k4me1_chow.bed \
      -b GSE77625/GSE77625_h3k27ac_chow.bed \
      > enhancer-like_chow.bed

If you haven't done so already, you should start a new file somewhere on your
laptop using Sublime Text (Mac) or Notepad++ (Windows). Paste these commands
into it to keep a record just like we did in R.

Let's do some spot-checks . . .

:Question: How many enhancer-like regions are there?
:Question: Is this more or less than we expect?
:Question: How do we know if we got the commands right?

:Exercise: Given the data I've provided and the files we've just created, how
           do we get intergenic enhancers in chow? (Check ``ls`` again for
           a reminder of what's available)

.. code-block:: bash

    bedtools intersect \
      -a enhancer-like_chow.bed \
      -b extra/transcripts.bed \
      -v \
      > intergenic_enhancer-like_chow.bed

:Question: Compared to our previous results, how many do we expect in the
           output (and why?)


Get TSSes
---------

We have transcripts, but not TSSes. Here's how to get the single 1-bp position
just upstream of TSSes. 

.. code-block:: bash

    bedtools flank \
      -l 1 \
      -r 0 \
      -s \
      -g extra/mm10.chromsizes \
      -i extra/transcripts.bed \
      > tsses.bed

:Question: How would you get 1kb upstream?

Find gained H3K27ac
-------------------
It's not the best way to do it, but a first-pass way of getting differential
regions is to do the intersection between conditions.

:Question: How would you find the H3K27ac present in HFD but not in chow?

.. code-block:: bash

    bedtools intersect \
      -a GSE77625/GSE77625_h3k27ac_hfd.bed \
      -b GSE77625/GSE77625_h3k27ac_chow.bed \
      -v \
      > gained_h3k27ac.bed

The paper performed the differential region calling in a more robust way, and
found very few differential regions. How many did we get?

:Exercise: How would you further restrict gained H3K4me1 sites to only keep
           those that *also* have gained H3K27ac sites?

:Exercise: How would you get *lost* H3K4me1 sites? And those that also lost
           H3K27ac?

TSSes with gained H3K27ac
-------------------------
We have TSSes. We have gained H3K27ac. Now we can figure out *which* TSSes have
gained H3K27ac:

.. code-block:: bash

    bedtools intersect \
      -a tsses.bed \
      -b gained_h3k27ac.bed \
      -u \
      > tsses_with_gained_h3k27ac.bed

Move to R
---------
To save time, we are skipping how to integrate these results with R. If you
want to know how to do that, see the `full version <05_chw_hfd_full.md>`_.
Otherwise, move on to the `minimal version <05_chow_hfd_minimal.Rmd>`_ to save
time.
