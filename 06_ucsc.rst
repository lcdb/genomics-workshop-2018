UCSC Genome Browser
===================

Kinds of data you can show
--------------------------

- BED files `BED format <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_
- bedGraph
- data that has chrom/start/stop. Can have a score.

BED 3::

    chr1  100  500
    chr1  150  300
    chr1  800  1000

BED 4::

    chr1  100  500   A
    chr1  150  300   B
    chr1  800  1000  C


BED 4 with track lines::

    track name="example track, dense" color="128,0,0" visibility=dense
    chr1  100  500   A
    chr1  150  300   B
    chr1  800  1000  C


    track name="example track, pack" color="0,255,0" visibility=pack
    chr1  100  500   A
    chr1  150  300   B
    chr1  800  1000  C

    track name="example track, squish" color="200,0,200," visibility=squish
    chr1  100  500   A
    chr1  150  300   B
    chr1  800  1000  C

    track name="example track, with score" visibility=pack useScore=1
    chr1  100  500   A  1000
    chr1  150  300   B  500
    chr1  800  1000  C  100

    track name="example track, with strand" visibility=pack

    chr1  100  500   A  0  -
    chr1  150  300   B  0  -
    chr1  800  1000  C  0  +


bedGraph: `bedGraph example <https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE69798&format=file>`_


Kinds of data you can't show
----------------------------
Raw reads (i.e. fastq). They need to be mapped first

Mapped reads. Visible, but they are hard to look at. They need to be aggregated
and usually scaled. BAM cannot be uploaded directly, but must instead be hosted
somewhere. http://helix.nih.gov/~dalerr/ex1.bam (this is hg19).

Data that are not in chrom-start-stop format, e.g.::

    Genes    Thy1
    Lypla1   1.21377
    Tcea1    0
    Atp6v1h  2.77439
    Oprk1    0
    Rb1cc1   2.96399
    Fam150a  0
    St18     0.721644
    Pcmtd1   3.34769
    ...
    ...


Things you can do with UCSC
---------------------------

BLAT
~~~~
Good hit example: TAGATTGTCTTTCGCTCTGA

Multi-hit example: ATGCACACACCACAAACACA

Multiple sequences example::

    >LNA-AsBGL3
    CGATCTCTGCTCTTAA
    >gRNA2
    CCAACACTATTAGATGTTC
    >gRNA3
    TAGATTGTCTTTCGCTCTGA
    >gRNA_LDB1
    TGGCAACTTTGACTATATGC
    >HbG1_3UTR_F
    atgcacacaccacaaacaca
    >BGL3_TSS
    ctaggctttttatagtttggggt

In-silco PCR
~~~~~~~~~~~~

Forward: GTGCATATTCTGAAACGGTAGTGACT
Reverse: TGCCTGGCTCCATCATATCA


Track search
------------
Sometimes what you need is only in a different assembly!

Track hubs
----------
Track hubs are sets of tracks that have been prepared by other people. Some are
public, some only work if you have the URL.

https://helix.nih.gov/~dalerr/adean/xiang-eto2/xiang-eto2.hub.txt


Caveats
-------
What is the y-axis? Rarely is it described. It may take
reverse-engineering the analysis to figure it out. See the
genomics section for an example of trying to figure out what the
scores in a BED file are for a published data set.

Example data
------------

Erythroid long noncoding RNA from Alvarez-Dominquez et al 2014, mm9. DOI:
10.1182/blood-2013-10-530683. `Table S2
<http://www.bloodjournal.org/highwire/filestream/319203/field_highwire_adjunct_files/2/TableS2.xls>`_.

Excel spreadsheet, separate worksheets for each class of ncRNA, mm9. Format::

    Name         Coordinates (mm9)
    lincRNA-EC1  chr5[+]23356389-23365750
    lincRNA-EC2  chr18[+]54581940-54621217
    lincRNA-EC3  chr4[-]109074924-109080172
    ...
    ...

Strategy: our goal is to convert::

    lincRNA-EC1  chr5[+]23356389-23365750

into::

    chr5    23356389   23365750   lincRNA-EC1  0   +



Physical domains from Sexton et al 2012, in fly embryro (dm3). PMID: 22265598.
`Table S1
<http://www.sciencedirect.com/science/MiamiMultiMediaURL/1-s2.0-S0092867412000165/1-s2.0-S0092867412000165-mmc1.xls/272196/html/S0092867412000165/b20f387f4a670f0bebb8bd53313b6424/mmc1.xls>`_

Excel spreadsheet, one worksheet. Format::

    Table S1: List of physical domains, related to Figure 3

    domain id	chrom	start	end	epigenetic class
    1	chr2L	67150	156849	Active
    2	chr2L	156850	250149	Null
    3	chr2L	250150	295749	Active
    4	chr2L	295750	421449	HP1 centromeric
    5	chr2L	421450	472049	HP1 centromeric
    6	chr2L	472050	488649	Active
    7	chr2L	488650	542149	PcG
    8	chr2L	542150	809749	PcG
    ...
    ...


GAF ChIP-seq: `GSE40646 <https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE40646>`_. Called peaks and WIG files. dm3.

Transcription in adult germline stem cells, PMID: 24835570.

FPKM, keyed by gene symbol. Format, for each sample::

    Genes    Thy1
    Lypla1   1.21377
    Tcea1    0
    Atp6v1h  2.77439
    Oprk1    0
    Rb1cc1   2.96399
    Fam150a  0
    St18     0.721644
    Pcmtd1   3.34769
    ...
    ...


