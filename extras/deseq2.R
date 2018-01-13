# This R script demonstrates the minimal steps to perform differential
# expression analysis with DESeq2. It expects a tab-delimited file of metadata
# and a file of counts data (rows=gene, cols=sample), typically as output from
# featureCounts.
#
# The output is a typical DESeq2 results file, in tab-delimited format.

require(DESeq2)

### load your gene count data into a data frame
counts.filename <- "replace_me_with_a_filename.tsv"
my.counts <- read.table(filename, header=FALSE, sep="\t")
my.counts <- as.matrix(my.counts)

### load your sample identifiers into a data frame
sample.filename <- "replace_me_with_another_filename.tsv"
my.samples <- read.table(filename, header=TRUE, sep="\t")

### the following assumes you have a factor variable "group"
###  in your sample table that contains levels "GFP" and "KD"

### you can see which columns you have with
###  colnames(my.samples)

dds <- DESeqDataSetFromMatrix(countData = my.counts,
                              colData = my.samples,
                              design = ~ group)
### run DESeq2
dds <- DESeq(dds, betaPrior=TRUE)

### Extract the results for the contrast comparing KD/GFP for the "group"
### factor
res <- results(dds, contrast = c("group","KD","GFP"))

### write your results to file
output.filename <- "replace_me_with_a_target_filename.tsv"
write.table(res, output.filename, row.names=FALSE, col.names=TRUE, quote=FALSE, sep="\t")
