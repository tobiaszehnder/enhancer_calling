#! /usr/bin/env Rscript-4

### This script extracts ATAC-seq peaks that are TSS-distal (1kb) and overlap an enhancer-cluster (maxgap 200bp) as predicted by CRUP.
### NEW: ATAC-seq peak regions are taken for overlap testing, but their according summits are written to file (+/- 250 bp)

args <- commandArgs(trailingOnly=T)
if (length(args) != 7) stop('Usage: ./atac_TSSdistal_over_crup.R <outfile.bed> <atac_peaks.bed> <atac_peaks.summits.bed> <crup_enhancer.bed> <ehmm_promoters.bed> <transcripts.bed> <genome_build>')

# load packages
packages <- c('rtracklayer', 'GenomicFeatures')
suppressMessages(null <- lapply(packages, library, character.only=T))

# load bed files
outfile <- args[1]
atac <- import.bed(args[2])
atac_summits <- import.bed(args[3])
crup <- import.bed(args[4])
promoters <- import.bed(args[5])
transcripts <- import.bed(args[6])
genome_build <- args[7]

# load transcripts from file (if passed) or pass
tss <- c(resize(transcripts, width=1, fix='start'), promoters)

# define overlaps
enhancers <- resize(atac_summits[suppressWarnings(which(!overlapsAny(atac, tss, maxgap=1000) & overlapsAny(atac, crup, maxgap=200)))], 500, 'center')

# write to file
export.bed(enhancers, outfile)
