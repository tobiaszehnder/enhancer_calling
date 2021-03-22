#! /usr/bin/env Rscript-4

### This script extracts ATAC-seq peaks that are TSS-distal (1kb) and overlap an enhancer-cluster (maxgap 200bp) as predicted by CRUP.

args <- commandArgs(trailingOnly=T)
if (length(args) != 6) stop('Usage: ./atac_TSSdistal_over_crup.R <outfile.bed> <atac_peaks.bed> <crup_enhancer.bed> <ehmm_promoters.bed> <transcripts.bed> <genome_build>')

# load packages
packages <- c('rtracklayer', 'GenomicFeatures')
suppressMessages(null <- lapply(packages, library, character.only=T))

# load bed files
outfile <- args[1]
atac <- import.bed(args[2])
crup <- import.bed(args[3])
promoters <- import.bed(args[4])
transcripts <- import.bed(args[5])
genome_build <- args[6]

# load transcripts from file (if passed) or pass
tss <- c(resize(transcripts, width=1, fix='start'), promoters)

# define overlaps
atac_tss_distal <- suppressWarnings(atac[!overlapsAny(atac, tss, maxgap=1000)])
enhancers <- suppressWarnings(atac_tss_distal[overlapsAny(atac_tss_distal, crup, maxgap=200)])

# write to file
export.bed(enhancers, outfile)
