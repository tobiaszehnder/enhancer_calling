#! /usr/bin/env Rscript-4

args <- commandArgs(trailingOnly=T)
if (length(args) != 3) stop('Usage: Rscript trim_regions_to_chromosome_boundaries.R regions.bed genome.sizes outfile.bed')

suppressMessages(library(rtracklayer))

gr <- import.bed(args[1])
sizes <- GRanges(apply(read.table(args[2]), 1, function(row) sprintf('%s:1-%s', row['V1'], row['V2'])))

gr_out <- suppressWarnings(intersect(gr, sizes))
export.bed(gr_out, args[3])
