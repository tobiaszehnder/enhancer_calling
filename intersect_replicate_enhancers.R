#! /usr/bin/env Rscript-4

### This script extracts the intersection of previously called replicate-specific enhancers (enhancers_per_replicate.R)
### It does so by taking the union of overlapping predictions, a more generous approach compared to the strict intersection, and thus accounting for the potential inaccuracy in the replicates' predictions.

args <- commandArgs(trailingOnly=T)
if (length(args) < 3) stop('Usage: ./enhancer_intersect_replicates.R enhancers_outfile.bed nrep enhancers_Rep1.bed ... enhancers_Rep[nrep].bed')
nrep <- args[2]
if (!(nrep == 1 | nrep == 2)) stop("Error: nrep must be 1 or 2")

# load packages
suppressMessages(library(rtracklayer))

# load bed files
e <- lapply(seq_len(nrep), function(i) {
  import.bed(args[i+2])
})

# combine: take the union of the regions that at least partially overlap between replicates
if (nrep == 1) {
  e_combined <- e[[1]]
} else {
  e_combined <- suppressMessages(GenomicRanges::union(e[[1]][overlapsAny(e[[1]],e[[2]])],e[[2]][overlapsAny(e[[2]],e[[1]])]))
}

# write to file
export.bed(e_combined, args[1])
