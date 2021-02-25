#! /usr/bin/env Rscript-4

### This script computes overlap of predicted enhancers in replicates.
# Usage: ./compute_stats.R > stats

args <- commandArgs(trailingOnly=T)
if (length(args) != 3) stop('Usage: ./enhancer_stats.R outfile n_replicates transcripts')
outfile <- args[1]
nreps <- as.integer(args[2])
transcripts_bed <- args[3]

packages <- c('GenomicRanges', 'rtracklayer', 'ggplot2', 'ggvenn', 'eulerr', 'polylabelr', 'GenSA')
null <- suppressMessages(lapply(packages, library, character.only=T, lib.loc='/project/process_seq_data/r_libs/'))

# arguments to be passed to the script
outdir <- dirname(outfile)
sample <- gsub('.enhancers.stats.pdf', '', basename(outfile))
e_final_bed <- gsub('.stats.pdf', '.bed', outfile)
e_filtered_bed <- sapply(seq_len(nreps), function(i) Sys.glob(sprintf('%s/%s*%s.enhancers.filtered.bed', outdir, sample, i)))
crup_bed <- sapply(seq_len(nreps), function(i) Sys.glob(sprintf('%s/crup/%s*%s.enhancers.crup.bed', outdir, sample, i)))
atac_bed <- sapply(seq_len(nreps), function(i) Sys.glob(sprintf('%s/atac_peaks/%s*%s.atac_peaks.bed', outdir, sample, i)))
promoters_bed <- sapply(seq_len(nreps), function(i) Sys.glob(sprintf('%s/promoter/%s*%s.promoters.bed', outdir, sample, i)))

# load bed files
e_final <- import.bed(e_final_bed)
e_filtered <- sapply(e_filtered_bed, import.bed)
crup <- sapply(crup_bed, import.bed)
atac <- sapply(atac_bed, import.bed)
promoters <- sapply(promoters_bed, import.bed)
tss <- resize(import.bed(transcripts_bed), width=1, fix='start')
tssProm <- sapply(seq_len(nreps), function(i) c(tss, promoters[[i]]))
names(e_filtered) <- names(crup) <- names(atac) <- names(promoters) <- c('Rep1', 'Rep2')
# replicate overlaps
# The sum of a replicate-specific and the intersection can be larger than all replicate-associated predictions because
# one replicate-associated prediction can overlap multiple predictions from the other replicate, thereby creating multiple overlap instances in `ov`
ov <- findOverlaps(e_filtered[[1]], e_filtered[[2]])
venn_reps <- euler(c('Rep1'=length(setdiff(seq_along(e_filtered[[1]]), queryHits(ov))),
                     'Rep2'=length(setdiff(seq_along(e_filtered[[2]]), subjectHits(ov))),
                     'Rep1&Rep2'=length(ov)))

# Enhancer Filtering
n_df <- sapply(seq_len(nreps), function(rep) {
  list(
    n_atac=length(atac[[rep]][!overlapsAny(atac[[rep]], crup[[rep]], maxgap=200) & !overlapsAny(atac[[rep]], tssProm[[rep]], maxgap=1000)]),
    n_atac_crup=length(atac[[rep]][overlapsAny(atac[[rep]], crup[[rep]], maxgap=200) & !overlapsAny(atac[[rep]], tssProm[[rep]], maxgap=1000)]),
    n_atac_tssProm=length(atac[[rep]][!overlapsAny(atac[[rep]], crup[[rep]], maxgap=200) & overlapsAny(atac[[rep]], tssProm[[rep]], maxgap=1000)]),
    n_crup=length(crup[[rep]][!overlapsAny(crup[[rep]], atac[[rep]], maxgap=200) & !overlapsAny(crup[[rep]], tssProm[[rep]], maxgap=1000)]),
    n_crup_tssProm=length(crup[[rep]][!overlapsAny(crup[[rep]], atac[[rep]], maxgap=200) & overlapsAny(crup[[rep]], tssProm[[rep]], maxgap=1000)]),
    n_tssProm=length(tssProm[[rep]][!overlapsAny(tssProm[[rep]], atac[[rep]], maxgap=1000) & !overlapsAny(tssProm[[rep]], crup[[rep]], maxgap=1000)]),
    n_atac_crup_tssProm=length(atac[[rep]][overlapsAny(atac[[rep]], crup[[rep]], maxgap=200) & overlapsAny(atac[[rep]], tssProm[[rep]], maxgap=1000)])
  )
})

venn_filtering <- list()
venn_tssProm <- list()
for (rep in seq_len(nreps)) {
  # Enhancer Filtering
  n_list <- n_df[,rep]
  venn_filtering[[rep]] <- euler(c('CRUP&ATAC'=n_list$n_atac_crup,
                                   'CRUP'=n_list$n_crup,
                                   'ATAC'=n_list$n_atac,
                                   'Promoter / TSS'=n_list$n_tssProm,
                                   'CRUP&Promoter / TSS'=n_list$n_crup_tssProm,
                                   'ATAC&Promoter / TSS'=n_list$n_atac_tssProm,
                                   'CRUP&ATAC&Promoter / TSS'=n_list$n_atac_crup_tssProm),
                                 shape='ellipse')  
  # Promoter / TSS: Rep1
  cols <- c("#9A504FAA", "#EF5649AA")
  ov_tss_prom <- findOverlaps(tss, promoters[[rep]], maxgap=1000)
  venn_tssProm[[rep]] <- euler(c('Promoter'=length(setdiff(seq_along(promoters[[rep]]), subjectHits(ov_tss_prom))),
                                 'TSS'=length(setdiff(seq_along(tss), queryHits(ov_tss_prom))),
                                 'Promoter&TSS'=length(ov_tss_prom)))
}

# plot to file
rep_cols <- c("#7FB8E0FF", "#f7df7fff")
filter_cols <- c("#7FB8E0FF", "#f7df7fff", "#CD534CFF")
tss_cols <-  c("#9A504FAA", "#EF5649AA")
pdf(outfile)
plot(venn_reps, quantities=list(type=c('counts')), fill=rep_cols, main='Replicate Overlap')
plot(venn_filtering[[1]], quantities=list(type=c('counts')), fill=filter_cols, main='Enhancer Filtering: Rep1')
plot(venn_tssProm[[1]], quantities=list(type=c('counts')), fill=tss_cols, main='Promoter / TSS: Rep1')
plot(venn_filtering[[2]], quantities=list(type=c('counts')), fill=filter_cols, main='Enhancer Filtering: Rep2')
plot(venn_tssProm[[2]], quantities=list(type=c('counts')), fill=tss_cols, main='Promoter / TSS: Rep2')
graphics.off()
