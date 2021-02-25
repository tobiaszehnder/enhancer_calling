#! /usr/bin/env Rscript-4

### This script extracts ATAC-seq peaks that are TSS-distal (1kb) and overlap an enhancer-cluster (maxgap 200bp) as predicted by CRUP.

args <- commandArgs(trailingOnly=T)
if (length(args) != 6) stop('Usage: ./atac_TSSdistal_over_crup.R <outfile.bed> <atac_peaks.bed> <crup_enhancer.bed> <ehmm_promoters.bed> <transcripts.bed/UCSC> <genome_build>')

# load packages
packages <- c('rtracklayer', 'GenomicFeatures')
suppressMessages(null <- lapply(packages, library, character.only=T))

# load bed files
atac <- import.bed(args[2])
crup <- import.bed(args[3])
promoters <- import.bed(args[4])
transcripts_file <- args[5]
genome_build <- args[6]

# load transcripts from file (if passed) or pass
tss <- promoters
if (transcripts_file=='UCSC') {
    if (startsWith(genome_build, 'mm') | startsWith(genome_build, 'hg')) {
	    tablename <- 'refGene'
	} else {
	  	tablename <- 'knownGene'
	}
	txdb <- tryCatch(makeTxDbFromUCSC(genome_build, tablename=tablename), error=function(cond) return(NA))
	if (!is.na(txdb)) {
		tss <- c(resize(transcripts(txdb), width=1, fix='start'), promoters)
		print('TSS filtering based on promoters predicted by eHMM and TSS from UCSC')
	} else {
	    print('No transcripts / TSS passed and no TxDb found on UCSC. TSS filtering will only be based on promoters predicted by eHMM.')
	}
} else if(file.exists(transcripts_file)) { # if transcripts were passed
    tss <- c(resize(import.bed(transcripts_file), width=1, fix='start'), promoters)
	print(sprintf('TSS filtering based on promoters predicted by eHMM and TSS from %s', transcripts_file))
} else {
    print('Path to transcripts file does not exist. TSS filtering will only be based on promoters predicted by eHMM.')
}

# define overlaps
atac_tss_distal <- atac[!overlapsAny(atac, tss, maxgap=1000)]
enhancers <- atac_tss_distal[overlapsAny(atac_tss_distal, crup, maxgap=200)]

# write to file
export.bed(enhancers, args[1])
