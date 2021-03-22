#! /usr/bin/env Rscript-4

args <- commandArgs(trailingOnly=T)
if (length(args) != 2) stop('Usage: prun mdl get_transcripts.R genome_build outdir')
build <- args[1]
outdir <- args[2]

# load packages
packages <- c('rtracklayer', 'GenomicFeatures')
suppressMessages(null <- lapply(packages, library, character.only=T))

if (startsWith(build, 'mm') | startsWith(build, 'hg')) {
    tablename <- 'refGene'
} else {
    tablename <- 'knownGene'
}   

txdb <- tryCatch(makeTxDbFromUCSC(build, tablename=tablename), error=function(cond) return(NA))
if (is.na(txdb)) {
    print('Path to transcripts file does not exist. TSS filtering will only be based on promoters predicted by eHMM.')
    trx <- GRanges()
} else {
    trx <- transcripts(txdb)
}

export.bed(trx, sprintf('%s/transcripts_%s_UCSC.bed', outdir, build))  
