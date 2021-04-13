#! /usr/bin/env Rscript-4

args <- commandArgs(trailingOnly=T)
if (length(args) != 2) stop('Usage: prun mdl get_transcripts.R genome_build outdir')
build <- args[1]
outdir <- args[2]

# load packages
packages <- c('rtracklayer', 'GenomicFeatures')
suppressMessages(null <- lapply(packages, library, character.only=T))

# fetch transcripts from four possible gene model sources (take the first that exists)
txdb <- {
    tryCatch(makeTxDbFromUCSC(build, tablename='refGene'), error=function(cond) return({
        tryCatch(makeTxDbFromUCSC(build, tablename='knownGene,'), error=function(cond) return({
            tryCatch(makeTxDbFromUCSC(build, tablename='ensGene'), error=function(cond) return({
                tryCatch(makeTxDbFromUCSC(build, tablename='ncbiRefSeq'), error=function(cond) return(NA))
            }))
        }))
    }))
 }

if (is.na(txdb)) {
    print(sprintf('Could not find %s TxDb on UCSC.', build))
    trx <- GRanges()
} else {
    trx <- transcripts(txdb)
}

export.bed(trx, sprintf('%s/transcripts_%s_UCSC.bed', outdir, build))  
