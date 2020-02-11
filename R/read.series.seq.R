read.series.seq <- function (fich="NULL")
# 
# Reads a sequence file. File must be saved from excel as TXT, with TAB
# separators
# returns the series of sequences that must be analysed by compseq.
{
if(fich=="NULL") stop("Missing file name\n")
seri<-scan(file=fich,what="character",sep="\n")
seri<-strsplit(seri,"\t")
seri
}

