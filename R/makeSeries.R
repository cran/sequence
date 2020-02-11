makeSeries <- function (x, fac) 
{
    chaine = NULL
    for (i in unique(levels(fac))) {
        u = paste(x[fac == i], collapse = "\t",sep="")
        u = paste(i,"\t",u,"\n", collapse = "",sep="")
        chaine = paste(chaine, u, collapse ="",sep="")
    }
    # chaine
    seri <- strsplit(chaine,"\n")
    m = length(seri[[1]])
    K = NULL
    for (i in 1:m) {
        k <- strsplit(as.character(seri[[1]][i]),"\t")
        K = c(K, k)
    }
    invisible(K)
}
