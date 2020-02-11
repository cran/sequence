geneseq <- function (nseq, lmin, lmax, order = 0, dico = NULL, mattrans = NULL) 
{
    if (is.null(dico) && is.null(mattrans)) 
        stop("***** ERROR ***** either mattrans or dico must be provided.\n")
    if (is.null(dico) && (order == 0)) 
        stop("***** ERROR ***** order=0 requires dico.\n")
    if (is.null(mattrans) && (order == 1)) 
        stop("***** ERROR ***** order=1 requires mattrans.\n")
    if (order == 0) {
        longlist <- round(runif(nseq, lmin, lmax))
        seri <- rep(list(NULL), nseq)
        for (i in 1:nseq) {
            n <- longlist[i]
            h <- c(paste("A", i, sep = ""), sample(dico, n, replace = TRUE))
            seri[[i]] <- h
        }
    }
    else if (order == 1) {
        if (!is.character(as.vector(mattrans[, 1]))) 
            stop("***** ERROR ***** The first column of mattrans must contain the character identifiers (dico).\n")
        if (dim(mattrans)[1] != dim(mattrans)[2] - 1) 
            stop("***** ERROR ***** mattrans must be a square matrix.\n")
        dico <- as.character(mattrans[, 1])
        p <- length(dico)
        longlist <- round(runif(nseq, lmin, lmax))
        seri <- rep(list(NULL), nseq)
        M <- as.matrix(mattrans[, 2:(p + 1)])
        for (i in 1:nseq) {
            n <- longlist[i]
            seri[[i]] <- paste("A", i, sep = "")
            init <- sample(1:p, 1)
            V <- rep(0, p)
            V[init] <- 1
            seri[[i]] <- c(seri[[i]], dico[init])
            for (j in 2:n) {
                Vp <- M %*% V
                CVp <- Vp
                for (k in 2:p) CVp[k] <- Vp[k] + CVp[k - 1]
                index <- min(which(runif(1) <= CVp))
                V <- rep(0, p)
                V[index] <- 1
                seri[[i]] <- c(seri[[i]], dico[index])
            }
        }
    }
    else stop("***** ERROR ***** only orders 0 and one are managed by this function.\n")
    seri
}
