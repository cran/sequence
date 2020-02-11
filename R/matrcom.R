matrcom <- function (x, printdata = FALSE, printdico = FALSE, printmat = FALSE, 
    printcom = TRUE) 
{
    if (!is.list(x)) 
        stop("argument is not a list")
    nseq <- length(x)
    for (i in 1:nseq) {
        if (!is.vector(x[[i]])) 
            stop("*** ERROR *** element ", i, " of the list is not a vector")
  
        if (length(x[[i]]) < 3) 
            stop("*** ERROR *** sequence ", i, " of the list has not enough elements")
    }
    text <- c(NULL)
    for (i in 1:nseq) {
        text <- c(text, c(x[[i]][2:length(x[[i]])]))
        dico <- levels(factor(text))
    }
    p <- length(dico)
    if (printdico) {
        message(" Consensus dictionnary:\n----------------------\n\n Number of items : ", 
            p, "\n")
        print(dico)
    }
    u <- matrix(data = 0, p, p)
    listmat <- list(u)
    ident <- c(NULL)
    for (i in 1:nseq) {
        listmat[[i]] <- mattrans(x[[i]][2:length(x[[i]])], 
            dico)
        ident[i] <- x[[i]][1]
    }
    message("\n LIST OF SEQUENCES IDENTIFIERS\n -----------------------------\n", 
        ident, "\n")
    if (printdata) {
        message("\n LISTING OF SEQUENCES\n --------------------\n\n")
        print(x)
    }
    if (printmat) {
        message("\n\n TRANSITION MATRIX\n -----------------\n\n")
        print(listmat)
    }
    mcom <- matrix(0, p, p)
    for (i in 1:nseq) mcom <- mcom + listmat[[i]]
    message("\n mcom \n\n")
    mcom <- data.frame(dico, mcom)
    names(mcom) <- c("id", dico)
    if (printcom) 
        print(mcom)
    invisible(mcom)
}
