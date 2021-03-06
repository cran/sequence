makeMatrans<-function (ser, printdata = FALSE, printdico = TRUE, printmat = FALSE) 
{
    if (!is.list(ser)) stop("The argument is not a list\n")
    nseq <- length(ser)
    for (i in 1:nseq) {
        if (!is.vector(ser[[i]])) 
            stop("*** ERROR *** the element ", i, " of the list is not a vector\n")
        if (length(ser[[i]]) < 3) 
            stop("*** ERROR *** The sequence ", i, " of the list has not enough elements\n")
    }
    text <- c(NULL)
    for (i in 1:nseq) {
        text <- c(text, c(ser[[i]][2:length(ser[[i]])]))
        dico <- levels(factor(text))
    }
    p <- length(dico)
    if (printdico) {
        message(" Consensus dictionary:\n----------------------\n\n Number of items : ", 
            p, "\n")
        print(dico)
    }
    u <- matrix(data = 0, p, p)
    listmat <- list(u)
    ident <- c(NULL)
    for (i in 1:nseq) {
        listmat[[i]] <- data.frame(mattrans(ser[[i]][2:length(ser[[i]])],dico),row.names=dico)
        names(listmat[[i]])=dico
        ident[i] <- ser[[i]][1]
    }
    mcom=0
    for(i in 1:nseq) mcom=mcom+listmat[[i]]
    if (printdico) 
        message("\n LIST OF IDENTIFIERS OF INDIVIDUALS\n ---------------------------------------\n", 
            ident, "\n")
    if (printdata) {
        message("\n PRINTING SEQUENCES\n ------------------------\n\n")
        print(ser)
    }
    if (printmat) {
        message("\n\n TRANSITION MATRIX\n ----------------------\n\n")
        print(listmat)
    }
    lmat=list() 
    lmat$listmat <- listmat
    lmat$dico <- dico
    lmat$mcom <- mcom
    class(lmat) = "listMat"
    invisible(lmat)

}
