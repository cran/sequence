makedico <-
function (x, printdata = FALSE, printdico = FALSE) 
{
    if (!is.list(x)) 
        stop("Argument not a list")
    nseq <- length(x)
    for (i in 1:nseq) {
        if (!is.vector(x[[i]])) 
            stop("*** ERROR *** element ", i, " of the list not a vector")
        if (length(x[[i]]) < 3) 
            warning("*** ERROR *** sequence ", i, " of the list has not enough elements")
    }
    text <- c(NULL)
    for (i in 1:nseq) {
        text <- c(text, c(x[[i]][2:length(x[[i]])]))
    }
    dico <- levels(factor(text))
    p <- length(dico)
    if(printdata){
	    for(i in 1:nseq) 
	    {print( x[[i]])
             message("\n")}
    }
    if (printdico) {
        message(" Dictionary consensus :\n----------------------\n\n Number of items : ", 
            p, "\n")
        print(dico)
    }
    invisible(dico)
}
