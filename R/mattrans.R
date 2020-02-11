mattrans <- function (x, dico = NULL, print = FALSE) 
{
    if (!is.character(x)) 
        stop("*** ERROR *** x is not a character vector")
    n <- length(x)
    y <- factor(x)
    if (is.null(dico)) 
        dico <- levels(y)
    p <- length(dico)
    if (print) 
        message("Dictionary : \n", dico, "\n\n")
    m <- matrix(data = 0, nrow = p, ncol = p)
    for (i in 1:(n - 1)) {
        ok <- FALSE
        source <- y[i]
        target <- y[i + 1]
        for (j in 1:p) {
            if (source == dico[j]) {
                k <- j
                ok <- TRUE
            }
            if (target == dico[j]) {
                l <- j
                ok <- TRUE
            }
        }
        if (!ok) 
            stop(" Item not in dictionary :", source, " or ", 
                target)
        m[k, l] <- m[k, l] + 1
    }
    m
}
