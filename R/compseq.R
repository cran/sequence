compseq<-function (ser, alpha = 0.05, meth = "ward.D", printdata = FALSE, 
    printdico = TRUE, printmat = FALSE, eps = 1e-07,clust=TRUE,pca=TRUE) 
{
    if (!is.list(ser)) 
        stop("The argument is not a list\n")
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
        listmat[[i]] <- mattrans(ser[[i]][2:length(ser[[i]])], 
            dico)
        ident[i] <- ser[[i]][1]
    }
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
    mdist <- matrix(0, nseq, nseq)
    msign <- mdist
    mstar <- mdist
    bonfer <- alpha/(nseq * (nseq - 1)/2)
    un <- rep(1, p)
    for (i in 1:(nseq - 1)) {
        for (j in (i + 1):nseq) {
            M <- listmat[[i]]
            N <- listmat[[j]]
            R <- M + N
            mM <- M %*% un
            mN <- N %*% un
            mR <- R %*% un
            dl <- p^2 - p
            mM[mM == 0] = 1
            mN[mN == 0] = 1
            mR[mR == 0] = 1
            dotR = t(t(R) %*% (diag(1/c(mR), p, p)))
            dotM = t(t(M) %*% (diag(1/c(mM), p, p)))
            dotN = t(t(N) %*% (diag(1/c(mN), p, p)))
            lR = log(dotR)
            lR[is.infinite(lR)] = 0
            lM = log(dotM)
            lM[is.infinite(lM)] = 0
            lN = log(dotN)
            lN[is.infinite(lN)] = 0
            lnV <- sum(R * lR) - sum(M * lM) - sum(N * lN)
            lnV <- abs(-2 * lnV)
            mdist[i, j] <- lnV
            if (i != j) 
                if (dl >= 1) 
                  msign[i, j] <- pchisq(lnV, df = dl)
                else msign[i, j] <- 9999
            if (j <= i) 
                mstar[i, j] <- " "
            else if (dl < 1) 
                mstar[i, j] <- "~"
            else if (msign[i, j] <= bonfer) 
                mstar[i, j] <- "*"
            else mstar[i, j] <- "NS"
        }
    }
    mdist <- mdist + t(mdist)
    message("\n Between sequences distance matrix\n\n\n")
    print(mdist)
    message("\n ChiSquare test ; exact probabilities\n\n\n")
    print(msign)
    message("\n ChiSquare test ; significance levels \n", " Degrees of freedom:", 
        dl, "\n\tNominal risk alpha = ", alpha, "\n\tBonferroni risk", 
        bonfer, "\n\n")
    print.table(mstar, quote = FALSE)
    rm(mstar)
    mcom=0
    for(i in 1:nseq) mcom=mcom+listmat[[i]]
    both=clust & pca
    one=(clust | pca) & !both
# Cluster analysis
if(both) {
	   opar=par(no.readonly=TRUE)
           on.exit(par(opar))
	   dev.new(width=14,height=7)
           par(mfrow=c(1,2))
         }  else if (one) dev.new()
if(clust)
{
    treeseq <- hclust(as.dist(mdist), method = meth)
    plot(treeseq, labels = ident, main = paste("Sequences clustering:\nmethod: ", meth),sub="",xlab="")
}
# Principal coordinates analysis (or Classical multidimensional scaling (MDS) of a data matrix.)
if(pca)
{
    d=dim(mdist)[1]
    message(" Principal coordinates analysis\n")
    message(" ------------------------------\n")
    corfac <- cmdscale(as.dist(mdist), k = d-1, eig = TRUE)
    message("\n Eigenvalues \n ")
    print(corfac$eig)
    message("\n Factors \n")
    print(corfac$GOF)
    message("\n")
    plot(corfac$points, xlab = "F1", ylab = "F2", main = "Principal coordinates analysis")
    for (i in 1:nseq) legend(corfac$points[i, 1], corfac$points[i, 
        2], ident[i], bty = "n", xjust = 0.5, yjust = 0.5)
    abline(h = 0)
    abline(v = 0)
}
    cseq <- NULL
    cseq$dico <- dico
    cseq$mdist <- mdist
    cseq$msign <- msign
    cseq$mcom <- mcom
    class(cseq) = "compseq"
    dev.new(width=7,height=7)
    dev.off()
    invisible(cseq)

}
