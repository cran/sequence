compmat <-
function (serMat, alpha = 0.05, meth = "ward.D", printdata = FALSE, 
    printdico = FALSE, printmat = FALSE, eps = 1e-07,clust=TRUE,pca=TRUE) 
{
    if (!is.list(serMat)) 
        stop("The argument is not a list\n")
    nmat <- length(serMat)
    for (i in 1:nmat) {
        if (!is.data.frame(serMat[[i]])) 
            stop("*** ERROR *** the element ", i, " of the list is not a data.frame\n")
        if (dim(serMat[[i]])[1]!=dim(serMat[[i]])[2] ) 
            stop("*** ERROR *** The element ", i, " of the list is not square matrix\n")
	if(any(dim(serMat[[i]])<2))
	    stop("*** ERROR *** The dimensions of matrix ",i," is too small. minimum = 2x2\n")
    }
    ldim=lapply(serMat,dim)
    if(length(unique(ldim))>1)  stop("*** ERROR *** Matrices must be of equal dimension\n")
    text <- c(NULL)
    ident=rep(" ",nmat)
    for (i in 1:nmat) {
        text <- c(text, c(names(serMat[[i]])))
        dico <- levels(factor(text))
	if(is.null(names(serMat[i]))) ident[i]=paste("M",i,sep="") else ident[i]=names(serMat)[i]  
    }
    p <- length(dico)
    if (printdico) {
        message(" Consensus dictionary:\n----------------------\n\n Number of items : ", 
            p, "\n")
        print(dico)
    }
    if (printdata) {
        message("\n PRINTING TRANSITION MATRICES\n ------------------------\n\n")
        print(serMat)
    }
    mdist <- matrix(0, nmat, nmat)
    msign <- mdist
    mstar <- mdist
    bonfer <- alpha/(nmat * (nmat - 1)/2)
    un <- rep(1, p)
    for (i in 1:(nmat - 1)) {
        for (j in (i + 1):nmat) {
            M <- as.matrix(serMat[[i]])
            N <- as.matrix(serMat[[j]])
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
    for(i in 1:nmat) mcom=mcom+serMat[[i]]
    both=clust & pca
    one=(clust | pca) & !both
# Cluster analysis
if(both) {
	        opar=par(no.readonly =TRUE)
                on.exit(par(opar))
	        dev.new(width=14,height=7)
                par(mfrow=c(1,2))
         }  
else if (one) dev.new()
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
    message("\n Valeurs propres ")
    print(corfac$eig)
    message("\n Facteurs \n")
    print(corfac$GOF)
    message("\n")
    plot(corfac$points, xlab = "F1", ylab = "F2", main = "Principal coordinates analysis")
    for (i in 1:nmat) legend(corfac$points[i, 1], corfac$points[i, 
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
