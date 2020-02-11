flux <-
function (tabcoord, tabtr, dic = NULL, fac = c(1, 2), Sort = NULL, 
    threshold = 0, scale = 0.1, cscale = 0.1,main="Transition graph",...) 
{
    n <- dim(tabtr)[1]
    p <- dim(tabtr)[2]
    if (is.character(as.vector(tabtr[, 1]))) {
        dic <- as.vector(tabtr[, 1])
        tabtr <- tabtr[, 2:p]
        p <- p - 1
    }
    nn <- sum(tabtr)
    nfac <- dim(tabcoord)[2] - 3
    idcoord <- as.character(tabcoord[, 1])
    index <- match(dic, idcoord)
    if (!all(idcoord[index] == dic)) 
        stop("***** ERROR ***** the dictionnaries of the coordinates table and of transition table are not identical *****\n")
    tabcoord <- tabcoord[index, 1:(nfac + 3)]
    w <- tabcoord[, 2]
    x <- tabcoord[, fac[1] + 2]
    y <- tabcoord[, fac[2] + 2]
    graph <- NULL
    for (i in 1:n) for (j in 1:n) {
        ni <- sum(tabtr[i, ])
        if (tabtr[i, j] != 0 && i != j && tabtr[i, j]/ni >= threshold) {
            fl <- smartArrow(A = c(x[i], y[i]), B = c(x[j], y[j]), 
                Ra = (w[i]/2) * cscale, Rb = (w[j]/2) * cscale, 
                ClegA = dic[i], ClegB = dic[j], width = tabtr[i, 
                  j] * w[i]/ni * scale, plot <- FALSE, trace = FALSE)
            flcorrect=!(class(fl)=="numeric")
            if(flcorrect) graph <- rbind(graph, as.vector(c(fl$left, fl$right, 
                fl$down, fl$up)))
        }
    }
    left <- min(graph[, 1])
    right <- max(graph[, 2])
    down <- min(graph[, 3])
    up <- max(graph[, 4])
    if (!is.null(Sort)) {
        if (Sort < 0) 
            decreasing = TRUE
        else decreasing = FALSE
        walk <- order(tabcoord[, abs(Sort) + 2], decreasing = decreasing)
        range <- walk
    }
    else range <- (1:n)
    plot(c(left, right), c(down, up), xlab = paste("factor", 
        as.character(fac[1])), ylab = paste("factor", as.character(fac[2])), 
        type = "n")
    for (i in range) for (j in range) {
        ni <- sum(tabtr[i, ])
        notNull=tabtr[i, j] != 0 
        offDiagonal= i != j
        aboveThreshold= tabtr[i, j]/sum(tabtr[i, ]) >= threshold
        if (notNull && offDiagonal && aboveThreshold) {
		A= c(x[i], y[i])
		B = c(x[j], y[j])
		Ra = w[i] * cscale/2
		Rb = w[j] * cscale/2
		ClegA = dic[i]
		ClegB = dic[j]
		width = tabtr[i, j] * w[i]/ni*scale
            fl <- smartArrow(A,B,Ra,Rb,ClegA,ClegB,width,...)
        }
    }
    abline(0, 0)
    abline(v = 0)
    title("Transition graph")
    invisible(graph)
}
