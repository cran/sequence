smartArrow <-
function (A = c(0, 2), B = c(2, 2), Ra = 0.2, Rb = 0.1, ClegA = "A", 
    ClegB = "B", width = 0.1, col = "lightgreen", ccol = "yellow", 
    density = NULL, angle = 0, marge = 1.3, debord = 1.4, plot = FALSE, 
    trace = TRUE) 
{
currentOptions=options()
on.exit(options(currentOptions))
    options(warn = 1)
    narcs = 100
    seg <- rotation(rbind(A, B), A, -pi/3)
    o <- seg[2, ]
    R <- sqrt(sum((B - A)^2))
    AG <- atan(Ra/R) * marge
    AD <- atan(Rb/R) * marge
    Adisp <- pi/3 - AG - AD
    if (Adisp <= 0) {
        warning("*** warning *** Too short distance between points ", 
            ClegA, " and ", ClegB, ". Arrow not drawn.\n")
        return(Adisp)
    }
    hauteur <- width * debord * sqrt(3/2)
    if (Adisp >= atan(hauteur/R)) 
        corps <- TRUE
    else corps <- FALSE
    basecenter <- rotation(rbind(o, A), o, -AG)
    if (corps) {
        arrow <- NULL
        a <- homoth(form = basecenter, centre = o, scale = (R - 
            (width/2))/R)
        b <- homoth(form = basecenter, centre = o, scale = (R + 
            (width/2))/R)
        arrow <- rbind(arrow, a[2, ], b[2, ])
        dalpha <- (Adisp - atan(hauteur/R))/narcs
        alpha <- 0
        u <- b
        for (i in 1:narcs) {
            u <- rotation(u, o, -dalpha)
            arrow <- rbind(arrow, u[2, ])
        }
        c <- u[2, ]
        d <- homoth(rbind(o, c), o, (R + (width * (debord - 1))/2)/R)
        arrow <- rbind(arrow, d[2, ])
        e <- rotation(basecenter, o, -Adisp)
        arrow <- rbind(arrow, e[2, ])
        f <- homoth(rbind(o, c), o, (R - (width * (debord + 1))/2)/R)
        arrow <- rbind(arrow, f[2, ])
        g <- homoth(rbind(o, c), o, (R - width)/R)
        arrow <- rbind(arrow, g[2, ])
        u <- g
        for (i in 1:narcs) {
            u <- rotation(u, o, dalpha)
            arrow <- rbind(arrow, u[2, ])
        }
    }
    if (!corps) {
        arrow <- NULL
        a <- homoth(form = basecenter, centre = o, scale = (R + 
            (width * debord/2))/R)
        arrow <- rbind(arrow, a[2, ])
        b <- homoth(form = basecenter, centre = o, scale = (R - 
            (width * debord/2))/R)
        arrow <- rbind(arrow, b[2, ])
        c <- rotation(basecenter, o, -Adisp)
        arrow <- rbind(arrow, c[2, ])
    }
    fl <- NULL
    fl$arrow <- arrow
    fl$bg <- col
    fl$A <- A
    fl$B <- B
    fl$Ra <- Ra
    fl$Rb <- Rb
    fl$drawcircles <- TRUE
    leftF <- min(fl$arrow[, 1])
    rightF <- max(fl$arrow[, 1])
    upF <- max(fl$arrow[, 2])
    downF <- min(fl$arrow[, 2])
    fl$left <- min(leftF, A[1] - Ra, B[1] - Rb, 0)
    fl$right <- max(rightF, A[1] + Ra, B[1] + Rb, 2)
    fl$up <- max(upF, A[2] + Ra, B[2] + Rb, 2)
    fl$down <- min(downF, A[2] - Ra, B[2] - Rb, 0)
    class(fl) <- "smartarrow"
    if (plot) 
        plot(x = c(fl$left, fl$right), y = c(fl$down, fl$up), 
            type = "n")
    if (trace) {
        symbols(rbind(A, B), circles = c(Ra, Rb), bg = ccol, 
            add = TRUE, inches = FALSE)
        # lines(fl$arrow, col = "black")
        legend(A[1], A[2], ClegA, bty = "n", xjust = 0.7, yjust = 0.5)
        legend(B[1], B[2], ClegB, bty = "n", xjust = 0.7, yjust = 0.5)
	if(is.null(density)) polygon(fl$arrow, col = col, border="black", angle = angle)
	else
	{
	polygon(fl$arrow, col = "white", border="black")
	polygon(fl$arrow, col = col, border="black", density=density, angle = angle)
	}
    }
    invisible(fl)
}
