add_legend <- function(...) {
        opar <- par(fig=c(0, 1, 0, 1), oma=c(1, 0, 0, 0),
                    mar=c(0, 0, 1, 0), new=TRUE)
        on.exit(par(opar))
        plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n', cex=1.3)
        legend(...)
}

#' 2 Parameters ####
#' plot2linear
#' Comparison of two epidemics (linear models) in the same plot
#'
#' @param y0_1 Epidemic 1: initial inoculum
#' @param y0_2 Epidemic 2: initial inoculum
#' @param r_1  Epidemic 1: rate of progress
#' @param r_2 Epidemic 2: rate of progress
#' @param maxt epidemic duration (days)
#' @return two rect lines
#' @export
#' @examples
#' plot2linear(0.01, 0.01, 0.03, 0.02, 30)

plot2linear <- function(y0_1, y0_2, r_1, r_2, maxt){
        par(mar=c(5,5,5,2))
        epi1 = paste0("y0 = ", y0_1, "   r = ", r_1 )
        epi2 = paste0("y0 = ", y0_2, "   r = ", r_2 )
        curve(
                y0_1 + r_1*x,
                from=0, to=maxt, lwd=3,
                xlab='Time (days)', ylab='Disease',
                xlim=c(0,maxt), ylim=c(0,1),
                cex.lab=2, cex.axis=1.5,
        )
        grid (NULL,NULL, lty = 6, col = "cornsilk2")
        text(0, 0.8, adj=0, cex=1.5,
             expression(y ~"="~ y[0] + x*r)
        )
        curve(
                y0_2 + r_2*x,
                from=0, to=maxt, col= "red",
                add=TRUE, lwd=3
        )
        add_legend("top",
                   legend=c(epi1, epi2),
                   col=c("black", "red"),
                   lty=c(1,1),
                   lwd=3,
                   horiz=TRUE, bty='n',
                   cex=1.2)
}

#' plot2expo
#'
#' Comparison of two epidemics (exponential models) in the same plot
#'
#' @param y0_1 Epidemic 1: initial inoculum
#' @param y0_2 Epidemic 2: initial inoculum
#' @param r_1  Epidemic 1: rate of progress
#' @param r_2 Epidemic 2: rate of progress
#' @param maxt epidemic duration (days)
#' @return two exponential lines
#' @export
#' @examples
#' plot2expo(0.01, 0.01, 0.05, 0.025, 100)

plot2expo <- function(y0_1, y0_2, r_1, r_2, maxt){
        par(mar=c(5,5,5,2))
        epi1 = paste0("y0 = ", y0_1, "   r = ", r_1 )
        epi2 = paste0("y0 = ", y0_2, "   r = ", r_2 )

        curve(
                y0_1 * exp(r_1*x),
                from=0, to=maxt, lwd=3,
                xlab='Time (days)', ylab='Disease',
                cex.lab=2, cex.axis=1.5,
                xlim=c(0,maxt), ylim=c(0,1)
        )
        text(0, 0.8, adj=0, cex=1.5,
             expression(y ~"="~ y[0] ^ ~ ~ r*x)
             )
        grid (NULL,NULL, lty = 6, col = "cornsilk2")

        curve(
                y0_2 * exp(r_2*x),
                from=0, to=maxt, col= "red",
                add=TRUE, lwd=3
        )
        add_legend("top",
                   legend=c(epi1, epi2),
                   col=c("black", "red"),
                   lty=c(1,1),
                   lwd=3,
                   horiz=TRUE, bty='n', cex=1.2)
}

#' plot2logis
#'
#' Comparison of two epidemics (logistic models) in the same plot
#'
#' @param y0_1 Epidemic 1: initial inoculum
#' @param y0_2 Epidemic 2: initial inoculum
#' @param r_1  Epidemic 1: rate of progress
#' @param r_2 Epidemic 2: rate of progress
#' @param maxt epidemic duration (days)
#' @return two logistic curves
#' @export
#' @examples
#' plot2logis(0.001, 0.001, 0.4,0.2, 60)

plot2logis <- function(y0_1, y0_2, r_1, r_2, maxt){
        par(mar=c(5,5,5,2))
        epi1 = paste0("y0 = ", y0_1, "   r = ", r_1 )
        epi2 = paste0("y0 = ", y0_2, "   r = ", r_2 )

        curve(
                1/(1+((1/y0_1)-1)*exp(-r_1*x)),
                from=0, to=maxt, lwd=3,
                xlab='Time (days)', ylab='Disease', cex.lab=2,
                xlim=c(0,maxt), ylim=c(0,1), cex.axis=1.5
        )
        grid (NULL,NULL, lty = 6, col = "cornsilk2")
        text(0, 0.8, adj=0, cex=1.5,
             expression(y ~"="~ frac(1, 1+((1/y[0])-1)^ ~ ~ -r*x))
        )
        curve(
                1/(1+((1/y0_2)-1)*exp(-r_2*x)),
                from=0, to=maxt, col= "red",
                add=TRUE, lwd=3
        )
        add_legend("top",
                   legend=c(epi1, epi2),
                   col=c("black", "red"),
                   lty=c(1,1),
                   lwd=3,
                   horiz=TRUE, bty='n', cex=1.2)
}

#' plot2mono
#'
#' Comparison of two epidemics (monomolecular models) in the same plot
#'
#' @param y0_1 Epidemic 1: initial inoculum
#' @param y0_2 Epidemic 2: initial inoculum
#' @param r_1  Epidemic 1: rate of progress
#' @param r_2 Epidemic 2: rate of progress
#' @param maxt epidemic duration (days)
#' @return two monomolecular curves
#' @export
#' @examples
#' plot2mono(0.01, 0.01, 0.03, 0.02, 100)

plot2mono <- function(y0_1, y0_2, r_1, r_2, maxt){
        par(mar=c(5,5,5,2))
        epi1 = paste0("y0 = ", y0_1, "   r = ", r_1 )
        epi2 = paste0("y0 = ", y0_2, "   r = ", r_2 )

        curve(
                1-(1-y0_1)*exp(-r_1*x),
                from=0, to=maxt, lwd=3,
                xlab='Time (days)', ylab='Disease', cex.lab=2,
                xlim=c(0,maxt), ylim=c(0,1), cex.axis=1.5
        )

        grid (NULL,NULL, lty = 6, col = "cornsilk2")

        text(0, 0.8, adj=0, cex=1.5,
             expression(y ~"="~ 1 - (1 - y[0]) ^ ~ ~ -r*x)
        )
        curve(
                1-(1-y0_2)*exp(-r_2*x),
                from=0, to=maxt, col= "red",
                add=TRUE, lwd=3
        )
        add_legend("top",
                   legend=c(epi1, epi2),
                   col=c("black", "red"),
                   lty=c(1,1),
                   lwd=3,
                   horiz=TRUE, bty='n', cex=1.2)
}

#' plot2gomp
#'
#' Comparison of two epidemics (gompertz models) in the same plot
#'
#' @param y0_1 Epidemic 1: initial inoculum
#' @param y0_2 Epidemic 2: initial inoculum
#' @param r_1  Epidemic 1: rate of progress
#' @param r_2 Epidemic 2: rate of progress
#' @param maxt epidemic duration (days)
#' @return two gompertz curves
#' @export
#' @examples
#' plot2gomp(0.01, 0.01, 0.05, 0.025, 100)

plot2gomp <- function(y0_1, y0_2, r_1, r_2, maxt){
        par(mar=c(5,5,5,2))
        epi1 = paste0("y0 = ", y0_1, "   r = ", r_1 )
        epi2 = paste0("y0 = ", y0_2, "   r = ", r_2 )

        curve(
                exp(-(-log(y0_1))*exp(-r_1*x)),
                from=0, to=maxt, lwd=3,
                xlab='Time (days)', ylab='Disease', cex.lab=2,
                xlim=c(0,maxt), ylim=c(0,1), cex.axis=1.5
        )

        grid (NULL,NULL, lty = 6, col = "cornsilk2")

        text(0, 0.8, adj=0, cex=1.5,
             expression(y ~"="~ exp(ln(y[0])* ~ ~  exp(-rx)))
        )

        curve(
                exp(-(-log(y0_2))*exp(-r_2*x)),
                from=0, to=maxt, col= "red",
                add=TRUE, lwd=3
        )
        add_legend("top",
                   legend=c(epi1, epi2),
                   col=c("black", "red"),
                   lty=c(1,1),
                   lwd=3,
                   horiz=TRUE, bty='n', cex=1.2)

}

#' plot3logis
#'
#' Comparison of three epidemics (logisic model) at the same plot
#'
#' @param y0_1 Epidemic 1: initial inoculum
#' @param y0_2 Epidemic 2: initial inoculum
#' @param y0_3 Epidemic 3: initial inoculum
#' @param r_1  Epidemic 1: rate of progress
#' @param r_2 Epidemic 2: rate of progress
#' @param r_1  Epidemic 3: rate of progress
#' @param maxt epidemic duration (days)
#' @return Three logistic curves
#' @export
#' @examples
#' plot3logis(0.01, 0.01, 0.01,  # Initial inoculum: Epid 1, Epid 2, Epid 3
#'            0.05, 0.025, 0.01, # Rate of progress: Epid 1, Epid 2, Epid 3
#'            100)               # Epidemic duration

plot3logis <- function(y0_1, y0_2, y0_3, r_1, r_2, r_3, maxt){
        par(mar=c(5,5,5,2))
        epi1 = paste0("y0 = ", y0_1, "   r = ", r_1 )
        epi2 = paste0("y0 = ", y0_2, "   r = ", r_2 )
        epi3 = paste0("y0 = ", y0_3, "   r = ", r_3 )

        curve(
                1/(1+(1-y0_1)/y0_1*exp(-r_1*x)),
                from=0, to=maxt, lwd=3,
                xlab='Time (days)', ylab='Disease', cex.lab=2,
                xlim=c(0,maxt), ylim=c(0,1), cex.axis=1.5
        )

        grid (NULL,NULL, lty = 6, col = "cornsilk2")

        text(0, 0.8, adj=0, cex=1.5,
             expression(y ~"="~ 1 - (1 - y[0]) ^ ~ ~ -r*x)
        )

        curve(
                1/(1+(1-y0_2)/y0_2*exp(-r_2*x)),
                from=0, to=maxt, col= "red",
                add=TRUE, lwd=3
        )
        curve(
                1/(1+(1-y0_3)/y0_3*exp(-r_3*x)),
                from=0, to=maxt, col= "green",
                add=TRUE, lwd=3
        )

               add_legend("top",
                   legend=c(epi1, epi2, epi3),
                   col=c("black", "red", "green"),
                   lty=c(1,1,1),
                   lwd=3,
                   horiz=TRUE, bty='n', cex=1.2)
}

#' plot2logis3M
#'
#' Comparison of three epidemics (logisic model) at the same plot.
#' from Madden, 2007: a/(1+((a/y0)-1)*exp(-r*x)),
#'
#' @param y0_1 Epidemic 1: initial inoculum
#' @param y0_2 Epidemic 2: initial inoculum
#' @param r_1  Epidemic 1: rate of progress
#' @param r_2 Epidemic 2: rate of progress
#' @param maxt epidemic duration (days)
#' @return Two logistic curves
#' @export
#' @examples
#' plot2logis3M(1, 0.7,       # Assintote: Epid 1, Epid 2
#'              0.001, 0.001, # Initial inoculum: Epid 1, Epid 2
#'              0.4, 0.4,     # Rate of progress: Epid 1, Epid 2
#'              50)           # Epidemic duration

plot2logis3M <- function(a_1, a_2, y0_1, y0_2, r_1, r_2, maxt){
        par(mar=c(5,5,5,2))
        epi1 = paste0("a = ", a_1, "  y0 = ", y0_1, "   r = ", r_1 )
        epi2 = paste0("a = ", a_2, "  y0 = ", y0_2, "   r = ", r_2 )

        curve(
                a_1/(1+((a_1/y0_1)-1)*exp(-r_1*x)),
                from=0, to=maxt, lwd=3,
                xlab='Time (days)', ylab='Disease', cex.lab=2,
                xlim=c(0,maxt), ylim=c(0,1), cex.axis=1.5
        )
        grid (NULL,NULL, lty = 6, col = "cornsilk2")
        text(0, 0.8, adj=0, cex=1.5,
             expression(y ~"="~ frac(a,(1+((a/y[0])-1)* ~ ~exp(-r*x))))
        )
        curve(
                a_2/(1+((a_2/y0_2)-1)*exp(-r_2*x)),
                from=0, to=maxt, col= "red",
                add=TRUE, lwd=3
        )
        add_legend("top",
                   legend=c(epi1, epi2),
                   col=c("black", "red"),
                   lty=c(1,1),
                   lwd=3,
                   horiz=TRUE, bty='n', cex=1.2)
}

#' plot2gomp3M
#'
#' Comparison of three epidemics (gompertz model) at the same plot. (from Madden, 2007)
#'
#' @param a_1 Epidemic 1: curve assintote
#' @param a_2 Epidemic 2: curve assintote
#' @param y0_1 Epidemic 1: initial inoculum
#' @param y0_2 Epidemic 2: initial inoculum
#' @param r_1  Epidemic 1: rate of progress
#' @param r_2 Epidemic 2: rate of progress
#' @param maxt epidemic duration (days)
#' @return Three gompertz curves
#' @export
#' @examples
#' plot2gomp3M(1, 0.7, 0.001, 0.001, 0.4, 0.4, 50)

plot2gomp3M <- function(a_1, a_2, y0_1, y0_2, r_1, r_2, maxt){
        par(mar=c(5,5,5,2))
        epi1 = paste0("a = ", a_1, " y0 = ", y0_1, "   r = ", r_1 )
        epi2 = paste0("a = ", a_2, " y0 = ", y0_2, "   r = ", r_2 )

        curve(
                a_1*(exp(-(-log(y0_1))*exp(-r_1*x))),
                from=0, to=maxt, lwd=3,
                xlab='Time (days)', ylab='Disease', cex.lab=2,
                xlim=c(0,maxt), ylim=c(0,1), cex.axis=1.5
        )
        grid (NULL,NULL, lty = 6, col = "cornsilk2")

        text(maxt*0.5, 0.1, adj=0, cex=1.5,
             expression(y ~"="~ a*(exp(-(-log(y[0]))*exp(-r*x))))
        )

        curve(
                a_2*(exp(-(-log(y0_2))*exp(-r_2*x))),
                from=0, to=maxt, col= "red",
                add=TRUE, lwd=3
        )
        add_legend("top",
                   legend=c(epi1, epi2),
                   col=c("black", "red"),
                   lty=c(1,1),
                   lwd=3,
                   horiz=TRUE, bty='n', cex=1.2)
}

#' plot2mono3M
#'
#' Comparison of three epidemics (monomolecular model) at the same plot. (from Madden, 2007)
#' Adaptado de Madden, 2007 -> y = 1 - (1-y0)*exp(-r * x)
#'
#' @param a_1 Epidemic 1: assintote
#' @param a_2 Epidemic 2: assintote
#' @param y0_1 Epidemic 1: initial inoculum
#' @param y0_2 Epidemic 2: initial inoculum
#' @param r_1  Epidemic 1: rate of progress
#' @param r_2 Epidemic 2: rate of progress
#' @param maxt epidemic duration (days)
#' @return Three gompertz curves
#' @export
#' @examples
#' plot2mono3M(
#' 0.8, 0.8,  # a:  Assintote
#' 0.4, 0.1,  # y0: Initial inoculum (Intercept)
#' 0.25, 0.25,# r:  progress rate
#' 20)        # maxt: epidemic duration

plot2mono3M <- function(a_1, a_2, y0_1, y0_2, r_1, r_2, maxt){
        par(mar=c(5,5,5,2))
        epi1 = paste0("a = ", a_1, "   y0 = ", y0_1, "   r = ", r_1 )
        epi2 = paste0("a = ", a_2, "   y0 = ", y0_2, "   r = ", r_2 )

        curve(
                a_1-(a_1-y0_1)*exp(-r_1*x),
                from=0, to=maxt, lwd=3,
                xlab='Time (days)', ylab='Disease', cex.lab=2,
                xlim=c(0,maxt), ylim=c(0,1), cex.axis=1.5
        )
        grid (NULL,NULL, lty = 6, col = "cornsilk2")

        text(maxt*0.5, 0.1, adj=0, cex=1.5,
             expression(y ~"="~ a-(a-y[0])* ~~exp(-r*x))
        )
        abline(v=0, h=c(y0_1, a_1), lty=2, col=1)
        text(x=0, y=c(y0_1, a_1), c("y0", "a"), adj= c(0.5,1), col = 1)

        curve(
                a_2-(a_2-y0_2)*exp(-r_2*x),
                from=0, to=maxt, col= "red",
                add=TRUE, lwd=3
        )
        abline(h=c(y0_2, a_2), lty=2, col=2)
        text(x=0, y=c(y0_2, a_2), c("y0", "a"), adj= c(0.5,1), col = 2)

        add_legend("top",
                   legend=c(epi1, epi2),
                   col=c("black", "red"),
                   lty=c(1,1),
                   lwd=3,
                   horiz=TRUE, bty='n', cex=1.2)
}

#' plot2gomp3J
#'
#' Comparison of three epidemics (gompertz model) at the same plot.
#' Alternative parametrization
#'
#' @param a_1 Epidemic 1: curve assintote
#' @param a_2 Epidemic 2: curve assintote
#' @param y0_1 Epidemic 1: time of start of epidemic
#' @param y0_2 Epidemic 2: time of start of epidemic
#' @param r_1  Epidemic 1: rate of progress
#' @param r_2 Epidemic 2: rate of progress
#' @param maxt epidemic duration (days)
#' @return Three gompertz curves
#' @export
#' @examples
#' plot2gomp3J(1, 1, 3, 6, 0.4, 0.4, 50)


plot2gomp3J <- function(a_1, a_2, b_1, b_2, r_1, r_2, maxt){
        par(mar=c(5,5,5,2))
        epi1 = paste0("a = ", a_1, "b = ", b_1, "   r = ", r_1 )
        epi2 = paste0("a = ", a_2, "b = ", b_2, "   r = ", r_2 )

        curve(
                a_1*(exp(-exp(b_1-r_1*x))),
                #a_1*(exp(-(-log(y0_1))*exp(-r_1*x))),
                from=0, to=maxt, lwd=3,
                xlab='Time (days)', ylab='Disease', cex.lab=2,
                xlim=c(0,maxt), ylim=c(0,1), cex.axis=1.5
        )
        grid (NULL,NULL, lty = 6, col = "cornsilk2")

         curve(
                a_2*(exp(-exp(b_2-r_2*x))),
                from=0, to=maxt, col= "red",
                add=TRUE, lwd=3
        )
        add_legend("top",
                   legend=c(epi1, epi2),
                   col=c("black", "red"),
                   lty=c(1,1),
                   lwd=3,
                   horiz=TRUE, bty='n', cex=1.2)
}

