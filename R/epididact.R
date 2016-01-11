#' Model HIR: Healthy - Infected - Removed
#'
#' Coupled differential equation models to study disease development
#' building on examples presented in madden et al. 2007, chapter 5
#'
#' @param x a numeric
#' @return something
#' @export
#'
HIR<-function(times,y,parms) {
        library(deSolve)
        dH<- -parms["beta"]*y["H"]*y["I"]
        dI<- parms["beta"]*y["H"]*y["I"] - parms["nu"]*y["I"]
        dR<- parms["nu"]*y["I"]
        list(c(dH,dI,dR), c(N=sum(y)))
}

#' Model HLIR: Healthy - Latent - Infected - Removed
#'
#' Coupled differential equation models to study disease development
#' building on examples presented in madden et al. 2007, chapter 5
#'
#' @param times
#' @param y
#' @param parms
#' @return something
#' @export
#'
HLIR<-function(times,y,parms) {
        library(deSolve)
        dH<- -parms["beta"]*y["H"]*y["I"]
        dL<- parms["beta"]*y["H"]*y["I"] - parms["omega"]*y["L"]
        dI<- parms["omega"]*y["L"] - parms["nu"]*y["I"]
        dR<- parms["nu"]*y["I"]
        list(c(dH,dL,dI,dR),c(N=sum(y)))

        }

add_legend <- function(...) {
        opar <- par(fig=c(0, 1, 0, 1), oma=c(1, 0, 0, 0),
                    mar=c(0, 0, 1, 0), new=TRUE)
        on.exit(par(opar))
        plot(0, 0, type='n', bty='n', xaxt='n', yaxt='n', cex=1.3)
        legend(...)
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
#'

plot2logis <- function(y0_1, y0_2, r_1, r_2, maxt){
        par(mar=c(5,5,5,2))
        epi1 = paste0("y0 = ", y0_1, "   r = ", r_1 )
        epi2 = paste0("y0 = ", y0_2, "   r = ", r_2 )

        curve(
                1/(1+(1-y0_1)/y0_1*exp(-r_1*x)),
                from=0, to=maxt, lwd=3,
                xlab='Time (days)', ylab='Disease', cex.lab=2,
                xlim=c(0,maxt), ylim=c(0,1), cex.axis=1.5
        )
        grid (NULL,NULL, lty = 6, col = "cornsilk2")
        curve(
                1/(1+(1-y0_2)/y0_2*exp(-r_2*x)),
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
#'

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

#' plot_2logis_logito
#'
#' Comparison two logits rects in the same plot
#'
#' @param y0_1 Epidemic 1: initial inoculum
#' @param y0_2 Epidemic 2: initial inoculum
#' @param r_1  Epidemic 1: rate of progress
#' @param r_2 Epidemic 2: rate of progress
#' @param maxt epidemic duration (days)
#' @return Dos logits rects
#' @export
#'

plot_2logis_logito <- function(y0_1, y0_2, r_1, r_2, maxt){
        par(mar=c(5,5,5,2))

        curve(qlogis(
                1/(1+(1-y0_1)/y0_1*exp(-r_1*x))),
              from=0, to=maxt, lwd=3,
              xlab='Time (days)', ylab='Disease', cex.lab=2,
              xlim=c(0,maxt),  cex.axis=1.5
        )
        abline(h=0,lty=2)
        curve(qlogis(
                1/(1+(1-y0_2)/y0_2*exp(-r_2*x))),
              from=0, to=maxt, col= "red",
              add=TRUE, lwd=3
        )
}


#' plot3logis
#'
#' Compara dos epidemias logisticas en el mismo grafico
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
#'

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
                   horiz=FALSE, bty='n', cex=1.2)
}

#' plot3model
#'
#' Compara modelos monomolecular, logistico y gompertz en el mismo grafico
#'
#' @param y0 intensidad de inóculo inicial
#' @param r  tasa de progreso
#' @param maxt tiempo de epidemia
#' @return tres curvas L, M y G en distintos gráficos
#' @export
#'

plot3model <- function(y0,r,maxt){

        par(mfrow=c(3,1), mar=c(5,5,0,0), oma=c(0,0,2,2))

        curve(1-(1-y0)*exp(-r*x),
              from=0, to=maxt, lwd=3,
              xlab='', ylab='', cex.lab=2,
              xlim=c(0,maxt), ylim=c(0,1), cex.axis=1.5
        )
        text(0.9*maxt, 0.1,"M", adj=c(0,0), cex = 2)

        curve(1/(1+(1-y0)/y0*exp(-r*x)),
              from=0, to=maxt, lwd=3,
              xlab='', ylab='Disease (0-1)', cex.lab=2,
              xlim=c(0,maxt), ylim=c(0,1), cex.axis=1.5
        )
        text(0.9*maxt,.1,"L", adj=c(0,0), cex = 2)

        curve(exp(-(-log(y0))*exp(-r*x)),
              from=0, to=100, lwd=3,
              xlab='Time', ylab='', cex.lab=2,
              xlim=c(0, maxt), cex.axis=1.5
        )
        text(0.9*maxt, 0.10,"G", adj=c(0,0), cex = 2)
}


