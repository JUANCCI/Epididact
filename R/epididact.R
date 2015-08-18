#' MODEL HIR, HEALTHY-INFECTED-REMOVED
#'
#'
#' COUPLED DIFFERENTIAL EQUATION MODELS TO STUDY DISEASE DEVELOPMENT
#' BUILDING ON EXAMPLES PRESENTED IN MADDEN ET AL. 2007, CHAPTER 5
#'
#' @param x a numeric
#' @return something
#' @export
#'
HIR<-function(times,y,parms) {
        dH<- -parms["beta"]*y["H"]*y["I"]
        dI<- parms["beta"]*y["H"]*y["I"] - parms["nu"]*y["I"]
        dR<- parms["nu"]*y["I"]
        list(c(dH,dI,dR), c(N=sum(y)))
}

HLIR<-function(times,y,parms) {
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

# Compara dos epidemias logisticas en el mismo grafico

plot2logis <- function(y0_1, y0_2, r_1, r_2, maxt){
        par(mar=c(5,5,5,2))
        epi1 = paste0("y0 = ", y0_1, "   r = ", r_1 )
        epi2 = paste0("y0 = ", y0_2, "   r = ", r_2 )

        curve(
                1/(1+(1-y0_1)/y0_1*exp(-r_1*x)),
                from=0, to=maxt, lwd=3,
                xlab='Tiempo (dias)', ylab='Enfermedad (0»1)', cex.lab=2,
                xlim=c(0,maxt), ylim=c(0,1), cex.axis=1.5,
        )

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

# plot2logis(0.001, 0.001, 0.2,0.4, 30)


# Compara dos epidemias monomolecular en el mismo grafico

plot2mono <- function(y0_1, y0_2, r_1, r_2, maxt){
        par(mar=c(5,5,5,2))
        epi1 = paste0("y0 = ", y0_1, "   r = ", r_1 )
        epi2 = paste0("y0 = ", y0_2, "   r = ", r_2 )

        curve(
                1-(1-y0_1)*exp(-r_1*x),
                from=0, to=maxt, lwd=3,
                xlab='Tiempo (dias)', ylab='Enfermedad (0»1)', cex.lab=2,
                xlim=c(0,maxt), ylim=c(0,1), cex.axis=1.5,
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

# plot2mono(0.1, 0.01, 0.03, 0.03, 100)

# Plot 2 rectas logito

plot_2logis_logito <- function(y0_1, y0_2, r_1, r_2, maxt){
        par(mar=c(5,5,5,2))

        curve(qlogis(
                1/(1+(1-y0_1)/y0_1*exp(-r_1*x))),
              from=0, to=maxt, lwd=3,
              xlab='Tiempo (dias)', ylab='Enfermedad (0»1)', cex.lab=2,
              xlim=c(0,maxt),  cex.axis=1.5,
        )
        abline(h=0,lty=2)
        curve(qlogis(
                1/(1+(1-y0_2)/y0_2*exp(-r_2*x))),
              from=0, to=maxt, col= "red",
              add=TRUE, lwd=3
        )
}


# plot_2logis_logito(0.001, 0.002, 0.4, 0.4, 30)

# Graficar monomolecular, logistico y gompertz

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
              xlim=c(0,maxt), ylim=c(0,1), cex.axis=1.5,
        )
#        abline(h=0.5, lty=2)
        text(0.9*maxt,.1,"L", adj=c(0,0), cex = 2)

        curve(exp(-(-log(y0))*exp(-r*x)),
              from=0, to=100, lwd=3,
              xlab='Time', ylab='', cex.lab=2,
              xlim=c(0, maxt), cex.axis=1.5
        )
        text(0.9*maxt, 0.10,"G", adj=c(0,0), cex = 2)
#        abline(h=0.36, lty=2)
}

# plot3model(0.01,0.2,30)
