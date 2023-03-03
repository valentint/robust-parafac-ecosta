##
##  sim_initall() and sim_initconv()
##
##  Present the results from the simulation for estimation of the
##  initial convergence in R-INT2.
##  - The simulation is conducted by sim_initall().
##  - The results are stored in a faile simx_initall_R3.rda or
##      simx_initall_R5.rda for R=3 or R=5
##  - The file contains a list of objects for each combination of
##      CONG=0.2, 0.3, ..., 0.7,
##      NOISE=0.15, 020 and 0.25 and
##      NOISE1=0.10, 0.15, 0.20
##  - The names of the objects are something like C2-0.15-0.1,
##      encoding the values of CONG, NOISE and NOISE1
##  - Each of this objects contains the following elements. These
##      elements are arrays of which the first dimension is the
##      simulation, i.e. i=1:nsim and the second dimension is the
##      tested initial convergence criteria: 0.1, 0.01, ..., 1e-8.
##      Thus most of the results are two-dimensional arrays 100 x 8,
##      some others have a third or higher dimension
##      - mse: nsim x 8, MSE comparing X and Xhat, removing the outliers
##      - mxB, mxC: nsim x 8, angle between the estimated subspace and the original
##          one, for loadings B and C respectively
##      - msev: nsim x 8 x 3, MSE between the loadings matrices A, B and C
##      - mseabc: nsim x 8, average MSE (mean of msev)
##      - time, iter: nsim x 8, time and number of iterations
##      - fit: nsim x 8 x 3, the three components of the fit: fit, ss and fp
##      - FC: nsim x 8 x F: factor congruence between the true and estimated
##          factors, based on the three loading matrices, FC <- FCa * FCb * FCc.
##          Good value is larger than 0.95
##      - FR: nsim x 8, Fault recoveries, min(FC) < 0.95, i.e. at least one
##          of the factor congruencies is less than 0.95
##      - TC: number of degenerate solutions/iterations based on tripple cosine,
##          i.e. solutions/iterations for which the tripple cosine
##          fell below -0.8 for more than 10 iterations.
##      - tottime: total time for execution of this cycle
##      - param: a list with all of the parameters with which the
##          simulation of this object was done
##
library(rrcov3way)
library(rrcov)

cong <- seq(0.2, 0.7, 0.1)
noise <- matrix(c(0.15, 0.10, 0.15, 0.15, 0.15, 0.20,
                  0.20, 0.10, 0.20, 0.15, 0.20, 0.20,
                  0.25, 0.10, 0.25, 0.15, 0.25, 0.20), ncol=2, byrow=TRUE)
lcong <- length(cong)
lnoise <- nrow(noise)

getInitconv <- function(simx_all, var="time", cong) {

    ## Find the number of CONG levels and noise levels
    ii <- which(!sapply(simx_all, is.null))
    cii <- names(simx_all)[ii]
    ccong <- unique(substr(cii, 1, 2))
    noise <- unique(substr(cii, 4, max(nchar(cii))))
    lcong <- length(ccong)
    lnoise <- length(noise)
    cat("\nlcong, lnoise: ", lcong, lnoise, "\n")

    incong <- 1:lcong
    if(!missing(cong))
        incong <- cong
    xout <- NULL
    for(i in incong) {
        for(j in 1:lnoise) {
            ix <- (i-1)*lnoise + j
            ##  print(ix)
            ##  print(names(simx_all)[ix])
            sx <- simx_all[[ix]]
            if(!is.null(sx)) {
                xout <- rbind(xout, sx[[var]])
            }
        }
    }

    xout
}

home <- "C:/users/valen/onedrive/myrepo/robust-parafac-ecosta/R"

nf <- 5
load(file.path(home, "Results", paste0("simx_initall_R", nf, ".rda")))
cat("\nLength: ", length(which(sapply(simx_all, FUN=function(x) !is.null(x)))), "simulations.\n")
cat("\nTotal time: ", tt <- round(sum(sapply(simx_all, FUN=function(x) ifelse(is.null(x), 0, x$tottime)))), "sec.,", round(tt/60/60, 1),"hours.\n")

oo <- getInitconv(simx_all, var="time")
windows(7,5)
boxplot(oo, ylab="Time", cex.axis=0.8, outline=FALSE, main=paste0("F=R=", nf))
colMedians(oo)
round(colMeans(oo), 4)

savePlot(file=paste0("initconv-time-R", nf, ".pdf"), type="pdf")

oo <- getInitconv(simx_all, var="iter")
boxplot(oo, ylab="Iter", cex.axis=0.8, outline=FALSE, main=paste0("F=R=", nf))
colMedians(oo)
round(colMeans(oo))

savePlot(file=paste0("initconv-iter-R", nf, ".pdf"), type="pdf")

oo <- getInitconv(simx_all, var="TC")
colSums(oo)

oo <- getInitconv(simx_all, var="FR")
colSums(oo)

## Percentage of FR in CONG=0.7, all levels of noise
oo <- getInitconv(simx_all, var="FR", cong=which(cong==0.7))
colSums(oo)
(aa <- round(100*colSums(oo)/nrow(oo),2))

library(xtable)
xtable(t(as.data.frame(aa)))
