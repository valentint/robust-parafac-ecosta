##  sim_all() and sim_RP()
##
##  sim_RP() conducts the simulation for
##      - a given number of factors (R=3 or R=5) and
##      - a given number of estimated factors (F=3, 4, 5 or 6) and
##      - a given set of CONG, NOISE and NOISE1 and
##      - a given type of outliers BL, GL or OG
##      - additionally the dimensions I, J and K can be controlled as
##          SIZEA: 100 x 20 x 20
##          SIZEB: 50  x 100 x 20
##          SIZEC: 100 x 100 x 20
##
##  sim_all() iterates through the values of CONG, NOISE and NOISE1:
##      CONG=0.3, ..., 0.7,
##      NOISE=0.15, 020 and 0.25 and
##      NOISE1=0.10, 0.15, 0.20
##
##  NOTE: the values of CONG, NOISE and NOISE1 are hard-coded in sim_all()!
##
##  - The results are stored in a file simx_all_size_RxFy_zz.rda where
##      - size = SIZEA (empty) or SIZEB or SIZEC
##      - x = 3 or 5
##      - y = 3, 4, 5 or 6
##      - zz = BL, GL or OG
##
##  - The file contains a list of objects for each combination of
##      CONG, NOISE and NOISE1
##  - The names of the objects are something like C3-0.15-0.1,
##      encoding the values of CONG, NOISE and NOISE1
##
##  - Each of this objects contains the following elements. These
##      elements are arrays of which the first dimension is the
##      simulation, i.e. i=1:nsim, the second dimension is
##      epsilon: 0.0, 0.10, 0.20 and the third fdimension is the
##      estimator: C, RALS or RINT2.
##      Thus most of the results are three-dimensional arrays 100 x 3 x 3,
##      some others have a third or higher dimension
##      - mse: nsim x 3 x 3, MSE comparing X and Xhat, removing the outliers
##      - mxB, mxC: nsim x 3 x 3, angle between the estimated subspace and the original
##          one, for loadings B and C respectively
##      - msev: nsim x 3 x 3 x 3, MSE between the loadings matrices A, B and C
##      - mseabc: nsim x 3 x 3, average MSE (mean of msev)
##      - time, iter: nsim x 3 x 3, time and number of iterations
##      - fit: nsim x 3 x 3 x 3, the three components of the fit: fit, ss and fp
##      - FC: nsim x 3 x 3 x F: factor congruence between the true and estimated
##          factors, based on the three loading matrices, FC <- FCa * FCb * FCc.
##          Good value is larger than 0.95
##      - FR: nsim x 3 x 3, Fault recoveries, min(FC) < 0.95, i.e. at least one
##          of the factor congruencies is less than 0.95
##      - TC: number of degenerate solutions/iterations based on tripple cosine,
##          i.e. solutions/iterations for which the tripple cosine
##          fell below -0.8 for more than 10 iterations.
##      - tottime: total time for execution of this cycle
##      - param: a list with all of the parameters with which the
##          simulation of this object was done
##
if(FALSE) {
    
    library(pracma)
    library(rrcov3way)
    library(rrcov)

    home <- "C:/VT/Work"

    source(file.path(home, "functions.R"))
    source(file.path(home, "sim_all.R"))

    sim_all(I=100, J=20, K=20, nsim=100, nf=3, nfe=3, conver=1e-8, outlier_type="BL")

    ##  simx=sim_RP(I=100, J=20, K=20, nsim=1, nf=3, nfe=3, conver=1e-08, start="svd", maxit=5000, outlier_type="BL", trace=FALSE)

}

###
##  Simulate for a grid of values for noise and CONG and fixed rank (nf),
##  number of factors to estimate (nfe) and outlier type
sim_all <- function(I=50, J=100, K=10, nsim=100, nf=3, nfe=nf, conver=1e-08, outlier_type=c("BL", "GL", "OG"), trace=FALSE) {

    outlier_type <- match.arg(outlier_type)

    cong <- seq(0.3, 0.7, 0.1)
    noise <- matrix(c(0.15, 0.10, 0.15, 0.15, 0.15, 0.20,
                      0.20, 0.10, 0.20, 0.15, 0.20, 0.20,
                      0.25, 0.10, 0.25, 0.15, 0.25, 0.20), ncol=2, byrow=TRUE)
    allmat <- NULL
    for(i in 1:length(cong)) {
        allmat <- rbind(allmat, cbind(rep(cong[i], nrow(noise)), noise))
    }
    lall <- nrow(allmat)

    simx_all <- vector("list", length=lall)
    names(simx_all) <- paste0("C", 10*allmat[,1], '-', allmat[,2], "-", allmat[,3])

    for(ix in 1:lall) {
        print(Sys.time())
        cat("\nCONG=", allmat[ix,1], ", NOISE=[", allmat[ix,2], ",",
            allmat[ix,3], "] - Simulation", ix, "out of", lall, "===\n")
        simx_all[[ix]] <- sim_RP(I=I, J=J, K=K, nsim=nsim, nf=nf, nfe=nfe, conver=conver,
            noise=allmat[ix,2], noise1=allmat[ix,3], cong=allmat[ix,1], outlier_type=outlier_type)
        save(simx_all, file=paste0("simx_all_R", nf, "F", nfe, "_", outlier_type, ".rda"))
    }
}

sim_RP <- function(I=50, J=100, K=10, nsim=100, nf=3, nfe=nf, conver=1e-08,
        start=c("svd", "random"), maxit=5000,
        noise=0.15, noise1=0.10, cong=0.5,
        outlier_type=c("BL", "GL", "OG", "bl", "gl", "og"),
        trace=FALSE)
{
    INITCONV <- 1e-2

    start <- match.arg(start)

    outlier_type <- tolower(match.arg(outlier_type))
    if(outlier_type == "bl") {
        c1 <- 10; c2 <- 0.1
    } else if(outlier_type == "gl") {
        c1 <- 10; c2 <- 0.0
    } else if(outlier_type == "og") {
        c1 <- 1; c2 <- 0.1
    }
    eps <- c(0, 0.1, 0.2); leps <- length(eps)
    cat("\nOutliers: ", toupper(outlier_type), "c1=", c1, "c2=", c2, "eps=", paste0(eps, collapse=","), "\n")

    methods <- c("C", "RALS", "RINT2"); lmeth <- length(methods)

    dnames <- list(1:nsim, eps, methods)
    mse <- array(NA, dim=c(nsim, leps, lmeth), dimnames=dnames)
    mxB <- array(NA, dim=c(nsim, leps, lmeth), dimnames=dnames)
    mxC <- array(NA, dim=c(nsim, leps, lmeth), dimnames=dnames)
    time <- array(NA, dim=c(nsim, leps, lmeth), dimnames=dnames)
    iter <- array(NA, dim=c(nsim, leps, lmeth), dimnames=dnames)
    fit <- array(NA, dim=c(nsim, leps, lmeth, 3), dimnames=list(1:nsim, eps, methods, c("fit", "fp", "ss")))

    FR <- array(NA, dim=c(nsim, leps, lmeth), dimnames=dnames)
    FC <- array(NA, dim=c(nsim, leps, lmeth, nf), dimnames=list(1:nsim, eps, methods, paste0("FC", 1:nf)))
    msev <- array(NA, dim=c(nsim, leps, lmeth, 3), dimnames=list(1:nsim, eps, methods, c("A", "B", "C")))
    mseabc <- array(NA, dim=c(nsim, leps, lmeth), dimnames=dnames)
    TC <- array(NA, dim=c(nsim, leps, lmeth), dimnames=dnames)

tic()

    for(j in 1:length(eps)) {

        xdat <- cp_gen(I=I, J=J, K=K, nsim=nsim, nf=nf,
            noise=noise, noise1=noise1, Acol=TRUE, Bcol=TRUE, Ccol=TRUE,
            congA=cong, congB=cong, congC=cong,
            eps=eps[j], type=outlier_type, c1=c1, c2=c2)

        for(i in 1:nsim)
        {
            cat("\nEPS=", eps[j], "Simulation", i, "out of", nsim, "\n")
            data <- list(X=xdat$Xs[[i]], A=xdat$As[[i]], B=xdat$Bs[[i]], C=xdat$Cs[[i]], iout=xdat$Os[[i]])

            ## Save the parameters to pass to the output
            if(j == 1 && i == 1) {
                param <- xdat$param
                param$nfe <- nfe
                param$start <- start
                param$eps <- eps
                param$initconv <- INITCONV
                param$conv <- conver
            }

            ## 1. Standard classical PARAFAC (ALS)
            par <- CPf(data$X, r=nfe, ort1=1, ort2=1, ort3=1, start=if(start=="svd") 0 else 1, conv=conver, maxit=maxit)
            imeth <- which(methods == "C")

            mse[i, j, imeth] <- mse2(data$X, par, iout=data$iout)
            mxB[i, j, imeth] <- mysubspace(data$B, par$B/norm(par$B, type="F"))
            mxC[i, j, imeth] <- mysubspace(data$C, par$C/norm(par$C, type="F"))
            time[i, j, imeth] <- par$cputime
            iter[i, j, imeth] <- par$iter
            fit[i, j, imeth, 1] <- par$fit
            fit[i, j, imeth, 2] <- par$fp
            fit[i, j, imeth, 3] <- par$ss

            TC[i, j, imeth] <- sum(par$ftiter[, 2] < -0.8)
            mx <- maxcong(ll0=list(A=data$A, B=data$B, C=data$C), ll=par, data$iout)
            mseabc[i, j, imeth] <- mx$MSE
            msev[i, j, imeth, ] <- mx$MSEV
            FC[i, j, imeth, ] <- mx$FC
            FR[i, j, imeth] <- min(mx$FC) < 0.95

            ## 2. Robust PARAFAC with ALS
            parr <- R_als(data$X, ncomp=nfe, type="als", start=start, conv=conver, maxit=maxit, ncomp.rpca=3, alpha=0.75)
            imeth <- which(methods == "RALS")

            mse[i, j, imeth] <- mse2(data$X, parr, iout=data$iout)
            mxB[i, j, imeth] <- mysubspace(data$B, parr$B/norm(parr$B, type="F"))
            mxC[i, j, imeth] <- mysubspace(data$C, parr$C/norm(parr$C, type="F"))
            time[i, j, imeth] <- parr$cputime
            iter[i, j, imeth] <- parr$robiter * parr$aveiter + parr$iter
            fit[i, j, imeth, 1] <- parr$fit
            fit[i, j, imeth, 2] <- parr$fp
            fit[i, j, imeth, 3] <- parr$ss

            TC[i, j, imeth] <- sum(parr$ftiter[,2] < -0.8)
            mx <- maxcong(ll0=list(A=data$A, B=data$B, C=data$C), ll=parr, data$iout)
            mseabc[i, j, imeth] <- mx$MSE
            msev[i, j, imeth, ] <- mx$MSEV
            FC[i, j, imeth, ] <- mx$FC
            FR[i, j, imeth] <- min(mx$FC) < 0.95

            ##  3. Robust PARAFAC with INT2
            parr2 <- R_als(data$X, ncomp=nfe, type="int2", start=start, initconv=INITCONV, conv=conver, maxit=maxit, ncomp.rpca=3, alpha=0.75)
            imeth <- which(methods == "RINT2")

            mse[i, j, imeth] <- mse2(data$X, parr2, iout=data$iout)
            mxB[i, j, imeth] <- mysubspace(data$B, parr2$B/norm(parr2$B, type="F"))
            mxC[i, j, imeth] <- mysubspace(data$C, parr2$C/norm(parr2$C, type="F"))
            time[i, j, imeth] <- parr2$cputime
            iter[i, j, imeth] <- parr2$robiter * parr2$aveiter + parr2$iter
            fit[i, j, imeth, 1] <- parr2$fit
            fit[i, j, imeth, 2] <- parr2$fp
            fit[i, j, imeth, 3] <- parr2$ss

            TC[i, j, imeth] <- sum(parr2$ftiter[,2] < -0.8)
            mx <- maxcong(ll0=list(A=data$A, B=data$B, C=data$C), ll=parr2, data$iout)
            mseabc[i, j, imeth] <- mx$MSE
            msev[i, j, imeth, ] <- mx$MSEV
            FC[i, j, imeth, ] <- mx$FC
            FR[i, j, imeth] <- min(mx$FC) < 0.95

       }
    }

    print(tottime <- toc())

    list(param=param, mse=mse, mxB=mxB, mxC=mxC, time=time, iter=iter, fit=fit, tottime=tottime,
        mseabc=mseabc, msev=msev, FC=FC, FR=FR, TC=TC)
}
