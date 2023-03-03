if(FALSE) {

    library(rrcov3way)
    library(rrcov)

    home <- "C:/VT/Work"
    home <- "C:/users/valen/onedrive/myrepo/robust-parafac-ecosta/R"
    source(file.path(home, "functions.R"))
    source(file.path(home, "sim-initconv.R"))

    set.seed(12345)
    sim_initall(I=50, J=50, K=50, nsim=100, nf=3, conver=1e-8)

    set.seed(567890)
    sim_initall(I=50, J=50, K=50, nsim=100, nf=5, conver=1e-8)

    simx=sim_initconv(I=20, J=20, K=20, nsim=10, nf=3, nfe=3, conver=1e-6, noise=0.15, noise1=0.10, cong=0.5)
}

## Simulate for a grid of values for noise and CONG
sim_initall <- function(I=50, J=50, K=50, nsim=100, nf=3, nfe=nf, conver=1e-8, trace=FALSE){

    cong <- seq(0.2, 0.7, 0.1)
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
        cat("\nCONG=", allmat[ix,1], ", NOISE=[", allmat[ix,2], ",", allmat[ix,3], " - Simulation", ix, "out of", lall, " ==============================\n")
        simx_all[[ix]] <- sim_initconv(I=I, J=J, K=K, nsim=nsim, nf=nf, nfe=nfe, conver=conver, noise=allmat[ix,2], noise1=allmat[ix,3], cong=allmat[ix,1])
        save(simx_all, file=paste0("simx_initall_R", nfe, ".rda"))
    }
}

sim_initconv <- function(I=50, J=100, K=10, nsim=100, nf=3, nfe=nf, conver=1e-08, start=c("svd", "random"), maxit=5000,
        noise=0.15, noise1=0.10, cong=0.5,
        trace=FALSE)
{

    start <- match.arg(start)

    eps <- c(1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8); leps <- length(eps)

    dnames <- list(1:nsim, eps)
    mse <- array(NA, dim=c(nsim, leps), dimnames=dnames)
    mxB <- array(NA, dim=c(nsim, leps), dimnames=dnames)
    mxC <- array(NA, dim=c(nsim, leps), dimnames=dnames)
    time <- array(NA, dim=c(nsim, leps), dimnames=dnames)
    iter <- array(NA, dim=c(nsim, leps), dimnames=dnames)
    fit <- array(NA, dim=c(nsim, leps, 3), dimnames=list(1:nsim, eps, c("fit", "fp", "ss")))

    FR <- array(NA, dim=c(nsim, leps), dimnames=dnames)
    FC <- array(NA, dim=c(nsim, leps, nf), dimnames=list(1:nsim, eps, paste0("FC", 1:nf)))
    msev <- array(NA, dim=c(nsim, leps, 3), dimnames=list(1:nsim, eps, c("A", "B", "C")))
    mseabc <- array(NA, dim=c(nsim, leps), dimnames=dnames)
    TC <- array(NA, dim=c(nsim, leps), dimnames=dnames)

tic()

    for(i in 1:nsim) {

        ## Fixed 20 percent outliers
        xdat <- cp_gen(I=I, J=J, K=K, nsim=nsim, nf=nf,
            noise=noise, noise1=noise1, Acol=TRUE, Bcol=TRUE, Ccol=TRUE,
            congA=cong, congB=cong, congC=cong,
            eps=0.2, type="bl", c1=10, c2=0.1)
        data <- list(X=xdat$Xs[[i]], A=xdat$As[[i]], B=xdat$Bs[[i]], C=xdat$Cs[[i]], iout=xdat$Os[[i]])

        for(j in 1:length(eps))
        {
            cat("\nINITCONV=", eps[j], "Simulation", i, "out of", nsim, "\n")

            ##  3. Robust PARAFAC with INT2
            parr2 <- R_als(data$X, ncomp=nfe, type="int2", start=start, initconv=eps[j], conv=conver, maxit=maxit, ncomp.rpca=3, alpha=0.5)

            mse[i, j] <- mse2(data$X, parr2, iout=data$iout)
            mxB[i, j] <- mysubspace(data$B, parr2$B/norm(parr2$B, type="F"))
            mxC[i, j] <- mysubspace(data$C, parr2$C/norm(parr2$C, type="F"))
            time[i, j] <- parr2$time
            iter[i, j] <- parr2$robiter * parr2$aveiter + parr2$iter
            fit[i, j, 1] <- parr2$fit
            fit[i, j, 2] <- parr2$fp
            fit[i, j, 3] <- parr2$ss

            TC[i, j] <- sum(parr2$ftiter[,2] < -0.8)
            mx <- maxcong(ll0=list(A=data$A, B=data$B, C=data$C), ll=parr2, data$iout)
            mseabc[i, j] <- mx$MSE
            msev[i, j, ] <- mx$MSEV
            FC[i, j, ] <- mx$FC
            FR[i, j] <- min(mx$FC) < 0.95
       }
    }

    print(tottime <- toc())

    list(mse=mse, mxB=mxB, mxC=mxC, time=time, iter=iter, fit=fit, tottime=tottime, mseabc=mseabc, msev=msev, FC=FC, FR=FR, TC=TC, param=xdat$param)
}
