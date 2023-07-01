##
##  Present the results from the simulation for comparison of the estimators.
##
##  - The simulation is conducted by sim_all().
##  - The results are stored in a file simx_all_size_RxFy_zz.rda where
##      - size = SIZEA (empty) or SIZEB or SIZEC
##      - x = 3 or 5
##      - y = 3, 4, 5 or 6
##      - zz = BL, GL or OG
##
##  - The file contains a list of objects for each combination of
##      CONG=0.3, ..., 0.7,
##      NOISE=0.15, 020 and 0.25 and
##      NOISE1=0.10, 0.15, 0.20
##  - The names of the objects are something like C3-0.15-0.1,
##      encoding the values of CONG, NOISE and NOISE1
##
##  - The content of each of the objects is described in sim_all.R.
##
library(rrcov3way)
library(rrcov)

cong <- seq(0.3, 0.7, 0.1)
noise <- matrix(c(0.15, 0.10, 0.15, 0.15, 0.15, 0.20,
                  0.20, 0.10, 0.20, 0.15, 0.20, 0.20,
                  0.25, 0.10, 0.25, 0.15, 0.25, 0.20), ncol=2, byrow=TRUE)
lcong <- length(cong)
lnoise <- nrow(noise)

## Returns string without leading or trailing white space
trim <- function (x) gsub("^\\s+|\\s+$", "", x)

##  Generate the name of a data set with specific simulation results
##  - Size normally is missing and this means the standard array size is uses: 100x20x20
##      It can be size="SIZEB", for 50x100x20 and size="SIZEC" for 100x100x20
##
get_simname <- function(home, R=3, F=3, type=c("BL", "GL", "OG"),
    size=c("SIZEA", "SIZEB", "SIZEC", "sizea", "sizeb", "sizec"),
    prefix="simx_all", subdir="Results") {

    type <- match.arg(type)
    size <- toupper(match.arg(size))
    if(size == "SIZEA")
        size <- NULL
    if(!missing(size) && !is.null(size) && nchar(size) > 0)
        prefix <- paste0(prefix, "_", size)
    stopifnot(F >= R)
    fname <- file.path(home, subdir, paste0(prefix, "_R", R, "F", F, ifelse(type == "BL", "", paste0("_", type)), ".rda"))
    print(fname)
    fname
}

##  Load a data set with specific simulation results
##  - Size normally is missing and this means the standard array size is uses: 100x20x20
##      It can be size="SIZEB", for 50x100x20 and size="SIZEC" for 100x100x20
##
sim_load <- function(home, R=3, F=3, type=c("BL", "GL", "OG"), size=NULL) {
    type <- match.arg(type)
    fname <- get_simname(home, R=R, F=F, type=type, size=size)
    load(fname)

    cat("\nLength: ", length(which(sapply(simx_all, FUN=function(x) !is.null(x)))), "simulations.")
    cat("\nTotal time: ", tt <- round(sum(sapply(simx_all, FUN=function(x) ifelse(is.null(x), 0, x$tottime)))), "sec.,", round(tt/60/60, 1),"hours.\n")

    simx_all
}

xunfold <- function(x) {

   dm <- dim(x)
   dnames <- dimnames(x)
    n <- dm[1]
    m <- dm[2]
    p <- dm[3]

    mat <- matrix(NA, nrow = n, ncol = m * p)
    for(j in 1:m)
        mat[, ((j - 1) * p + 1):(j * p)] <- x[, j, ]

    rownames(mat) <- dnames[[1]]
    colnames(mat) <- rep(dnames[[3]], length(dnames[[2]]))
    mat
}

## Extract the variable 'var' from the file simx_all, all combinations
##  of CONG, NOISE and NOISE1, nsim observations in each. With the given
##  CONG and noise, this is a 4500 x 9 matrix, where 9 = 3 values of
##  epsilon x 3 estimators C, RALS and RINT2. It is possible to restrict
##  to a given set of CONG and/or noise values, e.g. cong=which(cong==0.7) will return
##  a 900 x 9 matrix (where 900=3(noise) x 3(noise1) x 100 (nsim).
##
##  NOTE: it will not work correctly for fit, msev and FC, i.e.
##  when the variable has more dimensions. For 'fit' there is a separate function getFit().
##
getAggr <- function(simx_all, var="mse", cong, noise) {
    ii <- which(!sapply(simx_all, is.null))
    cii <- names(simx_all)[ii]
    ccong <- unique(substr(cii, 1, 2))
    cnoise <- unique(substr(cii, 4, max(nchar(cii))))
    lcong <- length(ccong)
    lnoise <- length(cnoise)
    cat("\nlcong, lnoise: ", lcong, lnoise, "\n")

    xout <- NULL
    xnames <- NULL

    incong <- 1:lcong
    innoise <- 1:lnoise
    if(!missing(cong))
        incong <- cong
    if(!missing(noise))
        innoise <- noise
    for(i in incong) {
        xi <- NULL
        for(j in innoise) {
            ix <- (i-1)*lnoise + j
            sx <- simx_all[[ix]]
            if(!is.null(sx)) {
                x <- xunfold(sx[[var]])
                if(is.null(xnames)) {
                    xnames <- colnames(x)
                }
                xi <- rbind(xi, x)
            }
        }
        print(dim(xi))
        xout <- rbind(xout, xi)
    }
    colnames(xout) <- xnames
    xout
}

## this is a variant of getAggr - but different: necessary because fit is
##  4-dimensional and cannot be handled by getAggr()
getFit <- function(simx_all, fit=c("fit", "fp", "ss"), cong) {

    fit <- match.arg(fit)

    ii <- which(!sapply(simx_all, is.null))
    cii <- names(simx_all)[ii]
    ccong <- unique(substr(cii, 1, 2))
    noise <- unique(substr(cii, 4, max(nchar(cii))))
    lcong <- length(ccong)
    lnoise <- length(noise)
    cat("\nlcong, lnoise: ", lcong, lnoise, "\n")

    xout <- NULL
    xnames <- NULL

    incong <- 1:lcong
    if(!missing(cong))
        incong <- cong
    for(i in incong) {
        xi <- NULL
        for(j in 1:lnoise) {
            ix <- (i-1)*lnoise + j
            sx <- simx_all[[ix]]
            if(!is.null(sx)) {
                xx <- sx[["fit"]][,,,fit]
                x <- xunfold(xx)
                if(is.null(xnames)) {
                    xnames <- colnames(x)
                }
                xi <- rbind(xi, x)
            }
        }
        print(dim(xi))
        xout <- rbind(xout, xi)
    }
    colnames(xout) <- xnames
    xout
}

##================================================================
home <- "C:/users/valen/onedrive/current/conferences-2022/EcoSta/Revision-1/R"

nf0   <- 3          # 3 or 5
nfe0  <- 3          # 3 or 4 or 6
type0 <- "BL"       # BL, GL, OG
##              R3F3, R3F4, R5F5, R5F6; R3F3_GL, R3F4_GL; R3F3_OG, R3F4_OG
size0 <- "SIZEB"    # NULL, SIZEB, SIZEC

simx_all <- sim_load(home, R=nf0, F=nfe0, type=type0, size=size0)

nf <- simx_all[[1]]$param$nf
nfe <- simx_all[[1]]$param$nfe
type <- toupper(simx_all[[1]]$param$type)
ctype <- ifelse(type == "BL", "", paste0("_", type))
size <- ifelse(simx_all[[1]]$param$J == 20, "SIZEA", ifelse(simx_all[[1]]$param$I==50 && simx_all[[1]]$param$J==100, "SIZEB", "SIZEC"))
csize <- ifelse(size=="SIZEA", "",paste0("_", size))
(main=paste0("F=", ifelse(nf==nfe, "R=", "R+1="), nfe))

##  Table 0
## FIT
oo <- getFit(simx_all, fit="fp")
fpx00 <- oo[,2]-oo[,3]        ## R-ALS - R-INT2 at 0% BL
fpx10 <- oo[,5]-oo[,6]        ## R-ALS - R-INT2 at 10% BL
fpx20 <- oo[,8]-oo[,9]        ## R-ALS - R-INT2 at 20% BL
nfpx00 <- length(which(abs(fpx00) >= 1e-4))
nfpx10 <- length(which(abs(fpx10) >= 1e-4))
nfpx20 <- length(which(abs(fpx20) >= 1e-4))
cat("\nWhere the difference in fit is larger than 1e-4?", nfpx20, "cases out of", nrow(oo),": ", round(100*nfpx20/nrow(oo), 2), " percent.\n")
cat("\nIn", round(100*length(which(fpx20 < 0))/nrow(oo), 1), "percent of the cases R-INT2 records better fit.\n")
fpx <- data.frame(V_0=c(100*nfpx00/nrow(oo), 100*length(which(fpx00 < 0))/nrow(oo)),
                  V_10=c(100*nfpx10/nrow(oo), 100*length(which(fpx10 < 0))/nrow(oo)),
                  V_20=c(100*nfpx20/nrow(oo), 100*length(which(fpx20 < 0))/nrow(oo)))
t(fpx)

if(FALSE) {
    ##  This is cumbersome - collect the FIT results from six files:
    ##  BL, GL and OG, by R3F3 and R3F4, and put them
    ##  into a data frame to create a Latex table.

    ##fpx_BL <- t(fpx)
    ##fpx_GL <- t(fpx)
    ##fpx_OG <- t(fpx)

    ##fpxR3F3 <- cbind(fpx_BL, fpx_GL)
    ##fpxR3F3 <- cbind(fpxR3F3, fpx_OG)

    ##fpxR3F4 <- cbind(fpx_BL, fpx_GL)
    ##fpxR3F4 <- cbind(fpxR3F4, fpx_OG)
    fpxR3F3
    fpxR3F4

    fpx <- rbind(fpxR3F3, fpxR3F4)
    fpx

    library(xtable)
    xtable(fpx, digits=1)
}

##  Figure 3 ===
## MSE
windows(7, 5)
oo <- getAggr(simx_all, var="mse")
boxplot(10000*oo, ylab="MSE", cex.axis=0.8, outline=FALSE, main=main)
mtext("x 1e-4", at=0.5)
mtext("no contamination", at=2, side=1, line=3)
mtext(paste0("10% ", type, " points"), at=5, side=1, line=3)
mtext(paste0("20% ", type, " points"), at=8, side=1, line=3)
round(colMedians(10000*oo), 3)
savePlot(file=paste0("MSE-R", nf, "F", nfe, ctype, csize, ".pdf"), type="pdf")

##  Figure 4 ===
##
##  Fix eps, noise (for the left panel) and cong (for the right panel)
eps <- 0.2
if(eps == 0.2) {
    ieps <- 3       # eps=0.2
    ceps <- "20"    # will be used in the file name, e.g. _BL20_
} else if(eps == 0.1) {
    ieps <- 2       # eps=0.1
    ceps <- "10"    # will be used in the file name, e.g. _BL10_
} else {
    ieps <- 1       # eps=0.0
    ceps <- "00"    # will be used in the file name, e.g. _BL00_ or _GL00_
}

acong <- seq(0.3, 0.7, 0.1)
anoise <- matrix(c(0.15, 0.10, 0.15, 0.15, 0.15, 0.20,
                  0.20, 0.10, 0.20, 0.15, 0.20, 0.20,
                  0.25, 0.10, 0.25, 0.15, 0.25, 0.20), ncol=2, byrow=TRUE)
allmat <- NULL
for(i in 1:length(acong)) {
    allmat <- rbind(allmat, cbind(rep(acong[i], nrow(anoise)), anoise))
}
lnoise <- nrow(anoise)
lcong <- length(acong)
lall <- nrow(allmat)

##  Figure 4, left panel
##  MSE: Aggregation by CONG (0.3, ..., 0.7)
param_noise=5       # Choose the noise combination: 20%-15
msecong <- NULL
for(i in 1:lcong) {
    msei <- NULL
    for(j in param_noise) {
        ix <- (i-1)*lnoise + j
        ##  print(names(simx_all)[ix])
        sx <- simx_all[[ix]]$mse[,ieps,]
        msei <- rbind(msei, sx)
    }
    ##  print(dim(msei))
    msecong <- cbind(msecong, msei)
}

nf <- simx_all[[1]]$param$nf
nfe <- simx_all[[1]]$param$nfe
type <- toupper(simx_all[[1]]$param$type)

windows(7,5)

mainFIG4_left <- paste0("F=", ifelse(nf==nfe, "R", "R+1"), "=", nfe, "; HO=", 100*anoise[param_noise, 1], "%", ", HE=",  100*anoise[param_noise, 2], "%")
xlab <- trim(paste0("CONG=", acong, "         ", collapse=""))
boxplot(msecong*10000, las=2, ylab="MSE", cex.axis=0.8, cex.lab=0.8, xlab=xlab, main=mainFIG4_left)
mtext("x 1e-4", at=1)
savePlot(file=paste0("MSE-R", nf, "F", nfe, "_CONG_", type, ceps, csize, ".pdf"), type="pdf")

##  Figure 4, right panel
##  MSE: Aggregation by Noise (15-10, ..., 25-20)
param_cong <- 3     ## Choose CONG=0.5
msenoise <- NULL
for(i in 1:lnoise) {
    msei <- NULL
    for(j in param_cong) {
        ix <- (j-1)*lnoise + i
        ##  print(ix)
        ##  print(names(simx_all)[ix])
        sx <- simx_all[[ix]]$mse[,ieps,]
        msei <- rbind(msei, sx)
    }
    ## print(dim(msei))
    msenoise <- cbind(msenoise, msei)
}

mainFIG4_right <- paste0("F=", ifelse(nf==nfe, "R", "R+1"), "=", nfe, "; CONG=", acong[param_cong])
xlab <- paste0("   ", 100*anoise[,1], "-", 100*anoise[,2], "   ", collapse="")
boxplot(msenoise*10000, las=2, ylab="MSE", cex.axis=0.8, cex.lab=0.8, xlab=xlab, main=mainFIG4_right)
mtext("x 1e-4", at=1)
mtext("HO-HE:", side=1, at=-1, line=3, cex=0.8)
savePlot(file=paste0("MSE-R", nf, "F", nfe, "_NOISE_", type, ceps, csize, ".pdf"), type="pdf")

##  Figure 5 ===
## Angle B subspace
windows(7, 5)
oo <- getAggr(simx_all, var="mxB")
boxplot(oo, ylab="B-sub", cex.axis=0.8, outline=FALSE, main=main)
mtext("no contamination", at=2, side=1, line=3)
mtext(paste0("10% ", type, " points"), at=5, side=1, line=3)
mtext(paste0("20% ", type, " points"), at=8, side=1, line=3)
round(colMedians(oo), 5)
savePlot(file=paste0("mxB-R", nf, "F", nfe, ctype, csize, ".pdf"), type="pdf")

##  Figure 6 ===
## Angle C subspace
windows(7, 5)
oo <- getAggr(simx_all, var="mxC")
boxplot(oo, ylab="C-sub", cex.axis=0.8, outline=FALSE, main=main)
mtext("no contamination", at=2, side=1, line=3)
mtext(paste0("10% ", type, " points"), at=5, side=1, line=3)
mtext(paste0("20% ", type, " points"), at=8, side=1, line=3)
round(colMedians(oo), 5)
savePlot(file=paste0("mxC-R", nf, "F", nfe, ctype, csize, ".pdf"), type="pdf")

##  Figure 7, upper panel ===
## Time
windows(5, 5)
oo <- getAggr(simx_all, var="time")
colMedians(oo)
oo <- oo[,c(8,9)]
boxplot(oo, ylab="Time", cex.axis=0.8, outline=FALSE, main=main)
savePlot(file=paste0("time-R", nf, "F", nfe, ctype, csize, ".pdf"), type="pdf")

##  Figure 7, lower panel ===
## Iter
windows(5, 5)
oo <- getAggr(simx_all, var="iter")
colMedians(oo)
oo <- oo[,c(8,9)]
boxplot(oo, ylab="Iter", cex.axis=0.8, outline=FALSE, main=main)
savePlot(file=paste0("iter-R", nf, "F", nfe, ctype, csize, ".pdf"), type="pdf")

##  Figure 8 ===
## Time by level of contamination
windows(7, 5)
oo <- getAggr(simx_all, var="time")
colMedians(oo)
boxplot(oo, ylab="Time", cex.axis=0.8, outline=FALSE, main=main)
mtext("no contamination", at=2, side=1, line=3)
mtext(paste0("10% ", type, " points"), at=5, side=1, line=3)
mtext(paste0("20% ", type, " points"), at=8, side=1, line=3)
savePlot(file=paste0("time-eps-R", nf, "F", nfe, ctype, csize, ".pdf"), type="pdf")

##  Table 7: elements

## FR
oo <- getAggr(simx_all, var="FR")
(mm <- colSums(oo))
fr <- round(matrix(100*colSums(oo)/nrow(oo),ncol=3, byrow=TRUE), 2)
colnames(fr) <- names(mm)[1:3]
rownames(fr) <- c("0%", "10%", "20%")
fr

## TC
oo <- getAggr(simx_all, var="TC")
(mm <- colSums(oo))
tc <- matrix(apply(oo, 2, FUN=function(x) length(which(x>0))), nrow=3, byrow=TRUE)
colnames(tc) <- names(mm)[1:3]
rownames(tc) <- c("0%", "10%", "20%")
tc

oo <- getAggr(simx_all, var="time")
(mm <- colMedians(oo))
tt <- round(matrix(colMedians(oo), ncol=3, byrow=TRUE), 2)
colnames(tt) <- names(mm)[1:3]
rownames(tt) <- c("0%", "10%", "20%")
tt

##================================================================
##  Table 7
##
##  FR and TC by rank and number of factors: for different levels of contamination
##
type <- "BL"
size <- "SIZEB"    # "SIZEB", "SIZEC"

## FR and TC - F=R=3
simx_all <- sim_load(home, R=3, F=3, type=type, size=size)

oo <- getAggr(simx_all, var="FR")
colSums(oo)
fr <- round(100*colSums(oo)/nrow(oo), 1)
f33 <- matrix(fr, ncol=3, byrow=TRUE)
colnames(f33) <- names(mm)[1:3]
rownames(f33) <- c("0%", "10%", "20%")

oo <- getAggr(simx_all, var="TC")
colSums(oo)
tc33 <- matrix(apply(oo, 2, FUN=function(x) length(which(x>0))), ncol=3, byrow=TRUE)
dimnames(tc33) <- dimnames(f33)

oo <- getAggr(simx_all, var="time")
colMedians(oo)
tt33 <- round(matrix(colMedians(oo), ncol=3, byrow=TRUE), 2)
dimnames(tt33) <- dimnames(f33)

## FR and TC - F=R+1=4
simx_all <- sim_load(home, R=3, F=4, type=type, size=size)

oo <- getAggr(simx_all, var="FR")
colSums(oo)
fr <- round(100*colSums(oo)/nrow(oo), 1)
f34 <- matrix(fr, ncol=3, byrow=TRUE)
dimnames(f34) <- dimnames(f33)

oo <- getAggr(simx_all, var="TC")
colSums(oo)
tc34 <- matrix(apply(oo, 2, FUN=function(x) length(which(x>0))), ncol=3, byrow=TRUE)
dimnames(tc34) <- dimnames(f33)

oo <- getAggr(simx_all, var="time")
colMedians(oo)
tt34 <- round(matrix(colMedians(oo), ncol=3, byrow=TRUE), 2)
dimnames(tt34) <- dimnames(f33)

## FR and TC - F=R=5
simx_all <- sim_load(home, R=5, F=5, type=type, size=size)

oo <- getAggr(simx_all, var="FR")
colSums(oo)
fr <- round(100*colSums(oo)/nrow(oo), 1)
f55 <- matrix(fr, ncol=3, byrow=TRUE)
dimnames(f55) <- dimnames(f33)

oo <- getAggr(simx_all, var="TC")
colSums(oo)
tc55 <- matrix(apply(oo, 2, FUN=function(x) length(which(x>0))), ncol=3, byrow=TRUE)
dimnames(tc55) <- dimnames(f33)

oo <- getAggr(simx_all, var="time")
colMedians(oo)
tt55 <- round(matrix(colMedians(oo), ncol=3, byrow=TRUE), 2)
dimnames(tt55) <- dimnames(f33)

## FR and TC - F=R+1=6
simx_all <- sim_load(home, R=5, F=6, type=type, size=size)

oo <- getAggr(simx_all, var="FR")
colSums(oo)
fr <- round(100*colSums(oo)/nrow(oo), 1)
f56 <- matrix(fr, ncol=3, byrow=TRUE)
dimnames(f56) <- dimnames(f33)

oo <- getAggr(simx_all, var="TC")
colSums(oo)
tc56 <- matrix(apply(oo, 2, FUN=function(x) length(which(x>0))), ncol=3, byrow=TRUE)
dimnames(tc56) <- dimnames(f33)

oo <- getAggr(simx_all, var="time")
colMedians(oo)
tt56 <- round(matrix(colMedians(oo), ncol=3, byrow=TRUE), 2)
dimnames(tt56) <- dimnames(f33)

frtab3 <- round(cbind(f33, f34), 1)
tctab3 <- round(cbind(tc33, tc34))
tttab3 <- round(cbind(tt33, tt34), 2)
frtab5 <- round(cbind(f55, f56), 1)
tctab5 <- round(cbind(tc55, tc56))
tttab5 <- round(cbind(tt55, tt56), 2)

library(xtable)
xtable(frtab3, digits=1)
xtable(tctab3, digits=0)
xtable(frtab5, digits=1)
xtable(tctab5, digits=0)

xtable(tttab3, digits=2)
xtable(tttab5, digits=2)

##================================================================
##
##  Table 8
##
##  Time and FR for good leverage points.
##  Actually, this is identical to Table 7 - just take the time and
##  FR for GL and R=3.
##

type <- "GL"
size <- NULL    # "SIZEB", "SIZEC"

simx_all <- sim_load(home, R=3, F=3, type=type, size=size)

## Time - F=R=3
oo <- getAggr(simx_all, var="time")
mm <- colMedians(oo)
t33 <- round(matrix(colMedians(oo), ncol=3, byrow=TRUE), 2)
colnames(t33) <- names(mm)[1:3]
rownames(t33) <- c("0%", "10%", "20%")

## FR - F=R=3
oo <- getAggr(simx_all, var="FR")
colSums(oo)
fr <- round(100*colSums(oo)/nrow(oo),2)
f33 <- matrix(fr, ncol=3, byrow=TRUE)
dimnames(f33) <- dimnames(t33)

simx_all <- sim_load(home, R=3, F=4, type=type, size=size)

## Time - F=R+1=4
oo <- getAggr(simx_all, var="time")
mm <- colMedians(oo)
t34 <- round(matrix(colMedians(oo), ncol=3, byrow=TRUE), 2)
dimnames(t34) <- dimnames(t33)

## FR - F=R+1=4
oo <- getAggr(simx_all, var="FR")
colSums(oo)
fr <- round(100*colSums(oo)/nrow(oo),2)
f34 <- matrix(fr, ncol=3, byrow=TRUE)
dimnames(f34) <- dimnames(t33)

xt3 <- cbind(t33, t34)
xf3 <- cbind(f33, f34)

gltab <- rbind(cbind(t33, t34), cbind(f33, f34))

library(xtable)
xtable(gltab)


##================================================================
##  Figure 9
##
##  This will not work for files which are not completed, like
##  simx_all_R3F4_GL.rda or simx_all_R5F6_GL.rda - due to the very
##  long computational time not all combinations of noise and CONG
##  were computed.
##
##  Also, it will not work for the larger tensors (50x100x20) because only
##  the cases 3/3 and 3/4 were computed.

##  Fix eps and type of outliers
type <- "GL"
size <- NULL
eps <- 0.2
if(eps == 0.2) {
    ieps <- 8       # eps=0.2, take columns 8 and 9
    ceps <- "20"    # will be used in the file name, e.g. _BL20_
} else if(eps == 0.1) {
    ieps <- 5       # eps=0.1, take columns 5 and 6
    ceps <- "10"    # will be used in the file name, e.g. _BL10_
} else {
    ieps <- 2       # eps=0.0, take columns 2 and 3
    ceps <- "00"    # will be used in the file name, e.g. _BL00_ or _GL00_
}

##  Time by CONG for R=3 and R=5
simx_all <- sim_load(home, R=3, F=3, type=type, size=size)
oo1 <- getAggr(simx_all, var="time", cong=1); oo1 <- oo1[, ieps:(ieps+1)]
oo2 <- getAggr(simx_all, var="time", cong=3); oo2 <- oo2[, ieps:(ieps+1)]
oo3 <- getAggr(simx_all, var="time", cong=5); oo3 <- oo3[, ieps:(ieps+1)]

simx_all <- sim_load(home, R=5, F=5, type=type, size=size)
bb1 <- getAggr(simx_all, var="time", cong=1); bb1 <- bb1[, ieps:(ieps+1)]
bb2 <- getAggr(simx_all, var="time", cong=3); bb2 <- bb2[, ieps:(ieps+1)]
bb3 <- getAggr(simx_all, var="time", cong=5); bb3 <- bb3[, ieps:(ieps+1)]

nf <- simx_all[[1]]$param$nf
nfe <- simx_all[[1]]$param$nfe

oldpar <- par(mfrow=c(3,1), mar=c(3, 4, 0, 0) + 0.1)
boxplot(cbind(oo1, bb1), ylab="Time, CONG=0.3", outline=FALSE)
boxplot(cbind(oo2, bb2), ylab="Time, CONG=0.5", outline=FALSE)
par(mar=c(4, 4, 0, 0) + 0.1)
boxplot(cbind(oo3, bb3), ylab="Time, CONG=0.7", outline=FALSE)
mtext("F=R=3", side=1, at=1.5, line=3)
mtext("F=R=5", side=1, at=3.5, line=3)
par(oldpar)

savePlot(file=paste0("time-CONG-", type, ceps, "_R.pdf"), type="pdf")

##  Time by CONG for R=3+1 and R=5+1
simx_all <- sim_load(home, R=3, F=4, type=type, size=size)
if(length(which(sapply(simx_all, FUN=function(x) !is.null(x)))) == length(simx_all)) {
    oo1 <- getAggr(simx_all, var="time", cong=1); oo1 <- oo1[, ieps:(ieps+1)]
    oo2 <- getAggr(simx_all, var="time", cong=3); oo2 <- oo2[, ieps:(ieps+1)]
    oo3 <- getAggr(simx_all, var="time", cong=5); oo3 <- oo3[, ieps:(ieps+1)]

    simx_all <- sim_load(home, R=5, F=6, type=type, size=size)
    if(length(which(sapply(simx_all, FUN=function(x) !is.null(x)))) == length(simx_all)) {
        bb1 <- getAggr(simx_all, var="time", cong=1); bb1 <- bb1[, ieps:(ieps+1)]
        bb2 <- getAggr(simx_all, var="time", cong=3); bb2 <- bb2[, ieps:(ieps+1)]
        bb3 <- getAggr(simx_all, var="time", cong=5); bb3 <- bb3[, ieps:(ieps+1)]

        oldpar <- par(mfrow=c(3,1), mar=c(3, 2, 0, 0) + 0.1)
        boxplot(cbind(oo1, bb1), outline=FALSE)
        boxplot(cbind(oo2, bb2), outline=FALSE)
        par(mar=c(4, 2, 0, 0) + 0.1)
        boxplot(cbind(oo3, bb3),  outline=FALSE)
        mtext("F=R+1=4", side=1, at=1.5, line=3)
        mtext("F=R+1=6", side=1, at=3.5, line=3)
        par(oldpar)

        savePlot(file=paste0("time-CONG-", type, ceps, "_Rplus1.pdf"), type="pdf")
    }
}
