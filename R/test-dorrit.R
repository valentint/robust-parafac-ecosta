##  Dorrit data
##
##  Dorrit fluorescence data
##  Data from http://www.models.life.ku.dk/dorrit and (probably corrected) by Mia Hibert
##  https://ucphchemometrics.com/datasets/
##
##  27 116  18

library(rrcov3way)
##  load("dorrit.rda")
data("dorrit")

dim(dorrit)

plot(pp <- Parafac(dorrit, ncomp=4, robust=TRUE))
savePlot(file="dorrit.pdf", type="pdf")

##  Compare ALS against INT2 (non-robust) ========================
##
##  1. Correct number of factors
set.seed(1234)
nf <- 4
par <- CPf(dorrit, r=nf, start=0, conv=1e-6, maxit=5000)
par2 <- int2(dorrit, r=nf, initconv=1e-02, conv=1e-06, start=0,  maxit=5000)

## OK, INT2 more than 30% faster

par$cputime
##        0.72

par2$cputime
##       0.5

round(100*(par$cputime-par2$cputime)/par$cputime, 2)

## Over-factoring
nf <- 5
par <- CPf(dorrit, r=nf, start=0, conv=1e-6, maxit=5000)
par2 <- int2(dorrit, r=nf, initconv=1e-02, conv=1e-06, start=0, maxit=5000)

## Here INT2 is more than 60% faster
par$cputime
##       3.06

par2$cputime
##       1.19

round(100*(par$cputime-par2$cputime)/par$cputime, 2)

## Compare R-ALS against R-INT2 ==================================
##
##  Correct number of factors
nf <- 4
parr <- R_als(dorrit, ncomp=nf, type="als", start="svd", conv=1e-8, maxit=5000, alpha=0.75)
parr2 <- R_als(dorrit, ncomp=nf, type="int2", start="svd", initconv=1e-2, conv=1e-8, maxit=5000, alpha=0.75)

parr$cputime
##      3.38

parr2$cputime
##      3.02

round(100*(parr$cputime-parr2$cputime)/parr$cputime, 2)

## Over-factoring
nf <- 5
parr <- R_als(dorrit, ncomp=nf, type="als", start="svd", conv=1e-8, maxit=5000, alpha=0.75)
parr2 <- R_als(dorrit, ncomp=nf, type="int2", start="svd", initconv=1e-2, conv=1e-8, maxit=5000, alpha=0.75)

parr$cputime
##      3.38

parr2$cputime
##      3.02

round(100*(parr$cputime-parr2$cputime)/parr$cputime, 2)

##================================================================
##
## Set or get the status of the random generator
myrng <- function(seed) {
    if(!missing(seed)) {
        set.seed(seed)
        return(seed)
    }

    eff_seed <- floor(2^15*runif(1)) + 1
    print(sprintf("Seed for session: %s", eff_seed))
    set.seed(eff_seed)
    eff_seed
}

##  Generate Table 9: for F=4 and F=5 ==========================================
##
##  Run R-ALS and R-INT2 several times and report the average number of iterations and times. 
##  Set nf to 4 or 5 for correct factor estimation and over-factoring respectively
##  This could be done also with microbenchmark...
##
##  Compare the resulting loadings matrices e.g. from the last repetition 
##  (with a given tolerance)
##  

nf <- 5
nrep <- 10
iter_inc <- rep(NA, nrep)
ainc <- rep(NA, nrep)
iter_rals <- rep(NA, nrep)
iter_rint2 <- rep(NA, nrep)
time_rals <- rep(NA, nrep)
time_rint2 <- rep(NA, nrep)
time_inc <- NA
for(i in 1:nrep) {

    ##  seed <- myrng()
    parr  <- R_als(dorrit, ncomp=nf, type="als", start="svd", conv=1e-8, maxit=5000, alpha=0.75)
    parr2 <- R_als(dorrit, ncomp=nf, type="int2", start="svd", initconv=1e-2, conv=1e-8, maxit=5000, alpha=0.75)

    iter_rals[i] <- parr$iter
    iter_rint2[i] <- parr2$iter
    time_rals[i] <- parr$cputime
    time_rint2[i] <- parr2$cputime
    iter_inc[i]  <- round(100*(parr$iter-parr2$iter)/parr$iter, 2)
    ainc[i] <- inc <- round(100*(parr$cputime-parr2$cputime)/parr$cputime, 2)

##    cat("\ni=", i, seed, parr$cputime, parr2$cputime, inc, fill=TRUE)
    cat("\ni=", i, parr$cputime, parr2$cputime, inc, fill=TRUE)
}

cat("\nITER: R-ALS, R-INT2, Improvement (percent)", mean(iter_rals), mean(iter_rint2), mean(iter_inc), "\n")
cat("\nTIME: R-ALS, R-INT2, Improvement (percent)", mean(time_rals), mean(time_rint2), mean(ainc), "\n")

##  Compare the three loading matrices: identical if rounded to
##  the third decimal position and reflected (adjusting the sign 
##  of the components as necessary). 

if(nf == 5) {
    all.equal(parr$A %*% diag(c(-1, -1, 1, -1, 1)), parr2$A, tolerance=1e-6, check.attributes=FALSE)
    all.equal(parr$B %*% diag(c(-1, -1, 1, -1, 1)), parr2$B, tolerance=1e-6, check.attributes=FALSE)
    all.equal(parr$C %*% diag(c(1, 1, 1, 1, 1)), parr2$C, tolerance=1e-5, check.attributes=FALSE)

    save(iter_rals, iter_rint2, time_rals, time_rint2, ainc, parr, parr2, file="dorrit_R5.rda")
}else {
    all.equal(parr$A %*% diag(c(-1, -1, 1, -1)), parr2$A, tolerance=1e-2, check.attributes=FALSE)
    all.equal(parr$B %*% diag(c(-1, -1, 1, 1)), parr2$B, tolerance=1e-2, check.attributes=FALSE)
    all.equal(parr$C %*% diag(c(1, 1, 1, -1)), parr2$C, tolerance=1e-3, check.attributes=FALSE)
    
    save(iter_rals, iter_rint2, time_rals, time_rint2, ainc, parr, parr2, file="dorrit_R4.rda")
}
