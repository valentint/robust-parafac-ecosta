##  do_scaclef():   Scale all three modes to have the same norm of each factor
##  do_inv():       Find the general inverse of a (rectangular) matrix a
##  maxcong():      Calculate the factor congruence between two PARAFAC solutions
##  mse2():         MSE between X and Xhat excluding outliers
##  mysubspace():   A variant of subspace() from pracma: If one of the matrices
##                  has more columns, remove the extra columns.
##
##  cp_gen():       Generate simulation data from PARAFAC model with homoscedastic
##                  and heteroscedastic noise, congruence and added outliers
##
##  CPf():          Standard ALS for CP - follows CPfunc from ThreeWay
##  CP_ALS():       Standard ALS for CP - calls rrcov3way:::cp_als()
##
##  atld():
##  int2:
##  cp_int2():
##  R-als():
##
##===============================================================================

library(rrcov)
library(rrcov3way)
library(pracma)

## Scale all three modes to have the same norm of each factor
if(FALSE) {
    library(rrcov3way)

    data(elind)
    cp <- Parafac(elind)
    sf1 <- do_scalef(cp$A, cp$B, cp$C)

    ##  The matrices A, B and C have the norm of all factors equal to 1
    norm_vec <- function(x) sqrt(sum(x^2))
    apply(sf1$A, 2, norm_vec)
    apply(sf1$B, 2, norm_vec)
    apply(sf1$C, 2, norm_vec)

    ##  Call the function do_scalef() once again on the already
    ##  scaled matrices - they remain the same of course.
    sf2 <- do_scalef(cp$A, cp$B, cp$C)
    all.equal(sf1, sf2)

}
do_scalef <- function(A, B, C) {
    q3 <- (diag(sqrt(diag(t(A) %*% A))) *
           diag(sqrt(diag(t(B) %*% B))) *
           diag(sqrt(diag(t(C) %*% C)))) ^ (1/3)

    C <- C %*% diag(1/sqrt(diag(t(C) %*% C))) %*% q3
    A <- A %*% diag(1/sqrt(diag(t(A) %*% A))) %*% q3
    B <- B %*% diag(1/sqrt(diag(t(B) %*% B))) %*% q3

    list(A=A, B=B, C=C)
}

##  Find the general inverse of a (rectangular) matrix a, i.e. find
##  a matrix G such that AGA == A.
##  First try to find   t(A %*% solve(crossprod(A))) and if it does not work,
##  calculate the general inverse using the function pinv() from package pracma.
##
if(FALSE) {
    library(rrcov3way)
    library(pracma)

    data(elind)
    cp <- Parafac(elind)
    A <- cp$A
    x1 <- t(A %*% solve(crossprod(A)))
    x2 <- pinv(A, 1e-12)
    x3 <- do_inv(A)
    all.equal(A %*% x1 %*% A, A)
    all.equal(A %*% x2 %*% A, A)
    all.equal(A %*% x3 %*% A, A)

}

do_inv <- function(a, tolerance=1e-12) {
    tryCatch(expr={t(a %*% solve(crossprod(a)))},
             error=function(msg) pracma::pinv(a, tol=tolerance))
}

##
##  Calculate the factor congruence between two PARAFAC solutions (two sets of
##  loading matrices A, B and C), by first finding the best permutation of the
##  target solution 'll' iterating through all posible permutations.
##
##  Then, rescale the loading matrices to equal norm of all factors, adjust the
##  sign of the factors in such way that the congruence of
##  each pair of factors is positive and calculate the MSE as the average of
##  the MSE for A, B and C.
##
if(FALSE) {                 # TESTING=============================

    set.seed(1234)
    xlist <- cp_gen(nsim=1, nf=2, noise=0.1, noise1=0)
    X <- xlist$Xs[[1]]

    cpint2 <- int2(X, r=3, initconv=1e-01, conv=1e-08, start=1, ort1=1, ort2=1, ort3=1, maxit=5000)
    cat("\nINT2 (f, fp, iter, tripcos, congruence[A,B,C]):", cpint2$f, cpint2$fp, cpint2$iter, cpint2$tripcos,
        round(diag(congruence(cpint2$A, xlist$As[[1]])), 2),
        round(diag(congruence(cpint2$B, xlist$Bs[[1]])), 2),
        round(diag(congruence(cpint2$C, xlist$Cs[[1]])), 2), "\n")

    (fc <- maxcong(ll0=list(A=xlist$As[[1]], B=xlist$Bs[[1]], C=xlist$Cs[[1]]),
                    ll=list(A=cpint2$A, B=cpint2$B, C=cpint2$C)))

    (fc <- maxcong(ll0=list(A=xlist$As[[1]], B=xlist$Bs[[1]], C=xlist$Cs[[1]]), cpint2))

}

maxcong <- function(ll0, ll, iout) {
    ## Generate all possible permutations of the elements of the vector v
    ##  and return a matrix each row of which contains  a permutation.
    perm <- function(v) {
        n <- length(v)
        if(n == 1) v
        else {
            X <- NULL
            for (i in 1:n)
                X <- rbind(X, cbind(v[i], perm(v[-i])))
            X
        }
    }

    # If ll contains not just the trhee matrices - probably a PARAFAC solution
    if(length(ll) > 3) {
        llx <- list(A=ll$A, B=ll$B, C=ll$C)
        ll <- llx
    }
    if(is.null(ll$A) || is.null(ll$B) || is.null(ll$C))
        stop("Invalid target: must contain matrices A, B and C")

    ## Remove the rows in A corresponding to the outliers in iout, if any
    if(!missing(iout) && !is.null(iout) && length(iout) > 0) {
        ll0$A <- ll0$A[-iout,]
        ll$A <- ll$A[-iout,]
    }

    n1 <- ncol(ll0$A)                   # n is the number of components
    n2 <- ncol(ll$A)                    # n is the number of components
    n <- max(n1, n2)
    P <- perm(1:n)                      # P is the permutation matrix
    perm_s <- perm_i <- 0
    for(i in 1:nrow(P))                 # all permutations
    {
        ## Calculate the sum of the (absolute values of the) congruences
        ##  between each pair of components in each of the trhee matrices
        ##  A, B and C.
        ##  Calculate the sum for the three matrices.
        s <- 0
        for(j in 1:3) {
            d <- abs(diag(congruence(ll0[[j]], ll[[j]][,P[i,]])))
            s <- s + sum(d)
            ##  cat("\n Perm, Fact:", i, j, d, sum(d), s, "\n")
        }

        ##  cat("\n --> Perm:", i, s, "\n")
        if(s > perm_s) {
            perm_s <- s         # the largest sum of congruences
            perm_i <- i         # the permuattion with the largest sum
        }
    }

    ## Choose the permutation maximizing the congruence
    ##  Calculate the factor congruence, based on the three loading matrices:
    ##      FC <- FCa * FCb * FCc
    ##  Correct the sign
    nmin <- min(n1, n2)
    perm_fc <- rep(1, nmin)
    for(i in 1:3) {
        ll[[i]] <- ll[[i]][,P[perm_i,]]
        d <- diag(congruence(ll0[[i]], ll[[i]]))
        perm_fc <- perm_fc * abs(d)

        rs <- sign(d)

##        cat("\ni, rs: ", i, rs,"\n")
##        print(d)

        ll[[i]][,1:nmin] <- ll[[i]][,1:nmin] %*% diag(rs)

##        print(ll[[i]])

    }
    names(perm_fc) <- paste0("F", 1:nmin)

    ## Rescale both solutions to equal norm of all factors
    ll0 <- do_scalef(ll0$A, ll0$B, ll0$C)
    ll <- do_scalef(ll$A, ll$B, ll$C)

    ##  Calculate the MSE for each pair of loadings matrices and find
    ##  the average
    msev <- rep(NA, 3)
    for(i in 1:3)
        msev[i] <- norm(ll[[i]][,1:nmin] - ll0[[i]][, 1:nmin], type="F")^2 / (nrow(ll[[i]])*nmin)

    mse <- mean(msev, na.rm=TRUE)

    list(FC=perm_fc, MSE=mse, MSEV=msev)
}

mse2 <- function(X, A, B, C, iout)
{
    if(is.list(A)) {
        ll <- A
        A <- ll$A
        B <- ll$B
        C <- ll$C
        rm(ll)
    }

    dimx <- dim(X)
    w <- rep(1, dimx[1])
    if(!missing(iout) && length(iout) > 0)
        w[iout] <- 0

    Xhat <- rrcov3way::toArray(A %*% t(krp(C, B)), dimx[1], dimx[2], dimx[3])
    a <- apply(X-Xhat, 1, function(x) {norm(x, type="F")^2})
    a <- sum(w*a)
    a/(sum(w) * dimx[2] * dimx[3])
}

##  If one of the matrices has more columns, remove the extra columns.
##  In case of classical PARAFAC with overfactoring, on data with
##  contamination, the one extra estimated column has extremal effect,
##  reducing drastically the angle between the subspaces. This effect
##  is not present with the robust methods, i.e. removing the extra
##  column doesnot change the angle in case of one of the robust methods.
mysubspace <- function(A, B) {
    if(ncol(A) > ncol(B))
        A <- A[, 1:ncol(B)]
    if(ncol(B) > ncol(A))
        B <- B[, 1:ncol(A)]

    subspace(A, B)
}

##================================================================
##  Generate simulation data from PARAFAC model with homoscedastic
##  and heteroscedastic noise, congruence and added outliers
##
##  Generate outliers:
##      - eps: contamination fraction
##      - type of the outliers:
##        none (eps=0) or
##        gl=good leverage points: c1 > 1 and c2 = 0, or
##        bl=bad leverage points: c1 > 1 and c2 > 0 or
##        og=residual outliers: c1=1 and c2 > 0
##      - In case of gl and bl the outliers are added to the "pure"
##          data and after that the noise is added. In case of og
##          the outliers are added to the noisy data.
##      - iout is a vector containing the outlier indexes
##

if(FALSE) {
    I <- 100
    J <- 20
    K <- 20
    nsim <- 2
    nf <- 3; noise <- 0.15; noise1 <- 0.1
    cong <- 0.9
    xdat <- cp_gen(I=I, J=J, K=K, nsim=nsim, nf=nf,
        noise=noise, noise1=noise1, Acol=TRUE, Bcol=TRUE, Ccol=TRUE,
        congA=cong, congB=cong, congC=cong, eps=0)

    cond(xdat$As[[1]])
    cond(xdat$As[[2]])
}

cp_gen <- function(I=20, J=20, K=20, nsim=200, nf=3,
    noise=0.05, noise1=0, Acol=TRUE, Bcol=TRUE, Ccol=TRUE,
    congA=0.5, congB=0.5, congC=0.5,
    eps=0, type=c("bl", "gl", "og"), c1=10, c2=0.1)
{

    ## Add factor collinearity c to the matrix A
    addcong <- function(A, c) {
        nf <- ncol(A)
        congr <- matrix(c, nf, nf)          # congruence matrix we want, e.g. t(A) %*% (A)
        diag(congr) <- 1
        R <- chol(congr)                    # upper triangular matrix

        A.qr <- qr(A)
        Q <- qr.Q(A.qr)
        R1 <- qr.R(A.qr)
        sgn <- sign(diag(R1))               # assicura elementi positivi
        R.new <- diag(sgn) %*% R            # assicura elementi positivi
        A <- Q %*% R.new
        A
    }

    type <- match.arg(type)         # type of outliers
    if(missing(c1) | missing(c2)) {
        if(type == "bl") {
            c1 <- 10; c2 <- 0.1
        } else if(type == "gl") {
            c1 <- 10; c2 <- 0.0
        } else if(type == "og") {
            c1 <- 1; c2 <- 0.1
        } else
            stop("Undefined outlier type")
    }

    param <- list(I=I, J=J, K=K, nsim=nsim, nf=nf, noise=noise, noise1=noise1,
        Acol=Acol, Bcol=Bcol, Ccol=Ccol, congA=congA, congB=congB, congC=congC,
        eps=eps, type=type, c1=c1, c2=c2)

    Xmat = array(NA, c(I, J, K))
    Xlist <- Alist <- Blist <- Clist <- Olist <- vector("list", nsim)

    for(i in 1:nsim)
    {
        Amat <- matrix(runif(I*nf), I, nf)
        Bmat <- matrix(runif(J*nf), J, nf)
        Cmat <- matrix(runif(K*nf), K, nf)

        ## Add factor collinearity as requested
        if(Acol)
            Amat <- addcong(Amat, congA)
        if(Bcol)
            Bmat <- addcong(Bmat, congB)
        if(Ccol)
            Cmat <- addcong(Cmat, congC)

        ##  Generate the "pure" (i.e. no noise) X array
        Xmat <- toArray(Amat %*% t(krp(Cmat, Bmat)), I, J, K)

        ## Generate homoscedastic noise
        # the coefficient sqrt(n/(1-n)) adds appropriate level of noise in terms of total variability
        # e.g. 5% -> 0.2294, 10% -> 0.3333, 20% -> 0.5
        E <- array(rnorm(I*J*K), c(I, J, K))
        for(k in 1:K) {
            E[,, k] <- sqrt(noise/(1-noise)) * E[,, k] %*% diag(sqrt(diag(t(Xmat[,,k]) %*% Xmat[,,k])) / sqrt(diag(t(E[,,k]) %*% E[,,k])))
        }

        ## Generate heteroschedastic noise
        E1 <- array(rnorm(I*J*K), c(I, J, K))
        for(k in 1:K) {
            E1[,,k] <- E1[,,k] * Xmat[,, k]
            E1[,,k] <- sqrt(noise1/(1-noise1)) * E1[,,k] %*% diag(sqrt(diag(t(Xmat[,,k]) %*% Xmat[,,k]))/sqrt(diag(t(E1[,,k]) %*% E1[,,k])))
        }

        ## Add outliers
        iout <- c()
        if(eps > 0) {
            iout <- sample(1:I, size=eps*I)

            if(type == "og") {
                if(c1 != 1)
                    warning("Bad leverage points are generated instead of residual outliers!")
                Xmat <- Xmat + E
                Xmat <- Xmat + E1
                Xmat[iout,,] <- Xmat[iout,,] + c2
            } else if(type == "gl") {
                if(c2 != 0)
                    warning("Bad leverage points are generated instead of good leverage points")
                Xmat[iout,,] <- Xmat[iout,,] * c1
                Xmat <- Xmat + E
                Xmat <- Xmat + E1
            } else if(type == "bl") {
                Xmat[iout,,] <- Xmat[iout,,] * c1 + c2
                Xmat <- Xmat + E
                Xmat <- Xmat + E1
            }
        } else {
            Xmat <- Xmat + E
            Xmat <- Xmat + E1
        }

        ## Prepare the output: A, B, C, X and O[utliers]
        Alist[[i]] <- Amat
        Blist[[i]] <- Bmat
        Clist[[i]] <- Cmat

        Xlist[[i]] <- Xmat

        if(length(iout) > 0)
            Olist[[i]] <- iout
    }

    list(As=Alist, Bs=Blist, Cs=Clist, Xs=Xlist, Os=Olist, param=param)
}

## Standard ALS for CP - follows CPfunc from ThreeWay
##
CPf <- function(X_, r, ort1=1, ort2=1, ort3=1, start=0, conv=1e-6, maxit=10000, A, B, C, mnor=FALSE, trace=FALSE)
{
    di <- dim(X_)

    n <- di[1]
    m <- di[2]
    p <- di[3]
    dn <- dimnames(X_)
    X <- ThreeWay::supermat(X_)$Xa

    ## matrice iterazioni e convergenza
    single_iter <- matrix(0, maxit)

    ####################################
    ftiter = matrix(0, maxit/10, 2)
    mintripcos = 0
    Xfit <- array(0, c(n,m,p))
    cputime = system.time({
        ssx = sum(X^2)
        if (start == 0) {
            if (n >= r) {
                AUT = eigen(X %*% t(X))
                A = AUT$vectors[, 1:r]
            }
            else {
                A = orth(matrix(runif(r * r, 0, 1), nrow = r) -
                             0.5)
                A = A[1:n, ]
            }
            Z = ThreeWay::permnew(X, n, m, p)
            if (m >= r) {
                AUT = eigen(Z %*% t(Z))
                B = AUT$vectors[, 1:r]
            }
            else {
                B = orth(matrix(runif(r * r, 0, 1), nrow = r) -
                             0.5)
                B = B[1:m, ]
            }
            Z = ThreeWay::permnew(Z, m, p, n)
            if (p >= r) {
                AUT = eigen(Z %*% t(Z))
                C = AUT$vectors[, 1:r]
            }
            else {
                C = orth(matrix(runif(r * r, 0, 1), nrow = r) -
                             0.5)
                C = C[1:p, ]
            }
        }
        if (start == 1) {
            if (n >= r) {
                A = orth(matrix(runif(n * r, 0, 1), nrow = n) -
                             0.5)
            }
            else {
                A = orth(matrix(runif(r * r, 0, 1), nrow = r) -
                             0.5)
                A = A[1:n, ]
            }
            if (m >= r) {
                B = orth(matrix(runif(m * r, 0, 1), nrow = m) -
                             0.5)
            }
            else {
                B = orth(matrix(runif(r * r, 0, 1), nrow = r) -
                             0.5)
                B = B[1:m, ]
            }
            if (p >= r) {
                C = orth(matrix(runif(p * r, 0, 1), nrow = p) -
                             0.5)
            }
            else {
                C = orth(matrix(runif(r * r, 0, 1), nrow = r) -
                             0.5)
                C = C[1:p, ]
            }
        }

        H = matrix(0, r, r^2)
        for (ii in 1:r) {
            H[ii, (ii - 1) * r + ii] = 1
        }
        H1 = ThreeWay::permnew(H, r, r, r)
        H1 = ThreeWay::permnew(B %*% H1, m, r, r)
        H1 = ThreeWay::permnew(C %*% H1, p, r, m)
        f = sum((X - A %*% H1)^2)

        fold = f + 2 * conv * f
        iter = 0
        BB = t(B) %*% B
        CC = t(C) %*% C
        while ((fold - f > conv * f | iter < 2) & f > conv^2 & iter < maxit) {
            fold = f
            Z1 = ThreeWay::permnew(X, n, m, p)
            Z1 = ThreeWay::permnew(t(B) %*% Z1, r, p, n)
            Z1 = ThreeWay::permnew(t(C) %*% Z1, r, n, r)
            XF = Z1 %*% t(H)
            if (ort1 == 1) {
                FF = BB * CC
                A = XF %*% solve(FF)
            }
            if (ort1 == 2) {
                SVD = svd(XF)
                A = SVD$u %*% t(SVD$v)
            }
            if (ort1 == 3) {
                FF = BB * CC
                SVD = svd(XF - matrix(1, n, 1) %*% apply(XF,
                                                         2, mean))
                A = SVD$u %*% t(SVD$v) + matrix(1, n, 1) %*%
                    apply(XF, 2, mean) %*% solve(FF)
            }
            AA = t(A) %*% A
            Z = ThreeWay::permnew(X, n, m, p)
            Z1 = ThreeWay::permnew(Z, m, p, n)
            Z1 = ThreeWay::permnew(t(C) %*% Z1, r, n, m)
            Z1 = ThreeWay::permnew(t(A) %*% Z1, r, m, r)
            XF = Z1 %*% t(H)
            if (ort2 == 1) {
                FF = AA * CC
                B = XF %*% solve(FF)
            }
            if (ort2 == 2) {
                SVD = svd(XF)
                B = SVD$u %*% t(SVD$v)
            }
            if (ort2 == 3) {
                FF = AA * CC
                SVD = svd(XF - matrix(1, m, 1) %*% apply(XF,
                                                         2, mean))
                B = SVD$u %*% t(SVD$v) + matrix(1, m, 1) %*%
                    apply(XF, 2, mean) %*% solve(FF)
            }
            BB = t(B) %*% B
            Z = ThreeWay::permnew(Z, m, p, n)
            Z1 = ThreeWay::permnew(Z, p, n, m)
            Z1 = ThreeWay::permnew(t(A) %*% Z1, r, m, p)
            Z1 = ThreeWay::permnew(t(B) %*% Z1, r, p, r)
            XF = Z1 %*% t(H)
            if (ort3 == 1) {
                FF = AA * BB
                C = XF %*% solve(FF)
            }
            if (ort3 == 2) {
                SVD = svd(XF)
                C = SVD$u %*% t(SVD$v)
            }
            if (ort3 == 3) {
                FF = AA * BB
                SVD = svd(XF - matrix(1, p, 1) %*% apply(XF,
                                                         2, mean))
                C = SVD$u %*% t(SVD$v) + matrix(1, p, 1) %*%
                    apply(XF, 2, mean) %*% solve(FF)
            }

            CC = t(C) %*% C
            if (ort3 == 1) {
                f = ssx - ThreeWay::tr(CC %*% FF)
            }
            else {
                H1 = ThreeWay::permnew(H, r, r, r)
                H1 = ThreeWay::permnew(B %*% H1, m, r, r)
                H1 = ThreeWay::permnew(C %*% H1, p, r, m)
                f = sum((X - A %*% H1)^2)
            }
            iter = iter + 1

            ### iterazioni e valori di convergenza
            single_iter[iter] <- abs(fold - f)
            ######################################

            if(iter%%10 == 0) {
                tripcos <- min(ThreeWay::phi(A, A) * ThreeWay::phi(B, B) * ThreeWay::phi(C, C))
                if(iter == 10)
                    mintripcos = tripcos
                if(tripcos < mintripcos)
                    mintripcos = tripcos
                if(trace && iter%%1000 == 0)
                    cat(paste("Minimal Triple cosine =", tripcos, sep = " "), fill = TRUE)

                ftiter[iter/10, ] = c(f, tripcos)
            }
            if(trace && iter %% 50 == 0)
                cat(paste("f=", f, "after", iter, "iters; diff.=", (fold - f), sep = " "), fill = TRUE)
        }
    })

    ftiter <- ftiter[1:(iter/10), , drop=FALSE]           # take only the first iter/10 rows

    ## Fit percentage
    fp = 100 - 100 * f/ssx

    ## Degeneracy problem if |tripcos|>.5
    tripcos = min(ThreeWay::phi(A, A) * ThreeWay::phi(B, B) * ThreeWay::phi(C, C))
    names(tripcos) = c("Minimal triple cosine")
    if (iter < 10) {
        mintripcos = tripcos
    }

    ## SCALE
    if(mnor) {
        qx <- do_scalef(A, B, C)
        A <- qx$A; B <- qx$B; C <- qx$C
    }

    ret <- list(fit=f, fp=fp, ss=ssx, A=A, B=B, C=C, iter=iter, cputime=cputime[1],
        tripcos=tripcos, mintripcos=mintripcos, ftiter=ftiter, single_iter=single_iter[1:iter],
        robust=FALSE, coda.transform="none")
    class(ret) <- "parafac"
    ret
}   # CPf

## Standard ALS for CP - calls rrcov3way:::cp_als()
##
CP_ALS <- function(X_, r, ort1=1, ort2=1, ort3=1, start=0, conv=1e-6, maxit=10000, A, B, C, mnor=FALSE, trace=FALSE)
{
    di <- dim(X_)

    n <- di[1]
    m <- di[2]
    p <- di[3]
    dn <- dimnames(X_)
    X <- unfold(X_)

    const <- c("none", "none", "none")
    const[1] <- if(ort1 == 1) "none" else if(ort1 == 2) "orth" else if(ort1 == 3) "zerocor" else stop("Invalid constraint, must be 1, 2 or 3!")
    const[2] <- if(ort2 == 1) "none" else if(ort2 == 2) "orth" else if(ort2 == 3) "zerocor" else stop("Invalid constraint, must be 1, 2 or 3!")
    const[3] <- if(ort3 == 1) "none" else if(ort3 == 2) "orth" else if(ort3 == 3) "zerocor" else stop("Invalid constraint, must be 1, 2 or 3!")

    start <- if(start==0) "svd" else if(start==1) "random" else stop("Invalid initialization method!")

    single_iter <- matrix(0, maxit)

    ret <- rrcov3way:::cp_als(X, n, m, p, ncomp=r, const=const, start=start, conv=conv, maxit=maxit, trace=trace)

    if(mnor) {
        qx <- do_scalef(ret$A, ret$B, ret$C)
        ret$A <- qx$A; ret$B <- qx$B; ret$C <- qx$C
    }

    ret
}

atld <- function(X, r, conv=1e-06, start=1, maxit=5000, mnor=FALSE)
{
    # Input
    # X is a threeway array (n x m x p)
    # r = number of factors
    # mnor is a flag: to obtain normalized loading matrices write TRUE

    # Output
    # out$A = loading matrix - first mode
    # out$B = loading matrix - second mode
    # out$C = loading matrix - third mode
    # out$fp = fit value expressed as a percentage
    # out$iter = number of iterations
    # out$cputime = time to convergence
    # out$tripcos is a measure of degeneracy

    #### ASSIGN INITIAL VALUES
    di <- dim(X)
    n <- di[1]
    m <- di[2]
    p <- di[3]
    dn <- dimnames(X)
    Xk <- rrcov3way::unfold(X)
    Xfit <- array(0,c(n,m,p))
    single_iter <- matrix(0, maxit)     # output
    ftiter <- matrix(0, maxit/10, 2)
    mintripcos <- 0
    epsilon <- 10 * .Machine$double.eps * max(di)

    cputime = system.time(
    {
        ss <- ssx <- sum(X^2)

        if(start == 0) {
            if(n >= r) {
                AUT = eigen(Xk %*% t(Xk))
                A = AUT$vectors[, 1:r]
            }
            else {
                A = pracma::orth(matrix(runif(r * r, 0, 1), nrow = r) - 0.5)
                A = A[1:n, ]
            }
            Z = rrcov3way::permute(Xk, n, m, p)
            if (m >= r) {
                AUT = eigen(Z %*% t(Z))
                B = AUT$vectors[, 1:r]
            }
            else {
                B = pracma::orth(matrix(runif(r * r, 0, 1), nrow = r) - 0.5)
                B = B[1:m, ]
            }
            Z = rrcov3way::permute(Z, m, p, n)
            if (p >= r) {
                AUT = eigen(Z %*% t(Z))
                C = AUT$vectors[, 1:r]
            }
            else {
                C = pracma::orth(matrix(runif(r * r, 0, 1), nrow = r) - 0.5)
                C = C[1:p, ]
            }
        }else if (start == 1) {
            if (n >= r) {
                A = pracma::orth(matrix(runif(n * r, 0, 1), nrow = n) - 0.5)
            }
            else {
                A = pracma::orth(matrix(runif(r * r, 0, 1), nrow = r) - 0.5)
                A = A[1:n, ]
            }
            if (m >= r) {
                B = pracma::orth(matrix(runif(m * r, 0, 1), nrow = m) - 0.5)
            }
            else {
                B = pracma::orth(matrix(runif(r * r, 0, 1), nrow = r) - 0.5)
                B = B[1:m, ]
            }
            if (p >= r) {
                C = pracma::orth(matrix(runif(p* r, 0, 1), nrow = p) - 0.5)
            }
            else {
                C = pracma::orth(matrix(runif(r * r, 0, 1), nrow = r) - 0.5)
                C = C[1:p, ]
            }
        }

        ###

        ## VT::16.02.2022 - fix the try()
        ##  PC <- try(t(C %*% solve(crossprod(C))))
        ##  if(inherits(PC, "try-error"))
        ##      PC <- pracma::pinv(C, tol=epsilon)
        PC <- do_inv(C, epsilon)

        f <- sum((Xk - A %*% t(rrcov3way::krp(C,B)))^2)
        LF = f + 2 * conv * f
        iter = 0
        #LF=0    #

        ## INITIALIZATION & OPTIMIZATION: ATLD Iterative steps
        while(abs((f - LF)/f) > conv && iter < maxit)
        {

            ##  cat("\n", iter, " f, LF, conv, (f-LF)/f: ", f, LF, conv, (f-LF)/f, "\n")

            iter <- iter + 1
            LF <- f

            Z1 <- rrcov3way::permute(Xk, n, m, p)
            A <- A %*% diag(1/sqrt(diag(t(A) %*% A))) %*% diag(sqrt(diag(t(C) %*% C)))
            PA <- do_inv(A, epsilon)

            B <- Z1 %*% rrcov3way::krp(t(PA), t(PC))
            Z1 <- rrcov3way::permute(Z1, m, p, n)
            B <- B %*% diag(1/sqrt(diag(t(B) %*% B))) %*% diag(sqrt(diag(t(A) %*% A)))
            PB <- do_inv(B, epsilon)


            C <- Z1 %*% rrcov3way::krp(t(PB),t(PA))
            C <- C %*% diag(1/sqrt(diag(t(C) %*% C))) %*% diag(sqrt(diag(t(B) %*% B)))
            PC <- do_inv(C, epsilon)

            A <- Xk %*% rrcov3way::krp(t(PC), t(PB))

            # Loss of fit
            f = sum((Xk - A %*% t(rrcov3way::krp(C,B)))^2)

            # Record Relative Fit
            single_iter[iter] <- abs((f - LF)/f)
        }
    })

    ## 100 - 100 *Lack-of-fit (LOF)
    fp = 100 - 100 * f/ss

    ## Degeneracy problem if |tripcos| > 0.5
    tripcos <- min(rrcov3way::congruence(A, A) * rrcov3way::congruence(B, B) * rrcov3way::congruence(C, C))
    names(tripcos) <- c("Minimal triple cosine")
    if(iter < 10) {
        mintripcos = tripcos
    }

    ## SCALE
    if(mnor) {
        qx <- do_scalef(A, B, C)
        A <- qx$A; B <- qx$B; C <- qx$C
    }

    pfac <- list(fit=f, fp=fp, ss=ssx, A=A, B=B, C=C, iter=iter, cputime=cputime[1],
        tripcos=tripcos, mintripcos=mintripcos, single_iter=single_iter[1:iter])
     pfac
}   # ATLD

int2 <- function (X, r, initconv=1e-01, conv=1e-06, start=0, ort1=1, ort2=1, ort3=1, maxit=5000, mnor=FALSE, trace=FALSE)
{
    # Input
    # X = Threeway array (n x m x p)
    # r = Number of factors
    # initconv = Convergence criterion for ATLD stage
    # conv = Convergence criterion for final ALS stage
    # maxit =	Maximal number of iterations
    # ort1 (A-mode constraints): 1 = No constrains, 2 = Orthonormality, 3 = Zero correlation
    # ort2 (B-mode constraints): 1 = No constrains, 2 = Orthonormality, 3 = Zero correlation
    # ort3 (C-mode constraints): 1 = No constrains, 2 = Orthonormality, 3 = Zero correlation
    # mnor is a flag: FALSE (default), TRUE (normalized loading matrices)

    ## Output
    # out$A = Score matrix - first mode
    # out$B = Loading matrix - second mode
    # out$C = Loading matrix - third mode
    # out$fp = Fit value expressed as percentage
    # out$iter = Number of iterations
    # out$cputime = Total Time to convergence (Stage n + Stage II)
    # out$tripcos = Minimal triple cosine product (measure of 2FDs)
    # out$Xfit = Reconstructed array, flattened with respect to the first mode,  dim= (n x JK)
    # out$single_iter = Relative fit for each iteration

    #### ASSIGN INITIAL VALUES
    di <- dim(X)
    n <- di[1]
    m <- di[2]
    p <- di[3]
    dn <- dimnames(X)
    Xk <- rrcov3way::unfold(X)

    single_iter <- matrix(0,maxit)# output
    ftiter <- matrix(0, maxit/10, 2)
    mintripcos <- 0
    eps <- .Machine$double.eps
    epsilon <- 10 * eps * norm(Xk, "1") * max(di)


    # constraints flags
    cputime <- system.time({
        ss <- ssx <- sum(X^2)

        if (start == 0) {
            if (n >= r) {
                AUT = eigen(Xk %*% t(Xk))
                A = AUT$vectors[, 1:r]
            }
            else {
                A = orth(matrix(runif(r * r, 0, 1), nrow = r) -
                             0.5)
                A = A[1:n, ]
            }
            Z = rrcov3way::permute(Xk, n, m, p)
            if (m >= r) {
                AUT = eigen(Z %*% t(Z))
                B = AUT$vectors[, 1:r]
            }
            else {
                B = orth(matrix(runif(r * r, 0, 1), nrow = r) -
                             0.5)
                B = B[1:m, ]
            }
            Z = rrcov3way::permute(Z, m, p, n)
            if (p >= r) {
                AUT = eigen(Z %*% t(Z))
                C = AUT$vectors[, 1:r]
            }
            else {
                C = orth(matrix(runif(r * r, 0, 1), nrow = r) -
                             0.5)
                C = C[1:p, ]
            }
        }
        if (start == 1) {
            if (n >= r) {
                A = orth(matrix(runif(n * r, 0, 1), nrow = n) -
                             0.5)
            }
            else {
                A = orth(matrix(runif(r * r, 0, 1), nrow = r) -
                             0.5)
                A = A[1:n, ]
            }
            if (m >= r) {
                B = orth(matrix(runif(m * r, 0, 1), nrow = m) -
                             0.5)
            }
            else {
                B = orth(matrix(runif(r * r, 0, 1), nrow = r) -
                             0.5)
                B = B[1:m, ]
            }
            if (p >= r) {
                C = orth(matrix(runif(p* r, 0, 1), nrow = p) -
                             0.5)
            }
            else {
                C = orth(matrix(runif(r * r, 0, 1), nrow = r) -
                             0.5)
                C = C[1:p, ]
            }
        }

        ## INITIALIZE ATLD
        if(trace)
            cat("\n Initialization of ATLD starting ...\n")

        PC <- do_inv(C, epsilon)
        f = sum((Xk - A %*% t(krp(C, B)))^2)
        LF = f + 2 * initconv * f
        iter = 0

        if(trace)
            cat("\n Initialization of ATLD ready. Starting ATLD iteration ...\n")

        ## ATLD Iterative steps
        while(abs((f - LF)/f) > initconv && iter < maxit) {
            iter <-  iter + 1
            LF <- f

            if(trace)
                cat("\n iter=", iter, "estimating A...\n")

            Z1 <-  rrcov3way::permute(Xk, n, m, p)
            A <- A %*% diag(1/sqrt(diag(t(A) %*% A))) %*% diag(sqrt(diag(t(C) %*% C)))
            PA <- do_inv(A, epsilon)
            B<-Z1 %*% krp(t(PA), t(PC))

            if(trace)
                cat("\n iter=", iter, "estimating B...\n")

            Z1 <- rrcov3way::permute(Z1, m, p, n)
            B <- B %*% diag(1/sqrt(diag(t(B) %*% B))) %*% diag(sqrt(diag(t(A) %*% A)))
            PB <- do_inv(B, epsilon)

            if(trace)
                cat("\n iter=", iter, "estimating C...\n")

            C <- Z1 %*% krp(t(PB), t(PA))
            C <- C %*% diag(1/sqrt(diag(t(C) %*% C))) %*% diag(sqrt(diag(t(B) %*% B)))
            PC <- do_inv(C, epsilon)

            A <- Xk %*% krp(t(PC), t(PB))

            ## Loss of fit
            f <- sum((Xk - A %*% t(krp(C, B)))^2)

            ## Record Relative Fit
            single_iter[iter] <- abs((f - LF)/f)

            if(trace)
                cat("\n iter=", iter, "f=", f, abs((f - LF)/f), "\n")

        }
        iter_opt <- iter

        ## REFINING SOLUTION
        ## INITIALIZE  ALS
        H <- matrix(0, r, r^2)
        for(ii in 1:r) {
            H[ii, (ii - 1) * r + ii] <- 1
        }
        H1 <- rrcov3way::permute(H, r, r, r)
        H1 <- rrcov3way::permute(B %*% H1, m, r, r)
        H1 <- rrcov3way::permute(C %*% H1, p, r, m)

        BB <- t(B) %*% B
        CC <- t(C) %*% C

        ## ALS iterations
        while((LF - f > conv * f | iter <= iter_opt+1) & f > conv^2 & iter < maxit) {
            iter = iter + 1
            LF = f
            Z1 = rrcov3way::permute(Xk, n, m, p)
            Z1 = rrcov3way::permute(t(B) %*% Z1, r, p, n)
            Z1 = rrcov3way::permute(t(C) %*% Z1, r, n, r)
            XF = Z1 %*% t(H)
            if (ort1 == 1) {
                FF = BB * CC
                A = XF %*% solve(FF)
            }
            if (ort1 == 2) {
                SVD = svd(XF)
                A = SVD$u %*% t(SVD$v)
            }
            if (ort1 == 3) {
                FF = BB * CC
                SVD = svd(XF - matrix(1, n, 1) %*% apply(XF,
                                                         2, mean))
                A = SVD$u %*% t(SVD$v) + matrix(1, n, 1) %*%
                    apply(XF, 2, mean) %*% solve(FF)
            }
            AA = t(A) %*% A
            Z = rrcov3way::permute(Xk, n, m, p)
            Z1 = rrcov3way::permute(Z, m, p, n)
            Z1 = rrcov3way::permute(t(C) %*% Z1, r, n, m)
            Z1 = rrcov3way::permute(t(A) %*% Z1, r, m, r)
            XF = Z1 %*% t(H)
            if (ort2 == 1) {
                FF = AA * CC
                B = XF %*% solve(FF)
            }
            if (ort2 == 2) {
                SVD = svd(XF)
                B = SVD$u %*% t(SVD$v)
            }
            if (ort2 == 3) {
                FF = AA * CC
                SVD = svd(XF - matrix(1, m, 1) %*% apply(XF,
                                                         2, mean))
                B = SVD$u %*% t(SVD$v) + matrix(1, m, 1) %*%
                    apply(XF, 2, mean) %*% solve(FF)
            }
            BB = t(B) %*% B
            Z = rrcov3way::permute(Z, m, p, n)
            Z1 = rrcov3way::permute(Z, p, n, m)
            Z1 = rrcov3way::permute(t(A) %*% Z1, r, m, p)
            Z1 = rrcov3way::permute(t(B) %*% Z1, r, p, r)
            XF = Z1 %*% t(H)
            if (ort3 == 1) {
                FF = AA * BB
                C = XF %*% solve(FF)
            }
            if (ort3 == 2) {
                SVD = svd(XF)
                C = SVD$u %*% t(SVD$v)
            }
            if (ort3 == 3) {
                FF = AA * BB
                SVD = svd(XF - matrix(1, p, 1) %*% apply(XF,
                                                         2, mean))
                C = SVD$u %*% t(SVD$v) + matrix(1, p, 1) %*%
                    apply(XF, 2, mean) %*% solve(FF)
            }
            CC = t(C) %*% C
            if (ort3 == 1) {
                f = ssx - rrcov3way::mtrace(CC %*% FF)
            }
            else {
                H1 = rrcov3way::permute(H, r, r, r)
                H1 = rrcov3way::permute(B %*% H1, m, r, r)
                H1 = rrcov3way::permute(C %*% H1, p, r, m)
                f = sum((Xk - A %*% H1)^2)
            }

            ##  Record Relative Fit
            single_iter[iter] <- abs((LF - f)/f)

            ##  TRIPLE COSINE
            if(iter%%10 == 0) {
                tripcos = min(rrcov3way::congruence(A, A) * rrcov3way::congruence(B, B) * rrcov3way::congruence(C, C))
                if (iter == 10)
                    mintripcos = tripcos
                if (tripcos < mintripcos)
                    mintripcos = tripcos
                if (trace && iter%%1000 == 0)
                    cat(paste("Minimal Triple cosine =", tripcos, sep = " "), fill = TRUE)

                ftiter[iter/10, ] = c(f, tripcos)
            }
            if(trace && iter %% 50 == 0)
                cat(paste("f=", f, "after", iter, "iters; diff.=", LF - f, sep = " "), fill = TRUE)
        }
    })

    ftiter <- ftiter[1:(iter/10), , drop=FALSE]           # take only the first iter/10 rows

    ## Fit percentage
    fp = 100 - 100 * f/ss

    ## Degeneracy problem if |tripcos|>.5
    tripcos = min(rrcov3way::congruence(A, A) * rrcov3way::congruence(B, B) * rrcov3way::congruence(C, C))
    names(tripcos) = c("Minimal triple cosine")
    if (iter < 10)
        mintripcos = tripcos

    ## SCALE
    if(mnor) {
        qx <- do_scalef(A, B, C)
        A <- qx$A; B <- qx$B; C <- qx$C
    }

    pfac <- list(fit=f, fp = fp, ss=ssx, A = A, B = B, C = C, iter = iter,
        iter_opt=iter_opt, cputime=cputime[1], tripcos=tripcos, mintripcos=mintripcos, ftiter=ftiter,
        single_iter=single_iter[1:iter], robust=FALSE, coda.transform="none")

    class(pfac) <- "parafac"

    pfac
}   # int2

cp_int2_test <- function(X, r, initconv=1e-01, conv=1e-06, start=0, ort1=1, ort2=1, ort3=1, maxit=5000, mnor=FALSE, trace=FALSE) {

    cpu <- system.time({
        atld <- atld(X, r=r, start=start, conv=initconv, maxit=maxit, mnor=FALSE)
        cp <- CPf(X, r=r, ort1=1, ort2=1, ort3=1, start=2, A=atld$A, B=atld$B, C=atld$C, conv=conv, maxit=maxit, mnor=mnor, trace=trace)
    })

    out <- cp
    out$iter <- atld$iter + cp$iter
    out$iter_opt <- atld$iter
    out$cputime <- atld$cputime + cp$cputime
    out
}

##===================================================================================

R_als <- function (X, ncomp=2, type=c("als", "int2"), initconv=1e-2,
    const="none", conv=1e-6, start="svd", maxit=10000,
    ncomp.rpca=0, alpha=0.75, robiter=100, crit=0.975,
    trace=FALSE)
{
    type <- match.arg(type)
    di <- dim(X)
    I <- di[1]
    J <- di[2]
    K <- di[3]
    dn <- dimnames(X)

    ssx <- sum(X^2)
    Xwide <- unfold(X)
    Ahat <- matrix(0, I, ncomp)

    ## define number of outliers the algorithm should resists
    h <- round(alpha*I)

cputime = system.time({

    ## Step 1 ROBPCA of XA
    if(trace)
        cat("\nStep 1. Perform robust PCA on the unfolded matrix.")

    outrobpca <- PcaHubert(Xwide, k=ncomp.rpca, kmax=ncol(Xwide), alpha=alpha, mcd=FALSE, trace=trace)
    Hset <- sort(sort(outrobpca@od, index.return=TRUE)$ix[1:h])
    Xhat <- Xwide[Hset,]
    fitprev <- 0
    changeFit <- 1 + conv
    aveiter <- iter <- 0

    if(trace)
        cat("\nStep 2. PARAFAC analysis. Start iteration, type=", toupper(type))

    start1 <- if(start == "svd") 0 else 1

    while (changeFit > conv & iter <= robiter)
    {
        iter <- iter + 1

        ##  Step 2 - PARAFAC analysis
        if(type == "als")
            ##ret <- rrcov3way:::cp_als(Xhat, h, J, K, ncomp=ncomp, const=const, conv=conv, start=start, maxit=maxit, trace=trace)
            ret <- CPf(toArray(Xhat, h, J, K), ncomp, conv=conv, start=start1, maxit=maxit, trace=trace)
        else
            ret <- int2(toArray(Xhat, h, J, K), ncomp, initconv=initconv, conv=conv, start=start1, maxit=maxit, trace=trace)

        aveiter <- aveiter + ret$iter

        Ah <- ret$A
        Bh <- ret$B
        Ch <- ret$C

        ## Step 3 - Fit the model
        KR <- krp(Ch, Bh)                   # Khatri-Rao product

##        for(i in 1:I) {
##            vJKx1 <- matrix(X[i,,], 1, J*K)
##            Ahat[i,] <- pracma::pinv(KR) %*% t(vJKx1)
##        }

        ## The above is equivalent to the following:
        Ahat <- Xwide %*% t(pracma::pinv(KR))

        Xfit <- Ahat %*% t(KR)

        ##  Step 4  - Computation of the residual distances
        rdsq <- apply((Xwide - Xfit)^2, 1, sum)
        rd <- sqrt(rdsq)

        Hset <- sort(sort(rdsq, index.return=TRUE)$ix[1:h])
        fit <- sum(rdsq[Hset])
        Xhat <- Xwide[Hset,]

        ##  Step 5  Fit of the model
        changeFit <- if(fitprev == 0) 1 + conv else abs(fit-fitprev)/fitprev

        if(trace)
            cat("\n---", toupper(type), iter, "Fit, Fitprev, changeFit, iter", fit, fitprev, changeFit, ret$iter)

        fitprev <- fit
     }
     if(trace)
        cat("\n--- ---\n")

    aveiter <- round(aveiter/iter)  # average number of ALS/INT2 iterations throughout the robust iteration

    ## Reweighting
    if(trace)
        cat("\nReweghting step, type=", toupper(type))

    cutoff.rd <- rrcov3way:::.cutoff.rd(rd, crit=crit, h)
    flag <- rd <= cutoff.rd
    Xflag <- X[flag,,]
    Xflag_wide <- unfold(Xflag)
    if(type == "als")
        ## ret <- rrcov3way:::cp_als(Xflag_wide, nrow(Xflag_wide), J, K, ncomp=ncomp, const=const, conv=conv, start=start, maxit=maxit, trace=trace)
        ret <- CPf(toArray(Xflag_wide, nrow(Xflag_wide), J, K), ncomp, conv=conv, start=start1, maxit=maxit, trace=trace)
    else
        ret <- int2(toArray(Xflag_wide, nrow(Xflag_wide), J, K), ncomp, initconv=initconv, conv=conv, start=start1, maxit=maxit, trace=trace)

    Arew <- ret$A
    Brew <- ret$B
    Crew <- ret$C

    KRrew <- krp(Crew, Brew)            # Khatri-Rao product

##    Arew <- matrix(0, I, ncomp)
##    for(i in 1:I) {
##        vJKx1 <- matrix(X[i,,], 1, J*K)
##        Arew[i,] <- pracma::pinv(KRrew) %*% t(vJKx1)
##    }

    ## The above is equivalent to the following:
    Arew <- Xwide %*% t(pracma::pinv(KRrew))
})

    Xfit <- Arew %*% t(KRrew)
    rdsq <- apply((Xwide - Xfit)^2, 1, sum)
    rd <- sqrt(rdsq)
    cutoff.rd <- rrcov3way:::.cutoff.rd(rd, crit=crit, h)
    fit <- sum(rdsq[rd <= cutoff.rd])
    fp <- 100*(1-fit/ssx)

    for(i in 1:ncomp) {
        Arew[,i] <- Arew[,i]*norm(as.matrix(Brew[,i]),type="F")*norm(as.matrix(Crew[,i]),type="F")
        Brew[,i] <- Brew[,i]/norm(as.matrix(Brew[,i]),type="F")
        Crew[,i] <- Crew[,i]/norm(as.matrix(Crew[,i]),type="F")
    }

    sd <- rrcov3way:::.cutoff.sd(Arew, alpha=alpha, crit=crit, robust=TRUE)
    flag <- rd <= cutoff.rd & sd$sd <= sd$cutoff.sd
    Xfit <- toArray(Xfit, I, J, K)

    ## dimnames back
    nfac <- paste0("F", 1:ncomp)
    dimnames(Arew) <- list(dn[[1]], nfac)
    dimnames(Brew) <- list(dn[[2]], nfac)
    dimnames(Crew) <- list(dn[[3]], nfac)
    dimnames(Xfit) <- dn
    names(rd) <- names(sd$sd) <- names(flag) <- dn[[1]]

    res <- list(fit=fit, fp=fp, ss=ssx, A=Arew, B=Brew, C=Crew, iter=ret$iter, const=ret$const,
                ftiter=ret$ftiter, aveiter=aveiter,
                flag=flag, Hset=Hset,
                Xhat=Xfit, rd=rd, cutoff.rd=cutoff.rd, sd=sd$sd, cutoff.sd=sd$cutoff.sd,
                alpha=alpha, pcaobj=outrobpca, robiter=iter, cputime=cputime[1], robust=TRUE, coda.transform="none")
    class(res) <- "parafac"
    res
}
