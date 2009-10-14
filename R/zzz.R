##########################################################################
##
##
## This software is distributed under the terms of the GNU GENERAL
## PUBLIC LICENSE Version 2, June 1991.  See the package LICENSE
## file for more information.
##
## Copyright (C) 2009 Jeff Gill, George Casella, and Jonathan Rapkin
## 
##########################################################################

#.onAttach <- function(...) {
 
  
#   date <- date()
#   x <- regexpr("[0-9]{4}", date)
#   this.year <- substr(date, x[1], x[1] + attr(x, "match.length") - 1)
   
   # echo output to screen
#   cat("##\n## glmdm pack \n")
#   cat("## Copyright (C) 2009-", this.year,
#      " Jeff Gill, George Casella, and Jonathan Rapkin\n", sep="")
  
#   require(compositions, quietly=TRUE)
#   require(MASS, quietly=TRUE)
#   require(msm, quietly=TRUE)
#   require(mvtnorm, quietly=TRUE)
#   require(e1071, quietly=TRUE)
#   require(coda, quietly=TRUE)
#   require(MCMCpack, quietly=TRUE)
#   require(class, quietly=TRUE)
#   require(lattice, quietly=TRUE)
#}

#.onUnload <- function(libpath) {
#    library.dynam.unload("glmdm", libpath)
#}

rmultinorm <- function(num.vals, mu.vec, vcmat, tol = 1e-08)  {
    k <- ncol(vcmat)

   if(length(mu.vec)!=k)
        stop(paste("rmultinorm error: rmultnorm: mu.vec vector wrong length:",length(mu.vec)))
    if(max(abs(vcmat - t(vcmat))) > tol)
        stop("rmultinorm error: variance-covariance matrix not symmetric")
    vs <- svd(vcmat)
    vcsqrt <- t(vs$v %*% (t(vs$u) * sqrt(vs$d)))
    ans.mat <- sweep(matrix(rnorm(num.vals*k), nrow = num.vals) %*% vcsqrt,2,mu.vec,"+")
    dimnames(ans.mat) <- list(NULL, dimnames(vcmat)[[2]])
    return(ans.mat)
}

logit <- function(Xb)  1/(1+exp(-Xb))
