\name{glmdm.linear}
\title{glmdm.linear}
\alias{glmdm.linear}
\description{R code for simulation of GLMDM}
\usage{glmdm.linear(Y, X, family=linear, num.reps=1000, a1=3, b1=2, d=0.25, MM=15,VV=30,...)}


\arguments{
    \item{Y}{Y}
    \item{X}{X}
    \item{family}{family}
    \item{num.reps}{num.reps}
    \item{a1}{a1}
    \item{b1}{b1}
    \item{d}{d}
    \item{MM}{MM}
    \item{VV}{VV}
    \item{...}{..}
}


\alias{C}
\alias{done.vec}
\alias{fm}
\alias{getm}
\alias{logit}
\alias{lower.m}
\alias{rmultinorm}
\alias{test.vec}
\alias{tol}
\alias{upper.m}
\alias{A}
\alias{A.K}
\alias{A.K.labels}
\alias{A.n}
\alias{a1}
\alias{aa}
\alias{B}
\alias{b1}
\alias{beta}
\alias{beta.post}
\alias{d}
\alias{glm.out}
\alias{hat.beta}
\alias{hat.var}
\alias{i}
\alias{initial.m}
\alias{K}
\alias{K.save}
\alias{like}
\alias{like.K}
\alias{like.K.save}
\alias{loglike}
\alias{loglike.s}
\alias{m}
\alias{m.save}
\alias{MM}
\alias{mu}
\alias{n}
\alias{n.i}
\alias{nk}
\alias{p}
\alias{psi}
\alias{psi.post}
\alias{q}
\alias{Sc}
\alias{Sh}
\alias{SIG}
\alias{SIG.inv}
\alias{sigma}
\alias{sigma.beta.sq}
\alias{sigma.sq}
\alias{ssas.df}
\alias{tau}
\alias{tau.sq}
\alias{tau.sq.post}
\alias{VV}
\alias{X}
\alias{X.beta1}
\alias{xi.post}
\alias{Y}
\alias{Z}
\alias{linear}
\alias{glmdm.linear}

\examples{
\dontrun{
{
data(scotvote)
attach(scotvote)
Y <- PerYesParl
X <- cbind(rep(1,length=length(Y)),PrivateHousingStarts,CouncilTax,Percentage5to15,
           PrimaryPTRatio,PerBirthsOut,PerClaimantFemale)
detach(scotvote)
glmdm.linear(Y,X,num.reps=5000,log=FALSE) 
}
}
}
