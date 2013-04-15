\name{glmdm}
\title{Generalized Linear Mixed Dirichlet Model}
\alias{glmdm}
\description{R code for simulation of GLMDM}
\usage{glmdm(formula, family=gaussian, data, num.reps=1000, a1=3, b1=2, 
	      d=0.25, MM=15, VV=30, ...)}

\arguments{
    \item{formula}{a symbolic description of the model to be fit.}
    \item{family}{a descreption of the error distribution and link function to be used in the model.}
    \item{data}{an optional data frame, list or environment containing the variables in the model.}
    \item{num.reps}{num.reps}
    \item{a1}{a1}
    \item{b1}{b1}
    \item{d}{d}
    \item{MM}{MM}
    \item{VV}{VV}
    \item{...}{..}
}

\alias{A.K.can}
\alias{A.K.labels.can}
\alias{A.n.can}
\alias{bb}
\alias{cand}
\alias{eta}
\alias{f.ratio}
\alias{f.y.can}
\alias{f.y.old}
\alias{glmdm}
\alias{j}
\alias{K.can}
\alias{L.m.hat}
\alias{L.m.s.hat}
\alias{like.K.can}
\alias{Lms.hat}
\alias{M}
\alias{m.hat}
\alias{m.hat.s}
\alias{m.hess.s}
\alias{m.hessian}
\alias{Mb}
\alias{mean}
\alias{mean.m}
\alias{meta}
\alias{mle.m}
\alias{mle.m.s}
\alias{mn}
\alias{mult.can}
\alias{mult.old}
\alias{mult.ratio}
\alias{new.q}
\alias{nu}
\alias{p.A.can}
\alias{p.A.old}
\alias{p.ratio}
\alias{pq}
\alias{psi.can}
\alias{rho}
\alias{Sca}
\alias{Sha}
\alias{test}
\alias{var.m}
\alias{X.betaM}
\alias{log}
\alias{...}

\examples{
\dontrun{

data(scotvote)
glmdm.linear.out <- glmdm(PerYesParl ~ PrivateHousingStarts + CouncilTax 
+ Percentage5to15 + PrimaryPTRatio + PerBirthsOut + PerClaimantFemale, 
data=scotvote, num.reps=5000)

data(ssas)
glmdm.probit.ssas <- glmdm(scotpar2 ~ househld + rsex + rage + relgsums + ptyallgs
+ idlosem + marrmus + ukintnat + natinnat + voiceuk3 + nhssat, data=ssas, 
family=binomial(link="probit"), num.reps=10000, log=TRUE)

data(asia)
glmdm.probit.asia <- glmdm(ATT ~ DEM + FED + SYS + AUT, data=asia, 
family=binomial(link="probit"), num.reps=10000, log=TRUE)

}
}

