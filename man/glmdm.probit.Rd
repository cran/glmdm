\name{glmdm.probit}
\title{glmdm.probit}
\alias{glmdm.probit}
\description{R code for simulation of GLMDM}
\usage{glmdm.probit(Y, X, num.reps=1000, a1=3, b1=2, d=0.25, MM=15,VV=30,...)}

\arguments{
    \item{Y}{Y}
    \item{X}{X}
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
\alias{glmdm.probit}
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
{
data(ssas)
attach(ssas)
Y <- ssas[,18]
X <- cbind(ssas$househld,ssas$rsex,ssas$rage,ssas$relgsums,ssas$ptyallgs,ssas$idlosem,
                 ssas$marrmus,ssas$ukintnat,ssas$natinnat,ssas$voiceuk3,ssas$nhssat,ssas$hincdif2,
                 ssas$unionsa,ssas$whrbrn,ssas$hedqual2)
detach(ssas)
glmdm.out <- glmdm.probit(Y,X,num.reps=10000)


data(asia)
Y <- asia[,1]
X <- as.matrix(asia[,-1]) 
glmdm.out <- glmdm.probit(Y,X,num.reps=10000,log=TRUE) 
}
}
}

