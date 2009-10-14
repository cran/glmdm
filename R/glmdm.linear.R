
####################################################################################################
# R CODE FOR SIMULATION OF GLMDM, (c) Kyung, Casella, and Gill, February 2008
####################################################################################################

####################################################################################################
# START glmdm.linear FUNCTION
####################################################################################################
glmdm.linear <- function(Y, X, family=linear, num.reps=1000, log=TRUE, a1=3, b1=2, d=0.25, MM=15,VV=30,...) {
# SET QUANTITIES FROM DATA
n <- nrow(X); K <- ncol(X)
mu <- 0; tau <- 1; 
# For the prior of m, we use Gamma(a,b)
# Here ab=MM and ab^2=VV. 
Sh <- (MM^2)/VV
Sc <- VV/MM

####################################################################################################
# RUN A LM TO GET hat.beta and hat.var
####################################################################################################

suppressWarnings(lm.out <- lm(Y ~ X))
hat.beta <- coefficients(summary(lm.out))[,1]
hat.var  <- (coefficients(summary(lm.out))[,2])^2

####################################################################################################
# FIX PRIOR DISTRIBUTIONS AND PRIOR PARAMETERS
####################################################################################################

# PRE-LOOP EFFICIENCIES
tau.sq        <- rinvgamma(1,shape=a1,scale=b1)
# If tau.sq is infinity; 1.7e+100 is reasonable big?!
if(is.finite(tau.sq)==FALSE) tau.sq <- 1.7e+100
tau           <- sqrt(tau.sq)
sigma.sq      <- 1
sigma         <- sqrt(sigma.sq)
sigma.beta.sq <- sigma.sq * d
p             <- length(hat.beta)     
aa            <- chol(solve(t(X) %*% X + diag(1/sigma.beta.sq, nrow=ncol(X))))
SIG           <- t(X)%*%X + diag(1/sigma.beta.sq, nrow=ncol(X))				# VAR/COVAR MATRIX
SIG.inv       <- solve(SIG)								# INVERT HERE FOR EFFICIENCY
Z             <- rep(NA,length=n)							# Z TO FILL-IN DURING LOOPING
bbeta         <- rmvnorm(1, mean=hat.beta, sigma=diag(hat.var))				# STARTING VALUES FROM ABOVE
bbeta         <- rbind( bbeta,matrix(rep(NA,num.reps*ncol(bbeta)),ncol=ncol(bbeta)) )	# MATRIX TO FILL IN DURING LOOP
n.i           <- rep(1,n);   								# LOOP NEEDS REP OF n'S
q             <- rDirichlet.acomp(1, alpha=n.i)						# PROBABILITY OF ASSIGNMENT
A.n 	      <- sample(1:n,n,replace=TRUE,q)				 		# CREATES THE ASSIGNMENTS.
A.K 	      <- table(A.n)						 		# CREATES THE LIST OF OCCUPIED
A.K.labels    <- as.numeric(names(A.K))  						# LOCATIONS OF OCCUPIED
K	      <- length(A.K)								# NUMBER OF OCCUPIED TABLES
B             <- rep(0,length=n)							# LENGTH n FOR ALL CASES
B[A.K.labels] <- 1									# 1 WHERE OCCUPIED >0
B             <- B*rnorm(n,0,tau)							# ASSIGN psi VALUES TO OCCUPIED
nk            <- rep(0,n); for (i in 1:n) nk[A.n[i]] <- nk[A.n[i]] + 1			# COUNTS OCCUPANTS AT TABLES
psi           <- B[A.n] 								# MAP TABLE psi'S TO CASES
like.K        <- 0									# LOG-LIKE STARTS AT ZERO
X.beta1       <- X%*%bbeta[1,]								# MULTIPLY X AND bbeta IN ADVANCE
for (i in 1:n)  
    like.K <- like.K + Y[i] * pnorm( X.beta1[i] + psi[i],log=TRUE) + 
		(1-Y[i]) * (1 - pnorm( X.beta1[i] + psi[i] )) + 
		dnorm(psi[i], mean=0, sd=tau,log=TRUE)
like.K <- exp(like.K + sum(lgamma(A.K)))

# SETUP TO SAVE GIBBS SAMPLER VALUES
tau.sq.post <- psi.post <- beta.post <- xi.post <- K.save <- like.K.save <- m.save <- NULL 

####################################################################################################
# NEW POSTERIOR FOR m FUNCTION
####################################################################################################
initial.m <- m <-10 ; #initial value of m
like      <- function(m.v){gamma(m.v)/gamma(m.v+n) * dgamma(m.v, shape=Sh, scale=Sc) * m.v^K}
loglike   <- function(m.v){lgamma(m.v)-lgamma(m.v+n) + dgamma(m.v, shape=Sh, scale=Sc, log=TRUE) + K*log(m.v)}
loglike.s <- function(m.v){lgamma(m.v)-lgamma(m.v+n) + dgamma(m.v, shape=Sh, scale=Sc, log=TRUE) + K*log(m.v) + nu*log(m.v)}

####################################################################################################
# CALL GIBBS FUNCTION
####################################################################################################
# MCMC LOOPING
for (M in 1:num.reps)  {
    #   if (M %% 100 == 0) print(paste("finished iteration",M))
        
    #print("M-H ON THE 'A' MATRIX")
    p.A.old  <- (m^K)
    f.y.old  <- like.K
    mult.old <- dmultinom(x=A.K, prob=q[A.K.labels])

    #print("CREATE NEW 'A' MATRIX, 'can' STANDS FOR CANDIDATE") 
    pq                <- nk +1   
    new.q             <- rDirichlet.acomp(1, alpha=pq)
    A.n.can           <- sample(1:n,n,replace=TRUE,new.q)                                           
    A.K.can           <- table(A.n.can)                                                           
    A.K.labels.can    <- as.numeric(names(A.K.can))                                             
    K.can             <- length(A.K.can)                                                       
    B                 <- rep(0,length=n)                                                  
    B[A.K.labels.can] <- 1                                                               
    B                 <- B*rnorm(n,0,tau)                                               
    psi.can           <- B[A.n.can]                                                        
    p.A.can           <- (m^K.can)  
    like.K.can        <- 0                                                                      
    X.betaM           <- X%*%bbeta[M,]                                                          
    for (i in 1:n)
	like.K.can    <- like.K.can - log(3^2)*n/2 - (1/(2*3^2))*(Y[i] - X.betaM[i] + psi.can[i])^2
    like.K.can        <- exp(like.K.can + sum(lgamma(A.K.can)))
    f.y.can           <- like.K.can
    mult.can<-dmultinom(x=A.K.can, prob=new.q[A.K.labels.can])

    #print("UPDATE 'A' and 'K', CREATE LIKELIHOOD")
    p.ratio <- p.A.can/p.A.old; f.ratio <- f.y.can/f.y.old; mult.ratio <- mult.can/mult.old
    if (is.finite(p.ratio) == FALSE) p.ratio       <- 1
    if (is.finite(f.ratio) == FALSE) f.ratio       <- 1
    if (is.finite(mult.ratio) == FALSE) mult.ratio <- 1

    # NOW UPDATE GIVEN WE ARE DONE WITH THE METROPOLIS STEP (ACCEPTING 'can' OR NOT)

    #print("UPDATE rho")
    rho <- p.ratio * f.ratio * mult.ratio 
    if (rho>runif(1))   {
	A.n        <- A.n.can 		# A.n <- 1:n for regular random effects model
   	A.K        <- A.K.can
    	A.K.labels <- A.K.labels.can    
    	K          <- K.can
    	psi        <- psi.can 
    	like.K     <- like.K.can
    }

    #print("UPDATE 'z': Truncated Normal")
    for (j in 1:n){
	mean = X[j,]%*%bbeta[M,] + psi[j]
        Z[j] <- rtnorm(1, mean=mean, sd=1, lower=-Inf, upper=0)
        if (Y[j]==1) Z[j] <- rtnorm(1, mean=mean, sd=1, lower=0, upper=Inf)
    }

    #print("UPDATE 'bbeta', 'tau' AND 'eta'")
    mn            <- SIG.inv %*% (t(X) %*% (Z-psi)) 
    Mb            <- t(aa) %*% array(rnorm(p), c(p, 1)) + mn 
    bbeta[M+1, ]  <- t(Mb) 
    bb            <- diag( sqrt( sigma.sq*tau.sq/(A.K + sigma.sq) ) )
    meta          <- (1/sigma.sq) * (bb^2) %*% (Z-X%*%bbeta[M+1,])[A.K] 
	eta           <- t(bb) %*% array(rnorm(K), c(K, 1)) + meta 
    tau.sq        <- rinvgamma(1,shape=(K)/2+a1,scale=sum(eta^2)/2+b1)
    # If tau.sq is infinity; 1.7e+100 is reasonable big?!
    if(is.finite(tau.sq)==FALSE) tau.sq <- 1.7e+100
    tau           <- sqrt(tau.sq)
    tau.sq.post   <- c(tau.sq.post,tau.sq)

    #print("UPDATE 'm'") 
    m.hat.s <- m.hess.s <- L.m.s.hat <- NULL
    mle.m         <- optim(par=initial.m, fn=loglike, method ="L-BFGS-B", lower = 0.01, upper = Inf, 
		          hessian = TRUE, control=list(fnscale=-1))
    m.hat         <- mle.m$par
    m.hessian     <- mle.m$hessian 
    L.m.hat       <- loglike(m.hat)
    for (nu in 1:2){ 
 	mle.m.s   <- optim(par=initial.m, fn=loglike.s, method ="L-BFGS-B", lower = 0.01, upper = Inf, 
		          hessian = TRUE, control=list(fnscale=-1))
	m.hat.s   <- c(m.hat.s, mle.m.s$par)
	m.hess.s  <- c(m.hess.s, mle.m.s$hessian) 
	Lms.hat   <- loglike.s(mle.m.s$par)
	L.m.s.hat <- c(L.m.s.hat, Lms.hat)
    }
        
    mean.m <- sqrt(m.hessian/m.hess.s[1])*exp(L.m.s.hat[1] - L.m.hat)
    var.m  <- sqrt(m.hessian/m.hess.s[2])*exp(L.m.s.hat[2] - L.m.hat) - mean.m^2
    Sha    <- (mean.m^2)/var.m
    Sca    <- var.m/mean.m  
    cand   <- rgamma(1, shape=Sha, scale=Sca) 
    m      <- cand*rho + m*(1-rho)
    m.save <- c(m.save,m) 
   if (M %% 100 == 0) print(paste("finished iteration",M))
}
return(bbeta)

} # END TO glmdm FUNCTION CALL

