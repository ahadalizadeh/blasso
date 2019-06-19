`blasso` <-
function(Y, X, iters, burn=0, thin=1, sampler = c("basic", "orthogonalized"), beta, sig2, tau, beta.prior = c("scaled", "classic"), fixsig2 = FALSE, sig2prior = NULL, fixtau = FALSE, tauprior = NULL, noisy = TRUE)#, norm = TRUE)
{
   call <- match.call()
   sampler <- match.arg(sampler)
   beta.prior <- match.arg(beta.prior)

   ## The prior on beta
   if(beta.prior=="scaled") {
      bprior <- 1
   } else {
      bprior <- 0
   }

   ## Do some error checking
   if(!is.matrix(X)) stop("X must be of class matrix")

   if(sig2 <= 0) stop("Parameter sig2 must be > 0")
   if(tau <= 0) stop("Parameter tau must be > 0")

   ## nonpositive values not allowed
   thin <- max(thin,1)

   n <- length(Y)
   p <- dim(X)[2]

   if(p>=n) stop("This Gibbs sampler is not implemented for the p >= n case. See function 'blasso.vs'")

   ## Demean the Y's
   Y <- Y - mean(Y)

   ## Rescale the Y's to have unit sample variance
   #Y <- Y/sqrt(var(Y))[1]

   ## Demean and rescale the X's to have unit sample variance
   X <- scale(X)[,]

   ## Demean the X's
   #X <- scale(X, TRUE, FALSE)[,]

   ## If requested (this is the default), scale the X's
   #if(norm==TRUE)
   #{
   #   if(sum(round(diag(t(X)%*%X))==1)!=p) 
   #   X <- scale(X, FALSE, TRUE)[,]/sqrt(n-1)
   #}

   ## Deal with sig2 prior
   if(fixsig2 == TRUE) {
      sig2prior <- c(NA, NA)
   } else {
      if(is.null(sig2prior)) sig2prior <- c(0,0)
   }

   ## Deal with tau prior
   if(fixtau == TRUE) {
      tauprior <- c(NA, NA)
   } else {
      if(is.null(tauprior)) stop("If not fixing tau, you must specify the prior")
   }

   if(sampler == "basic") {
      output <- .C("blassoGibbs", 
	 		    X = as.double(X), 
			    nn = as.integer(n),
			    pp = as.integer(p),
			    TT = as.integer(iters),
			    BB = as.integer(burn),
			    tthin = as.integer(thin),
			    ttau = as.double(tau),
			    ssig2 = as.double(sig2),
			    pphi = as.double(1.0),
                   sig2prior = as.double(sig2prior),
                   fits2 = as.integer(!fixsig2),
			    tauprior = as.double(tauprior),
			    fittau = as.integer(!fixtau),
			    modunc = as.integer(FALSE),
			    phiprior = as.double(1,1),
			    fitphi = as.integer(FALSE),
			    start = as.double(beta),
			    bdraws = as.double(rep(0,p*iters)),
			    sig2draws = as.double(rep(0,(!fixsig2)*iters + fixsig2*1)),
			    taudraws = as.double(rep(0,(!fixtau)*iters + fixtau*1)),
			    phidraws = as.double(0),
			    marginc = as.double(0),
			    rrb = as.integer(FALSE),
			    RB = as.double(rep(0,p)),
			    YtY = as.double(t(Y)%*%Y),
			    YtX = as.double(t(Y)%*%X),
			    NOISY = as.integer(noisy),
			    bprior = as.integer(bprior),
			    count = as.integer(0),
			    NAOK = TRUE,
			    package="blasso")

      ## Clean up the output
	 output[["X"]] <- matrix(output$X, ncol=p, byrow=F)
      output[["YtX"]] <- output[["RB"]] <- output[["rrb"]] <- output[["phidraws"]] <- output[["marginc"]] <- output[["fitphi"]] <- output[["phiprior"]] <- output[["modunc"]] <- output[["pphi"]] <- output[["tthin"]] <- NULL

   } else {
      XtXi <- solve(t(X)%*%X)
	 edecomp <- eigen(XtXi)
	 H <- edecomp$vectors
	 Lambda <- edecomp$values

	 if(any(H==0)) stop("An element of the H matrix is zero -- the sampler has not been implemented to accommodate this")
	 if(any(Lambda<=0)) stop("An eigenvalue is <= 0 -- the sampler has not been implemented to accommodate this")

	 yty <- t(Y)%*%Y
	 htxty <- t(H)%*%t(X)%*%Y

	 betahat <- XtXi%*%t(X)%*%Y

      output <- .C("blassoGibbsTrans",
	 		    nn = as.integer(n),
                   pp = as.integer(p),
			    TT = as.integer(iters),
			    BB = as.integer(burn),
			    tthin = as.integer(thin),
			    ttau = as.double(tau),
			    ssig2 = as.double(sig2),
			    sig2prior = as.double(sig2prior),
			    fits2 = as.integer(!fixsig2),
			    tauprior = as.double(tauprior),
			    fittau = as.integer(!fixtau),
			    start = as.double(beta),
			    bdraws = as.double(rep(0,p*iters)),
			    sig2draws = as.double(rep(0,(!fixsig2)*iters + fixsig2*1)),
			    taudraws = as.double(rep(0,(!fixtau)*iters + fixtau*1)),
			    bhat = as.double(betahat),
			    HH = as.double(H),
			    Lambda = as.double(Lambda),
			    YtY = as.double(yty),
			    HtXtY = as.double(htxty),
			    NOISY = as.integer(noisy),
			    bprior = as.integer(bprior),
			    count = as.double(NA),
			    NAOK = TRUE,
			    PACKAGE="blasso")

      ## Clean up the output
	 output[["bhat"]] <- output[["HH"]] <- output[["HtXtY"]] <- NULL

	 output[["H"]] <- H
	 output[["X"]] <- X
   }

   ## Clean up the output
   output[["nn"]] <- output[["pp"]] <- output[["YtY"]] <- output[["NOISY"]] <- output[["count"]] <- output[["start"]] <- output[["BB"]] <- output[["TT"]] <- output[["fits2"]] <- output[["fittau"]] <- output[["tthin"]] <- NULL

   output[["Y"]] <- Y

   if(fixsig2) { output[["sig2draws"]] <- output[["sig2prior"]] <- NULL }
   if(fixtau) { output[["taudraws"]] <- output[["tauprior"]] <- NULL }

   ## Organize some output
   output$bdraws <- matrix(output$bdraws, ncol=p, byrow=T)

   ## Rename some output (should be a better way to do this)
   names(output)[names(output) == "bdraws"] <- "beta"

   if(fixtau)
   {
      names(output)[names(output) == "ttau"] <- "tau"
   } else {
      names(output)[names(output) == "taudraws"] <- "tau"
	 output[["ttau"]] <- NULL
   }

   if(fixsig2)
   {
      names(output)[names(output) == "ssig2"] <- "sig2"
   } else {
      names(output)[names(output) == "sig2draws"] <- "sig2"
	 output[["ssig2"]] <- NULL
   }

   if(output$bprior == 1)
   {
      output$bprior <- "scaled"
   } else { 
      output$bprior <- "classic" 
   }
   names(output)[names(output) == "bprior"] <- "beta.prior"

   output[["sampler"]] <- sampler


   return(output)
}

