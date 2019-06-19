\name{blasso.vs}
\alias{blasso.vs}
\title{Variable selection Gibbs sampler for Bayesian Lasso regression model.}
\description{
  \code{blasso.vs} implements a variable selection Gibbs sampler for Bayesian lasso regression.
}
\usage{
blasso.vs(Y, X, iters, burn = 0, thin = 1, beta, sig2, tau, phi, 
   beta.prior = c("scaled", "classic"), fixsig2 = FALSE, sig2prior = NULL, fixtau = FALSE, 
   tauprior = NULL, fixphi = FALSE, phiprior = NULL, noisy = TRUE)
}
\arguments{
  \item{Y}{An \eqn{n \times 1}{n x 1} vector of response data.
     WARNING: \code{Y} is centered to have sample mean zero before the MCMC is run!}
  \item{X}{An \eqn{n \times p}{n x p} matrix of predictor variables.
     WARNING: the columns of \code{X} are centered to have sample mean zero before 
	the MCMC is run, and they are rescaled to have unit sample variance!}
  \item{iters}{Length of desired MCMC output}
  \item{burn}{Number of MCMC iterations before output is saved}
  \item{thin}{Thinning for chain (1 = save every iteration, 2 = save every other 
     iteration, etc.)}
  \item{beta}{A \eqn{p \times 1}{p x 1} vector of starting values for the regression 
     coefficients}
  \item{sig2}{Starting value for \eqn{\sigma^2}{sigma^2} if \code{fixsig2} is 
     \code{FALSE}, otherwise the value at which to fix \eqn{\sigma^2}{sigma^2}}
  \item{tau}{Starting value for \eqn{\tau}{tau} if \code{fixtau} is \code{FALSE}, 
     otherwise the value at which to fix \eqn{\tau}{tau}}
  \item{phi}{Starting values for \eqn{\phi}{phi} if \code{fixphi} is \code{FALSE}, 
     otherwise the value at which to fix \eqn{\phi}{phi}}
  \item{beta.prior}{The \dQuote{\code{classic}} prior for \eqn{\beta}{beta} does not 
     depend on \eqn{\sigma^2}{sigma^2}. The \dQuote{\code{scaled}} prior for 
	\eqn{\beta}{beta} contains the term \eqn{\tau/\sigma}{tau/sigma}}
  \item{fixsig2}{\code{TRUE} if \eqn{\sigma^2}{sigma^2} is to be fixed}
  \item{sig2prior}{Parameters for inverse-gamma prior on \eqn{\sigma^2}{sigma^2} (when 
     it is not fixed). The prior is parameterized so that \code{(a,b)} corresponds to a 
	prior expected value of \eqn{b/(a-1)}{b/(a-1)}. The default is (0,0), an improper 
	prior on \eqn{\sigma^2}{sigma^2}.}
  \item{fixtau}{\code{TRUE} if \eqn{\tau}{tau} is to be fixed}
  \item{tauprior}{Parameters for gamma prior on \eqn{\tau}{tau}. The prior is 
     parameterized so that \code{(r,s)} corresponds to a prior mean of \eqn{r/s}{r/s}.}
  \item{fixphi}{\code{TRUE} if \eqn{\phi}{phi} is to be fixed}
  \item{phiprior}{Parameters for beta distribution prior on \eqn{\phi}{phi}. The default 
     is (1,1), the uniform distribution.}
  \item{noisy}{If \code{TRUE} prints output describing sampler progress when R is run 
     from the command line}
}
\details{
  NOTE: \code{blasso.vs} mean centers both Y and the columns of X, and standardizes the 
  predictors X to have unit sample variance.  Please adjust your models and priors 
  accordingly!  No intercept term is included in the model.
}
\value{
  \item{X}{The mean-centered and rescaled design matrix}
  \item{Y}{The mean-centered response data}
  \item{beta}{A matrix with MCMC samples of \eqn{\beta}{beta}}
  \item{sig2}{If \code{fixsig2} is \code{FALSE}, a vector with MCMC samples of 
     \eqn{\sigma^2}{sigma^2}; otherwise the value at which \eqn{\sigma^2}{sigma^2} was fixed.}
  \item{tau}{If \code{fixtau} is \code{FALSE}, a vector with MCMC samples of 
     \eqn{\tau}{tau}; otherwise the value at which \eqn{\tau}{tau} was fixed.}
  \item{phi}{If \code{fixphi} if \code{FALSE}, a vector with MCMC samples of \eqn{\phi}{phi}}
  \item{marginc}{Rao-Blackwellized estimates of the marginal variable inclusion probabilities.}
  \item{sig2prior}{The parameters for the prior on \eqn{\sigma^2}{sigma^2} (if not 
     fixing \eqn{\sigma^2}{sigma^2})}
  \item{tauprior}{The parameters for the prior on \eqn{\tau}{tau} (if not fixing 
     \eqn{\tau}{tau})}
  \item{phiprior}{The parameters for the prior on \eqn{\phi}{phi} (if not fixing 
     \eqn{\phi}{phi})}
  \item{beta.prior}{The prior used for \eqn{\beta}{beta}.}
}
\references{Hans, C (2009) Model uncertainty and variable selection in Bayesian lasso regression. Statistics and Computing, 20, 221-229.\cr
	\url{http://www.stat.osu.edu/~hans}}
\author{Chris Hans <hans@stat.osu.edu>}
\note{Please contact Chris Hans <hans@stat.osu.edu> with comments or suggestions.
}
\seealso{\code{\link{blasso}}}
\examples{
   ## Generate some fake data
   X <- scale(matrix(rnorm(30),ncol=3))[,]
   Y <- matrix(rnorm(10, X\%*\%matrix(c(0,2,2.5),ncol=1),sd=0.5),ncol=1)

   # NOTE: Y and the columns of X are mean centered by the blasso.vs function,
   #       and the columns of X are rescaled to have unit sample variance
   #       (they already do in this example)
   
   mcmc.vs <- blasso.vs(Y, X, iters=10000, burn=100, thin=2, beta=c(1,1,1), sig2=1, tau=1, phi=0.5, tauprior=c(1,1))

}
\keyword{models}
\keyword{regression}
