\name{blasso}
\alias{blasso}
\title{Gibbs sampler for Bayesian Lasso regression model.}
\description{
  \code{blasso} implements MCMC for the Bayesian lasso regression model.
}
\usage{
blasso(Y, X, iters, burn = 0, thin = 1, sampler = c("basic", "orthogonalized"), beta, 
   sig2, tau, beta.prior = c("scaled", "classic"), fixsig2 = FALSE, sig2prior = NULL, 
   fixtau = FALSE, tauprior = NULL, noisy = TRUE)
}
\arguments{
  \item{Y}{An \eqn{n \times 1}{n x 1} vector of response data.
     WARNING: \code{blasso} centers \code{Y} to have sample mean zero before the MCMC is run!}
  \item{X}{An \eqn{n \times p}{n x p} matrix of predictor variables.
     WARNING: \code{blasso} centers the columns of \code{X} to have sample mean zero 
	before the MCMC is run, and they are rescaled to have unit sample variance!}
  \item{iters}{Length of desired MCMC output}
  \item{burn}{Number of MCMC iterations before output is saved}
  \item{thin}{Thinning for chain (1 = save every iteration, 2 = save every other 
     iteration, etc.)}
  \item{sampler}{The \dQuote{\code{basic}} sampler updates each component of 
     \eqn{\beta}{beta} one element at a time. The \dQuote{\code{orthogonalized}} 
	sampler transforms the regression coefficients to reduce autocorrelation.}
  \item{beta}{A \eqn{p \times 1}{p x 1} vector of starting values for the regression 
     coefficients}
  \item{sig2}{Starting value for \eqn{\sigma^2}{sigma^2} if \code{fixsig2} is 
     \code{FALSE}, otherwise the value at which to fix \eqn{\sigma^2}{sigma^2}}
  \item{tau}{Starting value for \eqn{\tau}{tau} if \code{fixtau} is \code{FALSE}, 
     otherwise the value at which to fix \eqn{\tau}{tau}}
  \item{beta.prior}{The \dQuote{\code{classic}} prior for \eqn{\beta}{beta} does not 
     depend on \eqn{\sigma^2}{sigma^2}. The \dQuote{\code{scaled}} prior for 
	\eqn{\beta}{beta} contains the term \eqn{\tau/\sigma}{tau/sigma}}
  \item{fixsig2}{\code{TRUE} if \eqn{\sigma^2}{sigma^2} is to be fixed}
  \item{sig2prior}{Parameters for inverse-gamma prior on \eqn{\sigma^2}{sigma^2} 
     (when it is not fixed). The prior is parameterized so that \code{(a,b)} corresponds 
	to a prior expected value of \eqn{b/(a-1)}{b/(a-1)}. The default is (0,0), an 
	improper prior on \eqn{\sigma^2}{sigma^2}.}
  \item{fixtau}{\code{TRUE} if \eqn{\tau}{tau} is to be fixed}
  \item{tauprior}{Parameters for gamma prior on \eqn{\tau}{tau}. The prior is 
     parameterized so that \code{(r,s)} corresponds to a prior mean of \eqn{r/s}{r/s}.}
  \item{noisy}{If \code{TRUE}, prints output describing sampler progress when R is run 
     from the command line}
}
\details{
  NOTE: \code{blasso} mean centers both Y and the columns of X, and standardizes the 
     predictors X to have unit sample variance.  Please adjust your models and priors 
	accordingly!  No intercept term is included in the model.  The function 
	\code{blasso} will only obtain samples from models with p <= n predictor variables.
	For the p > n case, see the function \code{blasso.vs}.
}
\value{
  \item{X}{The mean-centered and rescaled design matrix}
  \item{Y}{The mean-centered response data}
  \item{beta}{A matrix with MCMC samples of \eqn{\beta}{beta}}
  \item{sig2}{If \code{fixsig2} is \code{FALSE}, a vector with MCMC samples of 
     \eqn{\sigma^2}{sigma^2}; otherwise the value at which \eqn{\sigma^2}{sigma^2} was fixed.}
  \item{tau}{If \code{fixtau} is \code{FALSE}, a vector with MCMC samples of 
     \eqn{\tau}{tau}; otherwise the value at which \eqn{\tau}{tau} was fixed.}
  \item{sampler}{The type of sampler that was run}
  \item{sig2prior}{The parameters for the prior on \eqn{\sigma^2}{sigma^2} (if not 
     fixing \eqn{\sigma^2}{sigma^2})}
  \item{tauprior}{The parameters for the prior on \eqn{\tau}{tau} (if not fixing 
     \eqn{\tau}{tau})}
  \item{beta.prior}{The prior used for \eqn{\beta}{beta}.}
  \item{H}{If using the orthogonalized sampler, \code{H} is the matrix of eigenvectors 
     of \eqn{(X'X)^{-1}}{(X'X)^{-1}}}
  \item{Lambda}{If using the orthogonalized sampler, \code{Lambda} contains the 
     eigenvalues of \eqn{(X'X)^{-1}}{(X'X)^{-1}}}
}
\references{Hans, C (2009) Bayesian Lasso Regression. Biometrika, 96, 835-845. \cr
	\url{http://www.stat.osu.edu/~hans}
}
\author{Chris Hans <hans@stat.osu.edu>}
\note{Please contact Chris Hans <hans@stat.osu.edu> with comments or suggestions.
}
\seealso{\code{\link{blasso.vs}}}
\examples{
   ## Generate some fake data
   X <- scale(matrix(rnorm(30),ncol=3))[,]
   Y <- matrix(rnorm(10, X\%*\%matrix(c(-0.2,0.5,1.5),ncol=1), sd=0.8),ncol=1)

   # NOTE: Y and the columns of X are mean centered by the blasso function,
   #       and the columns of X and scaled to have unit variance
   #       (they already do in this example)
   # Use the basic sampler
   mcmc.basic <- blasso(Y, X, iters=2000, burn=100, thin=2, beta=c(0,0,0), sig2=1, tau=1, tauprior=c(1,1))


   # Use the orthogonalized sampler
   mcmc.orth <- blasso(Y, X, iters=2000, burn=100, thin=2, beta=c(0,0,0), sig2=1, tau=1, tauprior=c(1,1), sampler="orth")

}
\keyword{models}
\keyword{regression}
