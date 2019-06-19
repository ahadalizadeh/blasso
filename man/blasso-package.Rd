\name{blasso-package}
\alias{blasso-package}
\docType{package}
\title{
blasso: MCMC for Bayesian Lasso Regression
}
\description{
Three Gibbs samplers for the Bayesian Lasso regression model. Two of the Gibbs samplers - the basic and orthogonalized samplers - fit the \dQuote{full} model that uses all predictor variables. The variable selection Gibbs sampler accommodates model uncertainty by fitting a model that uses a point-mass/double-exponential mixture prior on the regression coefficients.
}
\details{
\tabular{ll}{
Package: \tab blasso\cr
Type: \tab Package\cr
Version: \tab 0.3\cr
Date: \tab 2010-07-08\cr
License: \tab GPL-2\cr
LazyLoad: \tab yes\cr
}
}
\author{
Chris Hans <hans@stat.osu.edu>

Maintainer: Chris Hans <hans@stat.osu.edu>
}
\references{
Hans, C. (2009) Bayesian lasso regression.  Biometrika, 96, 835-845.\cr
Hans, C. (2010) Model uncertainty and variable selection in Bayesian lasso regresion.  Statistics and Computing, 20, 221-229.\cr
\url{http://www.stat.osu.edu/~hans}
}
\keyword{models}
\keyword{regression}
