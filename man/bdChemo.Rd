\name{bdChemo}
\alias{bdChemo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Estimations of cell birth and death rates under drug interventions
}
\description{
	Generate posterior estimations (mean and credible intervals) of dose-specific cell birth and death rates, and dose-response curve summary statistics \eqn{GI_{50}}, \eqn{TGI}, \eqn{LC_{50}} and \eqn{IC_{50}}.
}
\usage{
bdChemo(x0, xtc, z, xt, bk,
Niter = c(1e5,2e5,3e5,1e6,1.5e6,2e6,3e6,4e6,5e6,6e6,7e6,8e6,9e6,1e7,3e7),
Nburn = 2e7, thin = 1e4, sa2l = .5, sa2m = 1, at = 10, bt = 1, al = 25, bl = 6,
ql = .025, qu = .975, plot.name = NULL, sample.return = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x0}{
     initial cell population size; a scalar (same for all concentrations)
}
  \item{xtc}{
     control group (without compound treatment) cell population size at follow-up time
}
  \item{z}{
     compound concentrations (in \eqn{\log_{10}} scale); a vector
}
  \item{xt}{
     cell population sizes at follow-up time corresponding to concentrations in z
}
  \item{bk}{
     independent background noise measurements
}
  \item{Niter}{
     number of iterations for different sampling stages (require at least 6 stages: only samples in the last stage will be used for birth and death estimations)
}
  \item{Nburn}{
     number of burn-in samples (must be greater than the second last element of Niter)
}
  \item{thin}{
     thinning parameter
}
  \item{sa2l}{
     prior normal distribution variance for \eqn{\alpha_\lambda}
}
  \item{sa2m}{
     prior normal distribution variance for \eqn{\alpha_\mu}
}
  \item{at}{
     prior inverse-gamma distribution shape parameter for \eqn{\tau_\lambda^2} and \eqn{\tau_\mu^2}
}
  \item{bt}{
     prior inverse-gamma distribution rate parameter for \eqn{\tau_\lambda^2} and \eqn{\tau_\mu^2}
}
  \item{al}{
     prior gamma distribution shape parameter for \eqn{l_\lambda^2} and \eqn{l_\mu^2}
}
  \item{bl}{
     prior gamma distribution rate parameter for \eqn{l_\lambda^2} and \eqn{l_\mu^2}
}
  \item{ql}{
     lower bound of credible interval intended to output
}
  \item{qu}{
     upper bound of credible interval intended to output
}
  \item{plot.name}{
     file name for cell growth, birth and death rates plots; if NULL (default), do not output any figure
}
  \item{sample.return}{
  	 whether posterior samples of \eqn{\lambda} and \eqn{\lambda} should be returned
  }
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
 A list with the following elements
	\item{hyperparameters}{
	posterior means and credible intervals of \eqn{\alpha_\lambda}, \eqn{\alpha_\mu}, \eqn{\tau_\lambda^2}, \eqn{\tau_\mu^2}, \eqn{\sigma^2}, \eqn{l_\lambda^2}, \eqn{l_\mu^2} and \eqn{\theta}; a matrix (columns for different variables)		
}
	\item{znew}{
	points where growth rates (means of Kendall Processes), birth and death rates are estimated
}
	\item{lambdas}{
	posterior means and credible intervals of birth rate \eqn{\lambda(z)}; a matrix (colums for different concentrations in \code{znew})
}
	\item{mus}{
	posterior means and credible intervals of death rate \eqn{\mu(z)}; a matrix (colums for different concentrations in \code{znew})
}
	\item{kmean}{
	posterior means and credible intervals of growth rate (Kendal Process mean) \eqn{m(z)}; a matrix (colums for different concentrations in \code{znew})
}
	\item{summary}{
	posterior means and credible intervals of \eqn{GI_{50}}, \eqn{TGI}, \eqn{LC_{50}}, \eqn{IC_{50}}; a matrix (columns for different variables)
}
	\item{post.lambda}{
	posterior samples of \eqn{\lambda} at \code{znew}, only returned if \code{sample.return = T}
}
	\item{post.mu}{
	posterior samples of \eqn{\mu} at \code{znew}, only returned if \code{sample.return = T}
}
 If "plot.name" is not NULL, also return plots of cell growth, birth and death estimations as shown in the main manuscript (curves and shades for mean and credible intervals, respectively)
}
\references{
 Liu Y. and Crawford F., Estimating dose-specific cell division and apoptosis rates from chemo-sensitivity experiments, submitted.

}
\author{
  Yiyi Liu
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.
library(bdChemo)
data(example)
bd.example = bdChemo(x0, xtc, z, xt, bk,
Niter = c(1e5,2e5,3e5,1e6,1.5e6,2e6,3e6,4e6,5e6,6e6,7e6,8e6,9e6,1e7,3e7),
Nburn = 2e7, thin = 1e4, sa2l = .5, sa2m = 1, at = 10, bt = 1, al = 25, bl = 6,
ql = .025, qu = .975, plot.name = "example.pdf", sample.return = F)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ chemosensitivity study }
\keyword{ Kendall Process }% __ONLY ONE__ keyword per line
