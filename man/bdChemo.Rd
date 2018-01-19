\name{bdChemo}
\alias{bdChemo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	Estimations of cell birth and death rates under drug interventions
}
\description{
	Generate posterior estimations (mean and credible intervals) of dose-specific cell birth and death rates, and dose-response curve summary statistics \eqn{GI_{50}}, \eqn{TGI} and \eqn{LC_{50}}}.
}
\usage{
bdChemo(z, x0, xt, x0c, xtc, bk, N = 3e7, sa2l = .5, sa2m = 1, at = 10, bt = 1, al = 25, bl = 6, ql = .025, qu = .975, plot.name = NULL, sample.return = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{z}{
     compound concentrations (in \eqn{\log_{10}} scale); a vector
}
  \item{x0}{
     initial cell population sizes corresponding to concentrations in z; a vector (same length as z)
}
  \item{xt}{
     cell population sizes at follow-up time corresponding to concentrations in z; a vector (same length as z)
}
  \item{xtc}{
     control group (without compound treatment) initial cell population size; a scaler (if there are replicates, take average); used only for computing GI50
}
  \item{xtc}{
     control group (without compound treatment) cell population size at follow-up time; a scaler (if there are replicates, take average); used only for computing GI50
}    
  \item{bk}{
     background noise measurements
}
  \item{N}{
     number of iterations;
}
  \item{sa2l}{
     prior normal distribution variance for \eqn{\alpha_\lambda}; sa2l controls the scale of \eqn{\phi_\lambda} in prior;
}
  \item{sa2m}{
     prior normal distribution variance for \eqn{\alpha_\mu}; controls the scale of \eqn{\phi_\mu} in prior;
}
  \item{at}{
     prior inverse-gamma distribution shape parameter for \eqn{\tau_\lambda^2} and \eqn{\tau_\mu^2};
}
  \item{bt}{
     prior inverse-gamma distribution rate parameter for \eqn{\tau_\lambda^2} and \eqn{\tau_\mu^2}; at and bt together control the scale of \eqn{\tau^2} in prior, which affects the fluctuation amplitude of the estimated curve (\eqn{\tau} is called amplitude parameter of the Gaussian kernel and the larger \eqn{\tau^2} is, the larger the amplitude tends to be)
}
  \item{al}{
     prior gamma distribution shape parameter for \eqn{l_\lambda^2} and \eqn{l_\mu^2}
}
  \item{bl}{
     prior gamma distribution rate parameter for \eqn{l_\lambda^2} and \eqn{l_\mu^2}; al and bl together control the scale of \eqn{l^2} in prior, which affects the smoothness of the estimated curve (\eqn{l} is called the length scale parameter of the Gaussian kernel, the larger \eqn{l^2} is, the smoother the curve tends to be) 
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
	posterior means and credible intervals of \eqn{GI_{50}}, \eqn{TGI}, \eqn{LC_{50}}; a matrix (columns for different variables)
}
	\itemn{post.lambda}{
	posterior samples of \eqn{\lambda} at \code{znew}, only returned if \code{sample.return = T}
}
	\itemn{post.mu}{
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
bd.example = bdChemo(z, x0, xt, x0c, xtc, bk, N = 1e5, plot.name = "example")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ chemosensitivity study }
\keyword{ Kendall Process }% __ONLY ONE__ keyword per line
