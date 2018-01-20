\name{example}
\docType{data}
\alias{example}
\title{Example chemosensitivity data}
\description{
	This dataset gives an example of inputs needed for bdChemo function.
}
\usage{data(example)}
\format{A list containing the following elements:

z: compound concentrations (in \eqn{\log_{10}} scale); a vector

x0: initial cell population sizes corresponding to concentrations in z; a vector (same length as z)

xt: cell population sizes at follow-up time corresponding to concentrations in z; a vector (same length as z)

x0c: control group (without compound treatment) initial cell population size; a scaler (if there are replicates, take average); used only for computing GI50

xtc: control group (without compound treatment) cell population size at follow-up time; a scaler (if there are replicates, take average); used only for computing GI50

bk: background noise measurements; a vector
}
\keyword{datasets}