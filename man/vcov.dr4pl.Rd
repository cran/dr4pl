% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dr4pl-vcov.R
\name{vcov.dr4pl}
\alias{vcov.dr4pl}
\alias{vcov.dr4pl_param}
\title{Obtain the variance-covariance matrix of the parameter estimators of a
4PL model.}
\usage{
\method{vcov}{dr4pl}(object, parm = NULL, use.Hessian = T, ...)

\method{vcov}{dr4pl_param}(object, dose, response, use.Hessian = T, ...)
}
\arguments{
\item{object}{An object of the dr4pl class.}

\item{parm}{parameters of the dr4pl object. Usually made with [dr4pl_theta].
The class of this object determines in which space the covariance is 
calculated for theta_2.}

\item{use.Hessian}{logical, if set to TRUE, the default, then
the Hessian matrix scaled by 1/2 is used as an approximation
to C.hat. Otherwise the First order Jacobian is used instead.}

\item{...}{dots for future extensions}

\item{dose}{dose levels}

\item{response}{response values}
}
\value{
The variance-covariance matrix of the parameter estimators of a 4PL
model whose columns are in the order of the upper asymptote, IC50, slope and lower
asymptote from left to right and whose rows are in the same order.

a covariance matrix. If the `parm` argument is of the class `dr4pl_log10`, then the
covariance of row/column 2 represents log10(theta_2). If theta is of `dr4pl_theta`,
then the covariance of row/column 2 represents theta_2 in linear space.
}
\description{
This function obtains the variance-covariance matrix of the parameter
estimators of a 4PL model. The variance-covariance matrix returned by this
function can be used to compute the standard errors and confidence intervals
for statistical inference.
}
\details{
This function obtains the variance-covariance matrix of the parameter
estimators of a 4PL model. The Hessian matrix is used to obtain the second order
approximation to the sum-of-squares loss function, and then the standard errors
are computed as the square roots of the half of the Hessian matrix. Please refer
to Subsection 5.2.2 of Seber and Wild (1989).
}
\examples{
obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_1)  # Fit a 4PL model to data
vcov(obj.dr4pl)  # Variance-covariance matrix of the parameters

obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_2)  # Fit a 4PL model to data
vcov(obj.dr4pl)  # Variance-covariance matrix of the parameters

}
\references{
\insertRef{Seber1989}{dr4pl}
}
