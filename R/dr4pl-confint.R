
#' @title Fit a 4 parameter logistic (4PL) model to dose-response data.
#'
#' @description Compute the approximate confidence intervals of the parameters of a
#' 4PL model based on the asymptotic normality of least squares estimators.
#'
#' @name confint.dr4pl
#'
#' @param object An object of the dr4pl class
#' @param parm parameters of the dr4pl object. Usually made with [dr4pl_theta]
#' @param level Confidence level
#' @param ... Other parameters to be passed to vcov
#'
#' @return A matrix of the confidence intervals in which each row represents a
#' parameter and each column represents the lower and upper bounds of the
#' confidence intervals of the corresponding parameters.
#'
#' @details This function computes the approximate confidence intervals of the
#' true parameters of a 4PL model based on the asymptotic normality of the least
#' squares estimators in nonlinear regression. The Hessian matrix is used to
#' obtain the second order approximation to the sum-of-squares loss function.
#' Please refer to Subsection 5.2.2 of Seber and Wild (1989).
#'
#' @examples
#' obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_1)  # Fit a 4PL model to data
#'
#' ## Use the data 'sample_data_1' to obtain confidence intervals.
#' confint(obj.dr4pl)  # 95% confidence intervals
#' confint(obj.dr4pl, level = 0.99)  # 99% confidence intervals
#'
#' theta <- FindInitialParms(x = sample_data_1$Dose, y = sample_data_1$Response)
#'
#' # Use the same data 'sample_data_1' but different parameter estimates to obtain
#' # confidence intervals.
#' confint(obj.dr4pl, parm = theta)
#'
#' @references
#' \insertRef{Seber1989}{dr4pl}
#'
#' @export
confint.dr4pl <- function(object, parm = NULL, level = 0.95, ...) {
  
  theta <- parm %theta% ParmToLog(coef(object))
  x <- object$data$Dose
  y <- object$data$Response
  calculate(theta, x, y, level = level, n = object$sample.size, ...)$ci.table
}






