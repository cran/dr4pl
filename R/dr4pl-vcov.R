#' @title Obtain the variance-covariance matrix of the parameter estimators of a
#' 4PL model.
#'
#' @description This function obtains the variance-covariance matrix of the parameter
#' estimators of a 4PL model. The variance-covariance matrix returned by this
#' function can be used to compute the standard errors and confidence intervals
#' for statistical inference.
#'
#' @name vcov.dr4pl
#'
#' @param object An object of the dr4pl class.
#' @param parm parameters of the dr4pl object. Usually made with [dr4pl_theta].
#' The class of this object determines in which space the covariance is 
#' calculated for theta_2.
#' @param ... dots for future extensions
#'
#' @return The variance-covariance matrix of the parameter estimators of a 4PL
#' model whose columns are in the order of the upper asymptote, IC50, slope and lower
#' asymptote from left to right and whose rows are in the same order.
#'
#' @details This function obtains the variance-covariance matrix of the parameter
#' estimators of a 4PL model. The Hessian matrix is used to obtain the second order
#' approximation to the sum-of-squares loss function, and then the standard errors
#' are computed as the square roots of the half of the Hessian matrix. Please refer
#' to Subsection 5.2.2 of Seber and Wild (1989).
#'
#' @examples
#' obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_1)  # Fit a 4PL model to data
#' vcov(obj.dr4pl)  # Variance-covariance matrix of the parameters
#'
#' obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_2)  # Fit a 4PL model to data
#' vcov(obj.dr4pl)  # Variance-covariance matrix of the parameters
#'
#' @references
#' \insertRef{Seber1989}{dr4pl}
#'
#' @return a covariance matrix. If the `parm` argument is of the class `dr4pl_log10`, then the
#' covariance of row/column 2 represents log10(theta_2). If theta is of `dr4pl_theta`,
#' then the covariance of row/column 2 represents theta_2 in linear space.
#'
#' @export
vcov.dr4pl <- function(object, parm = NULL, use.Hessian = T, ...) {
  
  theta <- parm %theta% ParmToLog(coef(object))
  x <- object$data$Dose  # Vector of dose levels
  y <- object$data$Response  # Vector of responses
  valid_dr4pl_param(theta)
  vcov.dr4pl_param(theta, x, y, use.Hessian = use.Hessian)
}



#' @rdname vcov.dr4pl
#' @param dose dose levels
#' @param response response values
#' @param use.Hessian logical, if set to TRUE, the default, then
#' the Hessian matrix scaled by 1/2 is used as an approximation
#' to C.hat. Otherwise the First order Jacobian is used instead.
#' @export
vcov.dr4pl_param <- function(object, dose, response, use.Hessian = T, ...) {
  C.hat <- Hessian(object, dose, response)/2

  if (!use.Hessian || !is_positive_semi_definite(C.hat)) {
    Jacobian <- DerivativeF(object, dose)
    C.hat <- t(Jacobian)%*%Jacobian
  }

  C.hat <- Matrix::nearPD(C.hat)$mat

  ind.mat.inv <- TRUE  # TRUE if matrix inversion is successful, FALSE otherwise
  vcov.mat <- try(solve(C.hat), silent = TRUE)  # Inverse matrix

  # If matrix inversion is infeasible, proceed with the Cholesky decomposition.
  if(inherits(vcov.mat, "try-error")) {

    vcov.Chol <- try(chol(C.hat, silent = TRUE))  # Cholesky decomposition

    # If the Cholesky decomposition is infeasible, use an approximated positve
    # definite Hessian matrix to obtain the variance-covariance matrix.
    if(inherits(vcov.Chol, "try-error")) {

      ind.mat.inv <- FALSE
      C.hat.pd <- Matrix::nearPD(C.hat)$mat/2
      vcov.mat <- solve(C.hat.pd)
    } else {

      vcov.mat <- chol2inv(vcov.Chol)
    }
  }

  if(!ind.mat.inv) {

    print("The hessian matrix is singular, so an approximated positive definite
          hessian matrix is used.")
  }

  colnames(vcov.mat) <- param_names(object)
  row.names(vcov.mat) <- param_names(object)
  return(vcov.mat)
}
