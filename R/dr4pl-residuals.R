
#' @title Compute dr4pl residuals.
#'
#' @rdname residuals-dr4pl
#' @param object A dr4pl or dr4pl_param object.
#' @param parm parameters of the dr4pl object. Usually made with [dr4pl_theta]
#' @param ... dots for future extensions
#'
#' @return Vector of residuals.
#' @export
residuals.dr4pl <- function(object, parm = NULL, ...) {
  theta <- parm %theta% coef(object)
  residuals(theta, object$data$Dose, object$data$Response)
}

#' @rdname residuals-dr4pl
#' @param dose Dose values
#' @param response Response values
#' @param ... dots for future extensions
#' @export
residuals.dr4pl_param <- function(object, dose, response, ...) {
  return(response - MeanResponse(object, dose))
}

