

#' Compute an estimated mean response.
#'
#' @name MeanResponse
#'
#' @param ... arguments to be passed to S3 methods
#' @param dr4pl dr4pl object
#' @param theta Parameters of the dr4pl object. Usually made with [dr4pl_theta]
#' @param x domain values for 4PL model. Values
#' should always be passed to this function on the
#' linear space.
#'
#' @return Predicted response values.
#'
#' @export
MeanResponse <- function(...) UseMethod("MeanResponse")

#' @rdname MeanResponse
#'
#' @export
MeanResponse.dr4pl <- function(dr4pl, theta = NULL, ...) {

  theta <- theta %theta% coef(dr4pl)
  MeanResponse(theta, dr4pl$data$Dose)
}

#' @rdname MeanResponse
#'
#' @export
MeanResponse.numeric <- function(theta, x, ...) {
  theta <- new_dr4pl_theta(theta)
  f <- theta[1] + (theta[4] - theta[1])/(1 + (x/theta[2])^theta[3])

  return(f)
}
#' @rdname MeanResponse
#'
#' @export
MeanResponse.dr4pl_theta <- function(theta, x, ...) {
  #x is always passed in linear space,
  #this function requires the linear form
  # of x to calculate the correct mean response
  f <- theta[1] + (theta[4] - theta[1])/(1 + (x/theta[2])^theta[3])

  return(f)
}

#' @rdname MeanResponse
#'
#' @export
MeanResponse.dr4pl_log10 <- function(theta, x, ...) {
  #x is always passed in linear space
  #this function requires the log10 form
  # of x to calculate the correct mean response
  f <- theta[1] + (theta[4] - theta[1])/
    (1 + 10^(theta[3]*(log10(x) - theta[2])))

  return(f)
}
