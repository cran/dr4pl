
#' @importFrom generics calculate
#' @export
generics::calculate

#' @title dr4pl-calculate
#' @description calculate various useful statistics.
#' @param x an object of class `dr4pl`
#' @param parm  parameters of the dr4pl object. Usually made with [dr4pl_theta]
#' @param level confidence level to calculate. Defaults to 0.95
#' @param ... extra arguments to be passed to [vcov.dr4pl]
#' @export
calculate.dr4pl <- function(x, parm = NULL, level = 0.95, ...) {
  theta <- parm %theta% ParmToLog(coef(x))
  dose <- x$data$Dose
  response <- x$data$Response
  calculate.dr4pl_param(theta, dose, response, n = x$sample.size, level = level, ...)
}



calculate.dr4pl_param <- function(x, dose, response, n, p = 4, level = 0.95, ...) {
  #Based on Saber Wiley section 5.1
  # Using 1/2H or Jacobean to Covariance matrix
  # this matrix, along with the sqrt of the error sum of square
  # can approximate the confidence intervals
  vcov.mat <- vcov(x, dose, response, ...)
  resid <- residuals(x, dose, response)
  ESS <- sqrt(sum(resid^2)/(n-p))
  q.t <- qt(1 - (1 - level)/2, df = n - p)
  std.err <- ESS*sqrt(diag(vcov.mat))
  ci.table <- cbind(x - q.t*std.err, x + q.t*std.err)

  colnames(ci.table) <- c(paste(100*(1 - level)/2, "%"),
                          paste(100*(1 - (1 - level)/2), "%"))

  row.names(ci.table) <- param_names(x)

  res <- list(
    est = x,
    vcov = vcov.mat,
    resid = resid,
    ESS = ESS,
    std.err = std.err,
    ci.table = ci.table
  )
  return(res)
}
