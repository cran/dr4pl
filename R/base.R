
# Squares of residuals.
#
# @param r Residuals.
#
# @return Squared residuals.
SquaredLoss <- function(r) {
  return(r^2)
}

# Absolute values of residuals.
#
# @param r Residuals.
#
# @return Absolute valued residuals.
AbsoluteLoss <- function(r) {

  return(abs(r))
}

# Values of Huber's loss function evaluated at residuals r.
#
# @param r Residuals.
# @return Huber's loss function values evaluated at residuals r.
HuberLoss <- function(r) {

  # The value 1.345 was suggested by Huber (1964).
  # See Huber, P. J. (1964). Robust Estimation of a Location Parameter. Annals of Statistics 53(1)
  const <- 1.345

  ret.val <- r^2  # Vector of results
  outer.term <- 2*const*abs(r) - const^2

  outer.idx <- (abs(r) > const)

  ret.val[outer.idx] <- outer.term[outer.idx]

  return(ret.val)
}

# Values of Tukey's biweight loss function evaluated at residuals r.
#
# @param r Residuals.
#
# @return result: Tukey's biweight loss function values evaluated at residuals r.
TukeyBiweightLoss <- function(r) {

  # The value 4.685 was suggested by Tukey.
  const <- 4.685

  ret.val <- (r^6)/(const^4) - 3*(r^4)/(const^2) + 3*r^2

  ret.val[abs(r) > const] <- const^2

  return(ret.val)
}

# Returns an loss function for given robust fitting method.
#
# @param method.robust NULL, absolute, Huber, or Tukey
#      - NULL: Sum of squares loss
#      - absolute: Absolute deviation loss
#      - Huber: Huber's loss
#      - Tukey: Tukey's biweight loss
#
# This will always expect a dr4pl_log10
# @return Value of the sum of squared residuals.
ErrFcn <- function(method.robust) {

  ### Check whether function arguments are appropriate.
  if(!(is.null(method.robust)||method.robust == "absolute"||method.robust == "Huber"||
     method.robust == "Tukey"||method.robust == "squared")) {

    abort(
      glue('The robust estimation method should be one of "squared",',
           '"absolute","Huber\" or "Tukey".'))
  }

  loss.fcn <- c()

  if(is.null(method.robust)||method.robust=="squared") {
    loss.fcn <- SquaredLoss
  } else if(method.robust == "absolute") {
    loss.fcn <- AbsoluteLoss
  } else if(method.robust == "Huber") {
    loss.fcn <- HuberLoss
  } else if(method.robust == "Tukey") {
    loss.fcn <- TukeyBiweightLoss
  }

  err.fcn <- function(retheta, x, y) {

    ### Check whether function arguments are appropriate.
    if(length(retheta) != 4) {

      abort("The number of parameters is not 4.")
    }
    if(length(x) != length(y)) {

      abort("The numbers of dose levels and responses should be the same.")
    }

    n <- length(y)
    f <- MeanResponse.dr4pl_log10(retheta, x)

    if(anyNA(f)) {

      abort("Some of the evaluated function values are NA's.")
    }

    return(sum(loss.fcn(y - f))/n)
  }

  return(err.fcn)
}



# Compute gradient values of the sum-of-squares loss function.
#
# @param retheta Parameters among which the IC50 parameter is logarithmically
#   transformed.
# @param x Doses.
# @param y Responses.
#
# @return Gradient values of the sum-of-squares loss function.
GradientSquaredLossLogIC50 <- function(retheta, x, y) {

  f <- MeanResponse.dr4pl_log10(retheta, x)  # Mean response values
  n <- length(x)  # Number of data observations

  return(-2*(y - f)%*%DerivativeF.dr4pl_log10(retheta, x)/n)
}




