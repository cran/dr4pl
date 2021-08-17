
# @title First Order Derivatives of the 4PL Model
# @name DerivativeF the derivative values of the mean response function.
# @param object Either a dr4pl object, or an
# object that can be coerced into a dr4pl_param object
# @param x Dose values
#
# @keywords internal
# @return Derivative values of the mean response function.
DerivativeF <- function(object, ...) UseMethod("DerivativeF")

# @rdname DerivativeF
#
DerivativeF.dr4pl <- function(object, theta = NULL){
  theta <- theta %theta% ParmToLog(coef(object))
  DerivativeF(theta, object$data$Dose)
}

# @rdname DerivativeF
DerivativeF.dr4pl_theta <- function(object, x) {
  theta <- object
  theta.1 <- theta[1]
  theta.2 <- theta[2]
  theta.3 <- theta[3]
  theta.4 <- theta[4]

  rho <- (x/theta.2)^theta.3

  ### Compute derivatives
  deriv.f.theta.1 <- 1 - 1/(1 + rho)
  deriv.f.theta.2 <- (theta.4 - theta.1)*theta.3/theta.2*rho/(1 + rho)^2
  deriv.f.theta.3 <- -(theta.4 - theta.1)*log(x/theta.2)*rho/(1 + rho)^2
  deriv.f.theta.4 <- 1/(1 + rho)

  # Handle limit cases
  deriv.f.theta.1[rho == Inf] <- 1
  deriv.f.theta.2[rho == Inf] <- 0
  deriv.f.theta.3[rho == Inf] <- 0
  deriv.f.theta.4[rho == Inf] <- 0

  deriv.f.theta.1[rho == 0] <- 0
  deriv.f.theta.2[rho == 0] <- 0
  deriv.f.theta.3[rho == 0] <- 0
  deriv.f.theta.4[rho == 0] <- 1

  deriv.f.theta <- cbind(deriv.f.theta.1, deriv.f.theta.2, deriv.f.theta.3, deriv.f.theta.4)

  # Check whether return values are appropriate
  if(anyNA(deriv.f.theta)) {
    abort(glue("Some of the derivative values are NA's.\n",
               "check values: {glue_collapse(x[Reduce(`|`, lapply(deriv.f.theta, is.na))], sep = ', ')}"))
  }

  return(deriv.f.theta)
}



# @rdname DerivativeF
DerivativeF.dr4pl_log10 <- function(object, x) {

  #cannot simply multiply by theta[2] because sometimes theta[2] in linear space
  # is too big and implies infinity - giving  NAN values
  #X is always given in linear space, convert to log space
  x <- log10(x)
  theta <- object
  theta.1 <- theta[1]
  theta.2 <- theta[2]
  theta.3 <- theta[3]
  theta.4 <- theta[4]

  rho <- 10^(theta.3*(x-theta.2))

  #Compute reparameterized derivatives
  deriv.f.theta.1 <- 1 - 1/(1 + rho)
  deriv.f.theta.2 <- (theta.4 - theta.1)*log(10)*theta.3*rho/((1 + rho)^2)
  deriv.f.theta.3 <- -(theta.4 - theta.1)*log(10)*rho*(x-theta.2)/((1 + rho)^2)
  deriv.f.theta.4 <-  1/(1 + rho)

  #Limit Cases
  deriv.f.theta.1[rho == Inf] <- 1
  deriv.f.theta.2[rho == Inf] <- 0
  deriv.f.theta.3[rho == Inf] <- 0
  deriv.f.theta.4[rho == Inf] <- 0

  deriv.f.theta.1[rho == 0] <- 0
  deriv.f.theta.2[rho == 0] <- 0
  deriv.f.theta.3[rho == 0] <- 0
  deriv.f.theta.4[rho == 0] <- 1

  deriv.f.theta <- cbind(deriv.f.theta.1,
                         deriv.f.theta.2,
                         deriv.f.theta.3,
                         deriv.f.theta.4)

  # Check whether return values are appropriate

  if(anyNA(deriv.f.theta))
    abort(glue("Some of the derivative values are NA's.\n",
               "check values: {glue_collapse(x[Reduce(`|`, lapply(deriv.f.theta, is.na))], sep = ', ')}"))

  return(deriv.f.theta)
}
