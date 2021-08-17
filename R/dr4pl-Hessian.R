# Compute the Hessian matrix of the sum-of-squares loss function.
#
# @param theta Parameters.
# @param x Doses.
# @param y Response.
#
# @return Hessian matrix of the sum-of-squares loss function.

Hessian <- function(theta, x, y) UseMethod("Hessian")

Hessian.dr4pl_theta <- function(theta, x, y) {
  n <- length(x)  # Number of observations
  p <- length(theta)  # Number of parameters

  theta.1 <- theta[1]
  theta.2 <- theta[2]
  theta.3 <- theta[3]
  theta.4 <- theta[4]

  # Second order derivatives of f
  second.deriv.f <- array(data = 0, dim = c(p, p, n))

  eta <- (x/theta.2)^theta.3  # Term needed in the Hessian matrix computation

  deriv.eta.2 <- -theta.3/theta.2*eta
  deriv.eta.3 <- eta*log(x/theta.2)

  second.deriv.f[1, 1, ] <- 0
  second.deriv.f[1, 2, ] <- deriv.eta.2/(1+eta)^2
  second.deriv.f[1, 3, ] <- deriv.eta.3/(1+eta)^2
  second.deriv.f[1, 4, ] <- 0

  second.deriv.f[2, 2, ] <- -((eta*theta.3*(theta.4-theta.1))/((theta.2^2)*(1+eta)^2))*
                              (1+theta.3-(2*theta.3*eta/(1+eta)))
  second.deriv.f[2, 3, ] <- (((theta.4-theta.1)*eta)/(theta.2*((1+eta)^2))) *
                              (1 + theta.3*log(x/theta.2) - ((2*theta.3*deriv.eta.3)/(1+eta)))
  second.deriv.f[2, 4, ] <- theta.3/theta.2*eta/(1 + eta)^2

  second.deriv.f[3, 3, ] <- - (((log(x/theta.2)^2)*(theta.4-theta.1)*eta)/((1+eta)^2))*(1-(2*eta/(1+eta)))
  second.deriv.f[3, 4, ] <- -log(x/theta.2)*eta/(1 + eta)^2

  second.deriv.f[4, 4, ] <- 0

  second.deriv.f <- (second.deriv.f + aperm(second.deriv.f, c(2, 1, 3)))/2

  # Substitute limits for second derivatives when eta are infinite
  second.deriv.f[, , eta == 0] <- 0
  second.deriv.f[, , eta == Inf] <- 0

  deriv.f <- DerivativeF(theta, x)
  residuals <- residuals(theta, x, y)

  hessian <- 2*t(deriv.f)%*%deriv.f -
    2*tensor(A = second.deriv.f, B = residuals, alongA = 3, alongB = 1)

  return(hessian)
}

# Compute the Hessian matrix of the sum-of-squares loss function with
# reparameterization.
#
# @param retheta Parameters of a 4PL model among which the EC50 parameter is
# in the log 10 dose scale.
# @param x Doses.
# @param y Responses.
#
# @return Hessian matrix of the sum-of-squares loss function in terms of
# reparameterized parameters.
Hessian.dr4pl_log10 <- function(theta, x, y) {


  n <- length(x)  # Number of observations
  p <- length(theta)  # Number of parameters


  theta.1 <- theta[1]
  theta.2 <- theta[2]
  theta.3 <- theta[3]
  theta.4 <- theta[4]

  # Second order derivatives of f
  second.deriv.f <- array(data = 0, dim = c(p, p, n))

  eta <- 10^(theta.3*(log10(x)-theta.2))  # Term needed in the Hessian matrix computation

  deriv.eta.2 <- -theta.3*log(10)*eta
  deriv.eta.3 <- (log10(x)-theta.2)*log(10)*eta


  second.deriv.f[1, 1, ] <- 0
  second.deriv.f[1, 2, ] <- deriv.eta.2/(1+eta)^2
  second.deriv.f[1, 3, ] <- deriv.eta.3/(1+eta)^2
  second.deriv.f[1, 4, ] <- 0


  second.deriv.f[2, 2, ] <- ((theta.4 - theta.1) * theta.3 * log(10) * deriv.eta.2 /
                               ((1+eta)^2)) * (1 - (2 * eta/(1 + eta)))
  second.deriv.f[2, 3, ] <- ((theta.4 - theta.1) * log(10) * eta /
                               ((1 + eta)^2)) * (theta.3 * log(10) * (log10(x) - theta.2) + 1 -
                                                   (2 * theta.3 * deriv.eta.3 /(1 + eta)))
  second.deriv.f[2, 4, ] <- (theta.3 * log(10) * eta)/((1 + eta)^2)

  second.deriv.f[3, 3, ] <- - (((theta.4 - theta.1) * (log10(x) - theta.2) * log(10) * deriv.eta.3) /
                                 ((1 + eta)^2)) * (1 - (2 * eta / (1 + eta)))
  second.deriv.f[3, 4, ] <- - (log10(x) - theta.2) * log(10) * eta / ((1 + eta)^2)

  second.deriv.f[4, 4, ] <- 0

  second.deriv.f <- (second.deriv.f + aperm(second.deriv.f, c(2, 1, 3)))/2

  # Substitute limits for second derivatives when eta are infinite
  second.deriv.f[, , eta == 0] <- 0
  second.deriv.f[, , eta == Inf] <- 0

  deriv.f <- DerivativeF(theta, x)
  residuals <- residuals(theta, x, y)

  hessian <- 2*t(deriv.f)%*%deriv.f -
    2*tensor(A = second.deriv.f, B = residuals, alongA = 3, alongB = 1)


  return(hessian)
}
