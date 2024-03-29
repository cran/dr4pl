#' @title FindHillBounds
#' @description Compute the Hill bounds based on initial parameter estimates and data.
#' @name FindHillBounds
#'
#' @param x Vector of doses.
#' @param y Vector of responses.
#' @param theta Parameters of a 4PL model.
#' @param level Confidence level to be used in computing the Hill bounds.
#' @param use.Hessian Indicator of whether the Hessian matrix (TRUE) or the
#' gradient vector is used in confidence interval computation.
#'
#' @return Data frame whose first column represents the bounds on the IC50 parameter
#' in log 10 scale and second column represents the bounds on the slope parameter.
#'
#' @details This function computes the Hill bounds based on initial parameter
#' estimates and data. It basically computes the confidence intervals of the
#' true parameters based on the variance-covariance matrix of a given initial
#' parameter estimates. The half of a hessian matrix is used as a
#' variance-covariance matrix. If matrix inversion of the variance-covariance matrix
#' is infeasible, a variation of the method in Wang et al. (2010) is used. The
#' parameter \code{level} is only for simulation.
#'
#' @author Hyowon An, \email{ahwbest@gmail.com}.
#'
#' @seealso \code{\link{FindInitialParms}}, \code{\link{FindLogisticGrids}}.
#'
#' @references \insertRef{Higham2002}{dr4pl} \insertRef{Wang2010}{dr4pl}
#'
#' @export
FindHillBounds <- function(x, y, theta,
                           use.Hessian = FALSE,
                           level = 0.9999) {

  # Check whether function arguments are appropriate.
  if(length(x) != length(y)) {

    abort("The numbers of dose levels and responses should be the same.")
  }

  n <- length(x)  # Number of observations in data
  retheta <- ParmToLog(theta)

  ## Compute confidence intervals of the true parameters.
  residuals <- residuals(retheta, x, y)

  ## When the Hessian matrix is used.
  if(use.Hessian) {

    Hessian <- Hessian(retheta, x, y)

    # Obtain a positie definite approximation of the Hessian matrix
    C.hat <- Matrix::nearPD(Hessian)$mat/2
  ## When the gradient vector is used.
  } else {

    deriv <- DerivativeF(retheta, x)

    C.hat <- Matrix::nearPD(t(deriv)%*%deriv)$mat
  }

  ind.mat.inv <- TRUE  # TRUE if matrix inversion is successful, FALSE otherwise
  vcov.mat <- try(solve(C.hat), silent = TRUE)  # Inverse matrix

  if(inherits(vcov.mat, "try-error")) {

    # Cholesky decomposition
    vcov.Chol <- try(chol(C.hat, silent = TRUE))

    if(inherits(vcov.Chol, "try-error")) {

      # Cholesky decomposition with a diagonal matrix as a lower bound
      vcov.Chol <- try(chol(0.99*C.hat + 0.01*diag(nrow(C.hat))))

      if(inherits(vcov.Chol, "try-error")) {

        ind.mat.inv <- FALSE
      } else {

        vcov.mat <- chol2inv(vcov.Chol)
      }

    } else {

      vcov.mat <- chol2inv(vcov.Chol)
    }
  }

  ## If matrix inversion is successful, compute the Hill bounds based on the
  ## confidence intervals of the true parameters based on the initial estimates.
  if(ind.mat.inv) {

    RSS <- sqrt(sum(residuals^2)/(n - 4))  # Residual sum of squares

    q.t <- qt(level, df = n - 4)  # Quantile of the t-distribution
    std.err <- RSS*sqrt(diag(vcov.mat))  # Standard error
    ci <- cbind(retheta - q.t*std.err, retheta + q.t*std.err)  # Confidence intervals

    # Perform constrained optimization
    bounds.retheta.2 <- ci[2, ]
    bounds.retheta.3 <- ci[3, ]
  ## If the variance-covariance matrix is unobtainable, mimic the method of
  ## Wang et al. (2010).
  } else {

    levels.log10.x <- unique(log10(x))  # Levels of log10 doses
    levels.log10.x <- levels.log10.x[levels.log10.x != -Inf]
    bounds.retheta.2 <- c(min(levels.log10.x), max(levels.log10.x))

    range.theta.3 <- c(0.01*retheta[3], 100*retheta[3])
    bounds.retheta.3 <- c(min(range.theta.3), max(range.theta.3))
  }

  Hill.bounds <- data.frame(LogTheta2 = bounds.retheta.2,
                            Theta3 = bounds.retheta.3)
  class(Hill.bounds) <- c("dr4pl_hill", class(Hill.bounds))
  return(Hill.bounds)
}

#' @title FindLogisticGrids
#' @description Compute the grids on the upper and lower asymptote parameters for the logistic
#' method based on initial parameter estimates and data.
#'
#' @name FindLogisticGrids
#'
#' @param x Vector of doses.
#' @param y Vector of responses.
#' @param retheta.init Parameters of a 4PL model among which the EC50 parameter is
#' in the log 10 dose scale.
#' @param use.Hessian Indicator of whether the Hessian matrix (TRUE) or the
#' gradient vector is used in confidence interval computation.
#'
#' @return Data frame whose first column represents the grid on the upper asymptote
#' parameter and second column represents the grid o the lower asymptote.
#'
#' @details This function computes the grids on the upper and lower asymptote
#' parameters based on initial parameter
#' estimates and data. It basically computes the confidence intervals of the
#' true parameters based on the variance-covariance matrix of the given initial
#' parameter estimates. If matrix inversion of the variance-covariance matrix
#' is infeasible, a variation of the method in Wang et al. (2010) is used.
#'
#' @author Hyowon An
#'
#' @seealso \code{FindHillBounds}, \code{FindInitialParms}
#'
#' @references
#' \insertRef{Wang2010}{dr4pl}
#'
#' @export
FindLogisticGrids <- function(x, y, retheta.init, use.Hessian = FALSE) {

  n <- length(x)  # Number of observations in data

  ## Compute confidence intervals of the true parameters.
  residuals <- residuals(retheta.init, x, y)

  ## When the Hessian matrix is used.
  if(use.Hessian) {

    Hessian <- Hessian(retheta.init, x, y)

    # Obtain a positie definite approximation of the Hessian matrix
    C.hat <- Matrix::nearPD(Hessian)$mat/2
    ## When the gradient vector is used.
  } else {

    deriv <- DerivativeF(retheta.init, x)

    C.hat <- t(deriv)%*%deriv
  }

  ind.mat.inv <- TRUE  # TRUE if matrix inversion is successful, FALSE otherwise
  vcov.mat <- try(solve(C.hat), silent = TRUE)  # Inverse matrix

  if(inherits(vcov.mat, "try-error")) {

    # Cholesky decomposition
    vcov.Chol <- try(chol(C.hat, silent = TRUE))

    if(inherits(vcov.Chol, "try-error")) {

      # Cholesky decomposition with a diagonal matrix as a lower bound
      vcov.Chol <- try(chol(0.99*C.hat + 0.01*diag(nrow(C.hat))))

      if(inherits(vcov.Chol, "try-error")) {

        ind.mat.inv <- FALSE
      } else {

        vcov.mat <- chol2inv(vcov.Chol)
      }

    } else {

      vcov.mat <- chol2inv(vcov.Chol)
    }
  }

  ## If matrix inversion is successful, compute the Hill bounds based on the
  ## confidence intervals of the true parameters based on the initial estimates.
  if(ind.mat.inv) {

    RSS <- sqrt(sum(residuals^2)/(n - 4))  # Residual sum of squares
    std.err <- RSS*sqrt(diag(vcov.mat))  # Standard error

    # Sizes of the grids on the upper and lower asymptote parameters
    grid.size.theta.1 <- std.err[1]
    grid.size.theta.4 <- std.err[4]

    # Grids on the upper and lower asymptote parameters
    grid.theta.1 <- c(-2*grid.size.theta.1, -grid.size.theta.1, 0,
                      grid.size.theta.1, 2*grid.size.theta.1)
    grid.theta.4 <- c(-2*grid.size.theta.4, -grid.size.theta.4, 0,
                      grid.size.theta.4, 2*grid.size.theta.4)

    grid.theta.1 <- retheta.init[1] + grid.theta.1
    grid.theta.4 <- retheta.init[4] + grid.theta.4

  ## If the variance-covariance matrix is unobtainable, use the 5% percentile
  ## of responses.
  } else {

    y.range <- diff(range(y))

    # Size of the grids on the upper and lower asymptote parameters
    grid.size <- 0.05*y.range

    grid.theta.1.4 <- c(-2*grid.size, -grid.size, 0, grid.size, 2*grid.size)
    grid.theta.1 <- retheta.init[1] + grid.theta.1.4
    grid.theta.4 <- retheta.init[4] + grid.theta.1.4
  }

  logistic.grids <- data.frame(Theta1 = grid.theta.1, Theta4 = grid.theta.4)

  return(logistic.grids)
}

#' @title FindInitialParms
#' @description Find initial parameter estimates for a 4PL model.
#'
#' @name FindInitialParms
#'
#' @param x Vector of dose levels
#' @param y Vector of responses
#' @param trend Indicator of whether the curve is a decreasing \eqn{\theta[3]<0}
#' or increasing curve \eqn{\theta[3]>0}. See \code{\link{dr4pl}} for detailed
#' explanation.
#' @param method.init Method of obtaining initial values of the parameters. See
#' \code{\link{dr4pl}} for detailed explanation.
#' @param method.robust Parameter to select loss function for the robust estimation
#' method to be used to fit a model. See \code{\link{dr4pl}} for detailed
#' explanation.
#'
#' @return Initial parameter estimates of a 4PL model in the order of the
#' upper asymptote, IC50, Slope and lower asymptote parameters.
#'
#' @export
FindInitialParms <- function(x, y,
                             trend = "auto",
                             method.init = "Mead",
                             method.robust = "squared") {

  n.grid.cells <- 25

  ### Types of available function arguments
  types.trend <- c("auto", "decreasing", "increasing")
  types.method.init <- c("logistic", "Mead")
  types.method.robust <- c("squared","absolute", "Huber", "Tukey")

  ### Check whether function arguments are appropriate.
  if(length(x) == 0 || length(y) == 0 || length(x) != length(y)) {

    abort("The same numbers of dose and response values should be supplied.")
  }
  if(!is.element(trend, types.trend)) {

    abort('The type of a dose-response curve should be one of "auto", "decreasing" and "increasing".')
  }
  if(!is.element(method.init, types.method.init)) {

    abort('The initialization method name should be one of "logistic" and "Mead".')
  }
  if(!(is.null(method.robust)||is.element(method.robust, types.method.robust))) { #Should add depricated warning for NULL argument?

    abort('The robust estimation method should be one of "squared", "absolute", "Huber" or "Tukey".')
  }
  ### Set the upper and lower asymptotes
  y.range <- diff(range(y))

  theta.1.4.zero <- 0.001*y.range  # This value will be added to the maximum and minimum

  y.max <- max(y) + theta.1.4.zero  # Maximum of responses
  y.min <- min(y) - theta.1.4.zero  # Minimum of responses

  ### Standard logistic method
  if(method.init == "logistic") {
    
    ## Obtain initial parameter estimates using the logistic method
    theta.Hill <- Logistic(x, y, y.max, y.min)

    retheta.Hill <- ParmToLog(theta.Hill)

    ## Obtain the logistic grids on the upper and lower asymptote parameters
    logistic.grids <- FindLogisticGrids(x, y, retheta.Hill)
    grid.theta.1 <- logistic.grids$Theta1
    grid.theta.4 <- logistic.grids$Theta4

    ## Matrix of initial parameter estimates
    theta.mat <- matrix(NA, nrow = nrow(logistic.grids)^2, ncol = 4)
    i.row <- 1

    for(theta.1.init in grid.theta.1) {

      for(theta.4.init in grid.theta.4) {

        # The upper asymptote should always be larger than the lower asymptote.
        if(theta.1.init<=theta.4.init) {

          next
        }

        theta.logistic <- Logistic(x, y, theta.1.init, theta.4.init)

        if(!is.null(theta.logistic)) {

          theta.mat[i.row, ] <- theta.logistic
        }

        i.row <- i.row + 1
      }
    }

    # This restricts approximation to decline curves only
    theta.mat <- theta.mat[!is.na(theta.mat[, 3]), ]
    if(trend == "decreasing") {

      theta.mat <- theta.mat[theta.mat[, 3]<0, ]
    } else if(trend == "increasing") {

      theta.mat <- theta.mat[theta.mat[, 3]>0, ]
    }

    if(nrow(theta.mat) == 0) {

      abort("The logistic method is not applicable for the input data. Please try Mead's method instead.")
    }

    loss.fcn <- ErrFcn(method.robust)
    losses <- rep(0, nrow(theta.mat))   # To save the error values

    for(i in 1:nrow(theta.mat)) {

      theta <- theta.mat[i, ]
      losses[i] <- loss.fcn(theta, x, y)
    }

    theta.init <- theta.mat[which.min(losses), ]
  ### Mead's method
  } else if(method.init == "Mead") {

    log.x <- log10(x)

    # Subtract the minimum response from responses to obtain positive values
    y.positive <- y - y.min

    # The grid on the slope parameter varies according to the trend option
    if(trend == "auto") {

      mu.3.vec <- 10^qcauchy(seq(from = 0, to = 1, length.out = n.grid.cells + 2))
    } else if(trend == "decreasing") {

      mu.3.vec <- 10^qcauchy(seq(from = 0.5, to = 1, length.out = n.grid.cells + 2))
    } else {

      mu.3.vec <- 10^qcauchy(seq(from = 0, to = 0.5, length.out = n.grid.cells + 2))
    }
    mu.3.vec <- mu.3.vec[-c(1, n.grid.cells + 2)]
    mu.3.vec <- mu.3.vec[mu.3.vec != 1]

    theta.mat <- matrix(0, nrow = length(mu.3.vec), ncol = 4)
 
    for(i in 1:length(mu.3.vec)) {

      mu.3 <- mu.3.vec[i]  # Current mu.3 value

      # Apply Mead's method
      data.Mead <- data.frame(Y = y.positive,
                              YMuX = y.positive*mu.3^log.x,
                              Response = 1)

      if(mu.3 < 1) {

        data.Mead <- data.Mead[data.Mead$YMuX != Inf, ]
      }

      lm.init <- stats::lm(Response ~ -1 + Y + YMuX, data = data.Mead)

      lm.init.coef <- lm.init$coefficients
      mu.1 <- lm.init.coef[1]
      mu.2 <- lm.init.coef[2]

      # The signs of mu.1 and mu.2 should be the same.
      if(sign(mu.1) != sign(mu.2)) {

        next
      }

      theta.mat[i, 1] <- 1/mu.1
      theta.mat[i, 2] <- (mu.1/mu.2)^(1/log10(mu.3))
      theta.mat[i, 3] <- -log10(mu.3)
    }

    theta.mat[, 1] <- theta.mat[, 1] + y.min
    theta.mat[, 4] <- theta.mat[, 4] + y.min

    theta.mat <- theta.mat[theta.mat[, 3] != 0, ]

    if(nrow(theta.mat) == 0) {

      abort("Mead's method is not applicable for the inupt data. Please try the logistic method instead.")
    }

    loss.fcn <- ErrFcn(method.robust)
    losses <- rep(0, nrow(theta.mat))   # To save the loss values

    for(i in 1:nrow(theta.mat)) {

      theta <- theta.mat[i, ]
      losses[i] <- loss.fcn(theta, x, y)
    }

    theta.init <- theta.mat[which.min(losses), ]
  }

  theta.init <- new_dr4pl_theta(theta.init)

  return(theta.init)
}

# Find initial parameter estimates by the logistic method.
#
# @param x Vector of dose levels
# @param y Vector of responses
# @param theta.1 Estimate of the upper asymptote parameter
# @param theta.4 Estimate of the lower asymptote parameter
#
# @return Parameter estimates of a 4PL model in the order from the upper asymptote,
# IC50, slope and lower asymptote.
Logistic <- function(x, y, theta.1, theta.4) {

  data.Hill <- data.frame(x = x, y = y)
  data.Hill <- subset(data.Hill, subset = (data.Hill$y<theta.1)&
                                          (data.Hill$y>theta.4)&
                                          (data.Hill$x>0))

  # There are not enough data to estimate paramter estimates
  if(nrow(data.Hill) <= 4) {

    return(NULL)
  }

  # Logit transformed responses
  data.Hill$LogitY <- log((data.Hill$y - theta.4)/(theta.1 - data.Hill$y))
  data.Hill$LogX <- log(data.Hill$x)  # Log transformed doses

  lm.Hill <- lm(LogitY ~ LogX, data = data.Hill)

  coef.Hill <- lm.Hill$coefficients
  theta.3 <- coef.Hill[2]
  theta.2 <- exp(-coef.Hill[1]/coef.Hill[2])
  theta <- c(theta.1,theta.2,theta.3,theta.4)
  if(anyNA(theta)|| (!is.na(theta[2])&&theta[2]<=0)) {
    return(NULL)
  }
  theta.init <- new_dr4pl_theta(theta)

  return(theta.init)
}
