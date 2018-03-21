
#' @title Fit a 4 parameter logistic (4PL) model to dose-response data.
#' 
#' @description Compute the approximate confidence intervals of the parameters of a 
#' 4PL model based on the asymptotic normality of least squares estimators.
#'   
#' @name confint.dr4pl
#'   
#' @param object An object of the dr4pl class
#' @param parm Parameters of a 4PL model
#' @param level Confidence level
#' @param ... Other parameters to be passed
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
  
  x <- object$data$Dose
  y <- object$data$Response
  
  # If parameter estimates are not provided by a user, then proceed with the
  # estimates stored in the input dr4pl object.
  if(is.null(parm)) {
    
    retheta <- ParmToLog(object$parameters)
  } else if(length(parm)!=4||parm[2]<=0) {
    
    stop("The input argument \"parm\" should be of length 4 and the IC50 estimate
         should be nonnegative.")
  } else {
    
    retheta <- ParmToLog(parm)
  }
  
  C.hat <- HessianLogIC50(retheta, x, y)/2  # Estimated variance-covariance matrix
  
  n <- object$sample.size  # Number of observations in data
  f <- MeanResponseLogIC50(x, retheta)
  
  ind.mat.inv <- TRUE  # TRUE if matrix inversion is successful, FALSE otherwise
  
  vcov.mat <- vcov(object, parm)  # Variance-covariance matrix
  
  if(!ind.mat.inv) {
    
    print("The hessian matrix is singular, so an approximated positive definite
          hessian matrix is used.")
  }
  
  RSS <- sqrt(sum((y - f)^2)/(n - 4))
  
  q.t <- qt(1 - (1 - level)/2, df = n - 4)
  std.err <- RSS*sqrt(diag(vcov.mat))  # Standard error
  ci.table <- cbind(retheta - q.t*std.err, retheta + q.t*std.err)
  ci.table[2, ] <- 10^ci.table[2, ]
  
  colnames(ci.table) <- c(paste(100*(1 - level)/2, "%"),
                          paste(100*(1 - (1 - level)/2), "%"))
  if(retheta[3]<=0) {
    
    row.names(ci.table) <- c("UpperLimit", "IC50", "Slope", "LowerLimit")
  } else {
    
    row.names(ci.table) <- c("UpperLimit", "EC50", "Slope", "LowerLimit")
  }
  
  return(ci.table)
}

#' @title Obtain coefficients of a 4PL model
#' 
#' @description This function obtains the coefficients of a 4PL model. Estimates
#' of the four parameters, the upper asymptote, IC50, slope and lower asymptote,
#' are returned.
#' 
#' @name coef.dr4pl
#' 
#' @param object A 'dr4pl' object
#' @param ... arguments passed to coef
#' 
#' @return A vector of parameters
#' 
#' @examples
#' obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_2)  # Fit a 4PL model to data
#' coef(obj.dr4pl)  # Print parameter estimates
#' 
#' obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_3)  # Fit a 4PL model to data
#' coef(obj.dr4pl)  # Print parameter estimates
#' 
#' @export
coef.dr4pl <- function(object, ...) {
  
  object$parameters
}

#' @title Perform the goodness-of-fit (gof) test for the 4PL model.
#' 
#' @description Perform the goodness-of-fit (gof) test for the 4PL model when there
#'   are at least two replicates for each dose level.
#'   
#' @name gof.dr4pl
#'   
#' @param object An object of the dr4pl class.
#' @param n.signif.digit Number of significant digits after the decimal point to be 
#' printed. The default value is 4, but users can change the value on their own.
#' 
#' @return A list of results in the order of a F-statistic value, p-value and a
#' degree of freedom.
#'   
#' @details Perform a goodness-of-fit (gof) test for the goodness of the 4PL models
#' for dose-response data. There should be at least two replicates at each dose
#' level. The test statistic follows the F distribution with degress of freedom
#' (n - 4) and (N - n) where N stands for the total number of observations and n
#' stands for the number of dose levels. For detailed explanation of the method,
#' please refer to Subsection 2.1.5 of Seber and Wild (1989).
#'
#' @examples
#' obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_4)  # Fit a 4PL model to data
#' gof.dr4pl(obj.dr4pl)  # Print the goodness-of-fit test results
#'
#' @references \insertRef{Seber1989}{dr4pl}
#'
#' @export
gof.dr4pl <- function(object, n.signif.digit = 4) {
  
  x <- object$data$Dose  # Dose levels
  y <- object$data$Response  # Responses
  
  J.i <- table(x)  # Numbers of observations at all dose levels
  n <- length(unique(x))  # Number of dose levels
  N <- object$sample.size  # Total number of observations
  p <- 4  # Number of parameters of the 4PL model is 4
  
  # Check whether function arguments are appropriate
  if(n <= 4) {
    
    stop("The number of dose levels should be larger than four to perform the
         goodness-of-fit test for the 4PL models.")
  }

  levels.x.sorted <- sort(unique(x))
  
  x.sorted <- sort(x)
  indices.x.sorted <- sort(x, index.return = TRUE)$ix
  y.sorted <- y[indices.x.sorted]
  
  # Average response values for different dose levels
  y.bar <- tapply(X = y.sorted, INDEX = x.sorted, FUN = mean)
  # Fitted response values
  y.fitted <- MeanResponse(levels.x.sorted, object$parameters)
  
  # Numbers of observations per dose level
  n.obs.per.dose <- tapply(X = y.sorted, INDEX = x.sorted, FUN = length)
  
  if(any(n.obs.per.dose <= 1)) {
    
    stop("There should be more than one observation for each dose level to perform
         the goodness-of-fit test.")
  }
  
  y.bar.rep <- rep(y.bar, times = n.obs.per.dose)
  
  ### Compute statistics values for the goodness-of-fit test.
  ss.lof <- J.i%*%(y.bar - y.fitted)^2  # Lack-of-fit sum of squares.
  ss.error <- sum((y.sorted - y.bar.rep)^2)  # Pure-error sum of squares.
  ss.vec <- c(ss.lof, ss.error)  # Vector of sums of squares.
  
  df.gof <- c(n - p, N - n)  #  Degrees of freedom.
  
  ms.vec <- ss.vec/df.gof
  
  gof.stat <- ms.vec[1]/ms.vec[2]
  gof.pval <- pf(gof.stat, df1 = n - p, df2 = N - n, lower.tail = FALSE)  # p-value

  result.tab <- cbind(round(df.gof, digits = n.signif.digit),
                      round(ss.vec, digits = n.signif.digit),
                      round(ms.vec, digits = n.signif.digit),
                      c(round(gof.stat, digits = n.signif.digit), ""),
                      c(round(gof.pval, digits = n.signif.digit), ""))
  row.names(result.tab) <- c("Lack-of-fit", "Pure-error")
  colnames(result.tab) <- c("d.f.", "Sum Sq", "Mean Sq", "F value", "Pr(>F)")
  
  cat("Goodnes-of-fit test\n\n")
  print(noquote(result.tab))
}

#' @title Obtain Inhibitory Concentrations (IC) of a dose-response curve
#' 
#' @description This function obtains estimates of the IC's of a dose-response 
#' curve. Typically the IC50 parameter is of interest, but sometimes IC10 or IC90
#' are important aspects of a dose-response curve. By controlling the function
#' argument, a user can obtain the IC's at various levels.
#' 
#' @name IC
#' 
#' @param object Object of the class `dr4pl` for which the IC values are obtained
#' @param inhib.percent Inhibited percentages at which thc IC values are obtained
#' @examples 
#' data.test <- data.frame(x = c(0.0001, 0.001, 0.01, 0.1, 1),
#'                         y = c(10, 9, 5, 1, 0))
#' obj.dr4pl <- dr4pl(y ~ x,
#'                    data = data.test)
#' IC(obj.dr4pl, inhib.percent = c(10, 90))
#' 
#' obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_4)  # Fit a 4PL model to data
#' IC(obj.dr4pl, inhib.percent = c(10, 50, 90))
#' 
#' @return IC values at the inhibited percentages provided by the argument 
#' \code{inhib.percent}
#' @export
IC <- function(object, inhib.percent) {
  
  ### Check whether function arguments are appropriate
  if(class(object) != "dr4pl") {
    
    stop("The object for which the IC values are obtained should be of the class
         \"dr4pl\".")
  }
  if(any(inhib.percent <= 0|inhib.percent >= 100)) {
    
    stop("Inhibited percentages should be between 0 and 100.")
  }
  
  theta <- object$parameters
  # Inhibited responses corresponding to inhibited percentages
  inhib.resp <- (inhib.percent*theta[1] + (100 - inhib.percent)*theta[4])/100
  IC.vec <- theta[2]*((theta[4] - inhib.resp)/(inhib.resp - theta[1]))^(1/theta[3])
  names(IC.vec) <- paste("InhibPercent:", inhib.percent, sep = "")
  
  return(IC.vec)
}

#' @title Make a plot of a 4PL model curve and data
#' 
#' @description This function displays a dose-response curve and data. As a default,
#' the x-axis represents dose levels in log 10 scale and the y-axis represents
#' responses. The black solid line represents a dose-response curve. The blue filled
#' circles represent data points and red triangles represent outliers.
#'   
#' @name plot.dr4pl
#' 
#' @param x `dr4pl' object whose data and mean response function will be plotted.
#' @param type.curve Indicator of the type of a dose-response curve. "all" indicates
#' that data and a curve will be plotted while "data" indicates that only data
#' will be plotted.
#' @param text.title Character string for the title of a plot with a default set to 
#'   "Dose response plot".
#' @param text.x Character string for the x-axis of the plot with a default set to 
#'   "Dose".
#' @param text.y Character string for the y-axis of the plot with a default set to 
#'   "Response".
#' @param indices.outlier Pass a vector indicating all indices which are outliers in 
#'   the data.
#' @param breaks.x Vector of desired break points for the x-axis
#' @param breaks.y Vector of desired break points for the y-axis
#' @param ... All arguments that can normally be passed to ggplot.
#' 
#' @examples
#' dr4pl.1 <- dr4pl(Response ~ Dose, data = sample_data_1)
#'
#' plot(dr4pl.1)
#' 
#' ## Able to further edit plots.
#' library(ggplot2) #needed to change color to green
#' dr4pl.1 <- dr4pl(Response ~ Dose,
#'                         data = sample_data_1,
#'                         text.title = "Sample Data Plot")
#'
#' a <- plot(dr4pl.1)
#' a + geom_point(color = "green", size = 5)
#' 
#' ## Bring attention to outliers using parameter indices.outlier.
#' dr4pl.3 <- dr4pl(Response ~ Dose,
#'                  data = drc_error_3,
#'                  method.init = "Mead",
#'                  method.robust = "absolute")
#' plot(dr4pl.3, indices.outlier = c(90, 101))
#' 
#' ## Change the plot title default with parameter text.title.
#' dr4pl.1 <- dr4pl::dr4pl(Response ~ Dose,
#'                         data = sample_data_1)
#' plot(dr4pl.1, text.title = "My New Dose Response plot")
#' 
#' ##Change the labels of the x and y axis to your need
#' library(drc)  # Needed to load 'decontaminants' data set
#' data.hpc <- subset(decontaminants, group %in% "hpc")
#' dr4pl.hpc <- dr4pl(count~conc, data = data.hpc)
#' plot(dr4pl.hpc,
#'      text.title = "hpc Decontaminants Plot",
#'      text.x = "Concentration",
#'      text.y = "Count")
#' 
#' @author Hyowon An, \email{ahwbest@gmail.com}
#' @author Justin T. Landis, \email{jtlandis314@gmail.com}
#' @author Aubrey G. Bailey, \email{aubreybailey@gmail.com}
#' 
#' @export
plot.dr4pl <- function(x,
                       type.curve = "all",
                       text.title = "Dose-response plot",
                       text.x = "Dose",
                       text.y = "Response",
                       indices.outlier = NULL,
                       breaks.x = NULL,
                       breaks.y = NULL,
                       ...) {
  
  ### Check whether function arguments are appropriate
  if(!is.character(text.title)) {
    
    stop("Title text should be characters.")
  }
  if(!is.character(text.x)) {
    
    stop("The x-axis label text should be characters.")
  }
  if(!is.character(text.y)) {
    
    stop("The y-axis label text should be characters.")
  }
  
  ### Draw a plot
  n <- x$sample.size
  color.vec <- rep("blue", n)
  shape.vec <- rep(19, n)
  
  if(!is.null(indices.outlier)) {
    
    color.vec[indices.outlier] <- "red"
    shape.vec[indices.outlier] <- 17
  }
  
  a <- ggplot2::ggplot(aes(x = x$data$Dose, y = x$data$Response), data = x$data)
  
  if(type.curve == "all") {
    
    a <- a + ggplot2::stat_function(fun = MeanResponse,
                                    args = list(theta = x$parameters),
                                    size = 1.2)
  }
  
  a <- a + ggplot2::geom_point(size = I(5), alpha = I(0.8), color = color.vec,
                               shape = shape.vec)
  
  a <- a + ggplot2::labs(title = text.title,
                         x = text.x,
                         y = text.y)
  
  # Set parameters for the grids
  a <- a + ggplot2::theme(strip.text.x = ggplot2::element_text(size = 16))
  a <- a + ggplot2::theme(panel.grid.minor = ggplot2::element_blank())
  a <- a + ggplot2::theme(panel.grid.major = ggplot2::element_blank())
  
  if(!is.null(breaks.x)) { 
    
    a <- a + ggplot2::scale_x_log10(breaks = breaks.x)
  } else { 
    
    a <- a + ggplot2::scale_x_log10()
  }
  if(!is.null(breaks.y)) {
    
    a <- a + ggplot2::scale_y_continuous(breaks = breaks.y)
  } else { 
    
    a <- a + ggplot2:: scale_y_continuous()
  }
  
  a <- a + ggplot2::theme_bw()
  # Test
  # Set parameters for the titles and text / margin(top, right, bottom, left)
  a <- a + ggplot2::theme(plot.title = ggplot2::element_text(size = 20, margin = ggplot2::margin(0, 0, 10, 0)))
  a <- a + ggplot2::theme(axis.title.x = ggplot2::element_text(size = 16, margin = ggplot2::margin(15, 0, 0, 0)))
  a <- a + ggplot2::theme(axis.title.y = ggplot2::element_text(size = 16, margin = ggplot2::margin(0, 15, 0, 0)))
  a <- a + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 16))
  a <- a + ggplot2::theme(axis.text.y = ggplot2::element_text(size = 16))
  
  return(a)
}

#' Print the dr4pl object to screen.
#'
#' @param object a dr4pl object to be printed
#' @param ... all normally printable arguments
#' 
#' @examples
#' ryegrass.dr4pl <- dr4pl(Response ~ Dose,
#'                         data = sample_data_1)
#' print(ryegrass.dr4pl)
#' 
#' obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_5)
#' print(obj.dr4pl)
print.dr4pl <- function(object, ...) {
  
  if(!inherits(object, "dr4pl")) {
    
    stop("The object to be printed should be of the class \'dr4pl\'.")
  }

  cat("Call:\n")
  print(object$call)
  
  cat("\nCoefficients:\n")
  print(object$parameters)
}

#' Print the dr4pl object summary to screen.
#' 
#' @param object a dr4pl object to be summarized
#' @param ... all normally printable arguments
#' 
#' @examples
#' library(drc)  # Needed for the data set 'ryegras'
#' dr4pl.ryegrass <- dr4pl(rootl ~ conc, data = ryegrass)
#' print(summary(dr4pl.ryegrass))
#' 
#' dr4pl.7 <- dr4pl(Response ~ Dose, data = sample_data_7)
#' print(summary(dr4pl.7))
print.summary.dr4pl <- function(object, ...) {
  
  cat("Call:\n")
  print(object$call)
  cat("\n")
  
  printCoefmat(object$coefficients, P.values = TRUE, has.Pvalue = TRUE)
}


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
#' @param ... Other function arguments to be passed to the default 'vcov' function.
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
#' @export
vcov.dr4pl <- function(object, ...) {
  
  x <- object$data$Dose  # Vector of dose levels
  y <- object$data$Response  # Vector of responses
  
  retheta <- ParmToLog(object$parameters)

  C.hat <- HessianLogIC50(retheta, x, y)/2  # Estimated variance-covariance matrix
  
  # If the Hessian matrix is not positive definite, use
  if(!matrixcalc::is.positive.semi.definite(C.hat)) {
    
    Jacobian <- DerivativeFLogIC50(retheta, x)
    C.hat <- t(Jacobian)%*%Jacobian
  }
  
  C.hat <- Matrix::nearPD(C.hat)$mat
  
  n <- object$sample.size  # Number of observations in data
  f <- MeanResponseLogIC50(x, retheta)
  
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
  
  colnames(vcov.mat) <- c("UpperLimit", "IC50", "Slope", "LowerLimit")
  if(retheta[3]<=0) {
    
    row.names(vcov.mat) <- c("UpperLimit", "IC50", "Slope", "LowerLimit")
  } else {
    
    row.names(vcov.mat) <- c("UpperLimit", "EC50", "Slope", "LowerLimit")
  }
  
  return(vcov.mat)
}


#' @title summary
#'
#' @description Print the dr4pl object summary.
#' 
#' @name summary.dr4pl
#' 
#' @param object a dr4pl object to be summarized
#' @param ... all normal summary arguments
#' 
#' @examples
#' obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_5)  # Fit a 4PL model to data
#' summary(obj.dr4pl)
#' 
#' obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_6)  # Fit a 4PL model to data
#' summary(obj.dr4pl)
#' 
#' @export
summary.dr4pl <- function(object, ...) {
  
  std.err <- sqrt(diag(vcov(object)))
  t.stats <- coef(object)/std.err
  ci <- confint(object)
  
  TAB <- data.frame(Estimate = object$parameters,
               StdErr = std.err,
               t.value = t.stats,
               p.value = 2*pt(-abs(t.stats), df = object$sample.size - 4))
  
  res <- list(call = object$call,
              coefficients = TAB)
  
  class(res) <- "summary.dr4pl"
  return(res)
}