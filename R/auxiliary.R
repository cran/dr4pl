

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


#' @name gof
#' @title Perform the goodness-of-fit (gof) test for a model.
#'
#' @description S3 method for a model object
#'
#' @param object A model object
#' @param ... dots for future extensions
#'
#' @return data.frame with goodness-of-fit results
#'
#' @export
gof <- function(object, ...) UseMethod('gof')

#' @name gof-dr4pl
#' @title Perform the goodness-of-fit (gof) test for the 4PL model.
#'
#' @description Perform the goodness-of-fit (gof) test for the 4PL model when there
#'   are at least two replicates for each dose level.
#'
#' @param object An object of the dr4pl class.
#' @param n.signif.digit Number of significant digits after the decimal point to be
#' printed. The default value is 4, but users can change the value on their own.
#' @param ... dots for future extensions
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
#' gof(obj.dr4pl)  # Print the goodness-of-fit test results
#'
#' @references \insertRef{Seber1989}{dr4pl}
#' @export
gof.dr4pl <- function(object, n.signif.digit = 4, ...) {

  x <- object$data$Dose  # Dose levels
  y <- object$data$Response  # Responses

  J.i <- table(x)  # Numbers of observations at all dose levels
  n <- length(unique(x))  # Number of dose levels
  N <- object$sample.size  # Total number of observations
  p <- 4  # Number of parameters of the 4PL model is 4

  # Check whether function arguments are appropriate
  if(n <= 4) {

    abort(
      glue("The number of dose levels should be larger than",
           "four to perform the goodness-of-fit test for the 4PL models.")
      )
  }

  levels.x.sorted <- sort(unique(x))

  x.sorted <- sort(x)
  indices.x.sorted <- sort(x, index.return = TRUE)$ix
  y.sorted <- y[indices.x.sorted]

  # Average response values for different dose levels
  y.bar <- tapply(X = y.sorted, INDEX = x.sorted, FUN = mean)
  # Fitted response values
  y.fitted <- MeanResponse(object$parameters, levels.x.sorted)

  # Numbers of observations per dose level
  n.obs.per.dose <- tapply(X = y.sorted, INDEX = x.sorted, FUN = length)

  if(any(n.obs.per.dose <= 1)) {

    abort(
      glue("There should be more than one observation for each",
           " dose level to performthe goodness-of-fit test.")
      )
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
  if(!inherits(object,"dr4pl")) {

    abort(
      glue('The object for which the IC values are obtained should',
           ' be of the class "dr4pl".'))
  }
  if(any(inhib.percent <= 0|inhib.percent >= 100)) {

    abort("Inhibited percentages should be between 0 and 100.")
  }

  theta <- object$parameters
  # Inhibited responses corresponding to inhibited percentages
  inhib.resp <- (inhib.percent*theta[1] + (100 - inhib.percent)*theta[4])/100
  IC.vec <- theta[2]*((theta[4] - inhib.resp)/(inhib.resp - theta[1]))^(1/theta[3])
  names(IC.vec) <- paste("InhibPercent:", inhib.percent, sep = "")

  return(IC.vec)
}


#' Print the dr4pl object to screen.
#'
#' @param x a dr4pl object to be printed
#' @param ... all normally printable arguments
#'
#' @examples
#' ryegrass.dr4pl <- dr4pl(Response ~ Dose,
#'                         data = sample_data_1)
#' print(ryegrass.dr4pl)
#'
#' obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_5)
#' print(obj.dr4pl)
#' @export
print.dr4pl <- function(x, ...) {

  cat("Call:\n")
  print(x$call)

  cat("\nCoefficients:\n")
  print(x$parameters)
}

#' Print the dr4pl object summary to screen.
#'
#' @param x a dr4pl object to be summarized
#' @param ... all normally printable arguments
#'
#' @examples
#' library(drc)  # Needed for the data set 'ryegras'
#' dr4pl.ryegrass <- dr4pl(rootl ~ conc, data = ryegrass)
#' print(summary(dr4pl.ryegrass))
#'
#' dr4pl.7 <- dr4pl(Response ~ Dose, data = sample_data_7)
#' print(summary(dr4pl.7))
#' @export
print.summary.dr4pl <- function(x, ...) {

  cat("Call:\n")
  print(x$call)
  cat("\n")

  printCoefmat(x$coefficients, P.values = FALSE, has.Pvalue = FALSE)
}



#' @title summary
#'
#' @description Print the summary of a dr4pl object.
#'
#' @name summary.dr4pl
#'
#' @param object a dr4pl object to be summarized
#' @param parm parameters of the dr4pl object. Usually made with [dr4pl_theta]
#' @param ... additional arguments to be passed to [calculate.dr4pl]
#'
#' @examples
#' obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_5)  # Fit a 4PL model to data
#' summary(obj.dr4pl)
#'
#' obj.dr4pl <- dr4pl(Response ~ Dose, data = sample_data_6)  # Fit a 4PL model to data
#' summary(obj.dr4pl)
#'
#' @export
summary.dr4pl <- function(object, parm = NULL, ...) {

  theta <- parm %theta% ParmToLog(coef(object))
  cal <- calculate(object, parm = theta, ...)
  TAB <- cbind(
    data.frame(Estimate = as.numeric(cal$est),
               StdErr = as.numeric(cal$std.err)),
    cal$ci.table)

  rownames(TAB) <- param_names(theta)
  res <- list(call = object$call,
              coefficients = TAB)

  class(res) <- "summary.dr4pl"
  return(res)
}

