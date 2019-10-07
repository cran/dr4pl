
#' @name dr4pl
#' 
#' @docType package
#' 
#' @import graphics 
#' @import stats
#' @import ggplot2
#' @import tensor
#' @importFrom Matrix nearPD
#' @importFrom Rdpack reprompt
#' @importFrom matrixcalc is.positive.definite
NULL
#' @title Fitting 4 Parameter Logistic (4PL) models to dose-response data.
#' 
#' @description This function fits a 4PL model to dose-response data. Users can
#' obtain fitted parameter estimates as return values. Using auxiliary functions
#' provided by this R package, users can plot a fitted dose-response curve and
#' obtain confidence intervals of true parameters. In addition, the goodness-of-fit
#' test for model adequacy of the 4PL models can be performed when replicates are
#' available for each dose level.
#' 
#' @export
dr4pl <- function(...)  UseMethod("dr4pl")

#' @describeIn dr4pl General 4PL model fitting function for analysis of
#'   dose-response relation.
#'
#' @param  formula Symbolic description of the model to be fit. Either of the
#' form 'response ~ dose' or as a data frame with response values in first
#' column and dose values in second column.
#' @param data Data frame containing variables in the model.
#' @param init.parm Either Null or a Vector of initial parameters to be optimized in the model.
#' \itemize{
#'   \item UpperLimit:  \eqn{\theta[1]}
#'   \item IC50/EC50: \eqn{\theta[2]}
#'   \item Slope: \eqn{\theta[3]}
#'   \item LowerLimit: \eqn{\theta[4]}
#' }
#' dr4pl assumes \eqn{\theta[1]>\theta[4]}. Note that when using upperl and
#' lowerl, the user may need to set this parameter because the estimated
#' parameters may not be within the feasible region. 
#' @param trend Indicator of whether a dose-response curve is a decreasing 
#' \eqn{\theta[3]<0} or increasing curve \eqn{\theta[3]>0}. The default is "auto" 
#' which indicates that the trend of the curve is automatically determined by
#' data. The option "decreasing" will impose a restriction \eqn{\theta[3]<=0} 
#' while the option "increasing" will impose a restriction \eqn{\theta[3]>=0} in an 
#' optimization process.
#' @param method.init Method of obtaining initial values of the parameters.
#' If this parameter is left unassigned, a default "Mead" method will be used.
#' Assign "logistic" to use the logistic method.
#' @param method.optim Method of optimization of the loss function specified by
#' \code{method.robust}. This function argument is directly passed to the function
#' \code{\link[stats]{constrOptim}} which is provided in the \pkg{base} package of R.
#' @param method.robust Parameter to select loss function for the robust estimation 
#' method to be used to fit a model. The argument NULL indicates the sum of squares
#' loss, "absolute" indicates the absolute deviation loss, "Huber" indicates Huber's
#' loss and "Tukey" indicates Tukey's biweight loss.
#' @param use.Hessian Indicator of whether the Hessian matrix (TRUE) or the
#' gradient vector is used in the Hill bounds.
#' @param level Confidence level to be used in Hill bounds computation.
#' @param failure.message Indicator of whether a message indicating attainment of
#' the Hill bounds and possible resolutions will be printed to the console (TRUE)
#' or hidden (FALSE).
#' @param upperl Either NULL or a numeric vector of length 4 that specifies the upper limit 
#' for the initial parameters of \eqn{c(\theta[1],\theta[2],\theta[3],\theta[4])} during the
#' optimization process. By default no upperl is assumed. If the user wants to constain only
#' some parameter values, set desired numeric bound in the appropriate position and fill Inf
#' to impose no upper bounds on other values. All upperl values must be greater than 
#' corresponding initalized parameters.
#' @param lowerl Either NULL or a numeric vector of length 4 that specifies the lower limit 
#' for the initial parameters of \eqn{c(\theta[1],\theta[2],\theta[3],\theta[4])} during the
#' optimization process. By default no lowerl is assumed. If the user wants to constain only
#' some parameter values, set desired numeric bound in the appropriate position and fill -Inf
#' to impose no lower bounds on other values. All lowerl values must be greater than 
#' corresponding initalized parameters.
#' @param ... Further arguments to be passed to \code{constrOptim}.
#' 
#' @return A 'dr4pl' object for which "confint", "gof", "print" and "summary"
#'   methods are implemented. For details, see the help page of each method.
#'   For example, type \code{?confint.dr4pl} to obtain the confidence intervals
#'   of parameters of the 'dr4pl' object.
#'   
#' @details This function fits a 4 parameter logistic (4PL) model to dose-response
#'   data. A formula of the model is
#'   \deqn{\theta[1]+(\theta[4]-\theta[1])/(1+(z/\theta[2])^\theta[3])}
#'
#'   \code{method.init} specifies an initialization method to get initial parameter
#'   estimates based on data. The currently supported initialization methods are
#'   "logistic" and 'Mead'. For further details, see the vignette.
#'
#'   \code{method.optim} specifies an optimization method to be used in
#'   "constrOptim" function. The currently supported optimization techniques
#'   include "Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN" and "Brent". For
#'   further details, see the help page of \code{\link[stats]{optim}}.
#'
#'   \code{method.robust} chooses a robust estimation method among 4 methods.
#'   The method of estimation is usually identified by the loss function of the
#'   method. This package supports 4 types of loss functions: sum-of-squares loss,
#'   absolute deviation loss, Huber's loss and Tukey's biweight loss. Each of
#'   loss function is explained in detail in the vignette.
#'   
#' @author Hyowon An, \email{ahwbest@gmail.com}
#' @author Justin T. Landis, \email{jtlandis314@gmail.com}
#' @author Aubrey G. Bailey, \email{aubreybailey@gmail.com}
#' 
#' @seealso \code{\link{confint.dr4pl}}, \code{\link{gof.dr4pl}},
#' \code{\link{print.dr4pl}}, \code{\link{summary.dr4pl}}
#' 
#' @export
dr4pl.formula <- function(formula,
                          data = list(),
                          init.parm = NULL,
                          trend = "auto",
                          method.init = "Mead",
                          method.robust = "squared",
                          method.optim = "Nelder-Mead",
                          use.Hessian = FALSE,
                          level = 0.9999,
                          failure.message = FALSE,
                          upperl = NULL,
                          lowerl = NULL,
                          ...) {
  
  mf <- model.frame(formula = formula, data = data)  # Model frame
  mm <- model.matrix(attr(mf, "terms"), data = mf)  # Model matrix
  
  # Check whether only one variable for doses was given in the formula.
  # Currently only one variable for doses is allowed.
  if(ncol(mm)>2) {
    
    stop("Only one indepedent variable should be specified for doses.")
  }
  
  # Dose-response models do not have intercepts.
  dose <- mm[, -1]
  response <- model.response(mf)
  
  obj <- dr4pl.default(dose = dose,
                       response = response,
                       init.parm = init.parm,
                       trend = trend,
                       method.init = method.init,
                       method.robust = method.robust,
                       method.optim = method.optim,
                       use.Hessian = use.Hessian,
                       level = level,
                       failure.message = failure.message,
                       upperl = upperl,
                       lowerl = lowerl,
                       ...)
  
  obj$call <- match.call()
  obj$formula <- formula

  return(obj)
}

#' @describeIn dr4pl Method for when formula argument is missing.
#'  dose and response arguments are necessary
#'   
#' @export
dr4pl.data.frame <- function(data,
                             dose,
                          response,
                          init.parm = NULL,
                          trend = "auto",
                          method.init = "Mead",
                          method.robust = "squared",
                          method.optim = "Nelder-Mead",
                          use.Hessian = FALSE,
                          level = 0.9999,
                          failure.message = FALSE,
                          upperl = NULL,
                          lowerl = NULL,
                          ...) {
  dose <- eval(substitute(dose),data)
  response <- eval(substitute(response),data)
  obj <- dr4pl.default(dose = dose,
                       response = response,
                       init.parm = init.parm,
                       trend = trend,
                       method.init = method.init,
                       method.robust = method.robust,
                       method.optim = method.optim,
                       use.Hessian = use.Hessian,
                       level = level,
                       failure.message = failure.message,
                       upperl = upperl,
                       lowerl = lowerl,
                       ...)
  
  obj$call <- match.call()
  
  return(obj)
}


#' @describeIn dr4pl Used in the default case, supplying a single dose and 
#'   response variable
#'   
#' @param dose Vector of dose levels
#' @param response Vector of responses
#'
#' @examples 
#'   ##Assign method.init = "logistic" to use logistic method of estimation.
#'   ##default method
#'   a <- dr4pl(dose = sample_data_1$Dose,
#'              response = sample_data_1$Response,
#'              method.init = "logistic")
#'   plot(a)
#'
#'   ##Use default or Assign method.init = "Mead" to use Mead's method of estimation.
#'   # Use method.robust to select desired loss function
#'   # formula method
#'   b <- dr4pl(formula = Response~Dose, 
#'              data = sample_data_4,
#'              method.init = "Mead", 
#'              method.robust = "Tukey" )
#'   plot(b)
#'   
#'   #data.frame method
#'   c <- dr4pl(data = sample_data_10,
#'              dose = Dose,
#'              response = Response)
#'   plot(c)
#' 
#'   ##compatable with ggplot
#'   library(ggplot2) #load ggplot2
#'   c <- dr4pl(Response~Dose, 
#'              data = drc_error_2, 
#'              method.optim = "CG", 
#'              trend = "decreasing" )
#'   d <- plot(c, x.breaks = c(.00135, .0135, .135, 1.35, 13.5))
#'   d + theme_grey()
#' @export
dr4pl.default <- function(dose,
                          response,
                          init.parm = NULL,
                          trend = "auto",
                          method.init = "Mead",
                          method.robust = "squared",
                          method.optim = "Nelder-Mead",
                          use.Hessian = FALSE,
                          level = 0.9999,
                          failure.message = FALSE,
                          upperl = NULL,
                          lowerl = NULL,
                          ...) {
  
  types.trend <- c("auto", "decreasing", "increasing")
  types.method.init <- c("logistic", "Mead")
  types.method.optim <- c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN")
  
  ### Check errors in functions arguments
  if(!is.numeric(dose)||!is.numeric(response)) {
    
    stop("Both doses and responses should be numeric.")
  }
  if(any(dose<0)) {
    
    stop("Dose levels should be nonnegative.")
  }
  if(length(dose) == 0 || length(response) == 0 || length(dose) != length(response)) {
    
    stop("The same numbers of dose and response values should be supplied.")
  }
  if(!is.element(method.init, types.method.init)) {
    
    stop("The initialization method name should be one of \"logistic\" and \"Mead\".")
  }
  if(!is.element(method.optim, types.method.optim)) {
    
    stop("The optimization method name should be one of \"Nelder-Mead\", \"BFGS\",
         \"CG\", \"L-BFGS-B\" and \"SANN\".")
  }
  if(!is.element(trend, types.trend)) {
    
    stop("The type of the \"trend\" parameter should be one of \"auto\", \"decreasing\" and \"increasing\".")
  }

  # Fit a 4PL model
  obj <- dr4plEst(dose = dose,
                  response = response,
                  init.parm = init.parm,
                  trend = trend,
                  method.init = method.init,
                  method.robust = method.robust,
                  method.optim = method.optim,
                  use.Hessian = use.Hessian,
                  level = level,
                  upperl = upperl,
                  lowerl = lowerl)

  obj$call <- match.call()
  class(obj) <- "dr4pl"
  
  # If any robust estimation method is indicated, report outliers to a user.
  if(!method.robust=="squared") {
    
    theta <- obj$parameters  # Robust parameter estimates
    residuals <- Residual(theta, dose, response)  # Residuals
    
    indices.outlier <- OutlierDetection(residuals)
    
    obj$idx.outlier <- indices.outlier
    obj$robust.plot <- plot(obj, indices.outlier = indices.outlier)
  }
  
  message.diagnosis <- NULL

  ### When convergence failure happens.
  if(obj$convergence == FALSE) {

    ## Decide the method of robust estimation which is more robust than the method
    ## input by a user.
    if(method.robust=="squared") {
      
      method.robust.new <- "absolute"
    } else if(is.element(method.robust, c("absolute", "Huber"))) {
      
      method.robust.new <- "Tukey"
    } else {
      
      stop("Convergence failure happened but no resolution could be found.")
    }
    
    n <- obj$sample.size  # Number of data points

    # Fit a 4PL model to data
    obj.robust <- dr4plEst(dose = dose,
                           response = response,
                           init.parm = init.parm,
                           trend = trend,
                           method.init = method.init,
                           method.robust = method.robust.new,
                           method.optim = method.optim,
                           use.Hessian = use.Hessian,
                           level = level,
                           upperl = upperl,
                           lowerl = lowerl)

    obj.robust$call <- match.call()
    class(obj.robust) <- "dr4pl"
    
    ## Detect outliers and report them.
    theta <- obj.robust$parameters  # Robust parameter estimates
    residuals <- Residual(theta, dose, response)  # Residuals
    
    indices.outlier <- OutlierDetection(residuals)
    
    obj.robust$idx.outlier <- indices.outlier
    obj.robust$robust.plot <- plot(obj.robust, indices.outlier = indices.outlier)

    ## Print different messages to the console depending on the convergence success
    ## of a robust fit.
    if(obj.robust$convergence) {
      
      message.diagnosis <- 
      paste("The Hill bounds have been hit during optimization, but other robust ",
            "estimation was succesful.\n",
            "Please refer to \"dr4pl.robust\" variable for diagnosis.\n",
            sep = "")
    } else {
      
      message.diagnosis <- 
      paste("The Hill bounds have been hit during optimization with ",
            obj$method.robust, " and ", obj.robust$method.robust, " methods.\n",                
            "Please try other initialization and robust estimation methods.\n", 
            sep = "")
    }
    
    obj$dr4pl.robust <- obj.robust
    obj$message.diagnosis <- message.diagnosis
  }
  
  if(failure.message&&!is.null(message.diagnosis)) {
    
    cat(message.diagnosis)
  }
  
  return(obj)
}

#' @title Private function to fit the 4PL model to dose-response data
#' 
#' @description Private function that actually fits the 4PL model to data. If the
#'   Hill bounds are attained at the end of optimization processes, then an
#'   indicator of convergence failure so that \code{\link{dr4pl.default}} can
#'   look for a remedy for convergence failure.
#' 
#' @name dr4plEst
#' 
#' @param dose Vector of dose levels
#' @param response Vector of responses
#' @param init.parm Vector of initial parameters of the 4PL model supplied by a
#'   user.
#' @param trend Indicator of whether a dose-response curve is a decreasing 
#' \eqn{\theta[3]<0} or increasing curve \eqn{\theta[3]>0}. The default is "auto" 
#' which indicates that the trend of the curve is automatically determined by
#' data. The option "decreasing" will impose a restriction \eqn{\theta[3]<=0} 
#' while the option "increasing" will impose a restriction \eqn{\theta[3]>=0} in an 
#' optimization process.
#' @param method.init Method of obtaining initial values of the parameters.
#' Should be one of "logistic" for the logistic method or "Mead" for the Mead
#' method. The default option is the Mead method.
#' @param method.robust Parameter to select loss function for the robust estimation 
#' method to be used to fit a model. The argument NULL indicates the sum of squares
#' loss, "absolute" indicates the absolute deviation loss, "Huber" indicates Huber's
#' loss and "Tukey" indicates Tukey's biweight loss.
#' @param method.optim Method of optimization of the parameters. This argument
#' is directly delivered to the \code{constrOptim} function provided in the
#' "base" package of R.
#' @param use.Hessian Indicator of whether the Hessian matrix (TRUE) or the
#' gradient vector is used in the Hill bounds.
#' @param level Confidence level to be used in Hill bounds computation.
#' @param upperl upper limit to init.parm
#' @param lowerl lower limit to init.parm
#' 
#' @return List of final parameter estimates, name of robust estimation, loss value
#' and so on.
dr4plEst <- function(dose, response,
                     init.parm,
                     trend,
                     method.init,
                     method.optim,
                     method.robust,
                     use.Hessian,
                     level,
                     upperl,
                     lowerl) {
  
  convergence <- TRUE
  x <- dose  # Vector of dose values
  y <- response  # Vector of responses
  n <- length(x)  # Number of observations
  
  # Choose the loss function depending on the robust estimation method
  loss.fcn <- ErrFcn(method.robust)
  # Currently only the gradient function for the squared loss is implemented
  grad <- GradientSquaredLossLogIC50
  
  tuning.barrier <- 1e-04  # Tuning parameter for the log Barrier method
  
  ### When initial parameter estimates are given
  if(!is.null(init.parm)) {

    # Check whether initial parameter estimates satisfy constraints
    if(init.parm[2] <= 0) {
      
      stop("The IC50 parameter should be positive.")
    }
    
    # Use given initial parameter estimates
    retheta.init <- init.parm
    retheta.init[2] <- log10(init.parm[2])
    
    names(retheta.init) <- c("Upper limit", "Log10(IC50)", "Slope", "Lower limit")
    if(init.parm[1]<init.parm[4]) {
      stop(paste("Choose init.parm such that \u03B8\u2081 is greater than \u03B8\u2084. \n Choose \u03B8\u2083 >0 for increaseing trend and \u03B8\u2083 <0 for decreasing trend. \n"))
    }

    constr.mat <- matrix(c(1, 0, 0, -1), nrow = 1, ncol = 4)
    constr.vec <- 0
    
    # Impose a constraint on the slope parameter based on the function argument
    # "trend".
    if(trend == "decreasing") {
      
      constr.mat <- rbind(constr.mat, matrix(c(0, 0, -1, 0), nrow = 1, ncol = 4))
      constr.vec <- c(constr.vec, 0)
    } else if(trend == "increasing") {
      
      constr.mat <- rbind(constr.mat, matrix(c(0, 0, 1, 0), nrow = 1, ncol = 4))
      constr.vec <- c(constr.vec, 0)
    }
    
    #Impose constraints on upper limit and or lower limit
    #"upperl" and "lowerl"
    if(!is.null(upperl)&!is.null(lowerl)){
      if(any(upperl<lowerl)) {
        stop("upperl must be greater than lowerl")
      }
    }
    if(!is.null(upperl)){
      if(is.numeric(upperl)&&(length(upperl)==4)) {
        if(!is.infinite(upperl[1])){
          constr.mat <- rbind(constr.mat, matrix(c(-1, 0, 0, 0), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, -1*upperl[1])
        }
        if(!is.infinite(upperl[2])){
          constr.mat <- rbind(constr.mat, matrix(c(0, -1, 0, 0), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, -1*log10(upperl[2]))
        }
        if(!is.infinite(upperl[3])){
          constr.mat <- rbind(constr.mat, matrix(c(0, 0, -1, 0), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, -1*upperl[3])
        }
        if(!is.infinite(upperl[4])){
          constr.mat <- rbind(constr.mat, matrix(c(0, 0, 0, -1), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, -1*upperl[4])
        }
      } else {
        stop("upperl must either be a numeric vector of length 4 or NULL")
      }
    } 
    if(!is.null(lowerl)){
      if(is.numeric(lowerl)&&(length(lowerl)==4)){
        if(!is.infinite(lowerl[1])){
          constr.mat <- rbind(constr.mat, matrix(c(1, 0, 0, 0), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, lowerl[1])
        }
        if(!is.infinite(lowerl[2])){
          constr.mat <- rbind(constr.mat, matrix(c(0, 1, 0, 0), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, log10(lowerl[2]))
        }
        if(!is.infinite(lowerl[3])){
          constr.mat <- rbind(constr.mat, matrix(c(0, 0, 1, 0), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, lowerl[3])
        }
        if(!is.infinite(lowerl[4])){
          constr.mat <- rbind(constr.mat, matrix(c(0, 0, 0, 1), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, lowerl[4])
        }
      } else {
        stop("lowerl must either be a numeric vector of length 4 or NULL")
      }
    }
    
    if(any(constr.mat%*%retheta.init<constr.vec)) {
      
      stop(paste("Initial parameter values are not in the interior of the feasible region.\n"))
    }
    
    # Fit a 4PL model to data
    optim.dr4pl <- constrOptim(theta = retheta.init,
                               f = loss.fcn,
                               grad = grad,
                               ui = constr.mat,
                               ci = constr.vec,
                               method = method.optim,
                               hessian = TRUE,
                               x = x,
                               y = y)

    loss <- optim.dr4pl$value
    hessian <- optim.dr4pl$hessian
    retheta <- optim.dr4pl$par
    
    theta <- retheta
    theta[2] <- 10^retheta[2]
    
  ### When initial parameter values are not given.
  } else {
    
    ## Obtain initial parameter estimates.
    theta.init <- FindInitialParms(x, y, trend, method.init, method.robust)
    retheta.init <- ParmToLog(theta.init)
    
    Hill.bounds <- FindHillBounds(x, y, theta.init, use.Hessian, level)
    
    constr.mat <- matrix(rbind(c(1, 0, 0, -1),
                               c(0, 1, 0, 0),
                               c(0, -1, 0, 0),
                               c(0, 0, 1, 0),
                               c(0, 0, -1, 0)),
                         nrow = 5,
                         ncol = 4)
    constr.vec <- c(0, Hill.bounds$LogTheta2[1], -Hill.bounds$LogTheta2[2],
                    Hill.bounds$Theta3[1], -Hill.bounds$Theta3[2])

    # Impose a constraint on the slope parameter based on the function argument
    # "trend".
    if(trend == "decreasing") {
      
      constr.mat <- rbind(constr.mat, matrix(c(0, 0, -1, 0), nrow = 1, ncol = 4))
      constr.vec <- c(constr.vec, 0)
    } else if(trend == "increasing") {
      
      constr.mat <- rbind(constr.mat, matrix(c(0, 0, 1, 0), nrow = 1, ncol = 4))
      constr.vec <- c(constr.vec, 0)
    }
    #Impose constraints on upper limit and or lower limit
    #"upperl" and "lowerl"
    if(!is.null(upperl)&!is.null(lowerl)){
      if(any(upperl<lowerl)) {
        stop("upperl must be greater than lowerl")
      }
    }
    if(!is.null(upperl)){
      if(is.numeric(upperl)&&(length(upperl)==4)) {
        if(!is.infinite(upperl[1])){
          constr.mat <- rbind(constr.mat, matrix(c(-1, 0, 0, 0), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, -1*upperl[1])
        }
        if(!is.infinite(upperl[2])){
          constr.mat <- rbind(constr.mat, matrix(c(0, -1, 0, 0), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, -1*log10(upperl[2]))
        }
        if(!is.infinite(upperl[3])){
          constr.mat <- rbind(constr.mat, matrix(c(0, 0, -1, 0), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, -1*upperl[3])
        }
        if(!is.infinite(upperl[4])){
          constr.mat <- rbind(constr.mat, matrix(c(0, 0, 0, -1), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, -1*upperl[4])
        }
      } else {
        stop("upperl must either be a numeric vector of length 4 or NULL")
      }
    } 
    if(!is.null(lowerl)){
      if(is.numeric(lowerl)&&(length(lowerl)==4)){
        if(!is.infinite(lowerl[1])){
          constr.mat <- rbind(constr.mat, matrix(c(1, 0, 0, 0), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, lowerl[1])
        }
        if(!is.infinite(lowerl[2])){
          constr.mat <- rbind(constr.mat, matrix(c(0, 1, 0, 0), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, log10(lowerl[2]))
        }
        if(!is.infinite(lowerl[3])){
          constr.mat <- rbind(constr.mat, matrix(c(0, 0, 1, 0), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, lowerl[3])
        }
        if(!is.infinite(lowerl[4])){
          constr.mat <- rbind(constr.mat, matrix(c(0, 0, 0, 1), nrow = 1, ncol = 4))
          constr.vec <- c(constr.vec, lowerl[4])
        }
      } else {
        stop("lowerl must either be a numeric vector of length 4 or NULL")
      }
    } 
    

    if(any(constr.mat%*%retheta.init<constr.vec)) {
      
      stop(paste("Initial parameter values are not in the interior of the feasible region.\n",
                 "Estimated Parameters:\n",
                 " UpperLimit: ",theta.init[1],
                 "\n ",ifelse(theta.init[3]<=0,"IC50: ","EC50: "), theta.init[2],
                 "\n Slope: ", theta.init[3],
                 "\n LowerLimit: ", theta.init[4],
                 "\nConsider setting init.parm argument to be within specified bounds of upperl, lowerl, and trend\n"))
    }

    # Fit the 4PL model
    optim.dr4pl <- constrOptim(theta = retheta.init,
                               f = loss.fcn,
                               grad = grad,
                               ui = constr.mat,
                               ci = constr.vec,
                               method = method.optim,
                               hessian = TRUE,
                               mu = tuning.barrier,
                               x = x,
                               y = y)
    
    loss <- optim.dr4pl$value
    hessian <- optim.dr4pl$hessian
    retheta <- optim.dr4pl$par
    
    theta <- LogToParm(retheta)
  }
  
  ### If the Hill bounds are hit.
  if(any(abs(constr.mat%*%retheta - constr.vec)<tuning.barrier)) {
    
    convergence <- FALSE
  } 
  
  # Data frame consisting of doses and responses
  data.dr4pl <- data.frame(Dose = dose, Response = response)

  name.robust <- method.robust
 
  
  list(convergence = convergence,
       data = data.dr4pl,
       hessian = hessian,
       loss.value = loss,
       method.robust = name.robust,
       parameters = theta,
       sample.size = n)
}
