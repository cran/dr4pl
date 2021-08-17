

#' @name dr4pl-param
#' @title Constructor for dr4pl theta parameter
#' @description As of version 2.0.0, dr4pl will
#' require the theta parameter to be made by this
#' function. This is to ensure the user is explicit
#' about what their parameter is, and what is being
#' optimized.
#' @param theta_1 Numeric. The upper asymptote of the 4pl model.
#' @param theta_2 Numeric. The value where Response is half way
#' between the upper and lower asymptotes. This parameter
#' is other wise known as IC50/EC50.
#' @param theta_3 Numeric. The slope of the 4pl model. When `theta_3>0`
#' the curve increases, and when `theta_3<0`, the curve
#' decreases.
#' @param theta_4 Numeric. The lower asymptote of the 4pl model.
#' @param isLog10 Logical value indicating if the second parameter
#' is on the log10 scale. The default, FALSE, informs dr4pl that
#' the value passed to the `theta_2` argument is on the linear
#' scale. TRUE will indicate that `theta_2` has already been
#' transformed with \link[base]{log10}. This ensures that the
#' parameters are handled properly later on.
#' @param x a dr4pl_param object to convert
#' @details
#' The function `dr4pl_theta` is a constructor function for the
#' end user. While the default values for `theta_1`, `theta_2`,
#' `theta_3`, `theta_4` are `NA`, certain functions of
#' `dr4pl` will not allow `NA` values. When `dr4pl_theta` is
#' used with the `init.parm` argument of [dr4pl], then
#' values with `NA` will be estimated, and those specified,
#' will be set as the initial parameter prior to optimization.
#' If `dr4pl_theta` is the object of S3 dispatch, then
#' no NA values are allowed. However if the object if
#' dispatch is a `dr4pl` object, and `dr4pl_theta` is passed
#' into the function as an additional object such as in [X], [Y],
#' then parameter non-NA values will be replaced in the appropriate
#' dr4pl parameter estimates for the purpose of said function.
#' @return an object of class `"dr4pl_param"`
#' @aliases dr4pl_param dr4pl_theta
#' @export
dr4pl_theta <- function(theta_1 = NA,
                        theta_2 = NA,
                        theta_3 = NA,
                        theta_4 = NA,
                        isLog10 = FALSE) {
  type <- if (isLog10) "dr4pl_log10" else "dr4pl_theta"
  new_dr4pl_param(x = c(theta_1, theta_2, theta_3, theta_4), type = type, allow.NA = T)
}

#' @rdname dr4pl-param
#' @param x a dr4pl_param object
#' @export
ParmToLog <- function(x) UseMethod('ParmToLog')
#' @rdname dr4pl-param
#' @export
ParmToLog.dr4pl_theta <- function(x){
  if(!is.na(x[2])){
    x[2] <- log10(x[2])
  }
  class(x) <- replace(class(x), class(x)=="dr4pl_theta", "dr4pl_log10")
  name_param(x)
}
#' @rdname dr4pl-param
#' @export
ParmToLog.dr4pl_log10 <- function(x) {
 # warning("Called ParmToLog when theta is already Log space.")
  x
}
#' @rdname dr4pl-param
#' @export
LogToParm <- function(x) UseMethod("LogToParm")
#' @rdname dr4pl-param
#' @export
LogToParm.dr4pl_theta <- function(x) {
  #warning("called LogToParm when theta is already Linear space.")
  x
}
#' @rdname dr4pl-param
#' @export
LogToParm.dr4pl_log10 <- function(x) {
  if (!is.na(x[2])){
    x[2] <- 10^(x[2])
  }
  class(x) <- replace(class(x), class(x)=="dr4pl_log10", "dr4pl_theta")
  name_param(x)
}

new_dr4pl_param <- function(x, type = NA, allow.NA = F) {
  obj_class <- if (!is.na(type)) c(type, "dr4pl_param") else "dr4pl_param"
  x <- structure(x, class = obj_class)
  valid_dr4pl_param(x, allow.NA = allow.NA)
  name_param(x)
}



new_dr4pl_theta <- function(x) new_dr4pl_param(x, type = "dr4pl_theta")
new_dr4pl_log10 <- function(x) new_dr4pl_param(x, type = "dr4pl_log10")

as_dr4pl_param <- function(x, ...) UseMethod("as_dr4pl_param")
as_dr4pl_param.numeric <- function(x, allow.NA = F, ...) {
  warn(glue('A numeric object is being coerced to a "dr4pl_param" object.',
            '`theta_2` is assumed to be in linear space. Please use ',
            '`dr4pl_theta()` to construct the theta parameter.'), 
       class = "dr4pl_param_warn",
       .frequency = "regularly",
       .frequency_id = "dr4pl_param_warn_numeric", x = x)
  x <- replace_theta(dr4pl_theta(), x)
  new_dr4pl_param(x, type = "dr4pl_theta", allow.NA = allow.NA)
}
as_dr4pl_param.dr4pl_param <- function(x, allow.NA = F, ...) {
  valid_dr4pl_param(x, allow.NA = allow.NA)
  x
}


# as_dr4pl_theta <- function(x, ...) UseMethod("as_dr4pl_theta")
# as_dr4pl_theta.numeric <- function(x, reference = dr4pl_theta(), allow.NA = F, ...) {
#   warning('A numeric object is being coerced to "dr4pl_param" object.',
#           '`theta_2` is assumed to be in linear space. Please use',
#           ' `dr4pl_theta()` to construct the theta parameter.', call. = F)
#   x <- replace_theta(reference, x)
#   new_dr4pl_param(x, type = "dr4pl_theta", allow.NA = allow.NA)
# }
# as_dr4pl_theta.dr4pl_theta  <- function(x, reference = dr4pl_theta(), allow.NA = F, ...) {
#   x <- replace_theta(reference, x)
#   new_dr4pl_param(x, type = "dr4pl_theta", allow.NA = allow.NA)
# }
# as_dr4pl_theta.dr4pl_log10 <- function(x, reference = dr4pl_theta(), allow.NA = F, ...) {
#   x <- replace_theta(reference, LogToParm(x))
#   new_dr4pl_param(x, type = "dr4pl_theta", allow.NA = allow.NA)
# }

name_param <- function(theta) UseMethod('name_param')
name_param.dr4pl_param <- function(theta) {
  names(theta) <- c("theta_1","theta_2","theta_3","theta_4")
  param_info(theta) <- param_names(theta)
  theta
}
name_param.numeric <- name_param.dr4pl_param
param_names <- function(theta) UseMethod('param_names')
param_names.dr4pl_param <- function(theta) {
  if(is.na(theta[3])||theta[3]<=0) {
    c("UpperLimit", "IC50", "Slope", "LowerLimit")
  } else {
    c("UpperLimit", "EC50", "Slope", "LowerLimit")
  }
}
param_names.numeric <- param_names.dr4pl_param
param_names.dr4pl_log10 <- function(theta) {
  if(is.na(theta[3])||theta[3]<=0) {
    c("UpperLimit", "Log10(IC50)", "Slope", "LowerLimit")
  } else {
    c("UpperLimit", "Log10(EC50)", "Slope", "LowerLimit")
  }
}

param_info <- function(theta) {
  attr(theta, "param_info")
}

`param_info<-` <- function(theta, value){
  attr(theta, "param_info") <- value
  theta
}

#' @export
print.dr4pl_param <- function(x, ...){
  txt <- mapply(FUN = function(nam, val){
    format(c(nam, as.character(val)), justify = "right")
  }, nam = param_info(x), val = x)
  txt <- apply(txt, 1,paste, collapse = "\t")
  txt <- paste(txt, collapse = "\n")
  cat(txt, "\n")
  invisible(x)
}



valid_dr4pl_param <- function(theta, allow.NA = FALSE) UseMethod("valid_dr4pl_param")
valid_dr4pl_param.numeric <- function(theta, allow.NA = FALSE) {
  abort("Numeric passed to `theta` argument. Please either construct theta with `dr4pl_theta` or extract it from a `dr4pl` object.")
}

valid_dr4pl_param.dr4pl_param <- function(theta, allow.NA = FALSE) {
  #require length 4
  abortifnot(
    length(theta)==4,
    "`dr4pl_param` should be of length 4.")
  # If NA's are not allowed - stop if NA's are detected
  if (!allow.NA) {
    abortifnot(
      !anyNA(theta),
      "`dr4pl_param` cannot have NA values.")
  }
  # enforce theta[1]>theta[4] iff both are already numeric and not NA
  if (!is.na(theta[1])&&!is.na(theta[4])) {
    abortifnot(
      theta[1]>theta[4],
      paste0("dr4pl theta_1 should be greater than theta_4,",
             " see `?dr4pl::dr4pl` for more details"))
  }
  invisible(theta)
}

valid_dr4pl_param.dr4pl_theta <- function(theta, allow.NA = FALSE) {
  NextMethod('valid_dr4pl_param')
  #when on linear scale and theta[2] is numeric
  # enforce a positive theta
  if (!is.na(theta[2])) {
    abortifnot(
      theta[2]>0,
      paste0("dr4pl theta_2 was assigned ", signif(theta[2], 3), ". dr4pl theta_2 should be greater than 0.")
    )
  }

}

