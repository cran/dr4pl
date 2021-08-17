
#' @importFrom generics augment
#' @export
generics::augment

#' @title Augment data with dr4pl
#' @name dr4pl-augment
#' @param x dr4pl object
#' @param data new data to use
#' @param ... For future extension. Should not be used.
#' @export
augment.dr4pl <- function(x,
                          data = NULL,
                          ...) {
  if (...length()>0) abort("... were used in `augment.dr4pl`. These dots are only present for future extension and should be empty.")
  data <- data %||% x$data
  mapping <- parse_dr4pl_mapping(x$call)
  if (inherits(mapping, "no_mapping")) abort("cannot use `augment.dr4pl` when dr4pl object was constructed with the `dr4pl.default` method.\nConstruct dr4pl with formula or data.frame method.")
  .dose <- eval(mapping$Dose, data)
  .resp <- eval(mapping$Response, data)
  .fitted <- MeanResponse(coef(x), x = .dose)
  .resid <- residuals(coef(x), dose = .dose, response = .resp)
  #.cooks <- cooks.distance(x)
  data[['.fitted']] <- .fitted
  data[['.resid']] <- .resid
 # data[['.cooksd']] <- .cooks
  if(requireNamespace("tibble", quietly = T)){
    return(tibble::as_tibble(data))
  } else {
    return(data)
  }

}


parse_dr4pl_mapping <- function(call){
 # browser()
  args <- call_args(call)
  args <- switch(names(args)[1],
         formula = mapping_parser("parse_formula", args),
         data = mapping_parser("parse_data", args),
         mapping_parser("parse_default", args))
  arg_parse(args)
}

mapping_parser <- function(.class, args) {
  structure(args, class = .class)
}

arg_parse <- function(args) UseMethod('arg_parse')
arg_parse.parse_formula <- function(args) {
  list(
    Response = args[[1]][[2]],
    Dose = args[[1]][[3]]
  )
}
arg_parse.parse_data <- function(args) {
  list(
    Response = args$response,
    Dose = args$dose
  )
}

arg_parse.parse_default <- function(args) {
  structure(list(), class = "no_mapping")
}

call_args <- function(call){
  args <- as.list(call[-1])
  stats::setNames(args, nm = names(args))
}
