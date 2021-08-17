


is_named <- function(x) !is.null(attr(x, 'names'))

dr4pl_theta_limits_ptype <- function(value = NA_real_) {
  structure(c(theta_1 = value,
              theta_2 = value,
              theta_3 = value,
              theta_4 = value))
}

dr4pl_theta_limits <- function(lowerl = NULL, upperl = NULL) {
  abortifnot(is.null(lowerl)||is.numeric(lowerl), "lowerl should be a numeric vector")
  abortifnot(is.null(upperl)||is.numeric(upperl), "upperl should be a numeric vector")
  structure(
    list(
      upperl = replace_theta(dr4pl_theta_limits_ptype(Inf), upperl),
      lowerl = replace_theta(dr4pl_theta_limits_ptype(-Inf), lowerl)
    ),
    class = "dr4pl_theta_limits"
  )
}

print.dr4pl_theta_limits <- function(x, ...){
  out <- data.frame(lowerl = x$lowerl,
                    upperl = x$upperl)
  print(out)
}

mold_constr <- function(object, ...) UseMethod("mold_constr")

mold_constr.dr4pl_theta_limits <- function(object, trend = "auto") {

  constr.mat <- matrix(c(1, 0, 0, -1), nrow = 1, ncol = 4)
  constr.vec <- 0
  constr.type <- ""
  upperl <- object$upperl
  lowerl <- object$lowerl
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
  if(any(upperl<lowerl)) {
    abort(glue("upperl must be greater than lowerl\n{print_output(const[upperl<lowerl,]}"))
  }
  if(!is.infinite(upperl[1])){
    constr.mat <- rbind(constr.mat, matrix(c(-1, 0, 0, 0), nrow = 1, ncol = 4))
    constr.vec <- c(constr.vec, -1*upperl[1])
    constr.type <- c(constr.type, "upperl")
  }
  if(!is.infinite(upperl[2])){
    constr.mat <- rbind(constr.mat, matrix(c(0, -1, 0, 0), nrow = 1, ncol = 4))
    constr.vec <- c(constr.vec, -1*log10(upperl[2]))
    constr.type <- c(constr.type, "upperl")
  }
  if(!is.infinite(upperl[3])){
    constr.mat <- rbind(constr.mat, matrix(c(0, 0, -1, 0), nrow = 1, ncol = 4))
    constr.vec <- c(constr.vec, -1*upperl[3])
    constr.type <- c(constr.type, "upperl")
  }
  if(!is.infinite(upperl[4])){
    constr.mat <- rbind(constr.mat, matrix(c(0, 0, 0, -1), nrow = 1, ncol = 4))
    constr.vec <- c(constr.vec, -1*upperl[4])
    constr.type <- c(constr.type, "upperl")
  }

  if(!is.infinite(lowerl[1])){
    constr.mat <- rbind(constr.mat, matrix(c(1, 0, 0, 0), nrow = 1, ncol = 4))
    constr.vec <- c(constr.vec, lowerl[1])
    constr.type <- c(constr.type, "lowerl")
  }
  if(!is.infinite(lowerl[2])){
    constr.mat <- rbind(constr.mat, matrix(c(0, 1, 0, 0), nrow = 1, ncol = 4))
    constr.vec <- c(constr.vec, log10(lowerl[2]))
    constr.type <- c(constr.type, "lowerl")
  }
  if(!is.infinite(lowerl[3])){
    constr.mat <- rbind(constr.mat, matrix(c(0, 0, 1, 0), nrow = 1, ncol = 4))
    constr.vec <- c(constr.vec, lowerl[3])
    constr.type <- c(constr.type, "lowerl")
  }
  if(!is.infinite(lowerl[4])){
    constr.mat <- rbind(constr.mat, matrix(c(0, 0, 0, 1), nrow = 1, ncol = 4))
    constr.vec <- c(constr.vec, lowerl[4])
    constr.type <- c(constr.type, "lowerl")
  }

  res <- structure(list(
    mat = constr.mat,
    vec = constr.vec,
    type = constr.type
  ), class = "dr4pl_cnstr")
  return(res)
}

mold_constr.dr4pl_hill <- function(object) {
  constr.mat <- matrix(rbind(c(0, 1, 0, 0),
                             c(0, -1, 0, 0),
                             c(0, 0, 1, 0),
                             c(0, 0, -1, 0)),
                       nrow = 4,
                       ncol = 4)
  constr.vec <- c(object$LogTheta2[1], -object$LogTheta2[2],
                  object$Theta3[1], -object$Theta3[2])
  constr.type <- rep("hill", 4)
  res <- structure(list(
    mat = constr.mat,
    vec = constr.vec,
    type = constr.type
  ), class = "dr4pl_cnstr")
  return(res)
}

dr4pl_cnstr <- function(mat = matrix(numeric(), ncol = 4),
                        vec = numeric(),
                        type = character()) {
  structure(
    list(
      mat = mat,
      vec = vec,
      type = type
    ),
    class = "dr4pl_cnstr"
  )
}

add_dr4pl_cnstr <- function(object = dr4pl_cnstr(), add, ...) {
  abortifnot(inherits(object, "dr4pl_cnstr"))
  add <- mold_constr(add, ...)
  structure(
    list(
      mat = rbind(object$mat, add$mat),
      vec = c(object$vec, add$vec),
      type = c(object$type, add$type)
    ),
    class = "dr4pl_cnstr"
  )
}

`[.dr4pl_cnstr` <- function(x, i) {
  structure(
    list(
      mat = x$mat[i,, drop = F],
      vec = x$vec[i],
      type = x$type[i]
    ),
    class = "dr4pl_cnstr"
  )

}

