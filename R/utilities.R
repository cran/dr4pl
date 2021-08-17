
# Utility functions

# if lhs is null, take rhs, if not null, coerce lhs to a dr4pl_theta with rhs as reference.
`%theta%` <- function(lhs, rhs) if(is.null(lhs)) rhs else as_dr4pl_param(lhs, allow.NA = F)


replace_theta <- function(x, replace) {
  repl_nm <- deparse(substitute(replace))
  if (is_named(x) && is_named(replace)){
    x_nms <- names(x)
    r_nms <- names(replace)
    abortifnot(
      all(r_nms%in%x_nms),
      glue("Cannot replace named argument ",
           glue_collapse(glue("`{setdiff(r_nms, x_nms)}`"), sep = ", "),
           " with any of the following names:\n\t",
           glue_collapse(glue("`{x_nms}`"), sep = ", "),
           "\nDo you have a typo?")
    )
    lgl <- x_nms %in% r_nms[!r_nms %in% names(replace)[is.na(replace)]]
    return(replace(x, lgl, replace[x_nms[lgl]]))
  } else if (length(replace)==4) {
    return(replace)
  } else if (!is.null(replace)){
    warn(glue("{repl_nm} is niether length 4 or named. Using defaults:\n{print_output(x)}"),)
  }
  return(x)
}

# @title check positive semi definite matrix
# @name is_positive_semi_definite
# @description drop and replace for matrixcalc::is.positive.semi.definite.
# This is the only function we use from said package, thus it does not
# seem right to import this with dr4pl
# @keywords internal
is_positive_semi_definite <- function(x, tol = 1e-08) {
  if (!is.matrix(x))
    stop("argument x is not a matrix")
  if (!is.numeric(x))
    stop("argument x is not a numeric matrix")
  if (nrow(x)!=ncol(x))
    stop("argument x is not a square matrix")
  if (!sum(x==t(x)) == (nrow(x)^2))
    stop("argument x is not a symmetric matrix")
  eigenvalues <- eigen(x, only.values = TRUE)$values
  n <- nrow(x)
  for (i in 1:n){
    if (abs(eigenvalues[i] < tol)) {
      eigenvalues[i] <- 0
    }
  }
  if (any(eigenvalues < 0))
    return(FALSE)
  return(TRUE)
}


abortifnot <- function(expr, msg = NULL) {
  expr_quo <- enquo(expr)
  msg <- msg %||% glue('expression `{as_label(expr_quo)}` was not TRUE')
  if(!eval_tidy(expr_quo))
    abort(msg)
}

print_output <- function(x, ...) {
  out <- capture.output(print(x, ...))
  paste0(out, collapse = "\n")
}

check_feasible_constraints <- function(theta, mat, vec) {
  res <- mat%*%theta<vec
  if(any(res)) {
    abort(
      glue(
        "Initial parameter values are not in the interior of the feasible region.\n",
        "Estimated Parameters:\n",
        print_output(theta),
        "\nFailed Constraints:\n",
        stop_constr_msg(theta, mat[res,, drop = F], vec[res])
      ), class = "dr4pl_infeasible_constraints"
    )
  }
}

stop_constr_msg <- function(theta, mat, vec){
  mat <- split(mat, seq_len(nrow(mat)))
  out <- mapply(constr_error_msg, theta = list(theta), mat = mat, vec = vec)
  paste0(out, "\n")
}
constr_error_msg <- function(theta, mat, vec){
  index <- mat!=0
  mat_vals <- mat[index]
  signs <- sign(mat_vals)*sign(theta[index])
  scalar <- ifelse(abs(mat_vals)==1, "", paste0(abs(mat_vals),"*"))
  msg_vals <- paste0(scalar,abs(theta[index]))
  msg_vals <- weave_signs(msg_vals, signs)
  wi <- which(index)
  theta_msg <- paste("theta",wi, sep = "_")
  if(any(wi==2)){
    theta_msg[wi==2] <- paste0("log10(",theta_msg[wi==2],")")
  }
  theta_msg <- paste0(scalar, theta_msg)
  theta_msg <- weave_signs(theta_msg, sign(mat_vals))


  paste(theta_msg, "=", msg_vals, " >= ", vec)

}

weave_signs <- function(x, signs) {
  n <- length(x)
  if(n>1){
    x[2:n] <- paste(ifelse(signs[2:n]>=0,"+","-"), x[2:n])
  }

  x[1] <- paste0(ifelse(signs[1]>=0,"","- "), x[1])
  paste(x ,collapse = " ")
}


report_class <- function(x) {
  glue("{as_label(enquo(x))}: <{glue_collapse(class(x), sep = '/')}>")
}
