library(dr4pl)
library(rlang)
obj <- dr4pl(Response ~ Dose, sample_data_4)
theta <- coef(obj)
theta_log <- ParmToLog(theta)

#a mix between rlang::quo and substitute - allows
# to sub for functions with env arg
sub2 <- function(expr, env) {
  browser()
  expr_sub <- substitute(expr)
  expr_sub <- do.call("substitute", list(expr_sub, env = env))
  quo_sub <- do.call("quo", list(expr_sub), envir = parent.frame())
  quo_sub
}


test_that("S3 Generic augment works", {
  .data <- augment(obj, sample_data_3)

  expect_identical(colnames(.data), c("Dose","Response",".fitted",".resid"))
  expect_identical(.data$.fitted, MeanResponse(theta, x = .data$Dose))
  expect_identical(.data$.resid, resid(theta, dose = .data$Dose, response = .data$Response))
})

test_that("S3 Generic calculate works", {
  expect_error(calculate(obj), NA)
  expect_error(calculate(obj, theta), NA)
  expect_error(calculate(obj, theta_log), NA)
})

test_that("S3 Generic confint works", {
  expect_error(confint(obj), NA)
  expect_error(confint(obj, theta), NA)
  expect_error(confint(obj, theta_log), NA)
})

test_that("S3 Generic resid/residual works", {

  expect_error(resid(obj), NA)
  expect_error(resid(obj, theta), NA)
  expect_error(resid(obj, theta_log), NA)
  expect_equal(resid(obj), resid(obj, theta_log))
  expect_error(suppressWarnings(resid(obj, c(theta_4 = 0))), "`dr4pl_param` cannot have NA values")
  expect_error(resid(obj, dr4pl_theta(theta_4 = 0)), "`dr4pl_param` cannot have NA values")

})

test_that("MeanResponse Identical Regardless of S3 class", {
  mr <- MeanResponse(obj)
  expect_equal(MeanResponse(theta, obj$data$Dose), mr)
  expect_equal(MeanResponse(theta_log, obj$data$Dose), mr)
  expect_equal(MeanResponse(as.numeric(theta), obj$data$Dose), mr)
  expect_error(suppressWarnings(MeanResponse(obj, c(theta_4 = 0))), "`dr4pl_param` cannot have NA values")
  expect_error(MeanResponse(obj, dr4pl_theta(theta_4 = 0)), "`dr4pl_param` cannot have NA values")
})

test_that("S3 vcov works", {
  expect_error(vcov(obj), NA)
  expect_error(vcov(obj, theta), NA)
  expect_error(vcov(obj, theta_log), NA)

})



test_that("DerivativeF() returns expected values", {
  deriv_linear <- dr4pl:::DerivativeF.dr4pl_theta
  deriv_log10 <- dr4pl:::DerivativeF.dr4pl_log10
  expect_equal_deriv <- function(data) {
    obj <- dr4pl(Response~Dose, data)
    d_linear <- deriv_linear(coef(obj), obj$data$Dose)
    d_log10 <- deriv_log10(ParmToLog(coef(obj)), obj$data$Dose)
    expect_equal(d_linear[,"deriv.f.theta.1"], d_log10[,"deriv.f.theta.1"])
    expect_equal(log(10)*coef(obj)[2]*d_linear[,"deriv.f.theta.2"], d_log10[,"deriv.f.theta.2"])
    expect_equal(d_linear[,"deriv.f.theta.3"], d_log10[,"deriv.f.theta.3"])
    expect_equal(d_linear[,"deriv.f.theta.4"], d_log10[,"deriv.f.theta.4"])
  }
  data_list <- lapply(paste0("sample_data_",1:13), get)

  lapply(data_list, expect_equal_deriv)
})

test_that("Hessian() returns expected values", {
  H_lin <- dr4pl:::Hessian.dr4pl_theta
  H_log <- dr4pl:::Hessian.dr4pl_log10
  expect_equal_Hessian <- function(data) {
    obj <- dr4pl(Response~Dose, data)
    theta <- coef(obj)
    d_lin <- H_lin(theta, obj$data$Dose, obj$data$Response)
    d_log <- H_log(ParmToLog(theta), obj$data$Dose, obj$data$Response)
    d_theta <- log(10)*theta[2]
    expect_equal(d_lin[,"deriv.f.theta.1"]*c(1,d_theta,1,1), d_log[,"deriv.f.theta.1"])
    #second deriv is hard to test....
    #expect_equal(d_lin[,"deriv.f.theta.2"]*d_theta*c(1,d_theta,1,1), d_log[,"deriv.f.theta.2"], tolerance = 1e-4)
    expect_equal(d_lin[,"deriv.f.theta.3"]*c(1,d_theta,1,1), d_log[,"deriv.f.theta.3"])
    expect_equal(d_lin[,"deriv.f.theta.4"]*c(1,d_theta,1,1), d_log[,"deriv.f.theta.4"])
  }
  data_list <- lapply(paste0("sample_data_",1:13), get)

  lapply(data_list, expect_equal_Hessian)
})
