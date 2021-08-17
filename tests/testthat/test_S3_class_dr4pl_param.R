library(dr4pl)

theta_lin <- dr4pl_theta(36000, 50, -.3, 0)
theta_log <- dr4pl_theta(36000, log10(50), -.3, 0, isLog10 = T)

test_that("dr4pl_theta constructs correct classes", {

  expect_s3_class(theta_lin, "dr4pl_param")
  expect_s3_class(theta_lin, "dr4pl_theta")

  expect_s3_class(theta_log, "dr4pl_param")
  expect_s3_class(theta_log, "dr4pl_log10")

})

test_that("dr4pl_theta() posses correct names/attributes", {
  #theta linear
  expect_equal(names(theta_lin), c("theta_1","theta_2","theta_3","theta_4"))
  expect_equal(attr(theta_lin, "param_info"), c("UpperLimit","IC50", "Slope", "LowerLimit"))
  expect_equal(attr(dr4pl_theta(theta_3 = .3), "param_info"), c("UpperLimit","EC50", "Slope", "LowerLimit"))
  #theta log10
  expect_equal(names(theta_log), c("theta_1","theta_2","theta_3","theta_4"))
  expect_equal(attr(theta_log, "param_info"), c("UpperLimit","Log10(IC50)", "Slope", "LowerLimit"))
  expect_equal(attr(dr4pl_theta(theta_3 = .3, isLog10 = T), "param_info"), c("UpperLimit","Log10(EC50)", "Slope", "LowerLimit"))
})

test_that("dr4pl_theta() throws expected errors", {
  expect_error(dr4pl_theta(100, -1, -1, 0), "should be greater than 0.")
  expect_error(dr4pl_theta(100, -1, -1, 0, isLog10 = T), NA)
  expect_error(dr4pl_theta(theta_1 = 0, theta_4 = 100), "dr4pl theta_1 should be greater than theta_4")
  expect_error(dr4pl_theta(), NA)
})

test_that("dr4pl_param values are identical", {
  expect_equal(ParmToLog(theta_lin), theta_log)
  expect_equal(LogToParm(theta_log), theta_lin)
})

test_that("dr4pl:::new_dr4pl_param throws expected errors", {
  new_param <- dr4pl:::new_dr4pl_param
  expect_error(new_param(100), "`dr4pl_param` should be of length 4")
  expect_error(new_param(c(100, 1, -1, NA), "dr4pl_theta"), "`dr4pl_param` cannot have NA values")
})


test_that("`init.parm` can take named numeric", {
  #ensure sentinel value is available to print message
  assign("dr4pl_param_warn_numeric", NULL, envir = rlang:::warning_freq_env)
  expect_warning(dr4pl(Response~Dose, sample_data_1, init.parm = c(theta_4 = 0)), regexp =  "A numeric object is being coerced to a \"dr4pl_param\" object")
  expect_error(dr4pl(Response~Dose, sample_data_1, init.parm = list(theta_4 = 0)), "no applicable method")
})

test_that("`upperl` and `lowerl` cantake named numeric", {
  expect_error(dr4pl(Response~Dose, sample_data_1, lowerl = c(theta_4 = 0)), NA)
  expect_error(dr4pl(Response~Dose, sample_data_1, lowerl = list(theta_4 = 0)), "lowerl should be a numeric vector")
  expect_error(dr4pl(Response~Dose, sample_data_1, upperl = list(theta_1 = 1e5)), "upperl should be a numeric vector")
})
