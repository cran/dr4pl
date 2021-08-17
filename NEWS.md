---
output: md_document
---

# dr4pl 2.0.0

## Main Reason for Major Update
* A second derivative used in the Hessian Matrix (with respects to `(theta_3)^2`) has been updated to the correct value. Confidence intervals for `theta_3` should be more accurate as a result of the correction. See 
* dr4pl's internals have changed to help improve clarity and consistency of the code. With this change, the API associated with the `dr4pl` call and its arguments have changed slightly. Many dr4pl internal functions now dispatch on this new `dr4pl_param` class.

## General Changes
* The dr4pl object has replaced 'parameters' with a new S3 object `dr4pl_param`. This class's main objective is to track of when the object is reparameterized internally. There is now a constructor function `dr4pl_theta` that will allow users to specify if their parameters are in the linear space or the log10 space when passing to `init.parm` argument. 
* `dr4pl_param` objects are converted to log10 space prior to optimization. To be consistent with previous version of dr4pl, reporting statistics will be on the reparameterized theta by default and will forgo any exponentiation to linear space. This is to increase clarity of which parameterization statistics like confidence intervals and standard errors refer to.
* The following functions are new S3 generics:
  * `augment`
  * `calculate`
  * `residuals`
* Many S3 methods that dispatch on a `dr4pl` object now have an optional `parm` argument. This can allow the user to specify a new `dr4pl_theta` object to evaluate in place of the theta parameter on the object. This applies to the following S3 generics:
  * `calculate`
  * `residuals`
  * `summary`
  * `vcov`
* `confint.dr4pl` and `vcov.dr4pl`'s `parm` argument now works as intended.
* You are no longer required to specify all 4 parameters when passing to arguments  `upperl` and `lowerl`. Specifying a named numeric vector in the form of `lowerl = c(theta_2=.5, theta_4 = 0)` will suffice.
* More informative error messages have been included if initial parameter estimates are outside the specified constraints.


# dr4pl 1.1.9

Bug fix update.
  dr4pl.data.frame will no longer fail when using tibbles as an input.
  "L-BFGS-B" is now a usable argument to the parameter 'method.optim' as originally intended.

General changes and updates 
  default value of method.robust is now set to "squared" to be more `conistent` with other values given.
  Added `upperl` and `lowerl` arguments to dr4pl call. User can now constrain optimization of theta parameters.
  This may require user to set `init.parm` as the estimated theta parameters will not always be within
  the feasible region specified by `upperl` and `lowerl`.
    
# dr4pl 1.1.8

Following proper semantic versioning.

# dr4pl 1.1.7

print.summary.dr4pl function no long provides t-statistics and p-values but now prints 95% confidence intervals instead. Consolidated Vignettes.

# dr4pl 1.1.6

dr4pl function now allows 'data' argument when using 'dose' and 'response' arguments. See dr4pl examples.






