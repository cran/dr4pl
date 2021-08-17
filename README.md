\#build status:

[![Build
Status](https://travis-ci.org/aubreybailey/dr4pl.svg?branch=master)](https://travis-ci.org/aubreybailey/dr4pl)

\#license:

[![Licence](https://img.shields.io/badge/licence-GPL--3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

\#cran status:

[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/dr4pl)](https://cran.r-project.org/package=dr4pl)

\#release version:

[![packageversion](https://img.shields.io/badge/GitHub%20Package%20version-2.0.0-orange.svg?style=flat-square)](https://github.com/DittmerLabUNC/dr4pl)

\#\#dr4pl

The package dr4pl (Dose Response 4 Parameter Logisitic model)
specializes in applying the 4 Parameter Logistic (4PL) model. The 4PL
model has been recognized as a major tool to analyze the relationship
between a dose and a response in pharmacological experiments. The
package dr4pl may be used to model increasing and decreasing curves. The
goal of dr4pl is to bring a statistical method which is capable of
handeling specific error cases of which other statistical packages
produce errors. Examples of Dose Response datasets that will produce
errors in other packages may be accessed by name once dr4pl is loaded
and these data sets are under the names of drc\_error\_1, drc\_error\_2,
drc\_error\_3, and drc\_error\_4. Along with these error data sets, this
package also supplies 13 standard example data sets for the 4PL model
under the name sample\_data\_1, sampel\_data\_2, etc. The package dr4pl
also alows for the user to decide how their theta variable is
approximated. The user may choose the default logistic model or use
Mead’s Method. Additionally, the user may decide between four loss
functions to minimize: Squared, Absolute, Huber, or Tukey’s biweight.
Please attempt each of the loss functions and choose the best fit from
plotting the dr4pl object.

Installation
------------

You can install dr4pl from github with:

``` r
# install.packages("devtools")
devtools::install_bitbucket("dittmerlab/dr4pl")
```

Example
-------

This is a basic example which shows you how to solve a common problem.
This example may be used with drc\_error\_1, drc\_error\_2,
drc\_error\_3, and drc\_error\_4:

``` r
## basic example code, datasets
## example requires the drc and dr4pl package to be loaded
library(dr4pl)
library(drc)
#> Loading required package: MASS
#> 
#> 'drc' has been loaded.
#> Please cite R and 'drc' if used for a publication,
#> for references type 'citation()' and 'citation('drc')'.
#> 
#> Attaching package: 'drc'
#> The following objects are masked from 'package:stats':
#> 
#>     gaussian, getInitial
a <- drc::drm(drc_error_1$Response~drc_error_1$Dose, fct = LL.4())
#> Error in optim(startVec, opfct, hessian = TRUE, method = optMethod, control = list(maxit = maxIt,  : 
#>   non-finite finite-difference value [4]
#> Error in drmOpt(opfct, opdfct1, startVecSc, optMethod, constrained, warnVal, : Convergence failed
plot(a)
#> Error in plot(a): object 'a' not found
```

``` r
## basic example code
## example requires the dr4pl package to be loaded
b <- dr4pl(drc_error_1$Response~drc_error_1$Dose, method.init = "logistic", method.robust = "Tukey") 
plot(b)
#> Warning: Transformation introduced infinite values in continuous x-axis

#> Warning: Transformation introduced infinite values in continuous x-axis
```

![](README_files/figure-markdown_github/example_solution-1.png)

``` r
summary(b)
#> Call:
#> dr4pl.formula(formula = drc_error_1$Response ~ drc_error_1$Dose, 
#>     method.init = "logistic", method.robust = "Tukey")
#> 
#>                Estimate      StdErr       2.5 %     97.5 %
#> UpperLimit   7.9134e+04  6.5396e+01  7.9004e+04 79262.7084
#> Log10(IC50) -1.2371e+01  2.5194e+00 -1.7347e+01    -7.3946
#> Slope       -7.3707e-02  2.2747e-02 -1.1864e-01    -0.0288
#> LowerLimit  -8.3931e+03  6.5427e+01 -8.5223e+03 -8263.8406
```

Updates
-------

With the release of `dr4pl (>= 2.0.0)`, `dr4pl` is getting a few quality
improvements as well as more accurate confidence intervals. First and
foremost, prior to `dr4pl (>= 2.0.0)`, `dr4pl` possessed an error in the
second derivative of the model function with respects to `(theta_3)^2`.
This has been corrected and a new vignette “dr4pl\_derivatives” exists
showing how each derivative was calculated for this package.

The `parameters` field of the `dr4pl` object has changed. This field
will now have an object inheriting from `dr4pl_param` S3 class. This is
likely the biggest change as a lot of internals now dispatch on this
object. The primary function of this object is to track if the `theta_2`
parameter is in the log10 space. Users were often confused as to which
parameter estimate they were looking at, to remedy this, when printing
the object, it will explicitly tell you the type.

``` r
obj <- dr4pl(Response~Dose, sample_data_3)
coef(obj)
#>       UpperLimit             IC50                 Slope        LowerLimit
#> 59858.0700194411 5.07949182697285    -0.611737121034004  980.820568608382
```

For backwards compatibility, `dr4pl` always stores the `dr4pl_theta`
version of the `dr4pl_param` object. The user can change the state of
the parameter with the `ParmToLog()` and `LogToParm()` functions.

``` r
ParmToLog(coef(obj))
#>       UpperLimit       Log10(IC50)                Slope        LowerLimit
#> 59858.0700194411 0.705820265870361   -0.611737121034004  980.820568608382
```

`dr4pl`’s convention was to always report the values in linear space
(even if the Confidence Intervals reflected the log10 parameter), but
this understandably so only increased confusion around how the values
were calculated. To remedy this, various S3 functions that dispatch with
`dr4pl` have an optional `theta` argument in which you can force a
calculation on a parameter set. To stay somewhat backwards compatible,
`dr4pl` will use the `log10` space by default since parameters are
always estimated in that space.

``` r
summary(obj) #uses log10 space by default
#> Call:
#> dr4pl.formula(formula = Response ~ Dose, data = sample_data_3)
#> 
#>                Estimate      StdErr       2.5 %     97.5 %
#> UpperLimit   5.9858e+04  6.1255e+02  5.8616e+04 61100.3693
#> Log10(IC50)  7.0582e-01  8.1390e-02  5.4075e-01     0.8709
#> Slope       -6.1174e-01  6.4982e-02 -7.4353e-01    -0.4799
#> LowerLimit   9.8082e+02  6.2518e+02 -2.8710e+02  2248.7443
summary(obj, theta = coef(obj)) #grab linear space of dr4pl_param
#> Call:
#> dr4pl.formula(formula = Response ~ Dose, data = sample_data_3)
#> 
#>                Estimate      StdErr       2.5 %     97.5 %
#> UpperLimit   5.9858e+04  6.1255e+02  5.8616e+04 61100.3693
#> Log10(IC50)  7.0582e-01  8.1390e-02  5.4075e-01     0.8709
#> Slope       -6.1174e-01  6.4982e-02 -7.4353e-01    -0.4799
#> LowerLimit   9.8082e+02  6.2518e+02 -2.8710e+02  2248.7443
```

Since the `dr4pl_param` object’s displayed names may change, the user is
able to programmatically select elements with theta\_1, theta\_2,
theta\_3, theta\_4 or by integer index.

``` r
ParmToLog(coef(obj))["theta_2"]
#>   theta_2 
#> 0.7058203
coef(obj)[2]
#>  theta_2 
#> 5.079492
```

`dr4pl` now imports the `rlang` package to allow for tidy evaluation
when using `dr4pl.data.frame`. This is not a major dependency as `dr4pl`
already imports `ggplot2` which imports `rlang`.

``` r
dr4pl(sample_data_3, dose = Dose, response = Response) 
#> Call:
#> dr4pl.data.frame(data = sample_data_3, dose = Dose, response = Response)
#> 
#> Coefficients:
#>       UpperLimit             IC50                 Slope        LowerLimit
#> 59858.0700194411 5.07949182697285    -0.611737121034004  980.820568608382
dr4pl(sample_data_3, dose = Dose/100000, response = Response) 
#> Call:
#> dr4pl.data.frame(data = sample_data_3, dose = Dose/1e+05, response = Response)
#> 
#> Coefficients:
#>       UpperLimit                 IC50                 Slope         LowerLimit
#> 56553.5168965252 0.000167802698462147    -0.880333688166774  -753.724715235523
```

Prior to `dr4pl (>= 2.0.0)`, users were required to specify a limit for
each parameter even if they wanted to only constrain one. Now the users
can supply a named numeric vector for ease of use.

``` r
dr4pl(sample_data_1, dose = Dose, response = Response, lowerl = c(theta_4 = 0)) #make lowerlimit positive 
#> Call:
#> dr4pl.data.frame(data = sample_data_1, dose = Dose, response = Response, 
#>     lowerl = c(theta_4 = 0))
#> 
#> Coefficients:
#>       UpperLimit             IC50                 Slope           LowerLimit
#> 112247.294407694 17.2973238818617    -0.384948923142149  0.00189050905753972
```

In addition, users will receive a more helpful error message if they
ever attempt to constrain the fit with `upperl` and `lowerl`.

``` r
dr4pl(Response~Dose, sample_data_4, lowerl = c(theta_4 = 0)) #make lowerlimit positive 
#> Error: Initial parameter values are not in the interior of the feasible region.
#> Estimated Parameters:
#>       UpperLimit      Log10(IC50)               Slope    LowerLimit
#> 31295.0216027415 4.67986724187866    -8.2357409544985       -57.218 
#> Failed Constraints:
#> theta_4 = - 57.218  >=  0
```

You may also pass a named numeric vector to `init.parm`, but you will
also receive a warning message as `init.parm` expects a `dr4pl_param`
object constructed from `dr4pl_theta()`. This warning will only appear
every 8 hours though.

``` r
dr4pl(Response~Dose, sample_data_4, init.parm = c(theta_4 = 0.1),lowerl = c(theta_4 = 0))
#> Warning: A numeric object is being coerced to a "dr4pl_param" object.`theta_2` is assumed to be in linear space. Please use `dr4pl_theta()` to construct the theta parameter.
#> This warning is displayed once every 8 hours.
#> Call:
#> dr4pl.formula(formula = Response ~ Dose, data = sample_data_4, 
#>     init.parm = c(theta_4 = 0.1), lowerl = c(theta_4 = 0))
#> 
#> Coefficients:
#>       UpperLimit             IC50                Slope         LowerLimit
#> 30599.6079444198 1020.20601140621    -8.03521868073897   3716.67380596223
```

`dr4pl_theta()` function allows the user to specify if the object is in
log10 space or not:

``` r
parm_lin <- dr4pl_theta(theta_2 = 2)
parm_log <- dr4pl_theta(theta_2 = 2, isLog10 = T)

parm_lin
#> UpperLimit   IC50    Slope   LowerLimit
#>         NA      2       NA           NA
parm_log
#> UpperLimit   Log10(IC50) Slope   LowerLimit
#>         NA             2    NA           NA

identical(parm_lin, LogToParm(parm_log))
#> [1] FALSE
```
