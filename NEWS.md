---
output: md_document
---

# dr4pl 1.1.6

dr4pl function now allows 'data' argument when using 'dose' and 'response' arguments. See dr4pl examples.

# dr4pl 1.1.7

print.summary.dr4pl function no long provides t-statistics and p-values but now prints 95% confidence intervals instead. Consolidated Vignettes.

# dr4pl 1.1.8

Following proper semantic versioning.

# dr4pl 1.1.9

Bug fix update.
  dr4pl.data.frame will no longer fail when using tibbles as an input.
  "L-BFGS-B" is now a usable argument to the parameter 'method.optim' as originally intended.

General changes and updates 
  default value of method.robust is now set to "squared" to be more conistent with other values given.
  Added upperl and lowerl arguments to dr4pl call. User can now constrain optimization of theta parameters.
    This may require user to set init.parm as the estimated theta parameters will not always be within
    the fesiable region specified by upperl and lowerl.