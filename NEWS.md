# spareg 1.1.2

* Change plot methods to display original scale for the threshold parameter when plot_along = "nu".
* Remove duplicated warning sign in parallel case.
* Make "name" argument optional in the constructor functions
* In addition to fast fitting in spar.cv where screening indicators and random projections are fixed before the cross-validation,
we allow for fast fitting where only the random projections are fixed in advance
and also allow for sampling screening and random projections in each fold of the  (full) cross-validation.
* only allow "val_measure" or "val_numactive" argument in plot.spar.cv. 

# spareg 1.1.1

* Improved man pages.

# spareg 1.1.0

* Enhance `print` method with type of measure being used for choosing $M$ and $\nu$
* Allow for cases where $p < n$.
* Allow specification of possible beta values in `simulate_spareg_data()`.
* Changed robust example in vignette to using a `poisson()` link.
* Reproducibility of parallel case is introduced using doRNG.
* Added extractor function `extract_spareg()` to extract fitted values,
residuals, and coefficients from a fitted model.
* Added extractor function `get_model()` to extract the best or the 1se model.
* Added extractor function `get_measure()` to extract the table of (cross-)
validation measures for the grid of `(nu, nummods)`.
* Added class `coefspar` with `print` and `summary` methods to improve the usability of
the `coef` method.
* Added extractor functions `get_intercept()` and `get_coef()`
to extract the intercept and coefficients (non-standardized) from the 
`coefspar` objects.
* Added `avg_type` i.e., type of averaging the marginal models to `spar()`. This
argument is used in computing the validation measure.
* Removed `coef` argument from `predict` method.
* Added `aggregate = c("mean", "median")` argument for `predict` method.

# spareg 1.0.0

* Initial CRAN submission.

