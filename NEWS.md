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

