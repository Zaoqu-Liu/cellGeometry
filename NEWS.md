News
=====

# cellGeometry 0.5.8
###### 03/02/2026
* **BUGFIX**: Fixed critical bug in `cellMarkers()` where `min_cells` filtering
was non-functional due to incorrect operator precedence.
* Improved numerical stability across all core functions:
  - Protected against division by zero in `dotprod()`, `equalweight()`, `scaleSphere()`
  - Added input clamping for `acos()` to prevent NaN from floating-point errors
  - Added singular matrix detection in `deconv()`
  - Fixed NA handling in outlier detection during multi-pass deconvolution
  - Added robust error handling in `rstudent()` and `cooks.distance()` methods
* Fixed edge case in `balance_cores()` when `cores = 1`.
* Unified code style: consistent spacing around operators throughout codebase.
* Added comprehensive test suite with 69+ tests covering:
  - Core mathematical functions
  - Scientific accuracy validation
  - Edge cases and numerical stability
  - Integration tests for full deconvolution workflow
* Updated maintainer and repository information.

# cellGeometry 0.5.7
###### 11/12/2025
* Change `log` argument in `deconvolute()` to `logged_bulk`. NB. this is a 
change of logic.

# cellGeometry 0.5.6
###### 07/10/2025
* Add `residuals.deconv()` to allow recalculation of full residuals matrix.

# cellGeometry 0.5.5
###### 10/09/2025
* Fix bug in `mergeMarkers()` if cellMarkers object has no cell grouping vector.

# cellGeometry 0.5.4
###### 09/09/2025
* Fix CRAN checks
* Switch from using `cat()` to `message()`

# cellGeometry 0.5.3
###### 26/08/2025
* Massive speed up (4-5x) in compensation optimisation

# cellGeometry 0.5.2
###### 30/07/2025

* Expanded to 3 methods for identifying outlier genes (var.e, Cook's distance,
Studentized residuals)
* Improved (nested) parallelisation of `tune_deconv()`
* Rewrite weights code (faster)

# cellGeometry 0.5.0
###### 14/07/2025

* Calculation of SE
* Detection and multipass removal of problematic genes with extreme residuals

# cellGeometry 0.2.0
###### 24/01/2025

* This is the initial build of cellGeometry
