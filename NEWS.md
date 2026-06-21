# CHANGES IN VERSION 1.0.0
* Initial release of the `scDiagnostics` package.

# CHANGES IN VERSION 1.4.0
* Add functionality to process SCE objects for PCA computation via the new `processPCA()` function. 
* Add functionality to downsample SCE objects for diagnostic functions.
* Add new diagnostic functions (`calculateTopLoadingGeneShifts()`, `compareMarkers()` and `calculateMMDPValue()`).
* Add graph integration diagnostic function in replacement of nearest neighbor diagnostic, `calculateGraphIntegration()`.
* Improve `regressPC()` function and plot method, which can now also regress against cell types and batches.
* Improve normalization for `plotMarkerExpression()` diagnostic function.
* Improve user control and options for plot methods.
* Update vignettes to reflect all new changes.

# CHANGES IN VERSION 1.6.0
* Renamed gene shift function for consistency (previously `calculateTopLoadingGeneShifts()`)
* Added gene specification parameter to `calculateGeneShifts()`
* Improved `calculateGeneShifts()` function, plot method, and color scheme

# CHANGES IN VERSION 1.8.0
* Added `calculateReconstructionError()` to detect out-of-distribution anomalies using cell-type-specific PCA reconstruction errors.
* Added `plot.calculateReconstructionErrorObject()` featuring robust visualization options (violin, boxplot, ridge, and ComplexHeatmap).
* Upgraded `detectAnomaly()` to resolve the curse of dimensionality by allowing Isolation Forests to run on the union of query and reference Highly Variable Genes (via `n_hvgs`) when `pc_subset = NULL`.
* Improved anomaly detection by switching default thresholding to a dynamic, data-driven Median Absolute Deviation method (`threshold_method = "MAD"`, `mad_multiplier = 2`) across relevant functions.
