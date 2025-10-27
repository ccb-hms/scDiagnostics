# CHANGES IN VERSION 1.0.0
* Initial release of the `scDiagnostics` package.

# CHANGES IN VERSION 1.4.0
* Add functionality to process SCE objects for PCA computation via the new `processPCA` function. 
* Add functionality to downsample SCE objects for diagnostic functions.
* Add new diagnostic functions (`calculateTopLoadingGeneShifts`, `compareMarkers` and `calculateMMDPValue`).
* Add graph integration diagnostic function in replacement of nearest neighbor diagnostic, `calculateGraphIntegration`.
* Improve `regressPC` function and plot method, which can now also regress against cell types and batches.
* Improve normalization for `plotMarkerExpression` diagnostic function.
* Improve user control and options for plot methods.
* Update vignettes to reflect all new changes.
