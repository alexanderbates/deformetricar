# deformetricar 0.2.0.9000 (development)

Modernisation for Deformetrica (>= 4.3) and current R package standards.

## New features

* `deformetrica_register()` and `deformetrica_shoot()`: fit and apply diffeomorphisms
  via the Deformetrica 4 `estimate` / `compute` CLI, replacing the defunct 2.1
  `ShootAndFlow3` path in `shootflow()`.
* `deformetrica_register_multi()`: fit ONE diffeomorphism to a whole set of matched
  objects at once (the flyconnectome/deformetricaLR recipe), with an optional shared
  landmark regulariser.
* Object-aware I/O: registration and shooting accept point matrices, `mesh3d`
  surfaces (SurfaceMesh + Current, faces preserved) and `nat` neurons
  (NonOrientedPolyLine, via the new `write_neuron_vtk()`); shooting returns the
  same class it was given.
* `find_deformetrica()` locates the executable on the `PATH`, via
  `options(deformetricar.exe=)`, or in a conda environment.
* Vignettes: *mosquito-to-fly* (Aedes → Drosophila JRC2018U) and *fafb-left-right*
  (matched-tract L-R bridging registration), each with a GIF.

## Infrastructure

* `Authors@R`, `Encoding: UTF-8`, `URL`/`BugReports`, roxygen markdown,
  `RoxygenNote` 7.3.2, testthat edition 3, and a `CITATION`.
* pkgdown site + GitHub Pages / R-CMD-check workflows.
* The legacy Deformetrica 2.1 `shootflow()` is retained but superseded by
  `deformetrica_shoot()`.
