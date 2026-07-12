# deformetricar 0.2.0.9000 (development)

A major modernisation: `deformetricar` is now a focused client for Deformetrica
(>= 4.3), meeting current R package standards.

## New features

* `deformetrica_register()` and `deformetrica_shoot()`: fit and apply diffeomorphisms
  via the Deformetrica 4 `estimate` / `compute` CLI.
* `deformetrica_register_multi()`: fit ONE diffeomorphism to a whole set of matched
  objects at once (the flyconnectome/deformetricaLR recipe), with an optional shared
  landmark regulariser. It now accepts a per-object `data_sigma` (scalar, named or
  positional vector) so objects can be weighted differently in one fit — e.g. strong
  homologous neuropils over a weak outer hull.
* `deformetrica_shoot(flow = TRUE)` returns every timepoint of the geodesic flow, for
  animating a warp.
* Object-aware I/O: registration and shooting accept point matrices, `mesh3d`
  surfaces (SurfaceMesh + Current, faces preserved) and `nat` neurons (PolyLine, via
  the new `write_neuron_vtk()`); shooting returns the same class it was given.
* `affine_prealign()`: centre + isotropic-scale + ICP pre-alignment returning an
  `apply()` closure that replays the identical transform on other objects.
* `split_mesh_lr()` and `mirror_lr_split()`: split a surface into left/right halves —
  by a midline plane, or (even for an unpaired midline object) by a self-mirror
  registration that labels each point by its travel toward its mirror partner.
* `ggplot_flow_gif()`: animate a geodesic flow as a `nat.ggplot` GIF, each object
  coloured over a translucent reference volume. It now also takes a `volume_col` (draw
  the brain hull as a soft coloured envelope, e.g. light pink) and a set of static
  `targets` drawn beneath the flow in a transparent greyscale ramp — so the coloured
  moving structures can be read landing on the fixed targets they should match.
* `install_deformetrica()`: one-call setup of the Deformetrica (>= 4.3) CLI in a
  reticulate-managed conda env / virtualenv; `find_deformetrica()` then resolves it
  automatically. `transform_vtk()` (renamed from the unexported `transform.vtk`, which
  R mistook for an S3 method) is now exported.
* `find_deformetrica()` locates the executable on the `PATH`, via
  `options(deformetricar.exe=)`, in an `install_deformetrica()` environment, or in a
  conda environment.
* Vignettes, each with GIFs: *mosquito-to-fly* (Aedes → Drosophila JRC2018F,
  whole-brain + central complex), *fafb-left-right* (FAFB brain mirror registration),
  and *l1-symmetrise-cognates* (symmetrise the L1 larval CNS from VFB CATMAID and find
  left-right cognates).

## Breaking changes

* `read.vtk()` / `write.vtk()` renamed to `read_vtk()` / `write_vtk()` so they no
  longer shadow `nat::write.vtk()` (an S3 generic). `nat`'s neuron VTK writes
  multi-point polylines, whereas Deformetrica 4 needs 2-point line segments, so the
  two cannot share a method — hence the distinct names.
* Removed the defunct Deformetrica 2.1 `shootflow()` / `ShootAndFlow3` path — use
  `deformetrica_shoot()`.
* Removed the L1 larval mirroring / cognate helpers (`symmetrisel1()`,
  `apply.mirror.affine()`, `otherside()`, `findcognate()` and friends) and their
  bundled transforms. That scientific workflow now lives in the *l1-symmetrise-cognates*
  vignette, driven by the live VFB dataset.

## Infrastructure

* `Authors@R`, `Encoding: UTF-8`, `URL`/`BugReports`, roxygen markdown, testthat
  edition 3, a `CITATION`, `NEWS.md`, pkgdown site + GitHub Pages / R-CMD-check
  workflows. Hard dependencies trimmed to `nat` and `Morpho`; `Rvcg`, `nat.ggplot`,
  `gifski` and the atlas/data packages are `Suggests`.
