# deformetricar modernisation (branch `dev`)

Goal: make `deformetricar` work with the **modern Deformetrica (>= 4.3)** library and meet current R
package standards. Informed directly by the BANC alignment warp work in `bancpipeline/deform`
(`deformetrica-align-warp.R`, `deformetrica-align-shoot.R`, `deform_shoot.sh`), which drove a full
Deformetrica 4.3.0 fit + shooting pipeline on O2.

## What Deformetrica 4.x changed (vs the 2.1 API this package was written for)

- **No more `ShootAndFlow3` / `sparseMatchings3` C++ binaries.** There is a single `deformetrica` CLI
  with subcommands `estimate` and `compute`.
- **Shooting is `deformetrica compute` + a `<model-type>Shooting` model.xml** referencing
  `<initial-control-points>` and `<initial-momenta>` (the estimator's
  `..._EstimatedParameters__ControlPoints.txt` / `...Momenta.txt`), plus `<deformation-parameters>`
  (kernel-width, kernel-type `torch`, number-of-timepoints). Output is
  `Shooting__GeodesicFlow__<obj>__tp_<i>__age_<t>.vtk`; take the max-tp file.
- **Kernel/device**: `<kernel-type>torch</kernel-type>` + `<kernel-device>auto|cpu|cuda</kernel-device>`
  in the model (NOT `use-cuda` in optimization_parameters.xml — Deformetrica 4.3 warns "Unknown entry"
  on that key and ignores it).
- **`Momenta.txt` has a leading `n_subjects n_cp dim` header row** that `ControlPoints.txt` does not —
  Deformetrica writes and reads it, but any external parser must account for it.
- **Object types** are `Landmark`, `PolyLine`, `SurfaceMesh` (not the old
  `PointCloud`/`OrientedPolyLine`/...). Attachment types: `Landmark` (ordered L2), `Varifold`, `Current`.

`R/deformetrica4.R` (new, this branch) wraps that path: `find_deformetrica()` +
`deformetrica_shoot(x, control_points, momenta, kernel_width, ...)`.

## Task list

### Done on `dev`
- [x] `R/deformetrica4.R`: `find_deformetrica()`, `deformetrica_shoot()` (Deformetrica 4 shooting apply),
      ported from the validated `bancpipeline/deform/deformetrica-align-shoot.R`.
- [x] DESCRIPTION: `Authors@R`, `Encoding: UTF-8`, `URL`/`BugReports`, `Roxygen: markdown`,
      `RoxygenNote: 7.3.2`, testthat edition 3, `Depends: R (>= 4.1)`.
- [x] `.Rbuildignore`: exclude the 37 MB `reg_neurons` blob, `fullmirror_nana.rds`, `R/playground.R`, dev docs.

### Done on `dev` (cont.)
- [x] `deformetrica_register()` (Deformetrica 4 estimate/fit path) + mesh support: `deformetrica_register`
      and `deformetrica_shoot` accept `mesh3d` (SurfaceMesh + Current, faces written to VTK; shoot returns a
      warped mesh/neuron/dotprops of the same class).
- [x] `tests/testthat/test-deformetrica4.R`: VTK point + mesh round-trip (nat/Rvcg data), find_deformetrica
      error path, register->shoot recovers a known translation (skip when the binary is absent).
- [x] `vignettes/mosquito-to-fly.Rmd`: Aedes mosquito brain -> JRC2018U fly, affine init (Morpho) then
      Deformetrica surface warp, homologous-neuropil validation + per-neuropil fallback, GIF frames from the
      geodesic flow. Heavy chunks `eval=FALSE` (need IBdb + GPU).
- [x] DESCRIPTION: VignetteBuilder + vignette Suggests (knitr, rmarkdown, nat.flybrains, insectbrainr, gifski).

### Concepts ported from flyconnectome/deformetricaLR (2026-07-10)
That repo (FAFB L-R bridging registration, Deformetrica 3.0.0.beta) contributed:
- [x] **Multi-object registration** — `deformetrica_register_multi()`: ONE diffeomorphism fit to a whole
      SET of matched objects at once (its `model.xml` had one `<object>` per cognate neuron tract). This is
      how a real bridging registration is built, vs a single source→target shape.
- [x] **Neuron tracts as `NonOrientedPolyLine`** — `write_neuron_vtk()` (VTK POLYDATA LINES) + object
      auto-classification (mesh→SurfaceMesh/Current, neuron→NonOrientedPolyLine/Varifold, points→Landmark).
- [x] **Landmark object as a global regulariser** — the `inc_tracts_lmarks` vs `inc_tracts_no_landmark`
      ablation: `deformetrica_register_multi(landmarks=list(source=,target=))` adds a shared Landmark object
      alongside the tracts. Per-object `data_sigma` supported.
- [x] **mirror → affine → deformable** L-R pipeline → `vignettes/fafb-left-right.Rmd` (with its own GIF).
- [ ] Still worth lifting: their `param_gen/` parameter-file generators (superseded by our XML writers) and
      the neuropil-landmark feathers (`regLandmarks/FCWB_unmirrored`) as ready-made anchor sets.

### Next (in rough priority)
- [ ] **insectbrainr may need a PR**: it is older and its Insect Brain Database mesh queries may be broken;
      reconcile the vignette's `insectbraindb_*` calls with the real API (preferred mosquito: **Aedes aegypti**),
      and if the mesh fetch is broken, fix it upstream (clone natverse/insectbrainr, dev branch).
- [ ] Regenerate `NAMESPACE` + `man/` with roxygen2 7.x (`devtools::document()`); export the new functions.
- [ ] Add a Deformetrica 4 **estimate** wrapper (`deformetrica_register()`) mirroring the
      `bancpipeline` fit: write ordered Landmark VTKs + Registration/DeterministicAtlas model, run
      `deformetrica estimate`, return control points + momenta. Fold in the object-id derivation fix
      (avoid `gsub`-ing a suffix that recurs in the object name).
- [ ] Supersede `shootflow()` (2.1) via lifecycle: keep it, mark `superseded`, point to
      `deformetrica_shoot()`. Do NOT delete — some users may still have Deformetrica 2.1.
- [ ] Relocate/shrink example data: move `reg_neurons` (37 MB) + `fullmirror_nana.rds` out of the repo
      root; ship a *small* `inst/extdata/` example registration + point set for tests/examples.
- [ ] Dependency trim: guard `Morpho`/`Rvcg`/`rgl`/`catmaid`/`nat.nblast` uses with
      `requireNamespace()` and move to `Suggests` (only `nat`, `xml2` are core).
- [ ] Tests (testthat 3e): a `deformetrica_shoot()` round-trip on a tiny fixture registration; skip on
      CI with `skip_if(is.null(find_deformetrica-safe))` when the binary is absent.
- [ ] `R CMD check` clean (0 errors/warnings); add `README.Rmd` + pkgdown.

## References
- Deformetrica: https://gitlab.com/icm-institute/aramislab/deformetrica (installed: 4.3.0)
- Working 4.x usage: `bancpipeline/deform/{deformetrica-align-warp,deformetrica-align-shoot}.R`,
  `deform/{deform_align_warp,deform_shoot}.sh`
- R package standards: posit-dev/skills `r-lib/` (testing-r-packages, cran-extrachecks, lifecycle).
