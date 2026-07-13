# Changelog

## deformetricar 1.0.0

A major modernisation: `deformetricar` is now a focused client for
Deformetrica (\>= 4.3), meeting current R package standards. This 1.0.0
release marks a clean break from the 10-year-old, L1-larva-specific
`ShootAndFlow3`/`sparseMatching3` era: the package is now template- and
species-agnostic. The old L1 symmetrisation helpers (`symmetrisel1()`,
`otherside()`, `apply.mirror.affine()`) are no longer exported — they
are rebuilt from scratch as a walkthrough in the L1 vignette, showing
how to compose a symmetrising warp from the general API.

### New features

- [`deformetrica_register()`](https://natverse.github.io/deformetricar/reference/deformetrica_register.md)
  and
  [`deformetrica_shoot()`](https://natverse.github.io/deformetricar/reference/deformetrica_shoot.md):
  fit and apply diffeomorphisms via the Deformetrica 4 `estimate` /
  `compute` CLI.
- [`deformetrica_register_multi()`](https://natverse.github.io/deformetricar/reference/deformetrica_register_multi.md):
  fit ONE diffeomorphism to a whole set of matched objects at once (the
  flyconnectome/deformetricaLR recipe), with an optional shared landmark
  regulariser. It accepts a per-object `data_sigma` and a per-object
  `object_kernel_width` (both scalar, named or positional vectors) so
  each object can be weighted and matched at its own scale in one fit —
  e.g. strong homologous neuropils over a weak outer hull.
- [`deformetrica_shoot()`](https://natverse.github.io/deformetricar/reference/deformetrica_shoot.md)
  guards the degenerate single-control-point case (which made
  Deformetrica’s torch kernel die with a cryptic `IndexError`) with a
  clear message.
- `deformetrica_shoot(flow = TRUE)` returns every timepoint of the
  geodesic flow, for animating a warp.
- Object-aware I/O: registration and shooting accept point matrices,
  `mesh3d` surfaces (SurfaceMesh + Current, faces preserved) and `nat`
  neurons (PolyLine, via the new
  [`write_neuron_vtk()`](https://natverse.github.io/deformetricar/reference/write_neuron_vtk.md));
  shooting returns the same class it was given.
- [`affine_prealign()`](https://natverse.github.io/deformetricar/reference/affine_prealign.md):
  centre + isotropic-scale + ICP pre-alignment returning an
  [`apply()`](https://rdrr.io/r/base/apply.html) closure that replays
  the identical transform on other objects.
- [`split_mesh_lr()`](https://natverse.github.io/deformetricar/reference/split_mesh_lr.md)
  and
  [`mirror_lr_split()`](https://natverse.github.io/deformetricar/reference/mirror_lr_split.md):
  split a surface into left/right halves — by a midline plane, or (even
  for an unpaired midline object) by a self-mirror registration that
  labels each point by its travel toward its mirror partner.
- [`ggplot_flow_gif()`](https://natverse.github.io/deformetricar/reference/ggplot_flow_gif.md):
  animate a geodesic flow as a GIF, each object coloured over a
  reference volume. It takes a `volume_col` (draw the brain hull as a
  soft coloured envelope, e.g. light pink), a set of static `targets`
  drawn beneath the flow in a transparent greyscale ramp (so the
  coloured moving structures can be read landing on the fixed targets
  they should match), and a `rotation_matrix` for the view. The general
  animation machinery now lives in
  [`nat.ggplot::ggneuron_gif()`](https://natverse.github.io/nat.ggplot/reference/ggneuron_gif.html)
  (\>= 1.1.0);
  [`ggplot_flow_gif()`](https://natverse.github.io/deformetricar/reference/ggplot_flow_gif.md)
  is a thin convenience wrapper around it, and the assembler falls back
  to `magick` when `gifski` (which needs a Rust toolchain) is
  unavailable.
- [`install_deformetrica()`](https://natverse.github.io/deformetricar/reference/install_deformetrica.md):
  one-call setup of the Deformetrica (\>= 4.3) CLI in a
  reticulate-managed conda env / virtualenv;
  [`find_deformetrica()`](https://natverse.github.io/deformetricar/reference/find_deformetrica.md)
  then resolves it automatically.
  [`transform_vtk()`](https://natverse.github.io/deformetricar/reference/transform_vtk.md)
  (renamed from the unexported `transform.vtk`, which R mistook for an
  S3 method) is now exported.
- [`find_deformetrica()`](https://natverse.github.io/deformetricar/reference/find_deformetrica.md)
  locates the executable on the `PATH`, via
  `options(deformetricar.exe=)`, in an
  [`install_deformetrica()`](https://natverse.github.io/deformetricar/reference/install_deformetrica.md)
  environment, or in a conda environment.
- Vignettes, each with GIFs: *mosquito-to-fly* (Aedes → Drosophila
  JRC2018F, whole-brain + central complex), *fafb-left-right* (FAFB
  brain mirror registration), and *l1-symmetrise-cognates* (symmetrise
  the L1 larval CNS from VFB CATMAID and find left-right cognates).

### Breaking changes

- `read.vtk()` /
  [`write.vtk()`](https://rdrr.io/pkg/nat/man/write.vtk.html) renamed to
  [`read_vtk()`](https://natverse.github.io/deformetricar/reference/read_vtk.md)
  /
  [`write_vtk()`](https://natverse.github.io/deformetricar/reference/write_vtk.md)
  so they no longer shadow
  [`nat::write.vtk()`](https://rdrr.io/pkg/nat/man/write.vtk.html) (an
  S3 generic). `nat`’s neuron VTK writes multi-point polylines, whereas
  Deformetrica 4 needs 2-point line segments, so the two cannot share a
  method — hence the distinct names.
- Removed the defunct Deformetrica 2.1 `shootflow()` / `ShootAndFlow3`
  path — use
  [`deformetrica_shoot()`](https://natverse.github.io/deformetricar/reference/deformetrica_shoot.md).
- Removed the L1 larval mirroring / cognate helpers (`symmetrisel1()`,
  `apply.mirror.affine()`, `otherside()`, `findcognate()` and friends)
  and their bundled transforms. That scientific workflow now lives in
  the *l1-symmetrise-cognates* vignette, driven by the live VFB dataset.

### Infrastructure

- `Authors@R`, `Encoding: UTF-8`, `URL`/`BugReports`, roxygen markdown,
  testthat edition 3, a `CITATION`, `NEWS.md`, pkgdown site + GitHub
  Pages / R-CMD-check workflows. Hard dependencies trimmed to `nat` and
  `Morpho`; `Rvcg`, `nat.ggplot`, `gifski` and the atlas/data packages
  are `Suggests`.
