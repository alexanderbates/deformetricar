# Fit a Deformetrica diffeomorphism between two point sets (registration)

The Deformetrica 4 estimate path: fit a diffeomorphism that deforms
`source` (the moving template) onto `target` (the fixed subject). With
ordered `source` \<-\> `target` row correspondences a point-to-point
Landmark attachment is used. Returns the estimated control points +
momenta, which
[`deformetrica_shoot()`](https://natverse.github.io/deformetricar/reference/deformetrica_shoot.md)
then applies to arbitrary new points.

## Usage

``` r
deformetrica_register(
  source,
  target,
  kernel_width,
  timepoints = 10L,
  object_type = "Landmark",
  attachment_type = c("Landmark", "Varifold", "Current"),
  noise_std = 0.01,
  max_iterations = 150L,
  device = c("auto", "cpu", "cuda"),
  deformetrica = NULL,
  workdir = tempfile("dfca_reg_"),
  verbose = FALSE
)
```

## Arguments

- source, target:

  N x 3 matrices of corresponding points (row *i* of `source` matches
  row *i* of `target`). Must have the same number of rows.

- kernel_width:

  Deformation kernel width (larger = stiffer / more global).

- timepoints:

  Geodesic integration steps.

- object_type:

  Deformetrica `deformable-object-type` (default "Landmark").

- attachment_type:

  "Landmark" (ordered L2, default), "Varifold" or "Current".

- noise_std:

  Data-attachment noise standard deviation.

- max_iterations:

  Optimiser iterations.

- device:

  `kernel-device`: "auto", "cpu" or "cuda".

- deformetrica:

  Optional path to the executable (auto-detected otherwise).

- workdir:

  Working directory for the run (a fresh tempdir by default).

- verbose:

  Show `deformetrica` output.

## Value

A `deformetricareg` object (a portable handle on the fitted
diffeomorphism, storing the control points + momenta inline). Apply it
with `nat::xform(x, reg)` or
[`deformetrica_shoot()`](https://natverse.github.io/deformetricar/reference/deformetrica_shoot.md);
it survives [`saveRDS()`](https://rdrr.io/r/base/readRDS.html). Its
elements `control_points`, `momenta`, `kernel_width`, `timepoints` and
`output_dir` are still accessible.

## See also

[`xformpoints.deformetricareg()`](https://natverse.github.io/deformetricar/reference/xformpoints.deformetricareg.md)
/
[`deformetrica_shoot()`](https://natverse.github.io/deformetricar/reference/deformetrica_shoot.md)
to apply the fit.
