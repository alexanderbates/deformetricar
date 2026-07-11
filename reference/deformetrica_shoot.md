# Apply a fitted Deformetrica diffeomorphism to 3D points (geodesic shooting)

Flows arbitrary points through the diffeomorphism defined by a fitted
registration's control points and initial momenta, by writing a
`<model-type>Shooting` model and calling `deformetrica compute`.

## Usage

``` r
deformetrica_shoot(
  x,
  control_points,
  momenta,
  kernel_width,
  timepoints = 10L,
  object_type = "Landmark",
  device = c("auto", "cpu", "cuda"),
  deformetrica = NULL,
  workdir = tempfile("dfca_shoot_"),
  flow = FALSE,
  verbose = FALSE
)
```

## Arguments

- x:

  An N x 3 matrix of coordinates, or anything
  [`nat::xyzmatrix()`](https://rdrr.io/pkg/nat/man/xyzmatrix.html)
  accepts.

- control_points, momenta:

  Path to the estimated Deformetrica `...ControlPoints.txt` /
  `...Momenta.txt`, or a numeric matrix. Passing the files Deformetrica
  itself wrote is strongly preferred: it writes momenta with a leading
  `n_subjects n_cp dim` header that its own reader expects, so a bare
  matrix is best-effort only.

- kernel_width:

  Deformation kernel width, in the coordinate units. Larger = stiffer /
  more global; must match the width the registration was fit with.

- timepoints:

  Number of geodesic integration steps; must match the fit.

- object_type:

  Deformetrica `deformable-object-type` (default "Landmark").

- device:

  `kernel-device`: "auto", "cpu" or "cuda".

- deformetrica:

  Optional path to the executable (auto-detected otherwise).

- workdir:

  Working directory for the run (a fresh tempdir by default).

- flow:

  If `TRUE`, return the *whole* geodesic flow (every timepoint) as a
  list of deformed objects (one per timepoint, named `tp_00`..`tp_NN`)
  rather than just the final shape. Useful for animating a warp.

- verbose:

  Show `deformetrica` output.

## Value

By default an N x 3 matrix of deformed coordinates (or a warped object
of the input class), in the input units and row order. With
`flow = TRUE`, a list of such objects, one per geodesic timepoint.
