# Affine pre-alignment of one object onto another

A diffeomorphism should not have to absorb a gross change of pose and
scale; that is this step's job. Both objects are centred, `source` is
isotropically scaled to `target`'s bounding-box diagonal, and a rigid
(or affine) ICP refines the orientation. The returned
[`apply()`](https://rdrr.io/r/base/apply.html) closure replays the
*same* transform on any other object, so a whole set of sub-objects
(e.g. neuropils) can be carried into the aligned frame consistently.

## Usage

``` r
affine_prealign(
  source,
  target,
  type = c("rigid", "similarity", "affine"),
  iterations = 30L
)
```

## Arguments

- source, target:

  Objects accepted by
  [`nat::xyzmatrix()`](https://rdrr.io/pkg/nat/man/xyzmatrix.html)
  (`mesh3d`, neuron, dotprops or an N x 3 matrix). `source` moves onto
  `target`.

- type:

  ICP type passed to
  [`Morpho::icpmat()`](https://rdrr.io/pkg/Morpho/man/icpmat.html):
  "rigid", "similarity" or "affine".

- iterations:

  ICP iterations.

## Value

A list with `aligned` (the transformed `source`), `apply` (a function
mapping any object through the identical transform) and the components
`scale`, `centre_source`, `centre_target` and `icp` (the 4x4 ICP
matrix).

## See also

[`deformetrica_register()`](https://natverse.github.io/deformetricar/reference/deformetrica_register.md)
for the non-rigid step that follows.
