# Assign left/right to an (even unpaired, midline) object by self-mirror warping

A mini use-case of the package. An object that straddles the midline (a
fly ellipsoid body, protocerebral bridge, ...) has no left/right label,
yet every point still belongs to a hemisphere. We recover that the way a
mirroring registration does (cf. `bancr::banc_lr`): reflect the object
across a nominal midline plane and, optionally, `refine` the reflection
with a Deformetrica warp of the object onto its own reflection so the
true (possibly curved) symmetry surface is found. Each point is then
labelled by *how far and in which direction it must travel to reach its
mirror partner* - points on the symmetry surface barely move; points to
one side travel toward the other.

## Usage

``` r
mirror_lr_split(
  x,
  axis = c("X", "Y", "Z"),
  mid = NULL,
  refine = FALSE,
  band = 0,
  positive_travel = c("L", "R"),
  ...
)
```

## Arguments

- x:

  A `mesh3d`.

- axis:

  Mirror axis, "X" (default), "Y" or "Z".

- mid:

  Midline coordinate on `axis` (defaults to the extent midpoint).

- refine:

  If `TRUE`, refine the flat reflection with a Deformetrica warp of `x`
  onto its reflection (finds a curved symmetry surface). Needs a working
  Deformetrica install; passes `...` (e.g. `kernel_width`, `device`)
  through to
  [`deformetrica_register()`](https://natverse.github.io/deformetricar/reference/deformetrica_register.md)/[`deformetrica_shoot()`](https://natverse.github.io/deformetricar/reference/deformetrica_shoot.md).

- band:

  Fraction of the axis extent within which \|travel\| is treated as
  midline (`"M"`) rather than L/R. Default 0 (every point is L or R).

- positive_travel:

  Which hemisphere a point that travels in the `+axis` direction belongs
  to ("L" default, or "R").

- ...:

  Passed to the Deformetrica calls when `refine = TRUE`.

## Value

A list with `L`/`R` sub-meshes, a per-vertex `side` factor and the
signed `travel` along `axis`.

## See also

[`split_mesh_lr()`](https://natverse.github.io/deformetricar/reference/split_mesh_lr.md)
for the plain-plane split.
