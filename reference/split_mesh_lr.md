# Split a surface mesh into left and right halves by a midline plane

Faces are assigned to a side by their centroid position along `axis`
relative to `mid`. For a cleaner split of a shape whose symmetry surface
is curved rather than a flat plane, see
[`mirror_lr_split()`](https://natverse.github.io/deformetricar/reference/mirror_lr_split.md).

## Usage

``` r
split_mesh_lr(x, mid = NULL, axis = c("X", "Y", "Z"), left = c("low", "high"))
```

## Arguments

- x:

  A `mesh3d`.

- mid:

  The midline coordinate on `axis`; defaults to the midpoint of the
  mesh's extent on that axis.

- axis:

  Split axis, "X" (default), "Y" or "Z".

- left:

  A face on the side with the *smaller* (`"low"`, default) or *larger*
  (`"high"`) `axis` coordinate is called the left half.

## Value

A list with `L` and `R` sub-meshes and the `mid` used.
