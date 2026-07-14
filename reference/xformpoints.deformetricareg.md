# Apply a Deformetrica registration as a natverse transform

[`deformetrica_register()`](https://natverse.github.io/deformetricar/reference/deformetrica_register.md)
and
[`deformetrica_register_multi()`](https://natverse.github.io/deformetricar/reference/deformetrica_register_multi.md)
return a `deformetricareg` — a portable handle on the fitted
diffeomorphism, and a
[`nat::xformpoints()`](https://rdrr.io/pkg/nat/man/xformpoints.html)
method. That makes a fitted registration a first-class natverse
transform: `nat::xform(x, reg)` warps any points / neuron / neuronlist /
mesh through it, and it composes with other registrations (e.g. in a
`nat.templatebrains` reglist). The control points and momenta are stored
*inline*, so a `deformetricareg` survives
[`saveRDS()`](https://rdrr.io/r/base/readRDS.html) and being moved to
another machine.

## Usage

``` r
# S3 method for class 'deformetricareg'
xformpoints(reg, points, ...)
```

## Arguments

- reg:

  A `deformetricareg`.

- points:

  An N x 3 matrix of coordinates.

- ...:

  Passed to
  [`deformetrica_shoot()`](https://natverse.github.io/deformetricar/reference/deformetrica_shoot.md)
  (e.g. `object_type`, `device`).

## Value

An N x 3 matrix of deformed coordinates.

## See also

[`deformetrica_shoot()`](https://natverse.github.io/deformetricar/reference/deformetrica_shoot.md),
which does the actual geodesic shooting.
