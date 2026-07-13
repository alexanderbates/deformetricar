# deformetricar: an R client for the Deformetrica shape-registration toolkit

**deformetricar** wraps the Deformetrica (\>= 4.3) command-line tools
for statistical shape analysis: it fits diffeomorphisms between point
clouds, neuron backbones and surface meshes, and applies them by
geodesic shooting. It also reads and writes VTK object formats.

## Package options

- `deformetricar.exe`:

  Path to (or name of) the `deformetrica` executable; consulted by
  [`find_deformetrica`](https://natverse.github.io/deformetricar/reference/find_deformetrica.md).
  Set it in your [`Rprofile`](https://rdrr.io/r/base/Startup.html) if
  Deformetrica is not on the `PATH` or in a conda environment.

## References

<https://www.deformetrica.org>

## See also

[`find_deformetrica`](https://natverse.github.io/deformetricar/reference/find_deformetrica.md),
[`deformetrica_register`](https://natverse.github.io/deformetricar/reference/deformetrica_register.md),
[`deformetrica_shoot`](https://natverse.github.io/deformetricar/reference/deformetrica_shoot.md)

## Author

**Maintainer**: Alexander Bates <alexander.shakeel.bates@gmail.com>

Authors:

- Gregory Jefferis <jefferis@gmail.com>
  ([ORCID](https://orcid.org/0000-0002-0587-9355))
