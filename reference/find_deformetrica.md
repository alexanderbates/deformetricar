# Locate the Deformetrica (\>= 4.3) command-line executable

Searches, in order: `options(deformetricar.exe=)`, the `deformetrica` on
the `PATH`, the reticulate-managed environments
[`install_deformetrica()`](https://natverse.github.io/deformetricar/reference/install_deformetrica.md)
creates (`"deformetrica"` then `"r-reticulate"`), then the common conda
location `~/.conda/envs/deformetrica/bin/deformetrica`.

## Usage

``` r
find_deformetrica(deformetrica = getOption("deformetricar.exe"))
```

## Arguments

- deformetrica:

  Optional explicit path to (or name of) the executable.

## Value

The resolved executable path.

## See also

[`install_deformetrica()`](https://natverse.github.io/deformetricar/reference/install_deformetrica.md)
to set one up from R.
