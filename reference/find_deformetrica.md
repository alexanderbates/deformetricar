# Locate the Deformetrica (\>= 4.3) command-line executable

Searches, in order: `options(deformetricar.exe=)`, the `deformetrica` on
the `PATH`, then the common conda location
`~/.conda/envs/deformetrica/bin/deformetrica`.

## Usage

``` r
find_deformetrica(deformetrica = getOption("deformetricar.exe"))
```

## Arguments

- deformetrica:

  Optional explicit path to (or name of) the executable.

## Value

The resolved executable path.
