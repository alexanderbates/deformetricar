<!-- badges: start -->
[![natverse](https://img.shields.io/badge/natverse-Part%20of%20the%20natverse-a241b6)](https://natverse.github.io)
[![Docs](https://img.shields.io/badge/docs-experimental-orange.svg)](https://alexanderbates.github.io/deformetricar/reference/)
[![R-CMD-check](https://github.com/alexanderbates/deformetricar/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/alexanderbates/deformetricar/actions/workflows/R-CMD-check.yaml)
[![Lifecycle: experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
<!-- badges: end -->

# deformetricar

**deformetricar** is an R client for the [Deformetrica](https://www.deformetrica.org)
shape-registration toolkit (`>= 4.3`). It lets you **fit** diffeomorphisms between
3D shapes — point clouds, neuron backbones and surface meshes — and **apply** them
(geodesic shooting) to arbitrary objects, with sensible defaults for registering
fly connectome neurons. It also reads and writes VTK object formats.

Deformetrica 4 replaced the old `ShootAndFlow3` C++ binaries with a single
`deformetrica` CLI (`estimate` / `compute`); `deformetricar` wraps that modern
interface.

<p align="center">
<img src="man/figures/mosquito_to_fly.gif" width="49%" alt="Aedes mosquito brain warping onto the Drosophila brain"/>
<img src="man/figures/mosquito_cx_to_fly.gif" width="49%" alt="Mosquito central complex warping onto the fly central complex"/>
</p>

*Left: the whole **Aedes aegypti** brain surface warped onto the **Drosophila**
JRC2018F template. Right: just the central complex — the mosquito's CBU/CBL onto the
fly's fan-shaped and ellipsoid bodies. Both are produced end-to-end by the
[mosquito-to-fly vignette](https://alexanderbates.github.io/deformetricar/articles/mosquito-to-fly.html).*

## Installation

```r
# install.packages("remotes")
remotes::install_github("alexanderbates/deformetricar")
```

You also need a working [Deformetrica (>= 4.3)](https://gitlab.com/icm-institute/aramislab/deformetrica)
install (e.g. in a conda environment). `find_deformetrica()` locates it on the
`PATH`, at `options(deformetricar.exe=)`, or in `~/.conda/envs/deformetrica/`.

## Quick start

```r
library(deformetricar)

# Fit a diffeomorphism between two corresponding point sets ...
fit <- deformetrica_register(source, target, kernel_width = 20)

# ... then apply it to any object (returns the same class you pass in)
warped <- deformetrica_shoot(new_points, fit$control_points, fit$momenta,
                             kernel_width = fit$kernel_width)

# Register a whole SET of matched objects at once (e.g. cognate neuron tracts)
fit <- deformetrica_register_multi(sources, targets, kernel_width = 20,
                                   landmarks = list(source = lm_s, target = lm_t))
```

## Articles

- **[Warping a mosquito brain onto the fly](https://alexanderbates.github.io/deformetricar/articles/mosquito-to-fly.html)** —
  register the *Aedes aegypti* brain onto *Drosophila* JRC2018U (affine init + surface warp), with a GIF.
- **[A FAFB left-right bridging registration](https://alexanderbates.github.io/deformetricar/articles/fafb-left-right.html)** —
  one diffeomorphism from matched cognate neuron tracts + neuropil landmarks.

## Citation

`citation("deformetricar")` cites the package, the
[natverse](https://doi.org/10.7554/eLife.53350) (Bates et al. 2020, *eLife*), and
the Deformetrica software (Bône et al. 2018; Durrleman et al. 2014). Please cite
all three when you use `deformetricar`.

## Acknowledgements

Deformetrica is developed by the [Aramis Lab](https://www.deformetrica.org) at the
Paris Brain Institute (ICM). `deformetricar` merely wraps it. Part of the
[natverse](https://natverse.github.io).
