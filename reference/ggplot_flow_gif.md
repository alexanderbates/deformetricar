# Animate a Deformetrica geodesic flow as a nat.ggplot GIF

Renders a set of warping objects across the timepoints of a geodesic
flow (as returned by
[`deformetrica_shoot()`](https://alexanderbates.github.io/deformetricar/reference/deformetrica_shoot.md)
with `flow = TRUE`) to a ping-pong GIF, each object in its own colour
over a translucent reference volume.

## Usage

``` r
ggplot_flow_gif(
  flows,
  cols = NULL,
  volume = NULL,
  volume_alpha = 0.12,
  alpha = 0.6,
  rotation_matrix = NULL,
  file = NULL,
  width = 900,
  height = 800,
  delay = 0.14,
  dpi = 96
)
```

## Arguments

- flows:

  A named list; each element is the per-timepoint list of objects for
  one structure (all lists the same length). Colours are taken per
  structure.

- cols:

  Named vector of colours (one per `flows` entry); recycled from a
  default palette if `NULL`.

- volume:

  Optional fixed reference object drawn translucent under every frame
  (e.g. the target brain surface).

- volume_alpha, alpha:

  Alphas for the reference volume and the flow objects.

- rotation_matrix:

  Optional 4x4 view matrix (as `nat.ggplot`/`rgl` use).

- file:

  Output GIF path. Frames are written next to it.

- width, height, delay, dpi:

  GIF frame size (px), per-frame delay (s) and dpi.

## Value

The GIF path if written (needs `gifski`), else the vector of frame PNGs.
