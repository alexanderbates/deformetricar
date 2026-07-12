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
  volume_col = "grey80",
  targets = NULL,
  target_cols = NULL,
  target_alpha = 0.18,
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
  (e.g. the target brain hull).

- volume_alpha, target_alpha, alpha:

  Alphas for the reference volume, the target objects and the flow
  objects.

- volume_col:

  Colour of the reference volume (default `"grey80"`; pass e.g. a light
  pink to render the brain hull as a soft envelope).

- targets:

  Optional named list of fixed reference objects drawn translucent
  *above* the volume but *below* the flow (e.g. the matched target
  neuropils each warping object should land on). Static across frames.

- target_cols:

  Named vector of colours for `targets` (one per entry); defaults to a
  transparent greyscale ramp so the coloured flow reads clearly on top.

- rotation_matrix:

  Optional 4x4 view matrix (as `nat.ggplot`/`rgl` use).

- file:

  Output GIF path. Frames are written next to it.

- width, height, delay, dpi:

  GIF frame size (px), per-frame delay (s) and dpi.

## Value

The GIF path if written (needs `gifski`, or falls back to `magick`),
else the vector of frame PNGs.
