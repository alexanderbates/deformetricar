# Write a neuron as a VTK NonOrientedPolyLine

Neuron backbones are the natural registration object for connectome work
(cf. the FAFB left-right bridging registration in
flyconnectome/deformetricaLR).

## Usage

``` r
write_neuron_vtk(x, file)
```

## Arguments

- x:

  A `nat` neuron, or a `neuronlist` (all backbones are written into one
  VTK, so a whole set of neurons can be registered as a single PolyLine
  object).

- file:

  Output path.

## Value

"complete" (invisibly via the file written).
