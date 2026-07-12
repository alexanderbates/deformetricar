# Package index

## Deformetrica 4 client

Fit and apply diffeomorphisms with a modern Deformetrica (\>= 4.3)
install.

- [`install_deformetrica()`](https://alexanderbates.github.io/deformetricar/reference/install_deformetrica.md)
  : Install Deformetrica into a managed Python environment
- [`find_deformetrica()`](https://alexanderbates.github.io/deformetricar/reference/find_deformetrica.md)
  : Locate the Deformetrica (\>= 4.3) command-line executable
- [`deformetrica_register()`](https://alexanderbates.github.io/deformetricar/reference/deformetrica_register.md)
  : Fit a Deformetrica diffeomorphism between two point sets
  (registration)
- [`deformetrica_register_multi()`](https://alexanderbates.github.io/deformetricar/reference/deformetrica_register_multi.md)
  : Fit ONE diffeomorphism to many matched objects (multi-object
  registration)
- [`deformetrica_shoot()`](https://alexanderbates.github.io/deformetricar/reference/deformetrica_shoot.md)
  : Apply a fitted Deformetrica diffeomorphism to 3D points (geodesic
  shooting)

## Alignment helpers

Prepare and post-process objects for a cross-brain registration.

- [`affine_prealign()`](https://alexanderbates.github.io/deformetricar/reference/affine_prealign.md)
  : Affine pre-alignment of one object onto another
- [`split_mesh_lr()`](https://alexanderbates.github.io/deformetricar/reference/split_mesh_lr.md)
  : Split a surface mesh into left and right halves by a midline plane
- [`mirror_lr_split()`](https://alexanderbates.github.io/deformetricar/reference/mirror_lr_split.md)
  : Assign left/right to an (even unpaired, midline) object by
  self-mirror warping

## Visualisation

Animate a geodesic flow with nat.ggplot.

- [`ggplot_flow_gif()`](https://alexanderbates.github.io/deformetricar/reference/ggplot_flow_gif.md)
  : Animate a Deformetrica geodesic flow as a nat.ggplot GIF

## Object I/O

Read and write neurons, meshes and point sets as VTK.

- [`write_neuron_vtk()`](https://alexanderbates.github.io/deformetricar/reference/write_neuron_vtk.md)
  : Write a neuron as a VTK NonOrientedPolyLine
- [`read_vtk()`](https://alexanderbates.github.io/deformetricar/reference/read_vtk.md)
  : Read a VTK format file
- [`write_vtk()`](https://alexanderbates.github.io/deformetricar/reference/write_vtk.md)
  : Write VTK file
- [`transform_vtk()`](https://alexanderbates.github.io/deformetricar/reference/transform_vtk.md)
  : Transform a .vtk object
