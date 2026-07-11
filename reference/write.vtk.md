# Write VTK file

Write VTK file

## Usage

``` r
write.vtk(
  points,
  filename,
  polygons = NULL,
  normals = NULL,
  datatype = c("float", "double"),
  title = filename
)
```

## Arguments

- points:

  3D coordinates (Nx3 matrix)

- filename:

  Path to output file

- polygons:

  Triangle indices (Mx3 matrix)

- normals:

  Normals (Nx3 matrix)

- datatype:

  .vtk datatype (defaults to float)

- title:

  Title of the .vtk file (defaults to filename)
