# Read a VTK format file

Read a VTK format file

## Usage

``` r
read_vtk(filename, item = c("points", "triangles", "normals"))
```

## Arguments

- filename:

  The path to the file on disk

- item:

  The element(s) within the file to read (defaults to points)

## Value

A matrix of points, indices (polygons) or normals
