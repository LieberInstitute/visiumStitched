# Fit spots to a new Visium-like array: LSAP approach

Given transformed pixel coordinates, modify the 'array_row' and
'array_col' columns to represent a larger Visium capture area containing
all capture areas in a common coordinate system. The number of array
rows/cols generally changes from the Visium standards of 78 and 128 (and
even may change in ratio between num rows and num cols).

## Usage

``` r
.fit_to_array_lsap(coords, inter_spot_dist_px)
```

## Arguments

- coords:

  A [`data.frame()`](https://rdrr.io/r/base/data.frame.html) containing
  capture areas of the same group, and containing columns 'key',
  'array_row', 'array_col', 'pxl_row_in_fullres', and
  'pxl_col_in_fullres'.

- inter_spot_dist_px:

  `numeric(1)` vector giving the pixel distance between any 2 spots in
  the new coordinates.

## Value

A [tibble](https://dplyr.tidyverse.org/reference/reexports.html) with
modified `array_row` + `array_col` columns, as well as new
`pxl_row_in_fullres_rounded` and `pxl_col_in_fullres_rounded` columns
representing the pixel coordinates rounded to the nearest exact array
coordinates.

## Details

Mapping to the proper array coordinates is framed as the linear sum
assignment problem, and solved using the Hungarian algorithm. This
approach is far slower than [`.fit_to_array()`](dot-fit_to_array.md),
running at O(n^3) with the number of spots, but guarantees a one-to-one
mapping of starting to target spots, at a small cost in the Euclidean
distance moved.

## Author

Nicholas J. Eagles
