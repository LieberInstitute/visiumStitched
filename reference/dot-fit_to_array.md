# Fit spots to a new Visium-like array: fast Euclidean approach

Given transformed pixel coordinates, modify the 'array_row' and
'array_col' columns to represent a larger Visium capture area containing
all capture areas in a common coordinate system. The number of array
rows/cols generally changes from the Visium standards of 78 and 128 (and
even may change in ratio between num rows and num cols).

## Usage

``` r
.fit_to_array(coords, inter_spot_dist_px)
```

## Arguments

- coords:

  A [`data.frame()`](https://rdrr.io/r/base/data.frame.html) whose rows
  represent capture areas of the same group, and containing columns
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

The mapping algorithm minimizes Euclidean distance of each source spot
to each target spot. Runtime is O(n) with the number of spots, making it
extremely fast. However, the Euclidean approach countintuitively may
result in duplicated mappings (one source to the same target) as well as
unexpected "holes" in the target array, which is often undesirable
downstream.

## Author

Nicholas J. Eagles
