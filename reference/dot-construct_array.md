# Construct a new Visium-like array encapsulating a set of spots

Given `coords` containing pixel coordinates of spots from potentially
multiple capture areas, return a new Visium-like array encapsulating all
such spots.

## Usage

``` r
.construct_array(coords, inter_spot_dist_px, buffer = 1)
```

## Arguments

- coords:

  A [`data.frame()`](https://rdrr.io/r/base/data.frame.html) with
  columns 'pxl_row_in_fullres' and 'pxl_col_in_fullres' whose rows
  contain spots from potentially multiple capture areas.

- inter_spot_dist_px:

  `numeric(1)` vector giving the pixel distance between any 2 spots in
  the new coordinates.

- buffer:

  `numeric(1)` vector giving the number of spot distances to pad the new
  array (on all sides) beyond the min/max pixel coordinates in `coords`.

## Value

A [tibble](https://dplyr.tidyverse.org/reference/reexports.html) with
columns 'array_row', 'array_col', 'pxl_row_in_fullres', and
'pxl_col_in_fullres', representing the new Visium-like array.

## Author

Nicholas J. Eagles
