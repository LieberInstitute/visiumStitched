# Map source spots to best target spots by solving the LSAP

Given `source_coords` and `target_coords`, both containing pixel
coordinates of spots, map each spot in `source_coords` to a unique spot
in `target_coords` such that the total squared Euclidean distance
between matched spots is minimized, with guaranteed one-to-one mapping.
This is done by solving the Linear Sum Assignment Problem (LSAP) using
the Hungarian algorithm. Return the `source_coords` with the newly
mapped `array_row` and `array_col` columns.

## Usage

``` r
.map_lsap(source_coords, target_coords)
```

## Arguments

- source_coords:

  A [`data.frame()`](https://rdrr.io/r/base/data.frame.html) containing
  the pixel coordinates (i.e. 'pxl_row_in_fullres' and
  'pxl_col_in_fullres') of starting spots from one capture area.

- target_coords:

  A [`data.frame()`](https://rdrr.io/r/base/data.frame.html) containing
  the pixel coordinates (i.e. 'pxl_row_in_fullres' and
  'pxl_col_in_fullres') of target spots which should just barely
  encompass the capture area in `source_coords`.

## Value

A [tibble](https://dplyr.tidyverse.org/reference/reexports.html) with
the same rows as `source_coords`, but with the `array_row` and
`array_col` columns (and rounded pixel coordinates) taken from the
best-matching spots in `target_coords`.

## Author

Nicholas J. Eagles
