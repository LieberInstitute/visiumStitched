# Add error metrics related to array-coordinate mapping

Given `tibble()`s before and after mapping to new array coordinates,
calculate metrics related to the suitability of the mapping.

## Usage

``` r
.add_error_metrics(coords, coords_new, inter_spot_dist_px)
```

## Arguments

- coords:

  A `tibble()` containing `array_row`, `array_col`, `key`,
  `pxl_col_in_fullres`, `pxl_row_in_fullres`,
  `pxl_col_in_fullres_rounded`, `pxl_row_in_fullres_rounded`, and
  `capture_area` columns, representing data before mapping to new array
  coordinates for one `group`.

- coords_new:

  A `tibble()` containing `array_row`, `array_col`, `key`,
  `pxl_col_in_fullres`, `pxl_row_in_fullres`,
  `pxl_col_in_fullres_rounded`, and `pxl_row_in_fullres_rounded`
  columns, representing data after mapping to new array coordinates for
  one `group`.

- inter_spot_dist_px:

  A `numeric(1)` giving the number of pixels between spots for the
  `group`.

## Value

A `tibble()` copy of `coords_new` with additional `shared_neighbors` and
`euclidean_error` columns.

## Details

Add column `shared_neighbors`, the fraction of neighbors a spot started
with that are retained after mapping; add column `euclidean_error`, the
number of multiples of the inter-spot distance a spot must move to be
placed in the new array coordinates.

## Author

Nicholas J. Eagles
