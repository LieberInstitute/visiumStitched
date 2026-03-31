# Calculate fraction of neighbors retained after mapping to new array coordinates

Given `tibble()`s before and after mapping to new array coordinates,
calculate for each spot the fraction of starting neighboring spots that
were retained in the new array-coordinate system. Add this metric and
return.

## Usage

``` r
.get_shared_neighbors(coords_new, coords)
```

## Arguments

- coords_new:

  A `tibble()` containing `array_row`, `array_col`, `key`, and
  `capture_area` columns, representing data after mapping to new array
  coordinates.

- coords:

  A `tibble()` containing `array_row`, `array_col`, `key`, and
  `capture_area` columns, representing data before mapping to new array
  coordinates.

## Value

A `tibble()` copy of `coords_new` with additional `shared_neighbors`
column.

## Author

Nicholas J. Eagles
