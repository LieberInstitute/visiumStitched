# Get keys of neighboring spots

For a given row of a `tibble()` containing array coordinates, find the
associated spot's neighbors (belonging to the same capture area) and
return their keys.

## Usage

``` r
.get_neighbors(i, coords)
```

## Arguments

- i:

  An `integer(1)` giving a row index in `coords`.

- coords:

  A `tibble()` containing `array_row`, `array_col`, `key`, and
  `capture_area` columns.

## Value

A [`character()`](https://rdrr.io/r/base/character.html) of neighboring
spot keys.

## Author

Nicholas J. Eagles
