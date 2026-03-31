# Check if coordinates are Visium-like

Sanity check designed to catch unforeseen bugs: halt if the tibble-like
`coords`, expected to contain columns 'array_row' and 'array_col',
represents an invalid Visium array.

## Usage

``` r
.validate_array(coords)
```

## Arguments

- coords:

  A [`data.frame()`](https://rdrr.io/r/base/data.frame.html) containing
  `'array_row'` and `'array_col'` columns calculated internally by
  [`add_array_coords()`](add_array_coords.md).

## Value

It returns `NULL` if all tests were correct.

## Author

Nicholas J. Eagles
