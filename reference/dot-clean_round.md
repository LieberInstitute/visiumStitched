# Round to the nearest integer, always rounding up at 0.5

This consistent behavior is favorable for our application, where we want
to minimize duplicate mappings of spots to new array coordinates.

## Usage

``` r
.clean_round(x)
```

## Arguments

- x:

  [`numeric()`](https://rdrr.io/r/base/numeric.html) vector.

## Value

A [`numeric()`](https://rdrr.io/r/base/numeric.html) vector rounded to
the nearest integer.

## Author

Nicholas J. Eagles
