# Return array coordinates fit to nearest spot with associated error

First, values of `x` are rounded to the nearest integer. Then, values of
`y` are rounded to the nearest valid integer under the constraint that
coordinates for x and y must be both odd or both even. These rounded
values are returned, along with the Euclidean distance needed to move x
and y from their original, non-integer values to their rounded values.

## Usage

``` r
.refine_fit(x, y, INTERVAL_X, INTERVAL_Y)
```

## Arguments

- x:

  [`numeric()`](https://rdrr.io/r/base/numeric.html) vector giving
  "ideal" array coordinates given every spot's transformed pixel
  coordinates.

- y:

  Same as x, though y must represent ideal array columns iff x
  represents array rows, and vice versa.

- INTERVAL_X:

  `numeric(1)` giving pixel distance between coordinate units used for
  `x` (e.g. if x represents ideal `array_col` values, `INTERVAL_X`
  represents pixel distance between spot columns).

- INTERVAL_Y:

  `numeric(1)` giving pixel distance between coordinate units used for
  `y`.

## Value

A `list` consisting of 3 unnamed
[`numeric()`](https://rdrr.io/r/base/numeric.html) vectors: rounded `x`,
rounded `y`, and the Euclidean distance in pixels from rounding both `x`
and `y`.

## Author

Nicholas J. Eagles
