# Merge overlapping spots

Given a stitched
[SpatialExperiment-class](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html),
merge overlapping (same array coordinates) spots by adding expression
(i.e. from `assays(spe)$counts`), returning a `SpatialExperiment` with
at most one spot per array location.

## Usage

``` r
merge_overlapping(spe)
```

## Arguments

- spe:

  A
  [SpatialExperiment-class](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
  with `colData(spe)` columns `array_row`, `array_col`, `key`, `group`,
  and `capture_area`.

## Value

A
[SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
with at most one spot per array location

## Details

`colData(spe)` and `spatialCoords(spe)` of the merged spots are taken
from the spots whose `exclude_overlapping` values are `TRUE`.

## Author

Nicholas J. Eagles

## Examples

``` r
if (!exists("spe")) {
    spe <- spatialLIBD::fetch_data(type = "visiumStitched_brain_spe")
}
#> 2026-04-01 15:04:37.644472 loading file /github/home/.cache/R/BiocFileCache/524411ae5eee_visiumStitched_brain_spe.rds%3Frlkey%3Dnq6a82u23xuu9hohr86oodwdi%26dl%3D1

#   Group colData by group and array coordinates
grouped_coldata <- colData(spe) |>
    dplyr::as_tibble() |>
    dplyr::group_by(group, array_row, array_col)

#   Find the first 100 keys that overlap other spots and don't, respectively
overlapping_keys <- grouped_coldata |>
    dplyr::filter(dplyr::n() > 1) |>
    dplyr::slice_head(n = 2) |>
    dplyr::ungroup() |>
    dplyr::slice_head(n = 100) |>
    dplyr::pull(key)
nonoverlapping_keys <- grouped_coldata |>
    dplyr::filter(dplyr::n() == 1) |>
    dplyr::ungroup() |>
    dplyr::slice_head(n = 100) |>
    dplyr::pull(key)

#   Built a small SPE containing some overlaps and some non-overlapping spots
small_spe <- spe[, c(overlapping_keys, nonoverlapping_keys)]

#   Merge overlapping spots
small_spe_merged <- merge_overlapping(small_spe)
#> Warning: Dropping assays other than 'counts' for merging.
#> Warning: Dropped reducedDims(spe) for merging
#> 'sample_id's are duplicated across 'SpatialExperiment' objects to cbind; appending sample indices.

#   All array coordinates have just one unique spot after merging
colData(small_spe_merged) |>
    dplyr::as_tibble() |>
    dplyr::group_by(group, array_row, array_col) |>
    dplyr::summarize(n = dplyr::n()) |>
    dplyr::pull(n) |>
    table()
#> `summarise()` has regrouped the output.
#> ℹ Summaries were computed grouped by group, array_row, and array_col.
#> ℹ Output is grouped by group and array_row.
#> ℹ Use `summarise(.groups = "drop_last")` to silence this message.
#> ℹ Use `summarise(.by = c(group, array_row, array_col))` for per-operation
#>   grouping (`?dplyr::dplyr_by`) instead.
#> 
#>   1 
#> 142 
```
