# Add info about how spots overlap among capture areas

Given a
[SpatialExperiment-class](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
and column name in its `colData`, return a modified copy of the
`SpatialExperiment` with additional `colData` columns:
`spe$exclude_overlapping` and `spe$overlap_key`.

## Usage

``` r
add_overlap_info(spe, metric_name)
```

## Arguments

- spe:

  A
  [SpatialExperiment-class](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
  with `colData(spe)` columns `array_row`, `array_col`, `key`, and
  `capture_area`.

- metric_name:

  `character(1)` in `colnames(colData(spe))`, where spots belonging to
  the capture area with highest average value for the metric take
  precedence over other spots.

## Value

A
[SpatialExperiment](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
object with additional `colData` columns `spe$exclude_overlapping` and
`spe$overlap_key`.

## Details

`spe$exclude_overlapping` is `TRUE` for spots with a higher-quality
overlapping capture area and `FALSE` otherwise.
[vis_clus](https://rdrr.io/pkg/spatialLIBD/man/vis_clus.html)
onlydisplays `FALSE` spots to prevent overplotting in regions of
overlap. `spe$overlap_key` gives comma-separated strings containing the
keys of any overlapping spots, and is the empty string otherwise.

## Author

Nicholas J. Eagles

## Examples

``` r
if (!exists("spe")) {
    spe <- spatialLIBD::fetch_data(type = "visiumStitched_brain_spe")
}
#> 2026-04-01 15:03:39.801849 loading file /github/home/.cache/R/BiocFileCache/524411ae5eee_visiumStitched_brain_spe.rds%3Frlkey%3Dnq6a82u23xuu9hohr86oodwdi%26dl%3D1

#    Find the mean of the 'sum_umi' metric by capture area to understand
#    which capture areas will be excluded in regions of overlap
SummarizedExperiment::colData(spe) |>
    dplyr::as_tibble() |>
    dplyr::group_by(capture_area) |>
    dplyr::summarize(mean_sum_umi = mean(sum_umi))
#> # A tibble: 3 × 2
#>   capture_area  mean_sum_umi
#>   <fct>                <dbl>
#> 1 V13B23-283_A1        2075.
#> 2 V13B23-283_C1        1188.
#> 3 V13B23-283_D1        2083.

spe <- add_overlap_info(spe, "sum_umi")

#    See how many spots were excluded by capture area
table(spe$exclude_overlapping, spe$capture_area)
#>        
#>         V13B23-283_A1 V13B23-283_C1 V13B23-283_D1
#>   FALSE          3920          4367          4596
#>   TRUE           1032            50             0

#    Examine how data about overlapping spots is stored (for the first
#    few spots with overlap)
head(spe$overlap_key[spe$overlap_key != ""])
#> [1] "TCAGAACGGCGGTAAT-1_V13B23-283_D1,TGCCCGATAGTTAGAA-1_V13B23-283_D1"
#> [2] "CGACTTGCCGGGAAAT-1_V13B23-283_D1,TCATGAAGCGCTGCAT-1_V13B23-283_D1"
#> [3] "CCTTAAGTACGCAATT-1_V13B23-283_D1"                                 
#> [4] "GAAGCCACTGATTATG-1_V13B23-283_D1"                                 
#> [5] "ACTCAAGTGCAAGGCT-1_V13B23-283_D1"                                 
#> [6] "GGCCTGCTTCTCCCGA-1_V13B23-283_D1"                                 
```
