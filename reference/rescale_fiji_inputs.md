# Write same-scale hires images for input to Fiji

Given a [`data.frame()`](https://rdrr.io/r/base/data.frame.html) of
sample information (`sample_info`) with columns `capture_area`, `group`,
and `spaceranger_dir`, Write new high-resolution images for use as input
to Fiji <https://imagej.net/software/fiji/>. Particularly when capture
areas come from different slides, there is a risk of significant scale
differences among SpaceRanger's `tissue_hires_image.png` images; that
is, the physical distance represented by a pixel from each capture area
may differ nontrivially, leading to a distance-distorted output image,
and inconsistent scaling when later transforming pixel coordinates. This
function writes approximately high-res images whose pixels are of equal
physical size within each `group`, then adds `intra_group_scalar` and
`group_hires_scalef` columns to `sample_info`. `intra_group_scalar`
gives the scalar by a which a given capture area's
`tissue_hires_image.png` image and pixel coordinates must be multiplied
to match the scale of other `group` members; `group_hires_scalef` gives
the new `tissue_hires_scalef` (as from SpaceRanger's
`scalefactors_json.json` file) appropriate for every capture area from
the group.

## Usage

``` r
rescale_fiji_inputs(sample_info, out_dir)
```

## Arguments

- sample_info:

  A [`data.frame()`](https://rdrr.io/r/base/data.frame.html) with
  columns `capture_area`, `group`, `fiji_xml_path`, `fiji_image_path`,
  `spaceranger_dir`, `intra_group_scalar`, and `group_hires_scalef`. The
  last two are made by `rescale_fiji_inputs()`.

- out_dir:

  A `character(1)` vector giving a path to a directory to place the
  output images, which must exist in advance.

## Value

A [tibble](https://dplyr.tidyverse.org/reference/reexports.html): a copy
of `sample_info` with additional columns `intra_group_scalar` and
`group_hires_scalef`.

## Author

Nicholas J. Eagles

## Examples

``` r
#    Define sample information for the example human brain data
sample_info <- dplyr::tibble(
    group = "Br2719",
    capture_area = c("V13B23-283_A1", "V13B23-283_C1", "V13B23-283_D1")
)
#   Add 'spaceranger_dir' column
sr_dir <- tempdir()
temp <- unzip(
    spatialLIBD::fetch_data("visiumStitched_brain_spaceranger"),
    exdir = sr_dir
)
#> 2026-04-01 15:05:00.732668 loading file /github/home/.cache/R/BiocFileCache/5244364f1fde_visiumStitched_brain_spaceranger.zip%3Frlkey%3Dbdgjc6mgy1ierdad6h6v5g29c%26dl%3D1
sample_info$spaceranger_dir <- file.path(
    sr_dir, sample_info$capture_area, "outs", "spatial"
)

#   Add Fiji-output-related columns
fiji_dir <- tempdir()
temp <- unzip(
    spatialLIBD::fetch_data("visiumStitched_brain_Fiji_out"),
    exdir = fiji_dir
)
#> 2026-04-01 15:05:02.920718 loading file /github/home/.cache/R/BiocFileCache/52444fc76ff3_visiumStitched_brain_fiji_out.zip%3Frlkey%3Dptwal8f5zxakzejwd0oqw0lhj%26dl%3D1
sample_info$fiji_xml_path <- temp[grep("xml$", temp)]
sample_info$fiji_image_path <- temp[grep("png$", temp)]

## Re-size images and add more information to the sample_info
out_dir <- tempdir()
sample_info_new <- rescale_fiji_inputs(sample_info, out_dir = out_dir)

#    Scale factors are computed that are necessary downstream (i.e. with
#    prep_fiji_*() functions)
sample_info_new[, setdiff(colnames(sample_info_new), colnames(sample_info))]
#> # A tibble: 3 × 2
#>   intra_group_scalar group_hires_scalef
#>                <dbl>              <dbl>
#> 1               1.00             0.0825
#> 2               1.00             0.0825
#> 3               1                0.0825

#    Image are produced that are ready for alignment in Fiji
list.files(out_dir)
#>  [1] "Br2719"                                   
#>  [2] "Br2719.png"                               
#>  [3] "Br2719.xml"                               
#>  [4] "HDF5Array_dataset_creation_global_counter"
#>  [5] "HDF5Array_dump"                           
#>  [6] "HDF5Array_dump_log"                       
#>  [7] "HDF5Array_dump_names_global_counter"      
#>  [8] "V13B23-283_A1"                            
#>  [9] "V13B23-283_A1.png"                        
#> [10] "V13B23-283_C1"                            
#> [11] "V13B23-283_C1.png"                        
#> [12] "V13B23-283_D1"                            
#> [13] "V13B23-283_D1.png"                        
#> [14] "bslib-33c615652fac282942ee0d98dcfba021"   
#> [15] "downlit"                                  
#> [16] "file89a412ef3c98"                         
#> [17] "file89a448cc94bc"                         
#> [18] "file89a454861104"                         
#> [19] "file89a4676af23f"                         
#> [20] "file89a46ed651cc"                         
#> [21] "file89a4e1da180"                          
```
