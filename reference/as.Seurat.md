# Convert a `SpatialExperiment` object to a `Seurat` object

Given a
[SpatialExperiment-class](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
object, first `as.Seurat()` is run, which operates on
[SingleCellExperiment-class](https://rdrr.io/pkg/SingleCellExperiment/man/SingleCellExperiment.html)
objects. The remaining components (images, spatial coordinates) are
added manually. The actual appearance of images are buggy for now.

## Usage

``` r
as.Seurat(
  spe,
  spatial_cols = c(tissue = "in_tissue", row = "array_row", col = "array_col", imagerow =
    "pxl_row_in_fullres", imagecol = "pxl_col_in_fullres"),
  verbose = TRUE
)
```

## Arguments

- spe:

  A
  [SpatialExperiment-class](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
  with `colData()` or `spatialCoords()` columns given by `spatial_cols`.
  This does not have to be a stitched `spe` object as this function
  should work with any type of `spe` objects.

- spatial_cols:

  A `character(5)` named vector mapping which `colData(spe)` or
  `spatialCoords(spe)` columns contain the `tissue`, `row`, `col`,
  `imagerow`, and `imagecol` information expected by Seurat.

- verbose:

  A `logical(1)` vector. If `TRUE`, print status update about the
  conversion process. This information can be useful for debugging.

## Value

A `Seurat` object.

## Details

Note that only the `lowres` images from `imgData(spe)` will be used.

## Author

Nicholas J. Eagles

## Examples

``` r
## Download some example data
spe_unstitched <- spatialLIBD::fetch_data(
    type = "spatialDLPFC_Visium_example_subset"
)[seq(100), seq(100)]
#> 2026-04-01 15:03:43.816041 loading file /github/home/.cache/R/BiocFileCache/656e2438ffe5_spatialDLPFC_spe_subset_example.rds%3Fdl%3D1

## Make the column names unique
colnames(spe_unstitched) <- spatialLIBD::add_key(spe_unstitched)$key
#> Overwriting 'spe$key'. Set 'overwrite = FALSE' if you do not want to overwrite it.

## Convert from a SpatialExperiment to a Seurat object
seur <- as.Seurat(spe_unstitched)
#> Running 'as.Seurat(spe)'...
#> Registered S3 method overwritten by 'spatstat.geom':
#>   method      from  
#>   plot.imlist imager
#> Warning: Keys should be one or more alphanumeric characters followed by an underscore, setting key from pc. to pc_
#> Warning: Keys should be one or more alphanumeric characters followed by an underscore, setting key from tsne. to tsne_
#> Warning: Keys should be one or more alphanumeric characters followed by an underscore, setting key from umap. to umap_
#> Warning: Keys should be one or more alphanumeric characters followed by an underscore, setting key from PC to PC_
#> Warning: Keys should be one or more alphanumeric characters followed by an underscore, setting key from UMAP to UMAP_
#> Warning: Keys should be one or more alphanumeric characters followed by an underscore, setting key from TSNE_perplexity80_ to TSNEperplexity80_
#> Warning: Keys should be one or more alphanumeric characters followed by an underscore, setting key from UMAP to UMAP_
#> Warning: Key ‘UMAP_’ taken, using ‘umapharmony_’ instead
#> Adding spot coordinates and images for sample Br6432_ant...
#> Returning converted object...
seur
#> An object of class Seurat 
#> 100 features across 100 samples within 1 assay 
#> Active assay: originalexp (100 features, 0 variable features)
#>  2 layers present: counts, data
#>  8 dimensional reductions calculated: X10x_pca, X10x_tsne, X10x_umap, PCA, UMAP, TSNE_perplexity80, HARMONY, UMAP.HARMONY
#>  1 image present: Br6432_ant

## Example with an stitched SPE object
if (!exists("spe")) {
    spe <- spatialLIBD::fetch_data(type = "visiumStitched_brain_spe")
}
#> 2026-04-01 15:03:49.312173 loading file /github/home/.cache/R/BiocFileCache/524411ae5eee_visiumStitched_brain_spe.rds%3Frlkey%3Dnq6a82u23xuu9hohr86oodwdi%26dl%3D1
seur_stitched <- as.Seurat(spe[seq(100), seq(100)])
#> Running 'as.Seurat(spe)'...
#> Warning: Keys should be one or more alphanumeric characters followed by an underscore, setting key from PC to PC_
#> Adding spot coordinates and images for sample Br2719...
#> Returning converted object...

## Let's look at our resulting Seurat object
seur_stitched
#> An object of class Seurat 
#> 100 features across 100 samples within 1 assay 
#> Active assay: originalexp (100 features, 0 variable features)
#>  2 layers present: counts, data
#>  1 dimensional reduction calculated: PCA
#>  1 image present: Br2719
```
