# Build stitched `SpatialExperiment`

First, read in capture-area-level `SpaceRanger`
<https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/running-pipelines/space-ranger-count>
outputs. Then, overwrite spatial coordinates and images to represent
group-level samples using `sample_info$group` (though keep original
coordinates in `colData` columns ending with the suffix `"_original"`).
Next, add info about overlaps (via `spe$exclude_overlapping` and
`spe$overlap_key`). Ultimately, return a
[SpatialExperiment-class](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
ready for visualization or downstream analysis.

## Usage

``` r
build_SpatialExperiment(
  sample_info,
  coords_dir,
  count_type = c("sparse", "HDF5"),
  data_type = c("raw", "filtered"),
  reference_gtf = NULL,
  gtf_cols = c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type"),
  calc_error_metrics = FALSE,
  algorithm = c("LSAP", "Euclidean")
)
```

## Arguments

- sample_info:

  A [`data.frame()`](https://rdrr.io/r/base/data.frame.html) with
  columns `capture_area`, `group`, `fiji_xml_path`, `fiji_image_path`,
  `spaceranger_dir`, `intra_group_scalar`, and `group_hires_scalef`. The
  last two are made by
  [`rescale_fiji_inputs()`](rescale_fiji_inputs.md).

- coords_dir:

  A `character(1)` vector giving the directory containing sample
  directories each with `tissue_positions.csv`,
  `scalefactors_json.json`, and `tissue_lowres_image.png` files produced
  from refinement with [prep_fiji_coords()](prep_fiji.md) and related
  functions.

- count_type:

  A `character(1)` vector passed to `type` from
  `SpatialExperiment::read10xVisiumWrapper`, defaulting to "sparse".

- data_type:

  A `character(1)` vector passed to `data` from
  `SpatialExperiment::read10xVisiumWrapper`, defaulting to "raw".

- reference_gtf:

  Passed to
  [`spatialLIBD::read10xVisiumWrapper()`](https://rdrr.io/pkg/spatialLIBD/man/read10xVisiumWrapper.html).
  If working on the same system where SpaceRanger was run, the GTF will
  be automatically found; otherwise a `character(1)` path may be
  supplied, pointing to a GTF file of gene annotation to populate
  `rowData()` with.

- gtf_cols:

  Passed to
  [`spatialLIBD::read10xVisiumWrapper()`](https://rdrr.io/pkg/spatialLIBD/man/read10xVisiumWrapper.html).
  Columns in the reference GTF to extract and populate `rowData()`.

- calc_error_metrics:

  A `logical(1)` vector indicating whether to calculate error metrics
  related to mapping spots to well-defined array coordinates. If `TRUE`,
  adds `euclidean_error` and `shared_neighbors` spot-level metrics to
  the `colData()`. The former indicates distance in number of inter-spot
  distances to "move" a spot to the new array position; the latter
  indicates the fraction of neighbors for the associated capture area
  that are retained after mapping, which can be quite time-consuming to
  compute.

- algorithm:

  A `character(1)` vector indicating which mapping algorithm to employ
  when computing group-wide array coordinates. The default of "LSAP" is
  generally recommended, as it guarantees one-to-one mappings at the
  cost of computational time and some Euclidean error. The faster
  alternative "Euclidean" minimizes Euclidean error but may produce
  duplicate mappings, which is generally undesirable downstream (for
  clustering, etc).

## Value

A
[SpatialExperiment-class](https://rdrr.io/pkg/SpatialExperiment/man/SpatialExperiment.html)
object with one sample per group specified in `sample_info` using
transformed pixel and array coordinates (including in the
`spatialCoords()`).

## Author

Nicholas J. Eagles

## Examples

``` r
########################################################################
#   Prepare sample_info
########################################################################

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
#> 2026-04-01 15:03:51.881521 loading file /github/home/.cache/R/BiocFileCache/5244364f1fde_visiumStitched_brain_spaceranger.zip%3Frlkey%3Dbdgjc6mgy1ierdad6h6v5g29c%26dl%3D1
sample_info$spaceranger_dir <- file.path(
    sr_dir, sample_info$capture_area, "outs", "spatial"
)

#   Add Fiji-output-related columns
fiji_dir <- tempdir()
temp <- unzip(
    spatialLIBD::fetch_data("visiumStitched_brain_Fiji_out"),
    exdir = fiji_dir
)
#> 2026-04-01 15:03:54.167904 loading file /github/home/.cache/R/BiocFileCache/52444fc76ff3_visiumStitched_brain_fiji_out.zip%3Frlkey%3Dptwal8f5zxakzejwd0oqw0lhj%26dl%3D1
sample_info$fiji_xml_path <- temp[grep("xml$", temp)]
sample_info$fiji_image_path <- temp[grep("png$", temp)]

## Re-size images and add more information to the sample_info
sample_info <- rescale_fiji_inputs(sample_info, out_dir = tempdir())

## Preparing Fiji coordinates and images for build_SpatialExperiment()
spe_input_dir <- tempdir()
prep_fiji_coords(sample_info, out_dir = spe_input_dir)
#> [1] "/tmp/Rtmp1ekYg3/Br2719/tissue_positions.csv"
prep_fiji_image(sample_info, out_dir = spe_input_dir)
#> [1] "/tmp/Rtmp1ekYg3/Br2719/tissue_lowres_image.png"
#> [2] "/tmp/Rtmp1ekYg3/Br2719/scalefactors_json.json" 

########################################################################
#   Build the SpatialExperiment
########################################################################

#    Since we don't have access to the original GTF used to run SpaceRanger,
#    we must explicitly supply our own GTF to build_SpatialExperiment(). We use
#    GENCODE release 32, intended to be quite close to the actual GTF used,
#    which is available from:
#    https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
bfc <- BiocFileCache::BiocFileCache()
gtf_cache <- BiocFileCache::bfcrpath(
    bfc,
    paste0(
        "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
        "release_32/gencode.v32.annotation.gtf.gz"
    )
)

## Now we can build the stitched SpatialExperiment object. Use the Euclidean
## algorithm for calculating new array coordinates because of speed in this
## example, but "LSAP" is generally recommended.
spe <- build_SpatialExperiment(
    sample_info, coords_dir = spe_input_dir, reference_gtf = gtf_cache,
    algorithm = "Euclidean"
)
#> Building SpatialExperiment using capture area as sample ID
#> 2026-04-01 15:04:03.363108 SpatialExperiment::read10xVisium: reading basic data from SpaceRanger
#> Warning: 'SpatialExperiment::read10xVisium' is deprecated.
#> Use 'VisiumIO::TENxVisium(List)' instead.
#> See help("Deprecated")
#> 2026-04-01 15:04:06.279749 read10xVisiumAnalysis: reading analysis output from SpaceRanger
#> 2026-04-01 15:04:06.60157 add10xVisiumAnalysis: adding analysis output from SpaceRanger
#> 2026-04-01 15:04:06.841388 rtracklayer::import: reading the reference GTF file
#> 2026-04-01 15:04:31.7866 adding gene information to the SPE object
#> Warning: Gene IDs did not match. This typically happens when you are not using the same GTF file as the one that was used by SpaceRanger. For example, one file uses GENCODE IDs and the other one ENSEMBL IDs. read10xVisiumWrapper() will try to convert them to ENSEMBL IDs.
#> Warning: Dropping 2226 out of 38606 genes for which we don't have information on the reference GTF file. This typically happens when you are not using the same GTF file as the one that was used by SpaceRanger.
#> 2026-04-01 15:04:33.262936 adding information used by spatialLIBD
#> Overwriting imgData(spe) with merged images (one per group)
#> Adding array coordinates and overlap info

## Let's explore the stitched SpatialExperiment object
spe
#> class: SpatialExperiment 
#> dim: 36380 14976 
#> metadata(0):
#> assays(1): counts
#> rownames(36380): ENSG00000243485.5 ENSG00000237613.2 ...
#>   ENSG00000198695.2 ENSG00000198727.2
#> rowData names(6): source type ... gene_type gene_search
#> colnames(14976): AAACAACGAATAGTTC-1_V13B23-283_A1
#>   AAACAAGTATCTCCCA-1_V13B23-283_A1 ... TTGTTTGTATTACACG-1_V13B23-283_D1
#>   TTGTTTGTGTAAATTC-1_V13B23-283_D1
#> colData names(31): sample_id in_tissue ... overlap_key
#>   exclude_overlapping
#> reducedDimNames(0):
#> mainExpName: NULL
#> altExpNames(0):
#> spatialCoords names(2) : pxl_col_in_fullres pxl_row_in_fullres
#> imgData names(4): sample_id image_id data scaleFactor
```
