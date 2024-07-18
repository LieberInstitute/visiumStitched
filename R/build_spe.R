#' Build stitched \code{SpatialExperiment}
#'
#' First, read in capture-area-level spaceranger outputs. Then, overwrite
#' spatial coordinates and images to represent group-level samples using
#' \code{sample_info$group} (though keep original coordinates in
#' \code{colData} columns ending in "_original"). Next, add info about
#' overlaps (via \code{spe$exclude_overlapping} and \code{spe$overlap_key}).
#' Ultimately, return a \code{SpatialExperiment} ready for visualization or
#' downstream analysis.
#'
#' @inheritParams add_array_coords
#' @param count_type A \code{character(1)} vector passed to \code{type} from
#' \code{SpatialExperiment::read10xVisium}, defaulting to "sparse"
#' @param reference_gtf Passed to [spatialLIBD::read10xVisiumWrapper()]
#' @param gtf_cols Passed to [spatialLIBD::read10xVisiumWrapper()]
#'
#' @return A \code{SpatialExperiment} object with one sample per group specified
#' in \code{sample_info} using transformed pixel and array coordinates (including
#' in the \code{spatialCoords}).
#'
#' @import SpatialExperiment DropletUtils
#' @importFrom SummarizedExperiment assays rowData
#' @importFrom readr read_csv
#' @importFrom spatialLIBD read10xVisiumWrapper
#' @export
#' @author Nicholas J. Eagles
#'
#' @examples
#' #   For internal testing
#' \dontrun{
#' sample_info <- readr::read_csv("dev/test_data/sample_info.csv")
#' coords_dir <- "dev/test_data"
#' spe <- build_spe(sample_info, coords_dir)
#' }
#'
#' ## TODO: add working examples
#' args(build_spe)
build_spe <- function(sample_info, coords_dir, count_type = "sparse", reference_gtf = NULL, gtf_cols = c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type")) {
    #   State assumptions about columns expected to be in sample_info
    expected_cols <- c("capture_area", "group", "spaceranger_dir")
    if (!all(expected_cols %in% colnames(sample_info))) {
        stop(
            sprintf(
                'Missing at least one of the following columns in "sample_info": "%s"',
                paste(expected_cols, collapse = '", "')
            )
        )
    }

    message("Building SpatialExperiment using capture area as sample ID")
    if(missing(reference_gtf)) {
        spe <- spatialLIBD::read10xVisiumWrapper(
            samples = dirname(sample_info$spaceranger_dir),
            sample_id = sample_info$capture_area,
            type = count_type,
            data = "raw",
            images = "lowres",
            load = FALSE,
            gtf_cols = gtf_cols
        )
    } else {
        spe <- spatialLIBD::read10xVisiumWrapper(
            samples = dirname(sample_info$spaceranger_dir),
            sample_id = sample_info$capture_area,
            type = count_type,
            data = "raw",
            images = "lowres",
            load = FALSE,
            reference_gtf = reference_gtf,
            gtf_cols = gtf_cols
        )
    }

    message("Overwriting imgData(spe) with merged images (one per group)")
    all_groups <- unique(sample_info$group)

    img_data <- readImgData(
        path = file.path(coords_dir, all_groups),
        sample_id = all_groups,
        imageSources = file.path(
            coords_dir, all_groups, "tissue_lowres_image.png"
        ),
        scaleFactors = file.path(
            coords_dir, all_groups, "scalefactors_json.json"
        ),
        load = TRUE
    )

    coldata_fixed <- colData(spe) |>
        dplyr::as_tibble() |>
        dplyr::mutate(
            capture_area = factor(sample_id),
            group = factor(
                sample_info$group[
                    match(capture_area, sample_info$capture_area)
                ]
            ),
            #   Not made a factor because of https://github.com/drighelli/SpatialExperiment/issues/151
            sample_id = as.character(group),
            barcode = colnames(spe),
            key = paste(barcode, capture_area, sep = "_")
        )

    spe <- SpatialExperiment(
        assays = SummarizedExperiment::assays(spe),
        rowData = SummarizedExperiment::rowData(spe),
        colData = coldata_fixed,
        spatialCoords = spatialCoords(spe),
        imgData = img_data
    )

    spe <- add_array_coords(spe, sample_info, coords_dir, overwrite = TRUE)
    spe <- add_overlap_info(spe, "sum_umi")

    colnames(spe) = spe$key

    return(spe)
}
