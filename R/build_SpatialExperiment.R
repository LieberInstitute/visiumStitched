#' Build stitched \code{SpatialExperiment}
#'
#' First, read in capture-area-level `SpaceRanger`
#' <https://www.10xgenomics.com/support/software/space-ranger/latest/analysis/running-pipelines/space-ranger-count>
#' outputs. Then, overwrite
#' spatial coordinates and images to represent group-level samples using
#' \code{sample_info$group} (though keep original coordinates in
#' \code{colData} columns ending with the suffix `"_original"`). Next, add info
#' about overlaps (via \code{spe$exclude_overlapping} and
#' \code{spe$overlap_key}).
#' Ultimately, return a
#' [SpatialExperiment-class][SpatialExperiment::SpatialExperiment-class] ready
#' for visualization or downstream analysis.
#'
#' @inheritParams add_array_coords
#' @param count_type A \code{character(1)} vector passed to \code{type} from
#' \code{SpatialExperiment::read10xVisium}, defaulting to "sparse".
#' @param reference_gtf Passed to [spatialLIBD::read10xVisiumWrapper()]. If
#' working on the same system where SpaceRanger was run, the GTF will be
#' automatically found; otherwise a `character(1)` path may be supplied,
#' pointing to a GTF file of gene annotation to populate \code{rowData()} with.
#' @param gtf_cols Passed to [spatialLIBD::read10xVisiumWrapper()]. Columns
#' in the reference GTF to extract and populate \code{rowData()}.
#'
#' @return A [SpatialExperiment-class][SpatialExperiment::SpatialExperiment-class]
#' object with one sample per group specified in \code{sample_info} using
#' transformed pixel and array coordinates (including  in the
#' \code{spatialCoords()}).
#'
#' @import SpatialExperiment DropletUtils
#' @importFrom SummarizedExperiment assays rowData
#' @importFrom readr read_csv
#' @importFrom spatialLIBD read10xVisiumWrapper
#' @export
#' @author Nicholas J. Eagles
#'
#' @examples
#' ########################################################################
#' #   Prepare sample_info
#' ########################################################################
#'
#' if (file.exists("sample_info.rds")) {
#'     sample_info <- readRDS("sample_info.rds")
#' } else {
#'     sample_info <- dplyr::tibble(
#'         group = "Br2719",
#'         capture_area = c("V13B23-283_A1", "V13B23-283_C1", "V13B23-283_D1")
#'     )
#'     #   Add 'spaceranger_dir' column
#'     sr_dir <- tempdir()
#'     temp <- unzip(
#'         spatialLIBD::fetch_data("visiumStitched_brain_spaceranger"),
#'         exdir = sr_dir
#'     )
#'     sample_info$spaceranger_dir <- file.path(
#'         sr_dir, sample_info$capture_area, "outs", "spatial"
#'     )
#'
#'     #   Add Fiji-output-related columns
#'     fiji_dir <- tempdir()
#'     temp <- unzip(
#'         spatialLIBD::fetch_data("visiumStitched_brain_Fiji_out"),
#'         exdir = fiji_dir
#'     )
#'     sample_info$fiji_xml_path <- temp[grep("xml$", temp)]
#'     sample_info$fiji_image_path <- temp[grep("png$", temp)]
#'
#'     ## Re-size images and add more information to the sample_info
#'     sample_info <- rescale_fiji_inputs(sample_info, out_dir = tempdir())
#'
#'     saveRDS(sample_info, "sample_info.rds")
#' }
#'
#' ## Preparing Fiji coordinates and images for build_SpatialExperiment()
#' spe_input_dir <- tempdir()
#' prep_fiji_coords(sample_info, out_dir = spe_input_dir)
#' prep_fiji_image(sample_info, out_dir = spe_input_dir)
#'
#' ########################################################################
#' #   Build the SpatialExperiment
#' ########################################################################
#'
#' #    Since we don't have access to the original GTF used to run SpaceRanger,
#' #    we must explicitly supply our own GTF to build_SpatialExperiment(). We use
#' #    GENCODE release 32, intended to be quite close to the actual GTF used,
#' #    which is available from:
#' #    https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-GRCh38-2024-A.tar.gz
#' bfc <- BiocFileCache::BiocFileCache()
#' gtf_cache <- BiocFileCache::bfcrpath(
#'     bfc,
#'     paste0(
#'         "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
#'         "release_32/gencode.v32.annotation.gtf.gz"
#'     )
#' )
#'
#' ## Now we can build the stitched SpatialExperiment object
#' spe <- build_SpatialExperiment(
#'     sample_info,
#'     coords_dir = spe_input_dir, reference_gtf = gtf_cache
#' )
#'
#' ## Let's explore the stitched SpatialExperiment object
#' spe
build_SpatialExperiment <- function(sample_info, coords_dir, count_type = "sparse", reference_gtf = NULL, gtf_cols = c("source", "type", "gene_id", "gene_version", "gene_name", "gene_type"), calc_error_metrics = FALSE) {
    ## For R CMD check
    sample_id <- capture_area <- group <- barcode <- NULL

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
    if (missing(reference_gtf)) {
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

    if (calc_error_metrics) {
        message("Adding array coordinates with error metrics and adding overlap info")
    } else {
        message("Adding array coordinates and overlap info")
    }

    spe <- add_array_coords(
        spe, sample_info, coords_dir,
        calc_error_metrics = calc_error_metrics
    )
    spe <- add_overlap_info(spe, "sum_umi")

    colnames(spe) <- spe$key

    return(spe)
}
