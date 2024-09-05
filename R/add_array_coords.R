#' Add transformed array and pixel coordinates to a \code{SpatialExperiment}
#'
#' Given a [SpatialExperiment-class][SpatialExperiment::SpatialExperiment-class],
#' sample information, and coordinates
#' produced from the refinement workflow, add array and pixel coordinates
#' appropriate for the linearly transformed capture areas making up each group
#' present in the [SpatialExperiment-class][SpatialExperiment::SpatialExperiment-class].
#'
#' Array coordinates are determined via an algorithm that fits each spot to
#' the nearest spot on a new, imaginary, Visium-like capture area. The imaginary
#' capture area differs from a real capture area only in its extent; array
#' coordinates still start at 0 but may extend arbitrarily beyond the normal
#' maximum indices of 77 and 127 to fit every capture area in each group
#' defined in the [SpatialExperiment-class][SpatialExperiment::SpatialExperiment-class].
#' The goal is to return well-defined
#' array coordinates in a consistent spatial orientation for each group, such
#' that downstream applications, such as clustering with `BayesSpace`, can
#' process each group as if it really were one capture area in the first place.
#' See
#' <https://research.libd.org/visiumStitched/articles/visiumStitched.html#defining-array-coordinates>
#' for more details.
#'
#' @param spe A
#' [SpatialExperiment-class][SpatialExperiment::SpatialExperiment-class] object.
#' @param sample_info A `data.frame()` with columns \code{capture_area},
#' \code{group}, \code{fiji_xml_path}, \code{fiji_image_path},
#' \code{spaceranger_dir}, \code{intra_group_scalar}, and
#' \code{group_hires_scalef}. The last two are made by `rescale_fiji_inputs()`.
#' @param coords_dir A \code{character(1)} vector giving the directory
#' containing sample directories each with \code{tissue_positions.csv},
#' \code{scalefactors_json.json}, and \code{tissue_lowres_image.png} files
#' produced from refinement with [prep_fiji_coords()][visiumStitched::prep_fiji_coords]
#' and related functions.
#' @param calc_error_metrics A \code{logical(1)} vector indicating whether to
#' calculate error metrics related to mapping spots to well-defined array
#' coordinates. If \code{TRUE}, adds \code{euclidean_error} and
#' \code{shared_neighbors} spot-level metrics to the \code{colData()}. The former
#' indicates distance in number of inter-spot distances to "move" a spot to the
#' new array position; the latter indicates the fraction of neighbors for the
#' associated capture area that are retained after mapping, which can be quite
#' time-consuming to compute.
#'
#' @return A [SpatialExperiment-class][SpatialExperiment::SpatialExperiment-class]
#' object with additional \code{colData}
#' columns \code{pxl_row_in_fullres_[suffix]} and \code{pxl_col_in_fullres_[suffix]}
#' with \code{[suffix]} values \code{original} and \code{rounded};
#' \code{array_row_original} and \code{array_col_original} columns; and
#' modified \code{colData()} columns \code{array_row} and
#' \code{array_col} and \code{spatialCoords()} with their transformed values.
#'
#' @import dplyr
#' @importFrom readr read_csv
#' @importFrom S4Vectors DataFrame
#' @importFrom rjson fromJSON
#' @export
#' @author Nicholas J. Eagles
#'
#' @examples
#' if (!exists("spe")) {
#'     spe <- spatialLIBD::fetch_data(type = "visiumStitched_brain_spe")
#' }
#'
#' ########################################################################
#' #   Prepare sample_info
#' ########################################################################
#'
#' if (file.exists("sample_info.rds")) {
#'     sample_info <- readRDS('sample_info.rds')
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
#' ## Preparing Fiji coordinates and images for build_spe()
#' spe_input_dir <- tempdir()
#' prep_fiji_coords(sample_info, out_dir = spe_input_dir)
#' prep_fiji_image(sample_info, out_dir = spe_input_dir)
#' 
#' ########################################################################
#' #   Add array coordinates
#' ########################################################################
#'
#' spe_new <- add_array_coords(spe, sample_info, tempdir())
#'
#' #    Several columns related to spatial coordinates were added
#' added_cols_regex <- "^(array|pxl)_(row|col)(_in_fullres)?_(original|rounded)$"
#' colnames(SummarizedExperiment::colData(spe_new))[
#'     grep(added_cols_regex, colnames(SummarizedExperiment::colData(spe_new)))
#' ]
#'
#' #    'array_row', 'array_col', and spatialCoords() were overwritten with
#' #    their transformed values
#' head(spe$array_row)
#' head(spe$array_col)
#' head(SpatialExperiment::spatialCoords(spe_new))
add_array_coords <- function(spe, sample_info, coords_dir, calc_error_metrics = FALSE) {
    ## For R CMD check
    key <- in_tissue <- NULL

    #   State assumptions about columns expected to be in sample_info
    expected_cols <- c("capture_area", "group")
    if (!all(expected_cols %in% colnames(sample_info))) {
        stop(
            sprintf(
                'Missing at least one of the following columns in "sample_info": "%s"',
                paste(expected_cols, collapse = '", "')
            )
        )
    }

    #   55-micrometer diameter for Visium spot; 100 micrometers between spots;
    #   65-micrometer spot diameter used in 'spot_diameter_fullres' calculation
    #   for spaceranger JSON. See documentation for respective quantities. The
    #   difference between 55 and 65 does indeed exist and is properly
    #   documented, but is likely a bug in the sense the choice was probably
    #   unintentional
    #   https://kb.10xgenomics.com/hc/en-us/articles/360035487812-What-is-the-size-of-the-spots-on-the-Visium-Gene-Expression-Slide-
    #   https://kb.10xgenomics.com/hc/en-us/articles/360035487892-How-much-space-is-there-between-spots-referred-to-as-white-space-
    #   https://support.10xgenomics.com/spatial-gene-expression/software/pipelines/latest/output/spatial
    SPOT_DIAMETER_JSON_M <- 65e-6
    INTER_SPOT_DIST_M <- 100e-6

    all_groups <- unique(sample_info$group)

    #   Read in tissue positions for all groups, track group and capture area,
    #   then compute adjusted array coordinates
    coords_list <- list()
    for (i in seq(length(all_groups))) {
        coords <- file.path(
            coords_dir, all_groups[i], "tissue_positions.csv"
        ) |>
            readr::read_csv(show_col_types = FALSE, progress = FALSE)

        #   From the spaceranger JSON, we have the spot diameter both in pixels
        #   and meters, and can therefore compute the image's pixel/m ratio.
        #   Then use that to compute the distance between spots in pixels
        sr_json <- rjson::fromJSON(
            file = file.path(
                coords_dir, all_groups[i], "scalefactors_json.json"
            )
        )
        px_per_m <- sr_json$spot_diameter_fullres / SPOT_DIAMETER_JSON_M
        inter_spot_dist_px <- INTER_SPOT_DIST_M * px_per_m

        #   Adjust 'array_row' and 'array_col' with values appropriate for the new
        #   coordinate system (a larger Visium grid with equal inter-spot distances)
        coords_list[[i]] <- .fit_to_array(coords, inter_spot_dist_px)

        if (calc_error_metrics) {
            coords_list[[i]] = coords |>
                mutate(
                    capture_area = stringr::str_split_i(key, '^[ACTG]+-1_', 2)
                ) |>
                .add_error_metrics(coords_list[[i]], inter_spot_dist_px)
        }
    }

    coords <- do.call(rbind, coords_list) |>
        #   Coordinates must be integers
        mutate(
            across(
                matches("^(array|pxl)_(row|col)(_in_fullres)?"),
                ~ as.integer(round(.x))
            )
        )

    #   Line up and potentially subset refined coords to those in 'spe'
    match_index <- match(spe$key, coords$key)
    if (any(is.na(match_index))) {
        stop("Unrecognized key(s) in refined coords.")
    }
    coords <- coords[match(spe$key, coords$key), ] |>
        select(-c(key, in_tissue))

    #   Retain "_original" copies of the array coordinates
    for (col_name in c("array_row", "array_col")) {
        spe[[paste0(col_name, "_original")]] <- spe[[col_name]]
    }

    #   Add transformed coordinates and rounded pixel coordinates as columns to
    #   colData. Don't add transformed pixel coordinates, which will become
    #   the spatialCoords()
    colData(spe) <- colData(spe) |>
        as.data.frame() |>
        #   Remove the original array coordinates
        select(-c("array_row", "array_col")) |>
        cbind(
            coords |>
                #   Don't add transformed pixel coordinates to colData
                select(-matches("^pxl_(row|col)_in_fullres$"))
        ) |>
        S4Vectors::DataFrame()

    #   Make transformed coordinates the default in the colData and
    #   spatialCoords. Retain "_original" copies of the coordinates
    for (col_name in c("pxl_row_in_fullres", "pxl_col_in_fullres")) {
        spe[[paste0(col_name, "_original")]] <- spatialCoords(spe)[, col_name]
        spatialCoords(spe)[, col_name] <- coords[[col_name]]
    }

    return(spe)
}
