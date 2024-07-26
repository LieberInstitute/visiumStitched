#' Apply transform info from Fiji XML output
#'
#' Given a `data.frame()` of sample information (\code{sample_info}) with
#' columns \code{capture_area}, \code{group}, and \code{fiji_xml_path},
#' expected to have one unique path to Fiji XML output per group, read in
#' the pixel coordinates from each capture area's \code{tissue_positions.csv}
#' file from SpaceRanger, and transform using the rotation matrix specified
#' by Fiji <https://imagej.net/software/fiji/>.
#' Write one new \code{tissue_positions.csv} file per group.
#'
#' @param out_dir A \code{character(1)} vector giving a path to a directory to
#' place the output pixel coordinates CSVs. Provided the parent directory
#' exists, \code{out_dir} will be created if necessary.
#' @inheritParams add_array_coords
#'
#' @return This function returns a `character()` with the file paths to the
#' `tissue_positions.csv` files it created.
#'
#' @import xml2
#' @importFrom stringr str_replace_all str_detect
#' @importFrom readr read_csv write_csv
#' @importFrom rjson fromJSON
#' @importFrom pkgcond suppress_warnings
#'
#' @family functions for parsing Fiji outputs
#'
#' @export
#' @author Nicholas J. Eagles
#'
#' @examples
#' #    Define sample information for the example human brain data
#' sample_info <- dplyr::tibble(
#'     group = "Br2719",
#'     capture_area = c("V13B23-283_A1", "V13B23-283_C1", "V13B23-283_D1")
#' )
#' #   Add 'spaceranger_dir' column
#' sr_dir <- tempdir()
#' temp <- unzip(
#'     spatialLIBD::fetch_data("visiumStitched_brain_spaceranger"),
#'     exdir = sr_dir
#' )
#' sample_info$spaceranger_dir <- file.path(
#'     sr_dir, sample_info$capture_area, "outs", "spatial"
#' )
#'
#' #   Add Fiji-output-related columns
#' fiji_dir <- tempdir()
#' temp <- unzip(
#'     spatialLIBD::fetch_data("visiumStitched_brain_Fiji_out"),
#'     exdir = fiji_dir
#' )
#' sample_info$fiji_xml_path <- temp[grep("xml$", temp)]
#' sample_info$fiji_image_path <- temp[grep("png$", temp)]
#'
#' ## Re-size images and add more information to the sample_info
#' sample_info <- rescale_fiji_inputs(sample_info, out_dir = tempdir())
#'
#' spe_input_dir <- tempdir()
#' out_file <- prep_fiji_coords(sample_info, out_dir = spe_input_dir)
#' out_file
#'
#' #    A file of spatial coordinates for the stitched Br2719 was produced
#' readr::read_csv(out_file)
prep_fiji_coords <- function(sample_info, out_dir) {
    ## For R CMD check
    group <- barcode <- key <- pxl_col_in_fullres <- pxl_row_in_fullres <- NULL

    TISSUE_COLNAMES <- c(
        "barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres",
        "pxl_col_in_fullres"
    )

    #   State assumptions about columns expected to be in sample_info
    expected_cols <- c(
        "capture_area", "group", "fiji_xml_path", "intra_group_scalar",
        "group_hires_scalef"
    )
    if (!all(expected_cols %in% colnames(sample_info))) {
        stop(
            sprintf(
                'Missing at least one of the following columns in "sample_info": "%s"',
                paste(expected_cols, collapse = '", "')
            )
        )
    }

    if (!all(file.exists(sample_info$fiji_xml_path))) {
        stop("All files in 'sample_info$fiji_xml_path' must exist.")
    }

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    out_paths <- list()
    for (this_group in unique(sample_info$group)) {
        this_sample_info <- sample_info |>
            dplyr::filter(group == this_group)

        if (length(unique(this_sample_info$fiji_xml_path)) > 1) {
            stop("Expected one unique path for 'fiji_xml_path' per group in 'sample_info'.")
        }

        #   Find all XML elements containing input image paths and
        #   transformation matrices
        transform_nodes <- this_sample_info$fiji_xml_path[1] |>
            read_xml() |>
            pkgcond::suppress_warnings(pattern = "Attribute o_width of element t2_patch") |>
            pkgcond::suppress_warnings(pattern = "Attribute o_height of element t2_patch") |>
            xml_find_all(".//t2_patch")

        #   Find paths to input images and the order the corresponding capture
        #   areas appear in these paths
        input_paths <- xml_attr(transform_nodes, "file_path")
        input_indices <- vapply(
            this_sample_info$capture_area, function(x) grep(x, input_paths), integer(1)
        )
        if (length(input_paths) != nrow(this_sample_info) || any(is.na(input_indices))) {
            stop("Expected each capture area to be present exactly once in the input filenames to Fiji for each group.")
        }

        coords_list <- list()

        #   Loop through all capture areas in this group
        for (i in seq(nrow(this_sample_info))) {
            #   Parse the rotation matrix from the Fiji XML, and scale
            #   translations from high to fullres
            rot <- transform_nodes[input_indices[i]] |>
                xml_attr("transform") |>
                stringr::str_replace_all("matrix|[\\(\\)]", "") |>
                strsplit(",") |>
                unlist() |>
                as.numeric() |>
                matrix(nrow = 2, ncol = 3)
            rot[, 3] <- rot[, 3] / this_sample_info$group_hires_scalef[1]

            #   Read in the raw tissue positions for this capture area, handling
            #   the old and new format for tissue positions
            coords_path <- list.files(
                this_sample_info$spaceranger_dir[i],
                "^tissue_positions(_list)?\\.csv$",
                full.names = TRUE
            )[1]
            if (stringr::str_detect(coords_path, "tissue_positions_list\\.csv$")) {
                coords <- read_csv(
                    coords_path,
                    col_names = FALSE, show_col_types = FALSE,
                    progress = FALSE
                )
                colnames(coords) <- TISSUE_COLNAMES
            } else {
                coords <- read_csv(
                    coords_path,
                    show_col_types = FALSE, progress = FALSE
                )
            }

            coords <- coords |>
                dplyr::rename(key = barcode) |>
                dplyr::mutate(
                    key = paste(
                        key, this_sample_info$capture_area[i],
                        sep = "_"
                    )
                )

            #   Take just the x and y coords, scale to match the rest of the
            #   group, and apply the rotation matrix
            coords_xy <- coords |>
                dplyr::select(pxl_col_in_fullres, pxl_row_in_fullres) |>
                dplyr::mutate(
                    pxl_col_in_fullres = this_sample_info$intra_group_scalar[i] *
                        pxl_col_in_fullres,
                    pxl_row_in_fullres = this_sample_info$intra_group_scalar[i] *
                        pxl_row_in_fullres,
                    ones = 1
                ) |>
                as.matrix()
            coords_xy <- t(rot %*% t(coords_xy))

            coords_list[[i]] <- coords |>
                dplyr::mutate(
                    pxl_col_in_fullres = coords_xy[, 1],
                    pxl_row_in_fullres = coords_xy[, 2]
                )
        }

        #   Merge coordinates for all capture areas in this group and write to
        #   CSV
        coords <- do.call(rbind, coords_list)
        this_out_dir <- file.path(out_dir, this_group)
        dir.create(this_out_dir, showWarnings = FALSE)
        out_paths[[i]] <- file.path(this_out_dir, "tissue_positions.csv")
        readr::write_csv(
            coords, out_paths[[i]],
            progress = FALSE
        )
    }

    return(unlist(out_paths))
}
