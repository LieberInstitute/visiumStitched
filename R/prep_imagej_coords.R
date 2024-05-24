#' Apply transform info from ImageJ XML output
#'
#' Given a \code{tibble} of sample information (\code{sample_info}) with
#' columns \code{capture_area}, \code{group}, and \code{imagej_xml_path},
#' expected to have one unique path to ImageJ XML output per group, read in
#' the pixel coordinates from each capture area's \code{tissue_positions.csv}
#' file from Spaceranger, and transform using the rotation matrix specified
#' by ImageJ. Write one new \code{tissue_positions.csv} file per group
#' a copy of \code{sample_info} with additional columns \code{transform_x},
#' \code{transform_y}, and \code{transform_theta} representing the rigid affine
#' transforms needed to be applied to capture areas within each group to
#' produce the proper relative arrangements.
#'
#' @param sample_info A \code{tibble} with columns \code{capture_area},
#' \code{group}, and \code{imagej_xml_path}
#' @param out_dir A character(1) vector giving a path to a directory to place
#' the output pixel coordinates CSVs. Provided the parent exists, \code{out_dir}
#' will be created if necessary.
#'
#' @return NULL
#'
#' @import xml2
#' @importFrom dplyr mutate filter select
#' @importFrom stringr str_replace_all
#' @importFrom readr write_csv
#' @importFrom tibble as_tibble
#' 
#' @export
#' @author Nicholas J. Eagles
#'
#' @examples
#' \dontrun{
#' #   For internal testing
#' sample_info <- readr::read_csv("dev/test_data/sample_info.csv")
#' prep_imagej_coords(sample_info) |>
#'     dplyr::select(tidyselect::starts_with('transform')) |>
#'     print()
#' }
#'
#' ## TODO: add working examples
#' args(prep_imagej_coords)
prep_imagej_coords <- function(sample_info, out_dir) {
    TISSUE_COLNAMES <- c(
        "barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres",
        "pxl_col_in_fullres"
    )

    if (!all(file.exists(sample_info$imagej_xml_path))) {
        stop("All files in 'sample_info$imagej_xml_path' must exist.")
    }
    
    new_sample_info_list = list()
    for (this_group in unique(sample_info$group)) {
        this_sample_info = sample_info |>
            dplyr::filter(group == this_group)

        if (length(unique(this_sample_info$imagej_xml_path)) > 1) {
            stop("Expected one unique path for 'imagej_xml_path' per group in 'sample_info'.")
        }

        #   Find all XML elements containing input image paths and
        #   transformation matrices
        transform_nodes = this_sample_info$imagej_xml_path[1] |>
            read_xml() |>
            xml_find_all('.//t2_patch')
        
        #   Find paths to input images and the order the corresponding capture
        #   areas appear in these paths
        input_paths = xml_attr(transform_nodes, 'file_path')
        input_indices = sapply(
            this_sample_info$capture_area, function(x) grep(x, input_paths)
        )
        if (length(input_paths) != nrow(this_sample_info) || any(is.na(input_indices))) {
            stop("Expected each capture area to be present exactly once in the input filenames to ImageJ for each group.")
        }

        coords_list = list()

        #   Loop through all capture areas in this group
        for (i in seq(nrow(this_sample_info))) {
            #   Parse the rotation matrix from the ImageJ XML
            rot = transform_nodes[input_indices[i]] |>
                xml_attr('transform') |>
                stringr::str_replace_all('matrix|[\\(\\)]', '') |>
                strsplit(',') |>
                unlist() |>
                as.numeric() |>
                matrix(nrow = 2, ncol = 3)
            
            #   Read in the raw tissue positions for this capture area
            coords = list.files(
                    this_sample_info$spaceranger_dir[i],
                    '^tissue_positions(_list)?\\.csv$',
                    full.names = TRUE
                )[1] |>
                read.csv(col.names = TISSUE_COLNAMES) |>
                tibble::as_tibble() |>
                dplyr::mutate(
                    key = paste(
                        barcode, this_sample_info$capture_area[i], sep = '_'
                    )
                )
            
            #   Take just the x and y coords, and apply the rotation matrix
            coords_xy = coords |>
                dplyr::select(pxl_row_in_fullres, pxl_col_in_fullres) |>
                dplyr::mutate(ones = 1) |>
                as.matrix()
            coords_xy = t(rot %*% t(coords_xy))

            coords_list[[i]] = coords |>
                dplyr::mutate(
                    pxl_row_in_fullres = coords_xy[,1],
                    pxl_col_in_fullres = coords_xy[,2]
                )
        }

        #   Merge coordinates for all capture areas in this group and write to
        #   CSV
        coords = do.call(rbind, coords_list)
        this_out_dir = file.path(out_dir, this_group)
        dir.create(this_out_dir, showWarnings = FALSE)
        readr::write_csv(coords, file.path(this_out_dir, 'tissue_positions.csv'))
    }

    return(NULL)
}
