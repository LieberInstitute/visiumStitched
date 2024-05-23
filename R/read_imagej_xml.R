#' Parse transform info from ImageJ XML output
#'
#' Given a \code{tibble} of sample information (\code{sample_info}) with
#' columns \code{capture_area}, \code{group}, and \code{imagej_xml_path},
#' expected to have one unique path to ImageJ XML output per group, return
#' a copy of \code{sample_info} with additional columns \code{transform_x},
#' \code{transform_y}, and \code{transform_theta} representing the rigid affine
#' transforms needed to be applied to capture areas within each group to
#' produce the proper relative arrangements.
#'
#' @param sample_info A \code{tibble} with columns \code{capture_area},
#' \code{group}, and \code{imagej_xml_path}
#'
#' @return A \code{tibble}: a copy of \code{sample_info} with additional columns
#' \code{transform_x}, \code{transform_y}, and \code{transform_theta}
#'
#' @import xml2
#' @importFrom dplyr mutate filter select
#' @importFrom stringr str_replace_all str_split_i
#' 
#' @export
#' @author Nicholas J. Eagles
#'
#' @examples
#' \dontrun{
#' #   For internal testing
#' sample_info <- readr::read_csv("dev/test_data/sample_info.csv")
#' read_imagej_xml(sample_info) |>
#'     dplyr::select(tidyselect::starts_with('transform')) |>
#'     print()
#' }
#'
#' ## TODO: add working examples
#' args(read_imagej_xml)
read_imagej_xml <- function(sample_info) {
    if (!all(file.exists(sample_info$imagej_xml_path))) {
        stop("All files in 'sample_info$imagej_xml_path' must exist.")
    }
    
    all_groups = unique(sample_info$group)
    new_sample_info_list = list()
    for (this_group in all_groups) {
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

        #   Specifically, the linear transform matrix ordered by capture area as
        #   it appears in 'this_sample_info'
        transforms = xml_attr(transform_nodes, 'transform')[input_indices]

        #   Add 'transform_x', 'transform_y' and 'transform_theta' to the sample
        #   info for this group
        new_sample_info_list[[this_group]] = this_sample_info |>
            dplyr::mutate(
                transform_x = transforms |>
                    stringr::str_replace_all('matrix|[\\(\\)]', '') |>
                    stringr::str_split_i(',', 5) |>
                    as.numeric(),
                transform_y = transforms |>
                    stringr::str_replace_all('matrix|[\\(\\)]', '') |>
                    stringr::str_split_i(',', 6) |>
                    as.numeric(),
                #   Take arccos of the [0, 0] element of the rotation matrix
                #   to determine theta
                transform_theta = transforms |>
                    stringr::str_replace_all('matrix|[\\(\\)]', '') |>
                    stringr::str_split_i(',', 1) |>
                    as.numeric() |>
                    acos(),
                #   Since arccos is positive by convention, multiply by -1 if
                #   element [1, 0] of the rotation matrix (i.e. sin(theta)) is
                #   negative
                multiplier = transforms |>
                    stringr::str_replace_all('matrix|[\\(\\)]', '') |>
                    stringr::str_split_i(',', 2) |>
                    as.numeric() |>
                    sign(),
                transform_theta = transform_theta * multiplier
            ) |>
            dplyr::select(-multiplier)
    }

    return(do.call(rbind, new_sample_info_list))
}
