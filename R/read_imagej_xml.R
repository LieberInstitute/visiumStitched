#' Parse transform info from ImageJ XML output
#'
#' Given a \code{tibble} of sample information (\code{sample_info}) with
#' columns \code{capture_area}, \code{group}, and \code{imagej_out_path},
#' expected to have one unique path to ImageJ XML output per group, return
#' a copy of \code{sample_info} with additional columns \code{transform_x},
#' \code{transform_y}, and \code{transform_theta} representing the rigid affine
#' transforms needed to be applied to capture areas within each group to
#' produce the proper relative arrangements.
#'
#' @param sample_info A \code{tibble} with columns \code{capture_area},
#' \code{group}, and \code{imagej_out_path}
#'
#' @return A \code{tibble}: a copy of \code{sample_info} with additional columns
#' \code{transform_x}, \code{transform_y}, and \code{transform_theta}
#'
#' @import xml2
#' @importFrom dplyr mutate filter
#' @importFrom stringr str_replace_all str_split_i
#' 
#' @export
#' @author Nicholas J. Eagles
#'
#' @examples
#' #   For internal testing
read_imagej_xml <- function(sample_info) {
    all_groups = unique(sample_info$group)
    new_sample_info_list = list()
    for (this_group in all_groups) {
        this_sample_info = sample_info |>
            dplyr::filter(group == this_group)

        if (length(unique(this_sample_info$imagej_out_path)) > 1) {
            stop("Expected one unique path for 'imagej_out_path' per group in 'sample_info'.")
        }

        transform_nodes = this_sample_info$imagej_out_path[1] |>
            read_xml() |>
            xml_find_all('.//t2_patch')
        
        input_paths = xml_attr(transform_nodes, 'file_path')
        input_indices = sapply(
            this_sample_info$capture_area, function(x) grep(x, input_paths)
        )
        if (length(input_paths) != nrow(this_sample_info) || any(is.na(input_indices))) {
            stop("Expected each capture area to be present exactly once in the input filenames to ImageJ for each group.")
        }
        transforms = xml_attr(transform_nodes, 'transform')[input_indices]

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
            )
    }

    return(do.call(rbind, new_sample_info_list))
}
