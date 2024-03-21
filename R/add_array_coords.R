#' Add transformed array and pixel coordinates to a \code{SpatialExperiment}
#'
#' Given a \code{SpatialExperiment}, sample information, and coordinates
#' produced from the Samui refinement workflow, add array and pixel coordinates
#' appropriate for the linearly transformed capture areas making up each group
#' present in the \code{SpatialExperiment}.
#' 
#' Array coordinates are determined via an algorithm that fits each spot to
#' the nearest spot on a new, imaginary, Visium-like capture area. The imaginary
#' capture area differs from a real capture area only in its extent; array
#' coordinates still start at 0 but may extend arbitrarily beyond the normal
#' maximum indices of 77 and 127 to fit every capture area in each group
#' defined in the \code{SpatialExperiment}. The goal is to return well-defined
#' array coordinates in a consistent spatial orientation for each group, such
#' that downstream applications, such as clustering with BayesSpace, can
#' process each group as if it really were one capture area in the first place.
#'
#' @param spe A \code{SpatialExperiment}
#' @param sample_info A \code{tibble} with columns \code{capture_area},
#' \code{group}, and \code{spaceranger_dir}
#' @param coords_dir A \code{character(1)} vector giving the directory
#' containing \code{tissue_positions_[group].csv} files produced from refinement
#' with Samui
#' @param overwrite A \code{logical(1)} vector indicating whether to overwrite
#' \code{spatialCoords(spe)}, and \code{colData(spe)} columns \code{array_row},
#' \code{array_col}, \code{pixel_row_in_fullres}, and
#' \code{pixel_col_in_fullres} with the transformed values. Note that the
#' original values are preserved even when TRUE, in versions of
#' \code{colData(spe)} columns ending in \code{_original}.
#'
#' @return A \code{SpatialExperiment} object with modified \code{colData}
#' columns \code{array_row}, \code{array_col}, \code{pixel_row_in_fullres}, and
#' \code{pixel_col_in_fullres}, and additional corresponding columns ending in
#' \code{_original}
#'
#' @export
#' @author Nicholas J. Eagles

add_array_coords = function(spe, sample_info, coords_dir, overwrite = TRUE) {
    all_groups = unique(sample_info$group)

    #   Read in tissue positions for all groups, track group and capture area,
    #   then merge into a single tibble
    tissue_list = list()
    for (i in seq(length(all_groups))) {
        tissue_list[[i]] = file.path(
                coords_dir,
                sprintf('tissue_positions_%s.csv', all_groups[i])
            ) |>
            read_csv(show_col_types = FALSE) |>
            mutate(
                group = all_groups[i],
                capture_area = str_extract(key, '_(.*)', group = 1)
            )
    }
    coords = do.call(rbind, tissue_list)
}
