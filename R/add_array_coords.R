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
#' @importFrom readr read_csv
#' @importFrom S4Vectors DataFrame
#' @importFrom rjson fromJSON
#' @export
#' @author Nicholas J. Eagles
#' 
#' @examples
#' #   For internal testing
#' \dontrun{
#'    library(HDF5Array)
#'    spe = loadHDF5SummarizedExperiment('dev/test_data/spe_filtered')
#'    sample_info = readr::read_csv('dev/test_data/sample_info.csv')
#'    coords_dir = 'dev/test_data'
#'    spe_new = add_array_coords(spe, sample_info, coords_dir, overwrite = TRUE)
#'}

add_array_coords = function(spe, sample_info, coords_dir, overwrite = TRUE) {
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

    all_groups = unique(sample_info$group)

    #   Read in tissue positions for all groups, track group and capture area,
    #   then compute adjusted array coordinates
    coords_list = list()
    for (i in seq(length(all_groups))) {
        coords = file.path(
                coords_dir,
                sprintf('tissue_positions_%s.csv', all_groups[i])
            ) |>
            readr::read_csv(show_col_types = FALSE) |>
            mutate(
                group = all_groups[i],
                capture_area = str_extract(key, '_(.*)', group = 1)
            )
        
        #   From the spaceranger JSON, we have the spot diameter both in pixels
        #   and meters, and can therefore compute the image's pixel/m ratio.
        #   Then use that to compute the distance between spots in pixels
        sr_json <- rjson::fromJSON(
            file = file.path(
                sample_info$spaceranger_dir[i], "scalefactors_json.json"
            )
        )
        px_per_m <- sr_json$spot_diameter_fullres / SPOT_DIAMETER_JSON_M
        inter_spot_dist_px <- INTER_SPOT_DIST_M * px_per_m

        #   Adjust 'array_row' and 'array_col' with values appropriate for the new
        #   coordinate system (a larger Visium grid with equal inter-spot distances)
        coords_list[[i]] <- .fit_to_array(coords, inter_spot_dist_px)
    }

    coord_cols = c(
        "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"
    )
    coords = do.call(rbind, coords_list) |>
        #   Coordinates must be integers
        mutate_at(coord_cols, ~ as.integer(round(.))) |>
        rename_with(coord_cols, ~ paste0(.x, '_transformed'))

    #   Add transformed coordinates as columns to colData
    temp = colnames(spe)
    colData(spe) = colData(spe) |>
        as_tibble() |>
        left_join(coords, by = "key") |>
        DataFrame()
    colnames(spe) = temp

    #   If 'overwrite', make transformed coordinates the default in the colData
    #   and spatialCoords
    if (overwrite) {
        spatialCoords(spe)$pxl_col_in_fullres = coords$pxl_col_in_fullres_transformed
        spatialCoords(spe)$pxl_row_in_fullres = coords$pxl_row_in_fullres_transformed

        #   Make transformed coordinates the default in the colData
        for(col_name in coord_cols) {
            spe[[col_name]] = coords[[paste0(col_name, "_transformed")]]
        }
    }

    return(spe)
}
