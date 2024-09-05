#' Get keys of neighboring spots
#'
#' For a given row of a `tibble()` containing array coordinates, find the
#' associated spot's neighbors (belonging to the same capture area) and return
#' their keys.
#'
#' @param i An `integer(1)` giving a row index in `coords`.
#' @param coords A `tibble()` containing `array_row`, `array_col`, `key`, and
#' `capture_area` columns.
#'
#' @return A `character()` of neighboring spot keys.
#'
#' @author Nicholas J. Eagles
#' @keywords internal
.get_neighbors <- function(i, coords) {
    this_array_row <- coords[[i, "array_row"]]
    this_array_col <- coords[[i, "array_col"]]
    this_capture_area <- coords[[i, "capture_area"]]

    #   Grab the keys of the neighboring spots and return
    neighbor_keys <- coords |>
        filter(
            capture_area == this_capture_area,
            (array_row == this_array_row & array_col == this_array_col - 2) |
                (array_row == this_array_row & array_col == this_array_col + 2) |
                (array_row == this_array_row - 1 & array_col == this_array_col - 1) |
                (array_row == this_array_row - 1 & array_col == this_array_col + 1) |
                (array_row == this_array_row + 1 & array_col == this_array_col - 1) |
                (array_row == this_array_row + 1 & array_col == this_array_col + 1)
        ) |>
        pull(key)
    return(neighbor_keys)
}

#' Calculate fraction of neighbors retained after mapping to new array
#' coordinates
#'
#' Given `tibble()`s before and after mapping to new array coordinates,
#' calculate for each spot the fraction of starting neighboring spots that
#' were retained in the new array-coordinate system. Add this metric and
#' return.
#'
#' @param coords_new A `tibble()` containing `array_row`, `array_col`, `key`,
#' and `capture_area` columns, representing data after mapping to new array
#' coordinates.
#' @param coords A `tibble()` containing `array_row`, `array_col`, `key`,
#' and `capture_area` columns, representing data before mapping to new array
#' coordinates.
#'
#' @return A `tibble()` copy of `coords_new` with additional `shared_neighbors`
#' column.
#'
#' @author Nicholas J. Eagles
#' @keywords internal
.get_shared_neighbors <- function(coords_new, coords) {
    coords_new$shared_neighbors <- vapply(
        seq_len(nrow(coords)),
        function(i) {
            n_before <- .get_neighbors(i, coords)
            n_after <- .get_neighbors(i, coords_new)
            return(mean(n_before %in% n_after))
        },
        numeric(1)
    )

    return(coords_new)
}

#' Add error metrics related to array-coordinate mapping
#'
#' Given `tibble()`s before and after mapping to new array coordinates,
#' calculate metrics related to the suitability of the mapping.
#'
#' Add column `shared_neighbors`, the fraction of neighbors a spot started
#' with that are retained after mapping; add column `euclidean_error`, the
#' number of multiples of the inter-spot distance a spot must move to be
#' placed in the new array coordinates.
#'
#' @param coords_new A `tibble()` containing `array_row`, `array_col`, `key`,
#' `pxl_col_in_fullres`, `pxl_row_in_fullres`, `pxl_col_in_fullres_rounded`,
#' and `pxl_row_in_fullres_rounded` columns, representing data after mapping
#' to new array coordinates for one `group`.
#' @param coords A `tibble()` containing `array_row`, `array_col`, `key`,
#' `pxl_col_in_fullres`, `pxl_row_in_fullres`, `pxl_col_in_fullres_rounded`,
#' `pxl_row_in_fullres_rounded`, and `capture_area` columns, representing data
#' before mapping to new array coordinates for one `group`.
#' @param inter_spot_dist_px A `numeric(1)` giving the number of pixels between
#' spots for the `group`.
#'
#' @return A `tibble()` copy of `coords_new` with additional `shared_neighbors`
#' and `euclidean_error` columns.
#'
#' @importFrom stringr str_split_i
#'
#' @author Nicholas J. Eagles
#' @keywords internal
.add_error_metrics <- function(coords, coords_new, inter_spot_dist_px) {
    coords_new <- coords_new |>
        mutate(
            euclidean_error = (
                (pxl_col_in_fullres - pxl_col_in_fullres_rounded)**2 +
                    (pxl_row_in_fullres - pxl_row_in_fullres_rounded)**2
            )**0.5 / inter_spot_dist_px,
            capture_area = stringr::str_split_i(key, "^[ACTG]+-1_", 2),
        ) |>
        .get_shared_neighbors(coords) |>
        select(-capture_area)

    return(coords_new)
}
