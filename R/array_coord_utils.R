#   Helper functions for add_array_coords()

#' Check if coordinates are Visium-like
#'
#' Sanity check designed to catch unforeseen bugs: halt if the tibble-like 
#' \code{coords}, expected to contain columns 'array_row' and 'array_col',
#' represents an invalid Visium array
#'
#' @param coords A \code{tibble} containing 'array_row' and 'array_col' columns
#' calculated internally by \code{add_array_coords()}
#'
#' @return NULL
#'
#' @importFrom dplyr filter summarize pull
#' @author Nicholas J. Eagles
validate_array <- function(coords) {
    #   Even array rows can only use even column indices
    all_even_cols = coords |>
        dplyr::filter(array_row %% 2 == 0) |>
        dplyr::summarize(a = all(array_col %% 2 == 0)) |>
        dplyr::pull(a)
    if (!all_even_cols) {
        stop("Internal bug: failed to produce an array with even-indexed columns for all even rows!")
    }

    #   Odd array rows can only use odd column indices
    all_odd_rows = coords |>
        dplyr::filter(array_row %% 2 == 1) |>
        dplyr::summarize(a = all(array_col %% 2 == 1)) |>
        dplyr::pull(a)
    if (!all_odd_rows) {
        stop("Internal bug: failed to produce an array with odd-indexed columns for all odd rows!")
    }

    #   Check lower bound of array row and col (note we're allowing arbitrary
    #   maximum values rather than the convention of 78 rows and 128 columns)
    if (!(min(coords$array_row) %in% c(0, 1)) || !(min(coords$array_col) %in% c(0, 1))) {
        stop("Internal bug: failed to produce an array starting at index 0 or 1!")
    }

    #   Check an eccentric detail of Visium arrays: (0, 0) cannot exist
    if (any((coords$array_row == 0) & (coords$array_col == 0))) {
        stop("Internal bug: the invalid array coordinate (0, 0) exists after fitting")
    }
}
