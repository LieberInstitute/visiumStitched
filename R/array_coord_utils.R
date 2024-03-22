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
.validate_array <- function(coords) {
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

#' Return array coordinates fit to nearest spot with associated error
#'
#' First, values of \code{x} are rounded to the nearest integer. Then, values
#' of \code{y} are rounded to the nearest valid integer under the constraint
#' that coordinates for x and y must be both odd or both even. These rounded
#' values are returned, along with the Euclidean distance needed to move x and
#' y from their original, non-integer values to their rounded values.
#'
#' @param x \code{numeric()} vector giving "ideal" array coordinates given every
#' spot's transformed pixel coordinates.
#' @param y Same as x, though y must represent ideal array columns iff x
#' represents array rows, and vice versa.
#' @param INTERVAL_X \code{numeric(1)} giving pixel distance between coordinate
#' units used for \code{x} (e.g. if x represents ideal \code{array_col} values,
#' \code{INTERVAL_X} represents pixel distance between spot columns).
#' @param INTERVAL_Y \code{numeric(1)} giving pixel distance between coordinate
#' units used for \code{y}.
#'
#' @return A \code{list} consisting of 3 unnamed \code{numeric()} vectors:
#' rounded \code{x}, rounded \code{y}, and the Euclidean distance in pixels from
#' rounding both \code{x} and \code{y}.
#'
#' @author Nicholas J. Eagles
.refine_fit <- function(x, y, INTERVAL_X, INTERVAL_Y) {
    #   Round x to the nearest integer, and track the error from doing so in the
    #   variable 'dx'
    dx <- x - clean_round(x)
    x <- x - dx

    #   Given x, round y to the nearest valid integer (y must be even iff x is),
    #   and track the error from doing so in the variable 'dy'
    dy <- rep(0, length(y))
    dy[x %% 2 == 0] <- y[x %% 2 == 0] - clean_round(y[x %% 2 == 0] / 2) * 2
    dy[x %% 2 == 1] <- y[x %% 2 == 1] - (clean_round(y[x %% 2 == 1] / 2 - 0.5) * 2 + 1)
    y <- y - dy

    #   Summarize error in Euclidean distance
    error <- sqrt((INTERVAL_X * dx)**2 + (INTERVAL_Y * dy)**2)
    return(list(x, y, error))
}

#' Round to the nearest integer, always rounding up at 0.5
#' 
#' This consistent behavior is favorable for our application, where we want to
#' minimize duplicate mappings of spots to new array coordinates
#' 
#' @param x \code{numeric()} vector.
#'
#' @return A \code{numeric()} vector rounded to the nearest integer.
#'
#' @author Nicholas J. Eagles
clean_round <- function(x) {
    return(floor(x) + ((x * 10) %% 10 >= 5))
}
