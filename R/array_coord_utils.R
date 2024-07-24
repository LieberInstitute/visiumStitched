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
#' @return It returns `NULL` if all tests were correct.
#'
#' @author Nicholas J. Eagles
#' @keywords internal
.validate_array <- function(coords) {
    ## For R CMD check
    array_row <- array_col <- tmp_array <- NULL

    #   Even array rows can only use even column indices
    all_even_cols <- coords |>
        dplyr::filter(array_row %% 2 == 0) |>
        dplyr::summarize(tmp_array = all(array_col %% 2 == 0)) |>
        dplyr::pull(tmp_array)
    if (!all_even_cols) {
        stop("Internal bug: failed to produce an array with even-indexed columns for all even rows!")
    }

    #   Odd array rows can only use odd column indices
    all_odd_rows <- coords |>
        dplyr::filter(array_row %% 2 == 1) |>
        dplyr::summarize(tmp_array = all(array_col %% 2 == 1)) |>
        dplyr::pull(tmp_array)
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
    return(invisible(NULL))
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
#' @keywords internal
.refine_fit <- function(x, y, INTERVAL_X, INTERVAL_Y) {
    #   Round x to the nearest integer, and track the error from doing so in the
    #   variable 'dx'
    dx <- x - .clean_round(x)
    x <- x - dx

    #   Given x, round y to the nearest valid integer (y must be even iff x is),
    #   and track the error from doing so in the variable 'dy'
    dy <- rep(0, length(y))
    dy[x %% 2 == 0] <- y[x %% 2 == 0] - .clean_round(y[x %% 2 == 0] / 2) * 2
    dy[x %% 2 == 1] <- y[x %% 2 == 1] - (.clean_round(y[x %% 2 == 1] / 2 - 0.5) * 2 + 1)
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
#' @keywords internal
.clean_round <- function(x) {
    return(floor(x) + ((x * 10) %% 10 >= 5))
}

#' Fit spots to a new Visium-like array
#'
#' Given transformed pixel coordinates, modify the 'array_row' and
#' 'array_col' columns to represent a larger Visium capture area containing
#' all capture areas in a common coordinate system. The number of
#' array rows/cols generally changes from the Visium standards of 78 and 128
#' (and even may change in ratio between num rows and num cols).
#'
#' Runtime is O(n) with the number of spots, making it much faster than say,
#' a distance-matrix-based approach running at O(n^2).
#'
#' @param coords A \code{tibble()} whose rows represent capture areas of the
#' same group, and containing columns 'array_row', 'array_col',
#' 'pxl_row_in_fullres', and 'pxl_col_in_fullres'.
#' @param inter_spot_dist_px \code{numeric(1)} vector giving the pixel distance
#' between any 2 spots in the new coordinates.
#'
#' @return A \code{tibble()} with modified 'array_row' + 'array_col' columns, as
#' well as new 'pxl_row_in_fullres_rounded' and 'pxl_col_in_fullres_rounded'
#' columns representing the pixel coordinates rounded to the nearest exact array
#' coordinates.
#'
#' @author Nicholas J. Eagles
#' @keywords internal
.fit_to_array <- function(coords, inter_spot_dist_px) {
    ## For R CMD check
    array_row <- array_col <- pxl_col_in_fullres <- pxl_row_in_fullres <- NULL

    MIN_ROW <- min(coords$pxl_col_in_fullres)
    INTERVAL_ROW <- inter_spot_dist_px * cos(pi / 6)

    MIN_COL <- min(coords$pxl_row_in_fullres)
    MAX_COL <- max(coords$pxl_row_in_fullres)
    INTERVAL_COL <- inter_spot_dist_px / 2

    #   Calculate what "ideal" array rows and cols should be (allowed to be any
    #   float). Don't round yet. Note array_row maps with pxl_col, while
    #   array_col maps backwards with pxl_row
    array_row_temp <- (coords$pxl_col_in_fullres - MIN_ROW) /
        INTERVAL_ROW

    array_col_temp <- (MAX_COL - coords$pxl_row_in_fullres) /
        INTERVAL_COL

    #   For now, find the nearest row first, then round to the nearest possible
    #   column given the row
    temp <- .refine_fit(array_row_temp, array_col_temp, INTERVAL_ROW, INTERVAL_COL)
    error_row_first <- temp[[3]]
    coords$array_row <- temp[[1]]
    coords$array_col <- temp[[2]]

    #   Perform the opposite order (column then row). When this ordering results
    #   in lower error, use it instead
    temp <- .refine_fit(array_col_temp, array_row_temp, INTERVAL_COL, INTERVAL_ROW)
    error_col_first <- temp[[3]]
    coords$array_row[error_row_first > error_col_first] <- temp[[2]][
        error_row_first > error_col_first
    ]
    coords$array_col[error_row_first > error_col_first] <- temp[[1]][
        error_row_first > error_col_first
    ]

    #   Now make new pixel columns based on just the array values (these columns
    #   give the coordinates for given array row/cols)
    coords$pxl_col_in_fullres_rounded <- MIN_ROW + coords$array_row * INTERVAL_ROW
    coords$pxl_row_in_fullres_rounded <- MAX_COL - coords$array_col * INTERVAL_COL

    #-------------------------------------------------------------------------------
    #   array (0, 0) does not exist on an ordinary Visium array. Move any such
    #   values to the nearest alternatives
    #-------------------------------------------------------------------------------

    #   Nearest points to (0, 0) are (0, 2) and (1, 1):
    array_02 <- c(MIN_ROW, MIN_COL + 2 * INTERVAL_COL)
    array_11 <- c(MIN_ROW + INTERVAL_ROW, MIN_COL + INTERVAL_COL)

    #   Determine the distances to those nearest points
    dist_coords <- coords |>
        dplyr::filter(array_row == 0, array_col == 0) |>
        dplyr::mutate(
            dist_02 = (pxl_col_in_fullres - array_02[1])**2 +
                (pxl_row_in_fullres - array_02[2])**2,
            dist_11 = (pxl_col_in_fullres - array_11[1])**2 +
                (pxl_row_in_fullres - array_11[2])**2,
        )

    #   Move any instances of (0, 0) to the nearest alternative
    indices <- (coords$array_row == 0) & (coords$array_col == 0)
    coords[indices, "array_row"] <- ifelse(
        dist_coords$dist_02 < dist_coords$dist_11, 0, 1
    )
    coords[indices, "array_col"] <- ifelse(
        dist_coords$dist_02 < dist_coords$dist_11, 2, 1
    )

    #   Verify the newly assigned array row and cols have reasonable values
    .validate_array(coords)

    return(coords)
}
