#   Helper functions for add_array_coords()

#' Check if coordinates are Visium-like
#'
#' Sanity check designed to catch unforeseen bugs: halt if the tibble-like
#' \code{coords}, expected to contain columns 'array_row' and 'array_col',
#' represents an invalid Visium array.
#'
#' @param coords A `data.frame()` containing `'array_row'` and
#' `'array_col'` columns
#' calculated internally by \code{add_array_coords()}.
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
#' minimize duplicate mappings of spots to new array coordinates.
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

#' Construct a new Visium-like array encapsulating a set of spots
#' 
#' Given \code{coords} containing pixel coordinates of spots from potentially
#' multiple capture areas, return a new Visium-like array encapsulating all
#' such spots.
#' 
#' @param coords A `data.frame()` with columns 'pxl_row_in_fullres' and
#' 'pxl_col_in_fullres' whose rows contain spots from potentially multiple
#' capture areas.
#' @param inter_spot_dist_px \code{numeric(1)} vector giving the pixel distance
#' between any 2 spots in the new coordinates.
#' 
#' @return A [tibble][dplyr::reexports] with columns 'array_row', 'array_col',
#' 'pxl_row_in_fullres', and 'pxl_col_in_fullres', representing the new
#' Visium-like array.
#' 
#' @author Nicholas J. Eagles
#' @keywords internal
.construct_array = function(coords, inter_spot_dist_px) {
    ## For R CMD check
    array_row <- array_col <- pxl_col_in_fullres <- pxl_row_in_fullres <- NULL

    MIN_ROW <- min(coords$pxl_col_in_fullres)
    MAX_ROW <- max(coords$pxl_col_in_fullres)
    INTERVAL_ROW <- inter_spot_dist_px * cos(pi / 6)
    NUM_ROWS = (MAX_ROW - MIN_ROW) / INTERVAL_ROW

    MIN_COL <- min(coords$pxl_row_in_fullres)
    MAX_COL <- max(coords$pxl_row_in_fullres)
    INTERVAL_COL <- inter_spot_dist_px / 2
    NUM_COLS = (MAX_COL - MIN_COL) / INTERVAL_COL

    #   First the even rows and even cols
    row_indices = 2 * seq(ceiling(NUM_ROWS / 2)) - 2
    col_indices = 2 * seq(ceiling(NUM_COLS / 2)) - 2
    row_coords = MIN_ROW + row_indices * INTERVAL_ROW
    col_coords = MAX_COL - col_indices * INTERVAL_COL
    new_array = dplyr::tibble(
        array_row = rep(row_indices, times = length(col_coords)),
        pxl_col_in_fullres = rep(row_coords, times = length(col_coords)),
        array_col = rep(col_indices, each = length(row_coords)),
        pxl_row_in_fullres = rep(col_coords, each = length(row_coords))
    )

    #   Next the odd rows and odd cols
    row_indices = 2 * seq(ceiling(NUM_ROWS / 2)) - 1
    col_indices = 2 * seq(ceiling(NUM_COLS / 2)) - 1
    row_coords = MIN_ROW + row_indices * INTERVAL_ROW
    col_coords = MAX_COL - col_indices * INTERVAL_COL
    new_array = rbind(
        new_array,
        dplyr::tibble(
            array_row = rep(row_indices, times = length(col_coords)),
            pxl_col_in_fullres = rep(row_coords, times = length(col_coords)),
            array_col = rep(col_indices, each = length(row_coords)),
            pxl_row_in_fullres = rep(col_coords, each = length(row_coords))
        )
    )
    
    .validate_array(new_array)

    return(new_array)
}

#' Map source spots to best target spots by solving the LSAP
#' 
#' Given \code{source_coords} and \code{target_coords}, both containing pixel
#' coordinates of spots, map each spot in \code{source_coords} to a unique
#' spot in \code{target_coords} such that the total squared Euclidean distance
#' between matched spots is minimized, with guaranteed one-to-one mapping. This
#' is done by solving the Linear Sum Assignment Problem (LSAP) using the
#' Hungarian algorithm. Return the \code{source_coords} with the newly mapped
#' \code{array_row} and \code{array_col} columns.
#' 
#' @param source_coords A `data.frame()` containing the pixel coordinates (i.e.
#' 'pxl_row_in_fullres' and 'pxl_col_in_fullres') of starting spots from one
#' capture area.
#' @param target_coords A `data.frame()` containing the pixel coordinates (i.e.
#' 'pxl_row_in_fullres' and 'pxl_col_in_fullres') of target spots which should
#' just barely encompass the capture area in \code{source_coords}.
#' 
#' @return A [tibble][dplyr::reexports] with the same rows as \code{source_coords},
#' but with the \code{array_row} and \code{array_col} columns (and rounded pixel
#' coordinates) taken from the best-matching spots in \code{target_coords}.
#' 
#' @importFrom clue solve_LSAP
#' 
#' @author Nicholas J. Eagles
#' @keywords internal
.map_lsap = function(source_coords, target_coords) {
    ## For R CMD check
    array_row <- array_col <- NULL

    if (nrow(source_coords) > nrow(target_coords)) {
        stop("Internal bug: cannot fit to a smaller array!")
    }

    x = as.matrix(
        source_coords[, c("pxl_col_in_fullres", "pxl_row_in_fullres")]
    )
    y = as.matrix(
        target_coords[, c("pxl_col_in_fullres", "pxl_row_in_fullres")]
    )

    #   Compute cost matrix using squared Euclidean distance
    x2 <- rowSums(x^2)
    y2 <- rowSums(y^2)
    C  <- outer(x2, rep(1, nrow(y))) + outer(rep(1, nrow(x)), y2) - 2*(x %*% t(y))

    #   We need a square matrix: pad with zero-cost rows
    if (nrow(y) > nrow(x)) {
        C_pad <- rbind(C, matrix(0, nrow = nrow(y) - nrow(x), ncol = nrow(y)))
    } else {
        C_pad <- C
    }

    #   Solve and grab just the real rows
    perm <- clue::solve_LSAP(C_pad)[seq_len(nrow(x))]

    fit_coords = cbind(
            source_coords |> dplyr::select(-c(array_row, array_col)),
            target_coords[perm, ] |>
                dplyr::rename(
                    pxl_col_in_fullres_rounded = pxl_col_in_fullres,
                    pxl_row_in_fullres_rounded = pxl_row_in_fullres
                ) |>
                dplyr::select(
                    array_row, array_col, pxl_col_in_fullres_rounded,
                    pxl_row_in_fullres_rounded
                )
        ) |>
        dplyr::as_tibble()

    return(fit_coords)
}

#' Fit spots to a new Visium-like array: LSAP approach
#'
#' Given transformed pixel coordinates, modify the 'array_row' and
#' 'array_col' columns to represent a larger Visium capture area containing
#' all capture areas in a common coordinate system. The number of
#' array rows/cols generally changes from the Visium standards of 78 and 128
#' (and even may change in ratio between num rows and num cols).
#'
#' Mapping to the proper array coordinates is framed as the linear sum
#' assignment problem, and solved using the Hungarian algorithm. This
#' approach is far slower than \code{.fit_to_array()}, running at O(n^3) with
#' the number of spots, but guarantees a one-to-one mapping of starting to 
#' target spots, at a small cost in the Euclidean distance moved.
#'
#' @param coords A `data.frame()` containing capture areas of the
#' same group, and containing columns 'key', 'array_row', 'array_col',
#' 'pxl_row_in_fullres', and 'pxl_col_in_fullres'.
#' @param inter_spot_dist_px \code{numeric(1)} vector giving the pixel distance
#' between any 2 spots in the new coordinates.
#'
#' @return A [tibble][dplyr::reexports] with modified \code{array_row} + \code{array_col}
#' columns, as well as new \code{pxl_row_in_fullres_rounded} and
#' \code{pxl_col_in_fullres_rounded} columns representing the pixel coordinates
#' rounded to the nearest exact array coordinates.
#'
#' @author Nicholas J. Eagles
#' @keywords internal
.fit_to_array_lsap = function(coords, inter_spot_dist_px) {
    #   Build a new Visium-like array encompassing all capture areas
    target_coords = .construct_array(coords, inter_spot_dist_px)
    
    coords = coords |>
        dplyr::mutate(
            capture_area = stringr::str_split_i(key, "^[ACTG]+-1_", 2)
        )
    
    coords_list = list()
    for (this_ca in unique(coords$capture_area)) {
        this_source_coords = coords |>
            dplyr::filter(capture_area == this_ca)
        this_target_coords = target_coords |>
            dplyr::filter(
                pxl_row_in_fullres >= min(this_source_coords$pxl_row_in_fullres),
                pxl_row_in_fullres <= max(this_source_coords$pxl_row_in_fullres),
                pxl_col_in_fullres >= min(this_source_coords$pxl_col_in_fullres),
                pxl_col_in_fullres <= max(this_source_coords$pxl_col_in_fullres)
            )
        
        coords_list[[this_ca]] = .map_lsap(this_source_coords, this_target_coords)
    }

    return(do.call(rbind, coords_list) |> dplyr::select(-capture_area))
}

#' Fit spots to a new Visium-like array: fast Euclidean approach
#'
#' Given transformed pixel coordinates, modify the 'array_row' and
#' 'array_col' columns to represent a larger Visium capture area containing
#' all capture areas in a common coordinate system. The number of
#' array rows/cols generally changes from the Visium standards of 78 and 128
#' (and even may change in ratio between num rows and num cols).
#'
#' The mapping algorithm minimizes Euclidean distance of each source spot to
#' each target spot. Runtime is O(n) with the number of spots, making it
#' extremely fast. However, the Euclidean approach countintuitively may
#' result in duplicated mappings (one source to the same target) as well as
#' unexpected "holes" in the target array, which is often undesirable
#' downstream.
#'
#' @param coords A `data.frame()` whose rows represent capture areas of the
#' same group, and containing columns 'array_row', 'array_col',
#' 'pxl_row_in_fullres', and 'pxl_col_in_fullres'.
#' @param inter_spot_dist_px \code{numeric(1)} vector giving the pixel distance
#' between any 2 spots in the new coordinates.
#'
#' @return A [tibble][dplyr::reexports] with modified \code{array_row} + \code{array_col}
#' columns, as well as new \code{pxl_row_in_fullres_rounded} and
#' \code{pxl_col_in_fullres_rounded} columns representing the pixel coordinates
#' rounded to the nearest exact array coordinates.
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

    #   Verify the newly assigned array row and cols have reasonable values
    .validate_array(coords)

    return(coords)
}
