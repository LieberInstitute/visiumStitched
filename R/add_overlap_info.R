#' Add info about how spots overlap among capture areas
#'
#' Given a [SpatialExperiment-class][SpatialExperiment::SpatialExperiment-class]
#' and column name in its \code{colData},
#' return a modified copy of the \code{SpatialExperiment} with additional \code{colData}
#' columns: \code{spe$exclude_overlapping} and \code{spe$overlap_key}.
#'
#' \code{spe$exclude_overlapping} is `TRUE` for spots with a higher-quality
#' overlapping capture area and `FALSE` otherwise.
#' [vis_clus][spatialLIBD::vis_clus] onlydisplays `FALSE` spots to
#' prevent overplotting in regions of overlap. \code{spe$overlap_key} gives
#' comma-separated strings containing the keys of any overlapping spots, and is
#' the empty string otherwise.
#'
#' @param spe A [SpatialExperiment-class][SpatialExperiment::SpatialExperiment-class]
#' with `colData(spe)` columns
#' \code{array_row}, \code{array_col}, \code{key}, and
#' \code{capture_area}.
#' @param metric_name \code{character(1)} in \code{colnames(colData(spe))}, where
#' spots belonging to the capture area with highest average value for the metric
#' take precedence over other spots.
#'
#' @return A [SpatialExperiment][SpatialExperiment::SpatialExperiment-class]
#' object with additional \code{colData} columns \code{spe$exclude_overlapping}
#' and \code{spe$overlap_key}.
#'
#' @importFrom S4Vectors DataFrame
#' @export
#' @author Nicholas J. Eagles
#'
#' @examples
#' spe <- spatialLIBD::fetch_data(type = "visiumStitched_brain_spe")
#'
#' #    Find the mean of the 'sum_umi' metric by capture area to understand
#' #    which capture areas will be excluded in regions of overlap
#' SummarizedExperiment::colData(spe) |>
#'     dplyr::as_tibble() |>
#'     dplyr::group_by(capture_area) |>
#'     dplyr::summarize(mean_sum_umi = mean(sum_umi))
#'
#' spe <- add_overlap_info(spe, "sum_umi")
#'
#' #    See how many spots were excluded by capture area
#' table(spe$exclude_overlapping, spe$capture_area)
#'
#' #    Examine how data about overlapping spots is stored (for the first
#' #    few spots with overlap)
#' head(spe$overlap_key[spe$overlap_key != ""])
add_overlap_info <- function(spe, metric_name) {
    ## For R CMD check
    sample_id <- array_row <- array_col <- key <- capture_area <- exclude_overlapping_mean_tmp <- NULL

    #   State assumptions about columns expected to be in the colData
    expected_cols <- c(
        "array_row", "array_col", "sample_id",
        "capture_area", "key", metric_name
    )
    if (!all(expected_cols %in% colnames(colData(spe)))) {
        stop(
            sprintf(
                'Missing at least one of the following colData columns: "%s"',
                paste(expected_cols, collapse = '", "')
            )
        )
    }

    ############################################################################
    #   Compute 'overlap_key'
    ############################################################################

    col_data <- colData(spe) |>
        as_tibble() |>
        group_by(sample_id, array_row, array_col) |>
        #   Comma-separated list of spot keys at these coordinates (including
        #   each spot itself!)
        mutate(overlap_key = paste(unique(key), collapse = ",")) |>
        ungroup()

    #   Since we only want overlapping spot keys, remove the "identity spot"
    col_data$overlap_key <- sapply(
        seq_len(nrow(col_data)),
        function(i) {
            sub(col_data$key[i], "", col_data$overlap_key[i], fixed = TRUE)
        }
    )

    #   Remove double commas (only possible in the middle) or leading/trailing
    #   commas
    col_data$overlap_key <- sub(
        ",,",
        ",",
        sub("^,|,$", "", col_data$overlap_key)
    )

    ############################################################################
    #   Compute 'exclude_overlapping'
    ############################################################################

    metrics_by_id <- col_data |>
        group_by(capture_area) |>
        summarize(exclude_overlapping_mean_tmp = mean(eval(sym(metric_name)))) |>
        ungroup() |>
        arrange(desc(exclude_overlapping_mean_tmp))

    col_data$exclude_overlapping <- sapply(
        seq_len(nrow(col_data)),
        function(i) {
            #   Don't exclude spots that don't overlap anything
            if (col_data$overlap_key[i] == "") {
                return(FALSE)
            }

            #   Get the capture areas this spot overlaps
            each_array <- col_data$capture_area[
                match(
                    strsplit(col_data$overlap_key[i], ",")[[1]],
                    col_data$key
                )
            ]

            #   Exclude this spot if the highest-quality overlapping capture
            #   area is higher quality than the source capture area
            exclude <- min(match(each_array, metrics_by_id$capture_area)) <
                match(col_data$capture_area[i], metrics_by_id$capture_area)

            return(exclude)
        }
    )

    #   Add colData back to the SpatialExperiment
    temp <- colnames(spe)
    colData(spe) <- S4Vectors::DataFrame(col_data)
    colnames(spe) <- temp

    return(spe)
}
