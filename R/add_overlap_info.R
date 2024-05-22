#' Add info about how spots overlap among capture areas
#'
#' Given a \code{SpatialExperiment} and column name in its \code{colData},
#' return a modified copy of the \code{SpatialExperiment} with additional \code{colData}
#' columns: \code{spe$exclude_overlapping} and \code{spe$overlap_key}
#'
#' \code{spe$exclude_overlapping} is TRUE for spots with a higher-quality overlapping
#' capture area and FALSE otherwise. \code{spot_plot} only displays FALSE spots to
#' prevent overplotting in regions of overlap. \code{spe$overlap_key} gives comma-
#' separated strings containing the keys of any overlapping spots, and is the empty
#' string otherwise.
#'
#' @param spe A \code{SpatialExperiment} with colData columns \code{exclude_overlapping},
#' \code{array_row_transformed}, \code{array_col_transformed}, \code{key}, and
#' \code{sample_id_original}
#' @param metric_name character(1) in \code{colnames(colData(spe))}, where
#' spots belonging to the capture area with highest average value for the metric
#' take precedence over other spots
#'
#' @return A \code{SpatialExperiment} object with additional \code{colData}
#' columns \code{spe$exclude_overlapping} and \code{spe$overlap_capture_area}
#'
#' @import dplyr tibble
#' @importFrom S4Vectors DataFrame
#' @export
#' @author Nicholas J. Eagles
#'
#' @examples
#' #   For internal testing
#' \dontrun{
#' library("HDF5Array")
#' spe <- loadHDF5SummarizedExperiment("dev/test_data/spe_filtered")
#' spe_new <- add_overlap_info(spe, metric_name = "sum_umi")
#' p_old <- spot_plot(
#'     spe,
#'     sample_id = "Br8492", var_name = "sample_id_original",
#'     is_discrete = TRUE
#' )
#' p_new <- spot_plot(
#'     spe_new,
#'     sample_id = "Br8492", var_name = "sample_id_original",
#'     is_discrete = TRUE
#' )
#' print(p_old)
#' print(p_new)
#' }
add_overlap_info <- function(spe, metric_name) {
    #   State assumptions about columns expected to be in the colData
    expected_cols <- c(
        "array_row_transformed", "array_col_transformed", "sample_id",
        "sample_id_original", "key", metric_name
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
        group_by(sample_id, array_row_transformed, array_col_transformed) |>
        #   Comma-separated list of spot keys at these coordinates (including
        #   each spot itself!)
        mutate(overlap_key = paste(unique(key), collapse = ",")) |>
        ungroup()

    #   Since we only want overlapping spot keys, remove the "identity spot"
    col_data$overlap_key <- sapply(
        1:nrow(col_data),
        function(i) {
            sub(col_data$key[i], "", col_data$overlap_key[i], fixed = TRUE)
        }
    )
    col_data$overlap_key <- sub("^,|,,", "", col_data$overlap_key)

    ############################################################################
    #   Compute 'exclude_overlapping'
    ############################################################################

    metrics_by_id <- col_data |>
        group_by(sample_id_original) |>
        summarize(mean_thing = mean(eval(sym(metric_name)))) |>
        ungroup() |>
        arrange(desc(mean_thing))

    col_data$exclude_overlapping <- sapply(
        1:nrow(col_data),
        function(i) {
            #   Don't exclude spots that don't overlap anything
            if (col_data$overlap_key[i] == "") {
                return(FALSE)
            }

            #   Get the capture areas this spot overlaps
            each_array <- col_data$sample_id_original[
                match(
                    strsplit(col_data$overlap_key[i], ",")[[1]],
                    col_data$key
                )
            ]

            #   Exclude this spot if the highest-quality overlapping capture
            #   area is higher quality than the source capture area
            exclude <- min(match(each_array, metrics_by_id$sample_id_original)) <
                match(col_data$sample_id_original[i], metrics_by_id$sample_id_original)

            return(exclude)
        }
    )

    #   Add colData back to the SpatialExperiment
    temp <- colnames(spe)
    colData(spe) <- S4Vectors::DataFrame(col_data)
    colnames(spe) <- temp

    return(spe)
}
