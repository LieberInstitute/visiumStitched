#' Merge overlapping spots
#'
#' Given a stitched [SpatialExperiment-class][SpatialExperiment::SpatialExperiment-class],
#' merge overlapping (same array coordinates) spots by adding
#' expression (i.e. from \code{assays(spe)$counts}), returning a
#' \code{SpatialExperiment} with at most one spot per array location.
#' 
#' \code{colData(spe)} and \code{spatialCoords(spe)} of the merged spots are
#' taken from the spots whose \code{exclude_overlapping} values are \code{TRUE}.
#'
#' @param spe A [SpatialExperiment-class][SpatialExperiment::SpatialExperiment-class]
#' with `colData(spe)` columns
#' \code{array_row}, \code{array_col}, \code{key}, \code{group}, and
#' \code{capture_area}.
#'
#' @return A [SpatialExperiment][SpatialExperiment::SpatialExperiment-class]
#' with at most one spot per array location
#'
#' @importFrom tidyr pivot_longer pivot_wider
#' @importFrom tibble rownames_to_column
#' @importFrom Matrix Matrix
#' @export
#' @author Nicholas J. Eagles
#'
#' @examples
merge_overlapping <- function(spe) {
    #   Find keys corresponding to spots that must be merged
    overlapping_keys = colData(spe) |>
        as_tibble() |>
        group_by(group, array_row, array_col) |>
        filter(n() > 1) |>
        pull(key)
    
    merged_data =
        #   For the spots to merge, form a tibble with spots as rows and genes
        #   and colData() columns as columns
        assays(spe)$counts[,match(overlapping_keys, spe$key)] |>
        as.matrix() |>
        t() |>
        cbind(
            colData(spe[,match(overlapping_keys, spe$key)]) |>
                as_tibble()
        ) |>
        as_tibble() |>
        #   Now pivot longer to have one gene per row
        tidyr::pivot_longer(
            cols = !colnames(colData(spe)), values_to = "expression",
            names_to = "gene_id"
        ) |>
        #   Average expression of a gene for a given group and array coordinates
        group_by(group, array_row, array_col, gene_id) |>
        mutate(expression = sum(expression)) |>
        ungroup() |>
        #   Now just take colData variables from the non-excluded spot
        filter(exclude_overlapping) |>
        tidyr::pivot_wider(names_from = 'gene_id', values_from = 'expression')
    
    #   Construct the SPE for just the merged spots
    spe_overlap = SpatialExperiment(
        assays = list(
            counts = merged_data |>
                select(!colnames(colData(spe))) |>
                t() |>
                Matrix::Matrix()
        ),
        colData = merged_data |>
            select(colnames(colData(spe))),
        rowData = rowData(spe),
        spatialCoords = spatialCoords(spe)[match(merged_data$key, spe$key),]
    )
    colnames(spe_overlap) = spe_overlap$key

    #   Combined the merged spots and non-overlapping spots
    spe = cbind(spe_overlap, spe[, !(spe$key %in% overlapping_keys)])
    
    return(spe)
}
