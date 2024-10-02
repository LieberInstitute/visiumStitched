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
#' @importFrom BiocGenerics cbind
#' @export
#' @author Nicholas J. Eagles
#'
#' @examples
#' if (!exists("spe")) {
#'     spe <- spatialLIBD::fetch_data(type = "visiumStitched_brain_spe")
#' }
#' 
#' #   Group colData by group and array coordinates
#' grouped_coldata = colData(spe) |>
#'     dplyr::as_tibble() |>
#'     dplyr::group_by(group, array_row, array_col)
#' 
#' #   Find the first 100 keys that overlap other spots and don't, respectively
#' overlapping_keys = grouped_coldata |>
#'     dplyr::filter(n() > 1) |>
#'     dplyr::slice_head(n = 2) |>
#'     dplyr::ungroup() |>
#'     dplyr::slice_head(n = 100) |>
#'     dplyr::pull(key)
#' nonoverlapping_keys = grouped_coldata |>
#'     dplyr::filter(n() == 1) |>
#'     dplyr::ungroup() |>
#'     dplyr::slice_head(n = 100) |>
#'     dplyr::pull(key)
#' 
#' #   Built a small SPE containing some overlaps and some non-overlapping spots
#' small_spe = spe[, c(overlapping_keys, nonoverlapping_keys)]
#' 
#' #   Merge overlapping spots
#' small_spe_merged = merge_overlapping(small_spe)
#' 
#' #   All array coordinates have just one unique spot after merging
#' colData(small_spe_merged) |>
#'     dplyr::as_tibble() |>
#'     dplyr::group_by(group, array_row, array_col) |>
#'     dplyr::summarize(n = n()) |>
#'     dplyr::pull(n) |>
#'     table()
#' 
merge_overlapping <- function(spe) {
    ## For R CMD CHECK
    array_row = array_col = exclude_overlapping = gene_id = group = key = NULL

    #   Find keys corresponding to spots that must be merged
    overlapping_keys = colData(spe) |>
        as_tibble() |>
        group_by(group, array_row, array_col) |>
        filter(n() > 1) |>
        pull(key)
    
    #   Nothing to merge
    if(length(overlapping_keys) == 0) {
        return(spe)
    }

    #   Check assays and halt or warn
    if (!('counts' %in% names(assays(spe)))) {
        stop("'counts' assay missing; unable to merge overlapping spots.")
    }
    if (length(assays(spe)) > 1) {
        warning("Dropping assays other than 'counts' for merging.")
    }
    assays(spe) = list(counts = assays(spe)$counts)

    if(length(reducedDims(spe)) > 0) {
        warning("Dropped reducedDims(spe) for merging")
    }
    reducedDims(spe) = list()
    
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
    spe = BiocGenerics::cbind(spe_overlap, spe[, !(spe$key %in% overlapping_keys)])
    
    return(spe)
}
