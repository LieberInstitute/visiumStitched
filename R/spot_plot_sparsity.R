#' Spatial plot of the proportion of select genes with nonzero expression by spot
#'
#' This function wraps around \code{spot_plot}, plotting a summary of expression
#' across multiple genes in a single spatial plot. In each spot, a proportion is
#' computed of how many \code{genes} have nonzero expression, and this quantity
#' is plotted across the entire capture area.
#'
#' @inheritParams spot_plot_z_score
#'
#' @return A \code{ggplot} object containing a "spot plot" of the specified
#' sample and genes
#'
#' @export
#' @author Nicholas J. Eagles
#' @import SpatialExperiment MatrixGenerics
#' @importFrom SummarizedExperiment colData assays
#' @family Spot plots summarizing expression of multiple genes simultaneously
#'
#' @examples
#'
#' #   Grab an example SpatialExperiment and suppose all of its spots should be
#' #   plotted (for spatialNAc, 'exclude_overlapping' will only have genuinely
#' #   overlapping spots be TRUE)
#' spe <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium_example_subset")
#' spe$exclude_overlapping <- FALSE
#'
#' white_matter_genes <- c(
#'     "ENSG00000197971", "ENSG00000131095", "ENSG00000123560",
#'     "ENSG00000171885"
#' )
#' spot_plot_sparsity(
#'     spe = spe,
#'     genes = white_matter_genes,
#'     sample_id = unique(spe)$sample_id[1]
#' )
spot_plot_sparsity <- function(spe, genes, sample_id, assayname = "counts", minCount = 0.1, ...) {
    #   Check validity of arguments
    spe <- multi_gene_validity_check(
        spe, genes, sample_id, assayname, minCount, ...
    )

    #   For each spot, compute proportion of marker genes with nonzero
    #   expression
    spe$prop_nonzero_marker <- colMeans(assays(spe)[[assayname]] > 0)

    #   Plot spatial distribution of this proportion
    p <- spot_plot(
        spe, sample_id,
        var_name = "prop_nonzero_marker",
        is_discrete = FALSE, minCount = minCount, assayname = assayname, ...
    )

    return(p)
}
