#' Spatial plot of multiple genes combined by averaging Z scores
#'
#' This function wraps around \code{spot_plot}, plotting a summary of expression
#' across multiple genes in a single spatial plot. The specified expression
#' assay is subsetted to the select genes and transposed before performing PCA.
#' Each PC thus represents a "spot profile", where larger elements represent
#' spots of greater expression variation within the gene set. The first such
#' PC is plotted spatially, negated if necessary so that the majority of
#' coefficients across genes are positive (this encourages higher values to
#' represent areas of higher expression of the genes).
#'
#' @inheritParams spot_plot_z_score
#'
#' @return A \code{ggplot} object containing a "spot plot" of the specified
#' sample and genes
#'
#' @export
#' @author Nicholas J. Eagles
#' @import SpatialExperiment
#' @importFrom SummarizedExperiment colData assays
#' @importFrom stats prcomp
#' @importFrom Matrix t
#' @family Spot plots summarizing expression of multiple genes simultaneously
#'
#' @examples
#'
#' #   Grab an example SpatialExperiment and suppose all of its spots should be
#' #   plotted (for spatialNAc, 'exclude_overlapping' will only have genuinely
#' #   overlapping spots be TRUE)
#' spe <- if (!exists("spe")) {
#'     spatialLIBD::fetch_data(type = "spatialDLPFC_Visium_example_subset")
#' }
#' spe$exclude_overlapping <- FALSE
#'
#' white_matter_genes <- c(
#'     "ENSG00000197971", "ENSG00000131095", "ENSG00000123560",
#'     "ENSG00000171885"
#' )
#' spot_plot_pca(
#'     spe = spe,
#'     genes = white_matter_genes,
#'     sample_id = unique(spe)$sample_id[1]
#' )
spot_plot_pca <- function(spe, genes, sample_id, assayname = "logcounts", minCount = 0, ...) {
    #   Check validity of arguments
    spe <- multi_gene_validity_check(
        spe, genes, sample_id, assayname, minCount, ...
    )

    pc_exp <- stats::prcomp(
        Matrix::t(assays(spe)[[assayname]]), center = TRUE, scale = TRUE
    )
    spe$pc_select_genes <- pc_exp$x[, "PC1"]

    #   Given that:
    #       - 'genes' is assumed to represent markers of the subregion (and
    #         thus their expression follows a similar pattern spatially)
    #       - the first PC captures the most variation
    #   Then each gene's coefficients to the first PC should tend to have
    #   the same sign. Next, the sign of each PC is arbitary, and we'd like
    #   plots to have positive values where expression is greater. If most
    #   genes have negative coefficients to the first PC, we reverse the
    #   sign of the coefficients to make visual intrepretation consistent
    if (mean(pc_exp$rotation[, 1] > 0) < 0.5) {
        spe$pc_select_genes <- -1 * spe$pc_select_genes
    }

    #   Plot spatial distribution of this proportion
    p <- spot_plot(
        spe, sample_id,
        var_name = "pc_select_genes",
        is_discrete = FALSE, minCount = minCount, assayname = assayname, ...
    )

    return(p)
}
