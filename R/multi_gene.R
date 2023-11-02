#' Spatial plot of multiple genes combined by averaging Z scores
#'
#' This function wraps around \code{spot_plot}, plotting a summary of expression
#' across multiple genes in a single spatial plot. A Z-score of each gene's
#' expression across all spots is taken (computed for each spot). To produce a
#' single value for each spot, Z-scores are averaged across all \code{genes}.
#' This value is plotted spatially with \code{spot_plot}.
#'
#' @param spe A \code{SpatialExperiment} with colData column \code{exclude_overlapping},
#' passed to \code{spatialLIBD::vis_gene} or \code{spatialLIBD::vis_clus}
#' @param genes character() of gene names to plot in combination, expected to be
#' in \code{rownames(spe)}
#' @param sample_id character(1) passed to \code{sampleid} in
#' \code{spatialLIBD::vis_gene} or \code{spatialLIBD::vis_clus}. Assumed to be a
#' donor, possibly consisting of several capture areas to plot at once
#' @param assayname character(1) passed to \code{spatialLIBD::vis_gene}
#' @param minCount numeric(1) passed to passed to \code{spatialLIBD::vis_gene}
#' @param ... Parameters accepted by \code{spot_plot}, excluding
#' \code{is_discrete} or \code{var_name}, which are handled internally
#' 
#' @return A \code{ggplot} object containing a "spot plot" of the specified
#' sample and genes
#'
#' @export
#' @author Nicholas J. Eagles
#' @import SpatialExperiment SummarizedExperiment rlang
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
#' white_matter_genes = c(
#'     "ENSG00000197971", "ENSG00000131095", "ENSG00000123560",
#'     "ENSG00000171885"
#' )
#' spot_plot_z_score(
#'    spe = spe,
#'    genes = white_matter_genes,
#'    sample_id = unique(spe)$sample_id[1],
#'    assayname = 'logcounts'
#' )
spot_plot_z_score = function(
        spe, genes, sample_id, assayname = "logcounts", minCount = 0, ...
    ) {
    #   Check validity of arguments
    if (!all(genes %in% rownames(spe))) {
        stop("The SpatialExperiment does not contain the selected genes in its rownames")
    }
    if (!('sample_id' %in% colnames(colData(spe))) || !(sample_id %in% spe$sample_id)) {
        stop(paste("'spe$sample_id' must exist and contain the ID", sample_id))
    }
    if (!(assayname %in% names(assays(spe)))) {
        stop(sprintf("'%s' is not an assay in 'spe'"))
    }
    if (!is.missing(var_name) || !is.missing(is_discrete)) {
        stop("The 'var_name' and 'is_discrete' parameters are internally handled and may not be specified through '...' arguments")
    }

    spe = spe[genes, spe$sample_id == sample_id]

    #   For each spot, average expression Z-scores across all selected genes
    gene_z = (assays(spe)[[assayname]] - rowMeans(assays(spe)[[assayname]])) /
        (rowSdDiffs(assays(spe)[[assayname]]))
    spe$temp_var = colMeans(gene_z, na.rm = TRUE)

    #   Plot spatial distribution of averaged expression Z-scores for this
    #   sample
    p = spot_plot(
        spe, sample_id, var_name = 'temp_var', is_discrete = FALSE,
        minCount = minCount, assayname = assayname, ...
    )

    return(p)
}
