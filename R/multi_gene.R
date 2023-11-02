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
#' @import SpatialExperiment SummarizedExperiment Matrix MatrixGenerics
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
#' spot_plot_z_score(
#'     spe = spe,
#'     genes = white_matter_genes,
#'     sample_id = unique(spe)$sample_id[1],
#'     assayname = "logcounts"
#' )
spot_plot_z_score <- function(spe, genes, sample_id, assayname = "logcounts", minCount = 0, ...) {
    #   Check validity of arguments
    .multi_gene_validity_check(
        spe, genes, sample_id, assayname, minCount, ...
    )

    spe <- spe[genes, spe$sample_id == sample_id]

    #   For each spot, average expression Z-scores across all selected genes
    a <- assays(spe)[[assayname]]
    gene_z <- (a - rowMeans(a)) / rowSds(a)
    spe$Z_score <- colMeans(gene_z, na.rm = TRUE)

    #   Plot spatial distribution of averaged expression Z-scores for this
    #   sample
    p <- spot_plot(
        spe, sample_id,
        var_name = "Z_score", is_discrete = FALSE,
        minCount = minCount, assayname = assayname, ...
    )

    return(p)
}

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
#' @import SpatialExperiment SummarizedExperiment Matrix MatrixGenerics
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
    .multi_gene_validity_check(
        spe, genes, sample_id, assayname, minCount, ...
    )

    spe <- spe[genes, spe$sample_id == sample_id]

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
#' @import SpatialExperiment SummarizedExperiment Matrix
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
#' spot_plot_pca(
#'     spe = spe,
#'     genes = white_matter_genes,
#'     sample_id = unique(spe)$sample_id[1]
#' )
spot_plot_pca <- function(spe, genes, sample_id, assayname = "logcounts", minCount = 0, ...) {
    #   Check validity of arguments
    .multi_gene_validity_check(
        spe, genes, sample_id, assayname, minCount, ...
    )

    spe <- spe[genes, spe$sample_id == sample_id]

    pc_exp = prcomp(t(assays(spe)[[assayname]]), center = TRUE, scale = TRUE)
    spe$pc_select_genes <- pc_exp$x[,'PC1']

    #   Given that:
    #       - 'genes' is assumed to represent markers of the subregion (and
    #         thus their expression follows a similar pattern spatially)
    #       - the first PC captures the most variation
    #   Then each gene's coefficients to the first PC should tend to have
    #   the same sign. Next, the sign of each PC is arbitary, and we'd like
    #   plots to have positive values where expression is greater. If most
    #   genes have negative coefficients to the first PC, we reverse the
    #   sign of the coefficients to make visual intrepretation consistent
    if (mean(pc_exp$rotation[,1] > 0) < 0.5) {
        spe$pc_select_genes = -1 * spe$pc_select_genes
    }

    #   Plot spatial distribution of this proportion
    p <- spot_plot(
        spe, sample_id,
        var_name = "pc_select_genes",
        is_discrete = FALSE, minCount = minCount, assayname = assayname, ...
    )

    return(p)
}

#'   Check the validity of arguments passed to \code{multi_gene.R} plotting functions
#'
#' @author Nicholas J. Eagles
#' @inheritParams spot_plot_z_score
#' @import SpatialExperiment SummarizedExperiment
#' @return NULL
.multi_gene_validity_check <- function(spe, genes, sample_id, assayname, minCount, ...) {
    #   'genes'
    if (!all(genes %in% rownames(spe))) {
        stop("The SpatialExperiment does not contain the selected genes in its rownames")
    }

    #   'sample_id'
    if (!("sample_id" %in% colnames(colData(spe))) || !(sample_id %in% spe$sample_id)) {
        stop(paste("'spe$sample_id' must exist and contain the ID", sample_id))
    }

    #   'assayname'
    if (!(assayname %in% names(assays(spe)))) {
        stop(sprintf("'%s' is not an assay in 'spe'", assayname))
    }

    #   Not-allowed '...' parameters
    if (any(c("var_name", "is_discrete") %in% names(list(...)))) {
        stop("The 'var_name' and 'is_discrete' parameters are internally handled and may not be specified through '...' arguments")
    }
}
