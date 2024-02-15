#' Check the validity of arguments passed to \code{multi_gene.R} plotting
#' functions
#'
#' Also subset \code{spe} to the selected sample and genes, dropping genes with
#' constant expression across spots
#'
#' @author Nicholas J. Eagles
#' @inheritParams spot_plot_z_score
#' @import SpatialExperiment
#' @importFrom SummarizedExperiment colData assays
#' @return \code{SpatialExperiment} subsetted to the specified sample and to
#' each of the non-constant-expression genes
multi_gene_validity_check <- function(spe, genes, sample_id, assayname, minCount, ...) {
    #   'genes'
    if (!all(genes %in% rownames(spe))) {
        stop("The SpatialExperiment does not contain the selected genes in its rownames")
    }
    if (length(genes) <= 1) {
        stop("The behavior of this plotting function is only well-defined for a vector of at least 2 genes")
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

    #   Each multi-gene plotting function expects at least 2 genes with
    #   non-constant expression across spots. Warn if some are dropped, but halt
    #   if less than 2 remain after dropping
    spe <- spe[genes, spe$sample_id == sample_id]

    good_indices <- which(rowSds(assays(spe)[[assayname]]) != 0)
    if (length(good_indices) < 2) {
        stop("After dropping genes with no expression variation, less than 2 genes were left")
    }
    if (length(genes) - length(good_indices) > 0) {
        warning(
            sprintf(
                "Dropping gene(s) '%s' which have no expression variation",
                paste(genes[-good_indices], collapse = "', '")
            )
        )
    }

    return(spe[good_indices, ])
}
