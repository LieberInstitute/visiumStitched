#' Spatial plots of discrete or continuous features for stitched-together
#' capture areas.
#'
#' This function is essentially a wrapper around spatialLIBD::vis_clus and
#' spatialLIBD::vis_gene, suitable for merged samples (each sample in the
#' SpatialExperiment 'spe' is a donor consisting of multiple capture areas, with
#' colData column 'exclude_overlapping' indicating overlapping spots to drop (to
#' prevent overplotting).
#'
#' Spot sizes are *almost* consistent among donors, regardless of full-
#' resolution image dimensions, when title is NULL, include_legend is FALSE,
#' and the plot is saved to a square output (e.g. PDF with 7in width and
#' height). However, ggplot does not seem to scale plots of different aspect
#' ratios exactly consistently when writing to PDF (untested for other formats)
#'
#' @param spe A \code{SpatialExperiment} with colData column \code{exclude_overlapping},
#' passed to \code{spatialLIBD::vis_gene} or \code{spatialLIBD::vis_clus}
#' @param sample_id character(1) passed to \code{sampleid} in
#' \code{spatialLIBD::vis_gene} or \code{spatialLIBD::vis_clus}. Assumed to be a
#' donor, possibly consisting of several capture areas to plot at once
#' @param image_id character(1) giving the name of the image (e.g. "lowres") to
#' plot, used both to determine an appropriate spot size and passed to
#' \code{spatialLIBD::vis_gene} or \code{spatialLIBD::vis_clus}
#' @param title character(1) giving the title of the plot
#' @param var_name character(1) passed to \code{geneid} for \code{spatialLIBD::vis_gene}
#' or \code{clustervar} for \code{spatialLIBD::vis_clus}
#' @param include_legend logical(1): if FALSE, remove the plot legend
#' @param is_discrete logical(1): if TRUE, use \code{spatialLIBD::vis_clus};
#' otherwise, use \code{spatialLIBD::vis_gene}
#' @param colors character() of colors passed to \code{cont_colors} in
#' \code{spatialLIBD::vis_gene} if not \code{is_discrete}
#' @param assayname character(1) passed to \code{spatialLIBD::vis_gene} if
#' not \code{is_discrete}
#' @param minCount numeric(1) passed to passed to \code{spatialLIBD::vis_gene} if
#' not \code{is_discrete}
#' @param spatial logical(1) passed to \code{sampleid} in
#' \code{spatialLIBD::vis_gene} or \code{spatialLIBD::vis_clus}
#'
#' @return A \code{ggplot} object containing a "spot plot" of the specified sample
#'
#' @export
#' @author Nicholas J. Eagles
#' @import viridisLite spatialLIBD ggplot2 SpatialExperiment SummarizedExperiment
#'
#' @examples
#'
#' #   Grab an example SpatialExperiment and suppose all of its spots should be
#' #   plotted (for spatialNAc, 'exclude_overlapping' will only have genuinely
#' #   overlapping spots be TRUE)
#' spe <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium_example_subset")
#' spe$exclude_overlapping <- FALSE
#'
#' #   Plot age spatially for the first sample
#' sample_id <- unique(spe$sample_id)[1]
#' p <- spot_plot(
#'     spe,
#'     sample_id = sample_id,
#'     title = sample_id, var_name = "age",
#'     include_legend = TRUE, is_discrete = FALSE, minCount = 0,
#'     assayname = "logcounts"
#' )
#' print(p)
spot_plot <- function(spe, sample_id, image_id = "lowres",
    title = sprintf("%s_%s", sample_id, var_name), var_name,
    include_legend = TRUE, is_discrete, colors = NULL,
    assayname = "logcounts", minCount = 0.5, spatial = FALSE) {
    #   This value was determined empirically, and results in good spot sizes.
    #   Note that it's sample-independent, and the final spot size to pass to
    #   'vis_gene' or 'vis_clus' uses this value along with the image
    #   dimensions, scale factors, and spatial coordinates for this particular
    #   sample
    IDEAL_POINT_SIZE <- 200

    ############################################################################
    #   Check validity of arguments
    ############################################################################

    # (Note that 'sample_id', 'var_name', 'assayname', 'minCount', and 'colors'
    # are not checked for validity here, since spatialLIBD functions handle
    # their validity)

    #   Check validity of spatial coordinates
    if (!all(c("pxl_col_in_fullres", "pxl_row_in_fullres") == sort(colnames(spatialCoords(spe))))) {
        stop("Abnormal spatial coordinates: should have 'pxl_row_in_fullres' and 'pxl_col_in_fullres' columns.")
    }

    #   State assumptions about columns expected to be in the colData
    expected_cols <- c("array_row", "array_col", "sample_id", "exclude_overlapping")
    if (!all(expected_cols %in% colnames(colData(spe)))) {
        stop(
            sprintf(
                'Missing at least one of the following colData columns: "%s"',
                paste(expected_cols, collapse = '", "')
            )
        )
    }

    #   Subset to specific sample ID, and ensure overlapping spots are dropped
    subset_cols <- (spe$sample_id == sample_id) &
        (is.na(spe$exclude_overlapping) | !spe$exclude_overlapping)
    if (length(which(subset_cols)) == 0) {
        stop("No non-excluded spots belong to this sample. Perhaps check spe$exclude_overlapping for issues.")
    }
    spe_small <- spe[, subset_cols]

    ############################################################################
    #   Compute an appropriate spot size for this sample
    ############################################################################

    #   Determine some pixel values for the horizontal bounds of the spots
    MIN_COL <- min(spatialCoords(spe_small)[, "pxl_row_in_fullres"])
    MAX_COL <- max(spatialCoords(spe_small)[, "pxl_row_in_fullres"])

    #   The distance between spots (in pixels) is double the average distance
    #   between array columns
    INTER_SPOT_DIST_PX <- 2 * (MAX_COL - MIN_COL) /
        (max(spe_small$array_col) - min(spe_small$array_col))

    #   Find the appropriate spot size for this donor. This can vary because
    #   ggplot downscales a plot the fit desired output dimensions (in this
    #   case presumably a square region on a PDF), and stitched images can vary
    #   in aspect ratio. Also, lowres images always have a larger image
    #   dimension of 1200, no matter how many spots fit in either dimension.
    small_image_data <- imgData(spe_small)[
        imgData(spe_small)$image_id == image_id,
    ]
    spot_size <- IDEAL_POINT_SIZE * INTER_SPOT_DIST_PX *
        small_image_data$scaleFactor / max(dim(small_image_data$data[[1]]))


    ############################################################################
    #   Produce the plot
    ############################################################################

    #   If the quantity to plot is discrete, use 'vis_clus'. Otherwise use
    #   'vis_gene'.
    if (is_discrete) {
        #   For 'vis_clus' only, supply a color scale if 'color' is not NULL
        if (is.null(colors)) {
            p <- vis_clus(
                spe_small,
                sampleid = sample_id, image_id = image_id, clustervar = var_name, auto_crop = FALSE,
                return_plots = TRUE, spatial = spatial, point_size = spot_size
            )
        } else {
            p <- vis_clus(
                spe_small,
                sampleid = sample_id, image_id = image_id, clustervar = var_name, auto_crop = FALSE,
                return_plots = TRUE, spatial = spatial, colors = colors,
                point_size = spot_size
            )
        }
    } else {
        p <- vis_gene(
            spe_small,
            sampleid = sample_id, image_id = image_id, geneid = var_name, return_plots = TRUE,
            spatial = spatial, point_size = spot_size, assayname = assayname,
            cont_colors = viridisLite::plasma(21), alpha = 1, auto_crop = FALSE,
            minCount = minCount
        )
    }

    #   Remove the legend if requested
    if (!include_legend) {
        p <- p + theme(legend.position = "none")
    }

    #   Overwrite the title
    p <- p + labs(title = title)

    return(p)
}
