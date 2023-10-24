#   Wrapper around spatialLIBD::vis_clus and spatialLIBD::vis_gene, suitable
#   for merged samples (each sample in the SpatialExperiment 'spe' is a donor
#   consisting of multiple capture areas, with colData column
#   'exclude_overlapping' indicating overlapping spots to drop (to prevent
#   overplotting).
#
#   Spot sizes are *almost* consistent among donors, regardless of full-
#   resolution image dimensions, when title is NULL, include_legend is FALSE,
#   and the plot is saved to a square output (e.g. PDF with 7in width and
#   height). However, ggplot does not seem to scale plots of different aspect
#   ratios exactly consistently when writing to PDF (untested for other formats)
#
#   Return a spot plot of sample 'sampleid', assumed to be a donor. 'coldatavar'
#   (character(1)) must be a valid colname in colData(spe).
#
#   spe:            passed to 'spe' in either 'vis_gene' or 'vis_clus'
#   sample_id:      passed to 'sampleid'
#   title:          title for the plot, expected to be one line (avoid use of
#                   "\n")
#   var_name:       passed to 'geneid' for 'vis_gene' and 'clustervar' for
#                   'vis_clus'
#   include_legend: (logical) if FALSE, remove the legend
#   is_discrete:    (logical) if TRUE, use 'vis_clus'; otherwise, use 'vis_gene'
#   colors:         passed to 'colors' for 'vis_gene' if not [is_discrete]
#   assayname:      passed to 'assayname' for 'vis_gene' if not [is_discrete]
#   minCount:       passed to 'minCount' for 'vis_gene' if not [is_discrete]
#
#   Returns a ggplot object
spot_plot <- function(
        spe, sample_id, title, var_name, include_legend, is_discrete,
        colors = NULL, assayname = "logcounts", minCount = 0.5
    ) {
    IDEAL_POINT_SIZE <- 200
    IMAGE_ID <- "lowres"

    #   Subset to specific sample ID, and ensure overlapping spots are dropped
    spe_small = spe[
        ,
        (spe$sample_id == sample_id) &
        (is.na(spe$exclude_overlapping) | !spe$exclude_overlapping)
    ]
    
    #   Determine some pixel values for the horizontal bounds of the spots
    MIN_COL = min(spatialCoords(spe_small)[, 'pxl_row_in_fullres'])
    MAX_COL = max(spatialCoords(spe_small)[, 'pxl_row_in_fullres'])

    #   The distance between spots (in pixels) is double the average distance
    #   between array columns 
    INTER_SPOT_DIST_PX = 2 * (MAX_COL - MIN_COL) /
        (max(spe_small$array_col) - min(spe_small$array_col))
    
    #   Find the appropriate spot size for this donor. This can vary because
    #   ggplot downscales a plot the fit desired output dimensions (in this
    #   case presumably a square region on a PDF), and stitched images can vary
    #   in aspect ratio. Also, lowres images always have a larger image
    #   dimension of 1200, no matter how many spots fit in either dimension.
    small_image_data = imgData(spe_small)[
        imgData(spe_small)$image_id == IMAGE_ID,
    ]
    spot_size = IDEAL_POINT_SIZE * INTER_SPOT_DIST_PX *
       small_image_data$scaleFactor / max(dim(small_image_data$data[[1]]))

    #   If the quantity to plot is discrete, use 'vis_clus'. Otherwise use
    #   'vis_gene'.
    if (is_discrete) {
        #   For 'vis_clus' only, supply a color scale if 'color' is not NULL
        if (is.null(colors)) {
            p <- vis_clus(
                spe_small,
                sampleid = sample_id, clustervar = var_name, auto_crop = FALSE,
                return_plots = TRUE, spatial = FALSE, point_size = spot_size
            )
        } else {
            p <- vis_clus(
                spe_small,
                sampleid = sample_id, clustervar = var_name, auto_crop = FALSE,
                return_plots = TRUE, spatial = FALSE, colors = colors,
                point_size = spot_size
            )
        }
    } else {
        p <- vis_gene(
            spe_small,
            sampleid = sample_id, geneid = var_name, return_plots = TRUE,
            spatial = FALSE, point_size = spot_size, assayname = assayname,
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
