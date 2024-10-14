#' Convert a \code{SpatialExperiment} object to a \code{Seurat} object
#'
#' Given a [SpatialExperiment-class][SpatialExperiment::SpatialExperiment-class]
#' object, first \code{as.Seurat()} is run, which operates on
#' [SingleCellExperiment-class][SingleCellExperiment::SingleCellExperiment-class]
#' objects. The remaining components (images, spatial coordinates) are added
#' manually. The actual appearance of images are buggy for now.
#'
#' Note that only the `lowres` images from `imgData(spe)` will be used.
#'
#' @param spe A
#' [SpatialExperiment-class][SpatialExperiment::SpatialExperiment-class] with
#' \code{colData()} or \code{spatialCoords()}
#' columns given by \code{spatial_cols}. This does not have to be a stitched
#' `spe` object as this function should work with any type of `spe` objects.
#' @param spatial_cols A `character(5)` named vector mapping which `colData(spe)`
#' or `spatialCoords(spe)` columns contain the `tissue`, `row`, `col`,
#' `imagerow`, and `imagecol` information expected by Seurat.
#' @param verbose A \code{logical(1)} vector. If `TRUE`, print status update
#' about the conversion process. This information can be useful for debugging.
#'
#' @return A \code{Seurat} object.
#'
#' @export
#' @author Nicholas J. Eagles
#' @import SpatialExperiment spatialLIBD Seurat
#' @importFrom SummarizedExperiment colData
#' @importFrom grDevices col2rgb
#' @importFrom methods new
#'
#' @examples
#' ## Download some example data
#' spe_unstitched <- spatialLIBD::fetch_data(
#'     type = "spatialDLPFC_Visium_example_subset"
#' )[seq(100), seq(100)]
#'
#' ## Make the column names unique
#' colnames(spe_unstitched) <- spatialLIBD::add_key(spe_unstitched)$key
#'
#' ## Convert from a SpatialExperiment to a Seurat object
#' seur <- as.Seurat(spe_unstitched)
#' seur
#'
#' ## Example with an stitched SPE object
#' if (!exists("spe")) {
#'     spe <- spatialLIBD::fetch_data(type = "visiumStitched_brain_spe")
#' }
#' seur_stitched <- as.Seurat(spe[seq(100), seq(100)])
#'
#' ## Let's look at our resulting Seurat object
#' seur_stitched
as.Seurat <- function(
        spe,
        spatial_cols = c(
            "tissue" = "in_tissue",
            "row" = "array_row",
            "col" = "array_col",
            "imagerow" = "pxl_row_in_fullres",
            "imagecol" = "pxl_col_in_fullres"
        ),
        verbose = TRUE) {
    SPOT_DIAMETER <- 55e-6

    #   Ensure all necessary columns are present in colData
    required_col_names <- c("tissue", "row", "col", "imagerow", "imagecol")
    if (!all(required_col_names %in% names(spatial_cols))) {
        missing_col_names <- required_col_names[!(required_col_names %in% names(spatial_cols))]
        stop(
            sprintf(
                "Expected the following named elements spatial_cols: '%s'",
                paste(missing_col_names, collapse = "', '")
            )
        )
    }

    required_cols <- spatial_cols[required_col_names]
    col_info <- cbind(colData(spe), SpatialExperiment::spatialCoords(spe))
    if (!all(required_cols %in% colnames(col_info))) {
        missing_cols <- required_cols[!(required_cols %in% colnames(col_info))]
        stop(
            sprintf(
                "Expected the following columns in colData(spe) or spatialCoords(spe): '%s'",
                paste(missing_cols, collapse = "', '")
            )
        )
    }

    #   Uniqueness of spot names
    if (any(duplicated(colnames(spe)))) {
        stop("Seurat requires colnames(spe) to be unique")
    }

    #   Low-res images must exist for each sample ID
    if (sum(imgData(spe)$image_id == "lowres") < length(unique(spe$sample_id))) {
        stop("Each sample ID must have a low-resolution image for conversion")
    }

    if (verbose) message("Running 'as.Seurat(spe)'...")
    seur <- Seurat::as.Seurat(spe)

    for (sample_id in unique(spe$sample_id)) {
        if (verbose) {
            message(
                sprintf(
                    "Adding spot coordinates and images for sample %s...",
                    sample_id
                )
            )
        }
        spe_small <- spe[, spe$sample_id == sample_id]

        coords <- col_info[, required_cols, drop = FALSE]
        colnames(coords) <- required_col_names
        coords$tissue <- as.integer(coords$tissue)
        coords <- as.data.frame(coords)

        this_img <- array(
            t(col2rgb(imgRaster(spe_small))),
            dim = c(dim(imgRaster(spe_small)), 3)
        ) / 256

        seur@images[[sample_id]] <- new(
            Class = "VisiumV1",
            image = this_img,
            scale.factors = scalefactors(
                spot = NA, fiducial = NA, hires = NA,
                lowres = imgData(spe)[
                    imgData(spe_small)$image_id == "lowres", "scaleFactor"
                ]
            ),
            coordinates = coords,
            spot.radius = SPOT_DIAMETER / scaleFactors(spe_small),
            assay = "originalexp",
            key = paste0(sample_id, "_")
        )
    }

    if (verbose) message("Returning converted object...")
    return(seur)
}
