#' Convert a \code{SpatialExperiment} object to a \code{Seurat} object
#'
#' Given a \code{SpatialExperiment} object, first \code{as.Seurat()} is run,
#' which operates on \code{SingleCellExperiment} objects. The remaining
#' components (images, spatial coordinates) are added manually. The actual
#' appearance of images are buggy for now.
#'
#' @param spe A \code{SpatialExperiment} with colData columns \code{in_tissue},
#' \code{array_row_transformed}, \code{array_col_transformed},
#' \code{pxl_row_in_fullres_transformed}, and \code{pxl_col_in_fullres_transformed}
#' @param verbose A logical(1) vector. If true, print status updates about the
#' conversion process
#'
#' @return A \code{Seurat} object
#'
#' @export
#' @author Nicholas J. Eagles
#' @import SpatialExperiment spatialLIBD Seurat
#' @importFrom SummarizedExperiment colData
#' @importFrom grDevices col2rgb
#'
#' @examples
#'
#' spe <- spatialLIBD::fetch_data(type = "spatialDLPFC_Visium_example_subset")
#'
#' ## Convert from a SpatialExperiment to a Seurat object
#' seur <- spe_to_seurat(spe)
#' seur
spe_to_seurat <- function(spe, verbose = TRUE) {
    SPOT_DIAMETER <- 55e-6

    #   Ensure all necessary columns are present in colData
    required_cols <- c(
        "sample_id", "in_tissue", "array_row_transformed", "array_col_transformed",
        "pxl_row_in_fullres_transformed", "pxl_col_in_fullres_transformed"
    )
    if (!all(required_cols %in% colnames(colData(spe)))) {
        missing_cols <- required_cols[!(required_cols %in% colnames(colData(spe)))]
        stop(
            sprintf(
                "Expected the following columns in colData(spe): '%s'",
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

        coords <- colData(spe_small)[
            ,
            c(
                "in_tissue", "array_row_transformed", "array_col_transformed",
                "pxl_row_in_fullres_transformed", "pxl_col_in_fullres_transformed"
            )
        ]
        colnames(coords) <- c("tissue", "row", "col", "imagerow", "imagecol")
        coords$tissue <- as.integer(coords$tissue)
        coords <- as.data.frame(coords)

        this_img <- array(
            t(col2rgb(imgRaster(spe_small))),
            dim = c(dim(imgRaster(spe_small)), 3)
        ) / 256

        seur@images[[sample_id]] <- Seurat:::VisiumV1(
            image = this_img,
            scale.factors = scalefactors(
                spot = NA, fiducial = NA, hires = NA,
                lowres = imgData(spe)[imgData(spe_small)$image_id == "lowres", "scaleFactor"]
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
