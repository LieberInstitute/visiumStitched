#' Write same-scale hires images for input to ImageJ
#'
#' Given a \code{tibble} of sample information (\code{sample_info}) with
#' columns \code{capture_area}, \code{group}, and \code{spaceranger_dir},
#' Write new high-resolution images for use as input to ImageJ. Particularly
#' when capture areas come from different slides, there is a risk of significant
#' scale differences among Spaceranger's \code{tissue_hires_image.png} images;
#' that is, the physical distance represented by a pixel from each capture area
#' may differ nontrivially, leading to a distance-distorted output image, and
#' inconsistent scaling when later transforming pixel coordinates. This function
#' writes approximately high-res images whose pixels are of equal physical size
#' within each \code{group}, then adds \code{intra_group_scalar} and
#' \code{group_hires_scalef} columns to \code{sample_info}. \code{intra_group_scalar}
#' gives the scalar by a which a given capture area's hires image and pixel
#' coordinates must be multiplied to match the scale of other \code{group}
#' members; \code{group_hires_scalef} gives the new \code{tissue_hires_scalef}
#' (as from Spaceranger's \code{scalefactors_json.json} file) appropriate for
#' every capture area from the group.
#'
#' @param sample_info A \code{tibble} with columns \code{capture_area},
#' \code{group}, and \code{spaceranger_dir}
#' @param out_dir A character(1) vector giving a path to a directory to place
#' the output images. Provided the parent exists, \code{out_dir}
#' will be created if necessary.
#'
#' @return A \code{tibble}: a copy of \code{sample_info} with additional columns
#' \code{intra_group_scalar} and \code{group_hires_scalef}
#'
#' @import imager
#' @importFrom dplyr mutate select group_by ungroup
#' @importFrom rjson fromJSON
#' 
#' @export
#' @author Nicholas J. Eagles
#'
#' @examples
#' \dontrun{
#' #   For internal testing
#' sample_info <- readr::read_csv("dev/test_data/sample_info.csv")
#' prep_imagej_coords(sample_info, tempdir())
#' }
#' 
#' ## TODO: add working examples
#' args(rescale_imagej_inputs)
rescale_imagej_inputs <- function(sample_info, out_dir) {
    dir.create(out_dir, showWarnings = FALSE)

    #   Read in high-res scalefactors and spot diameters for all samples
    sample_info$tissue_hires_scalef = vapply(
        sample_info$spaceranger_dir,
        function(x) {
            rjson::fromJSON(
                file = file.path(x, "scalefactors_json.json")
            )[['tissue_hires_scalef']]
        },
        numeric(1)
    )

    sample_info$spot_diameter_fullres = vapply(
        sample_info$spaceranger_dir,
        function(x) {
            rjson::fromJSON(
                file = file.path(x, "scalefactors_json.json")
            )[['spot_diameter_fullres']]
        },
        numeric(1)
    )

    sample_info = sample_info |>
        dplyr::group_by(group) |>
        dplyr::mutate(
            #   For spot coordinates: scale up so all represent the same
            #   distance per unit (pixel)
            intra_group_scalar = max(spot_diameter_fullres) /
                spot_diameter_fullres,
            #   Scale for high-res images from Spaceranger to also ensure
            #   consistent distance per pixel
            intra_group_scalar_image = intra_group_scalar *
                max(tissue_hires_scalef) / tissue_hires_scalef,
            #   Ratio of pixel/distance at high-res over at full-res
            group_hires_scalef = max(tissue_hires_scalef)
        ) |>
        dplyr::ungroup()
    
    for (i in seq(nrow(sample_info))) {
        this_image = load.image(
            file.path(sample_info$spaceranger_dir[i], 'tissue_hires_image.png')
        )

        #   Rescale according to the previously calculated factor to ensure
        #   consistent distance per pixel within the group
        this_image = resize(
            this_image,
            as.integer(
                sample_info$intra_group_scalar_image[i] * dim(this_image)[1]
            ),
            as.integer(
                sample_info$intra_group_scalar_image[i] * dim(this_image)[2]
            )
        )

        save.image(
            this_image,
            file.path(out_dir, sprintf('%s.png', sample_info$capture_area[i]))
        )
    }

    sample_info = sample_info |>
        dplyr::select(
            -c(
                spot_diameter_fullres, tissue_hires_scalef,
                intra_group_scalar_image
            )
        )
    
    #   Return 'sample_info' with additional columns 'intra_group_scalar' and
    #   'group_hires_scalef'
    return(sample_info)
}
