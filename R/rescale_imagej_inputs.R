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
#' @importFrom imager load.image resize save.image
#' @importFrom rjson fromJSON
#'
#' @export
#' @author Nicholas J. Eagles
#'
#' @examples
#' #    Define sample information for the example human brain data
#' sample_info = dplyr::tibble(
#'     group = "Br2719",
#'     capture_area = c("V13B23-283_A1", "V13B23-283_C1", "V13B23-283_D1")
#' )
#' #   Add 'spaceranger_dir' column
#' sr_dir = tempdir()
#' temp = unzip(
#'     spatialLIBD::fetch_data("visiumStitched_brain_spaceranger"), exdir = sr_dir
#' )
#' sample_info$spaceranger_dir = file.path(
#'     sr_dir, sample_info$capture_area, 'outs', 'spatial'
#' )
#'
#' #   Add ImageJ-output-related columns
#' imagej_dir = tempdir()
#' temp = unzip(
#'     spatialLIBD::fetch_data("visiumStitched_brain_ImageJ_out"), exdir = imagej_dir
#' )
#' sample_info$imagej_xml_path = temp[grep('xml$', temp)]
#' sample_info$imagej_image_path = temp[grep('png$', temp)]
#'
#' out_dir = tempdir()
#' sample_info_new = rescale_imagej_inputs(sample_info, out_dir = out_dir)
#'
#' #    Scale factors are computed that are necessary downstream (i.e. with
#' #    prep_imagej_*() functions)
#' sample_info_new[, setdiff(colnames(sample_info_new), colnames(sample_info))]
#'
#' #    Image are produced that are ready for alignment in Fiji
#' list.files(out_dir)
rescale_imagej_inputs <- function(sample_info, out_dir) {
    ## For R CMD check
    group <- spot_diameter_fullres <- intra_group_scalar <- tissue_hires_scalef <- intra_group_scalar_image <- NULL

    #   State assumptions about columns expected to be in sample_info
    expected_cols <- c("capture_area", "group", "spaceranger_dir")
    if (!all(expected_cols %in% colnames(sample_info))) {
        stop(
            sprintf(
                'Missing at least one of the following columns in "sample_info": "%s"',
                paste(expected_cols, collapse = '", "')
            )
        )
    }

    dir.create(out_dir, showWarnings = FALSE)

    #   Read in high-res scalefactors and spot diameters for all samples
    sample_info$tissue_hires_scalef <- vapply(
        sample_info$spaceranger_dir,
        function(x) {
            rjson::fromJSON(
                file = file.path(x, "scalefactors_json.json")
            )[["tissue_hires_scalef"]]
        },
        numeric(1)
    ) |> unname()

    sample_info$spot_diameter_fullres <- vapply(
        sample_info$spaceranger_dir,
        function(x) {
            rjson::fromJSON(
                file = file.path(x, "scalefactors_json.json")
            )[["spot_diameter_fullres"]]
        },
        numeric(1)
    ) |> unname()

    sample_info <- sample_info |>
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
        this_image <- imager::load.image(
            file.path(sample_info$spaceranger_dir[i], "tissue_hires_image.png")
        )

        #   Rescale according to the previously calculated factor to ensure
        #   consistent distance per pixel within the group
        this_image <- imager::resize(
            this_image,
            as.integer(
                sample_info$intra_group_scalar_image[i] * dim(this_image)[1]
            ),
            as.integer(
                sample_info$intra_group_scalar_image[i] * dim(this_image)[2]
            )
        )

        imager::save.image(
            this_image,
            file.path(out_dir, sprintf("%s.png", sample_info$capture_area[i]))
        )
    }

    sample_info <- sample_info |>
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
