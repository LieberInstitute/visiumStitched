#' Create low-res images and scale factors from high-res ImageJ output images
#'
#' After stitching all groups in \code{sample_info} with ImageJ, images of
#' various resolutions (pixel dimensions) are left. This function creates copies
#' of each image whose largest dimension is \code{lowres_max_size} pixels. It
#' also creates a corresponding \code{scalefactors_json.json} file much like
#' Spaceranger's. In conjunction with \code{prep_imagej_image()}, this function
#' prepares for building the \code{SpatialExperiment} with \code{build_spe()}.
#'
#' @param sample_info A \code{tibble} with columns \code{capture_area},
#' \code{group}, \code{imagej_image_path}, \code{spaceranger_dir},
#' \code{intra_group_scalar}, and \code{group_hires_scalef}
#' @param out_dir A character(1) vector giving a path to a directory to place
#' the output image(s) and scale factors. Provided the parent exists, \code{out_dir}
#' will be created if necessary.
#' @param lowres_max_size An integer(1) vector: the resolution (number of
#' pixels) of the larger dimension of the output image(s), considered to be "low
#' resolution".
#'
#' @return NULL
#'
#' @importFrom imager load.image resize save.image
#' @importFrom rjson fromJSON toJSON
#'
#' @export
#' @author Nicholas J. Eagles
#'
#' @examples
#' #    Define sample information for the example LS data 
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
#' sample_info = rescale_imagej_inputs(sample_info, out_dir = tempdir())
#' 
#' spe_input_dir = tempdir()
#' prep_imagej_image(
#'     sample_info, out_dir = spe_input_dir, lowres_max_size = 1000
#' )
#' 
#' #    A "low resolution" stitched image was produced, which has 1000
#' #    pixels in its largest dimension
#' this_image <- imager::load.image(
#'     file.path(spe_input_dir, "Br2719", 'tissue_lowres_image.png')
#' )
#' dim(this_image)
prep_imagej_image <- function(sample_info, out_dir, lowres_max_size = 1200) {
    #   State assumptions about columns expected to be in sample_info
    expected_cols <- c(
        "capture_area", "group", "imagej_image_path", "intra_group_scalar",
        "group_hires_scalef", "spaceranger_dir"
    )
    if (!all(expected_cols %in% colnames(sample_info))) {
        stop(
            sprintf(
                'Missing at least one of the following columns in "sample_info": "%s"',
                paste(expected_cols, collapse = '", "')
            )
        )
    }

    if (!all(file.exists(sample_info$imagej_image_path))) {
        stop("All files in 'sample_info$imagej_image_path' must exist.")
    }

    dir.create(out_dir, showWarnings = FALSE)

    for (this_group in unique(sample_info$group)) {
        this_sample_info <- sample_info |>
            dplyr::filter(group == this_group)

        if (length(unique(this_sample_info$imagej_image_path)) > 1) {
            stop("Expected one unique path for 'imagej_image_path' per group in 'sample_info'.")
        }

        this_image <- imager::load.image(this_sample_info$imagej_image_path[1])

        #   Combine info about the original scalefactors of the first capture
        #   area with group-related scalars to form a new scalefactors JSON
        #   for the whole stitched group
        sr_json <- rjson::fromJSON(
            file = file.path(
                this_sample_info$spaceranger_dir[1], "scalefactors_json.json"
            )
        )

        low_over_hi <- lowres_max_size / max(dim(this_image)[seq(2)])
        sr_json <- list(
            tissue_hires_scalef = this_sample_info$group_hires_scalef[[1]],
            tissue_lowres_scalef = this_sample_info$group_hires_scalef[[1]] *
                low_over_hi,
            spot_diameter_fullres = sr_json$spot_diameter_fullres *
                this_sample_info$intra_group_scalar[[1]]
        )

        this_image <- imager::resize(
            this_image,
            as.integer(low_over_hi * dim(this_image)[1]),
            as.integer(low_over_hi * dim(this_image)[2])
        )

        #   Save the lowres image and scalefactors JSON in a subdirectory of
        #   'out_dir' named with the current group
        this_out_dir <- file.path(out_dir, this_group)
        dir.create(this_out_dir, showWarnings = FALSE)
        imager::save.image(
            this_image, file.path(this_out_dir, "tissue_lowres_image.png")
        )
        write(
            rjson::toJSON(sr_json),
            file.path(this_out_dir, "scalefactors_json.json")
        )
    }

    return(invisible(NULL))
}
