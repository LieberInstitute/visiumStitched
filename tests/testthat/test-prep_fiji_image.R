test_that(
    "prep_fiji_image",
    {
        ########################################################################
        #   Prepare sample_info
        ########################################################################

        sample_info <- tibble(
            group = "Br2719",
            capture_area = c("V13B23-283_A1", "V13B23-283_C1", "V13B23-283_D1")
        )
        #   Add 'spaceranger_dir' column
        sr_dir <- tempdir()
        temp <- unzip(fetch_data("visiumStitched_brain_spaceranger"), exdir = sr_dir)
        sample_info$spaceranger_dir <- file.path(
            sr_dir, sample_info$capture_area, "outs", "spatial"
        )

        #   Add Fiji-output-related columns
        fiji_dir <- tempdir()
        temp <- unzip(fetch_data("visiumStitched_brain_Fiji_out"), exdir = fiji_dir)
        sample_info$fiji_xml_path <- temp[grep("xml$", temp)]
        sample_info$fiji_image_path <- temp[grep("png$", temp)]

        sample_info <- rescale_fiji_inputs(sample_info, out_dir = tempdir())

        ########################################################################
        #   Tests
        ########################################################################

        spe_input_dir <- tempdir()
        out_paths_actual = prep_fiji_image(
            sample_info,
            out_dir = spe_input_dir, lowres_max_size = 900
        )

        #   The expected output files should be produced
        out_paths_expected = c(
            file.path(spe_input_dir, "Br2719", "tissue_lowres_image.png"),
            file.path(spe_input_dir, "Br2719", "scalefactors_json.json")
        )
        expect_equal(all(file.exists(out_paths_expected)), TRUE)
        expect_equal(setequal(out_paths_actual, out_paths_expected), TRUE)

        #   The image should have the correct maximal dimension size
        this_image <- imager::load.image(out_paths_expected[1])
        expect_equal(max(dim(this_image)), 900)
    }
)
