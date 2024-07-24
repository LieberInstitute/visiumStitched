test_that(
    "prep_fiji_coords",
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
        out_file_actual <- prep_fiji_coords(sample_info, out_dir = spe_input_dir)

        #   The expected output file should be produced
        out_file_expected <- file.path(
            spe_input_dir, "Br2719", "tissue_positions.csv"
        )
        expect_equal(out_file_actual, out_file_expected)
        expect_equal(file.exists(out_file_actual), TRUE)

        #   The tissue positions should have the expected columns
        coords <- readr::read_csv(out_file_expected, show_col_types = FALSE)
        expect_equal(
            colnames(coords),
            c(
                "key", "in_tissue", "array_row", "array_col",
                "pxl_row_in_fullres", "pxl_col_in_fullres"
            )
        )
    }
)
