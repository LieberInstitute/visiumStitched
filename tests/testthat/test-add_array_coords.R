test_that(
    "add_array_coords",
    {
        if (!exists("spe")) {
            spe <- spatialLIBD::fetch_data(type = "visiumStitched_brain_spe")
        }

        ########################################################################
        #   Prepare sample_info
        ########################################################################

        if (file.exists("sample_info.rds")) {
            sample_info <- readRDS('sample_info.rds')
        } else {
            sample_info <- dplyr::tibble(
                group = "Br2719",
                capture_area = c("V13B23-283_A1", "V13B23-283_C1", "V13B23-283_D1")
            )
            #   Add 'spaceranger_dir' column
            sr_dir <- tempdir()
            temp <- unzip(
                spatialLIBD::fetch_data("visiumStitched_brain_spaceranger"),
                exdir = sr_dir
            )
            sample_info$spaceranger_dir <- file.path(
                sr_dir, sample_info$capture_area, "outs", "spatial"
            )

            #   Add Fiji-output-related columns
            fiji_dir <- tempdir()
            temp <- unzip(
                spatialLIBD::fetch_data("visiumStitched_brain_Fiji_out"),
                exdir = fiji_dir
            )
            sample_info$fiji_xml_path <- temp[grep("xml$", temp)]
            sample_info$fiji_image_path <- temp[grep("png$", temp)]

            ## Re-size images and add more information to the sample_info
            sample_info <- rescale_fiji_inputs(sample_info, out_dir = tempdir())

            saveRDS(sample_info, "sample_info.rds")
        }

        spe_input_dir <- tempdir()
        prep_fiji_coords(sample_info, out_dir = spe_input_dir)
        prep_fiji_image(sample_info, out_dir = spe_input_dir)

        ########################################################################
        #   Tests
        ########################################################################

        #   Remove any colData columns that should be added by add_array_coords()
        added_cols_regex <- "^(array|pxl)_(row|col)(_in_fullres)?_(original|rounded)$"
        temp <- colnames(spe)
        colData(spe) <- colData(spe) |>
            as_tibble() |>
            mutate(across(matches(added_cols_regex), ~NULL)) |>
            DataFrame()
        colnames(spe) <- temp

        spe_new <- add_array_coords(spe, sample_info, spe_input_dir)

        #   6 columns should've been added, matching the specific naming
        #   pattern
        expect_equal(
            length(grep(added_cols_regex, colnames(colData(spe_new)))), 6
        )

        #   "Original" columns should actually have their original values
        expect_identical(spe$array_row, spe_new$array_row_original)
        expect_identical(spe$array_col, spe_new$array_col_original)
        expect_identical(
            spatialCoords(spe)[, "pxl_row_in_fullres"],
            spe_new$pxl_row_in_fullres_original
        )
        expect_identical(
            spatialCoords(spe)[, "pxl_col_in_fullres"],
            spe_new$pxl_col_in_fullres_original
        )
    }
)
