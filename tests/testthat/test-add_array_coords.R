test_that(
    "add_array_coords",
    {
        spe <- fetch_data(type = "Visium_LS_spe")

        ########################################################################
        #   Prepare sample_info
        ########################################################################

        sample_info = tibble(
            group = "Br2719",
            capture_area = c("V13B23-283_A1", "V13B23-283_C1", "V13B23-283_D1")
        )
        #   Add 'spaceranger_dir' column
        sr_dir = tempdir()
        temp = unzip(fetch_data("Visium_LS_spaceranger"), exdir = sr_dir)
        sample_info$spaceranger_dir = file.path(
            sr_dir, sample_info$capture_area, 'outs', 'spatial'
        )

        #   Add ImageJ-output-related columns
        imagej_dir = tempdir()
        temp = unzip(fetch_data("Visium_LS_ImageJ_out"), exdir = imagej_dir)
        sample_info$imagej_xml_path = temp[grep('xml$', temp)]
        sample_info$imagej_image_path = temp[grep('png$', temp)]

        sample_info = rescale_imagej_inputs(sample_info, out_dir = tempdir())

        spe_input_dir = tempdir()
        prep_imagej_coords(sample_info, out_dir = spe_input_dir)
        prep_imagej_image(sample_info, out_dir = spe_input_dir)

        ########################################################################
        #   Tests
        ########################################################################

        #   Remove any colData columns that should be added by add_array_coords()
        added_cols_regex = "^(array|pxl)_(row|col)_(in_fullres)?_(transformed|original|rounded)$"
        temp = colnames(spe)
        colData(spe) = colData(spe) |>
            as_tibble() |>
            mutate(across(matches(added_cols_regex), ~ NULL)) |>
            DataFrame()
        colnames(spe) = temp

        spe_new = add_array_coords(
            spe, sample_info, spe_input_dir, overwrite = FALSE
        )

        #   12 columns should've been added, matching the specific naming
        #   pattern
        expect_equal(
            length(grep(added_cols_regex, colnames(colData(spe_new)))), 12
        )

        #   Array and spatial coords shouldn't be overwritten
        expect_identical(spe$array_row, spe_new$array_row)
        expect_identical(spe$array_col, spe_new$array_col)
        expect_identical(spatialCoords(spe), spatialCoords(spe_new))

        spe_new = add_array_coords(
            spe, sample_info, spe_input_dir, overwrite = TRUE
        )

        #   12 columns should've been added, matching the specific naming
        #   pattern
        expect_equal(
            length(grep(added_cols_regex, colnames(colData(spe_new)))), 12
        )

        #   spatialCoords should be updated with transformed coordinates
        expect_equal(
            spe_new$pxl_row_in_fullres_transformed,
            spatialCoords(spe)[['pxl_row_in_fullres']]
        )
        expect_equal(
            spe_new$pxl_col_in_fullres_transformed,
            spatialCoords(spe)[['pxl_col_in_fullres']]
        )
    }
)
