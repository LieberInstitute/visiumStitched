test_that(
    "rescale_imagej_inputs",
    {
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

        ########################################################################
        #   Tests
        ########################################################################

        sample_info_new = rescale_imagej_inputs(
            sample_info, out_dir = tempdir()
        )

        #   Exactly two columns should've been added
        expect_equal(
            setdiff(colnames(sample_info_new), colnames(sample_info)),
            c('intra_group_scalar', 'group_hires_scalef')
        )

        #   Exactly one capture area should have an intra_group_scalar of 1
        expect_equal(length(which(sample_info_new$intra_group_scalar == 1)), 1)
    }
)
