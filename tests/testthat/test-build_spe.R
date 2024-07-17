test_that(
    "build_spe",
    {
        bfc <- BiocFileCache()
        gtf_cache <- bfcrpath(
            bfc,
            paste0(
                "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
                "release_32/gencode.v32.annotation.gtf.gz"
            )
        )

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

        spe = build_spe(
            sample_info, coords_dir = spe_input_dir, reference_gtf = gtf_cache
        )

        ########################################################################
        #   Tests
        ########################################################################

        #   A SpatialExperiment with groups as samples is returned. All capture
        #   areas are present in the 'capture_area' colData() column
        expect_equal(is(spe, "SpatialExperiment"), TRUE)
        expect_equal(setequal(spe$sample_id, sample_info$group), TRUE)
        expect_equal(setequal(spe$capture_area, sample_info$capture_area), TRUE)
    }
)
