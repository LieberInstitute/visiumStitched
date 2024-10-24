test_that(
    "build_SpatialExperiment",
    {
        bfc <- BiocFileCache::BiocFileCache()
        gtf_cache <- BiocFileCache::bfcrpath(
            bfc,
            paste0(
                "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/",
                "release_32/gencode.v32.annotation.gtf.gz"
            )
        )

        ########################################################################
        #   Prepare sample_info
        ########################################################################

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

        spe_input_dir <- tempdir()
        prep_fiji_coords(sample_info, out_dir = spe_input_dir)
        prep_fiji_image(sample_info, out_dir = spe_input_dir)

        spe <- pkgcond::suppress_warnings(
            build_SpatialExperiment(sample_info, coords_dir = spe_input_dir, reference_gtf = gtf_cache),
            pattern = "GTF file as the one that was used by SpaceRanger"
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
