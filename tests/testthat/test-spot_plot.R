test_that(
    "spot_plot",
    {
        spe <- fetch_data(type = "spatialDLPFC_Visium_example_subset")
        sample_id <- unique(spe$sample_id)[1]

        #   'exclude_overlapping' isn't defined yet
        expect_error(
            spot_plot(
                spe, sample_id,
                title = sample_id, var_name = "age",
                include_legend = TRUE, is_discrete = FALSE, minCount = 0,
                assayname = "logcounts"
            ),
            "exclude_overlapping"
        )

        #   All spots are excluded
        spe$exclude_overlapping <- TRUE
        expect_error(
            spot_plot(
                spe, sample_id,
                title = sample_id, var_name = "age",
                include_legend = TRUE, is_discrete = FALSE, minCount = 0,
                assayname = "logcounts"
            ),
            "No non-excluded spots"
        )

        #   Here we essentially just check that the function runs without errors
        spe$exclude_overlapping <- TRUE
        expect_equal(
            class(
                spot_plot(
                    spe, sample_id,
                    title = sample_id, var_name = "age",
                    include_legend = TRUE, is_discrete = FALSE, minCount = 0,
                    assayname = "logcounts"
                )
            ),
            c("gg", "ggplot")
        )
    }
)
