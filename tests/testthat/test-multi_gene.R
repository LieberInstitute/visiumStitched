library(SpatialExperiment)

test_that(
    "spot_plot_z_score",
    {
        spe <- fetch_data(type = "spatialDLPFC_Visium_example_subset")

        #   Bad gene names (gene symbols are not in rownames)
        sample_id <- unique(spe$sample_id)[1]
        genes = 'MOBP'
        expect_error(
            spot_plot_z_score(spe, genes, sample_id),
            "does not contain the selected genes"
        )
        
        #   Bad sample ID
        sample_id = 'asdasdasdasd'
        genes = rownames(spe)[1:10]
        expect_error(
            spot_plot_z_score(spe, genes, sample_id),
            "must exist and contain the ID"
        )

        #   Bad assay name
        sample_id <- unique(spe$sample_id)[1]
        genes = rownames(spe)[1:10]
        expect_error(
            spot_plot_z_score(spe, genes, sample_id, assayname = "asfasdad"),
            "is not an assay"
        )

        #   Specifying arguments technically accepted by "..." but internally
        #   handled and nonsensical
        expect_error(
            spot_plot_z_score(spe, genes, sample_id, is_discrete = TRUE),
            "internally handled and may not be specified"
        )
        expect_error(
            spot_plot_z_score(spe, genes, sample_id, var_name = TRUE),
            "internally handled and may not be specified"
        )

        #   Here we essentially just check that the function runs without errors
        expect_equal(
            class(spot_plot_z_score(spe, genes, sample_id)),
            c("gg", "ggplot")
        )
    }
)
