library(SpatialExperiment)

#   Given a SpatialExperiment and multi-gene-plotting function 'plot_fun'
#   from this package, run tests and return NULL
test_plot_fun <- function(spe, plot_fun) {
    #   Bad gene names (gene symbols are not in rownames)
    sample_id <- unique(spe$sample_id)[1]
    genes <- "MOBP"
    expect_error(
        plot_fun(spe, genes, sample_id),
        "does not contain the selected genes"
    )

    #   Bad sample ID
    sample_id <- "asdasdasdasd"
    genes <- rownames(spe)[1:10]
    expect_error(
        plot_fun(spe, genes, sample_id),
        "must exist and contain the ID"
    )

    #   Bad assay name
    sample_id <- unique(spe$sample_id)[1]
    genes <- rownames(spe)[1:10]
    expect_error(
        plot_fun(spe, genes, sample_id, assayname = "asfasdad"),
        "is not an assay"
    )

    #   Specifying arguments technically accepted by "..." but internally
    #   handled and nonsensical
    expect_error(
        plot_fun(spe, genes, sample_id, is_discrete = TRUE),
        "internally handled and may not be specified"
    )
    expect_error(
        plot_fun(spe, genes, sample_id, var_name = TRUE),
        "internally handled and may not be specified"
    )

    #   Here we essentially just check that the function runs without errors
    expect_equal(
        class(plot_fun(spe, genes, sample_id)),
        c("gg", "ggplot")
    )
}

spe <- fetch_data(type = "spatialDLPFC_Visium_example_subset")
spe$exclude_overlapping <- FALSE

test_that("spot_plot_z_score", {
    test_plot_fun(spe, spot_plot_z_score)
})
test_that("spot_plot_sparsity", {
    test_plot_fun(spe, spot_plot_sparsity)
})
