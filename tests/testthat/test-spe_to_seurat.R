test_that(
    "spe_to_seurat",
    {
        spe <- fetch_data(type = "spatialDLPFC_Visium_example_subset")

        #   Should be missing several "transformed" colData columns
        expect_error(
            spe_to_seurat(spe, verbose = FALSE),
            "Expected the following columns"
        )

        #   Add the required transformed columns
        spe$array_row_transformed <- spe$array_row
        spe$array_col_transformed <- spe$array_col
        spe$pxl_row_in_fullres_transformed <- spatialCoords(spe)[, "pxl_row_in_fullres"]
        spe$pxl_col_in_fullres_transformed <- spatialCoords(spe)[, "pxl_col_in_fullres"]

        #   Should be missing several "transformed" colData columns
        expect_error(
            spe_to_seurat(spe, verbose = FALSE),
            "Seurat requires colnames\\(spe\\) to be unique"
        )

        #   Now (apparently) remove low-res images
        colnames(spe) <- spe$key
        temp <- imgData(spe)
        imgData(spe)$image_id <- "hires"
        expect_error(
            spe_to_seurat(spe, verbose = FALSE),
            "Each sample ID must have a low-resolution image for conversion"
        )
        imgData(spe) <- temp

        #   Remove most reducedDims, since too many names fail Seurat
        #   conventions and produce warnings that don't signify problems with
        #   'spe_to_seurat'
        colnames(SingleCellExperiment::reducedDims(spe)[[4]]) <- paste0(
            colnames(SingleCellExperiment::reducedDims(spe)[[4]]), "_"
        )
        SingleCellExperiment::reducedDims(spe) <- list(
            "first_rd" = SingleCellExperiment::reducedDims(spe)[[4]]
        )

        #   Now mostly just check that nothing fails during conversion
        seur <- spe_to_seurat(spe, verbose = FALSE)
        expect_equal(as.character(class(seur)), "Seurat")
        expect_equal(imgData(spe)$sample_id, names(seur@images))
    }
)
