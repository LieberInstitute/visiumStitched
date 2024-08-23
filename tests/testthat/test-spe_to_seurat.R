test_that(
    "spe_to_seurat",
    {
        spe <- fetch_data(type = "spatialDLPFC_Visium_example_subset")[seq(100), seq(100)]

        ## Make the column names unique
        colnames(spe) <- spatialLIBD::add_key(spe)$key

        #   Should be missing several "transformed" colData columns
        expect_error(
            spe_to_seurat(
                spe,
                spatial_cols = c(
                    "tissue" = "in_tissue",
                    "row" = "array_row_transformed",
                    "col" = "array_col_transformed",
                    "imagerow" = "pxl_row_in_fullres_transformed",
                    "imagecol" = "pxl_col_in_fullres_transformed"
                ),
                verbose = FALSE
            ),
            "Expected the following columns"
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
