test_that(
    ".add_error_metrics",
    {
        if (!exists("spe")) {
            spe <- spatialLIBD::fetch_data(type = "visiumStitched_brain_spe")
        }

        #   Use the existing object's colData to define coords before and after
        #   aligning to the nearest new array coordinates (as well as define the
        #   distance between spot centroids in pixels)
        coords_new <- colData(spe) |>
            cbind(spatialCoords(spe)) |>
            as_tibble() |>
            select(
                array_row, array_col, key, pxl_col_in_fullres,
                pxl_row_in_fullres, pxl_col_in_fullres_rounded,
                pxl_row_in_fullres_rounded, capture_area,
                array_row_original, array_col_original
            ) |>
            #   Fix the fact that SpatialExperiment::mirrorObject was applied,
            #   but it doesn't properly adjust the rounded pixel coordinates
            mutate(
                pxl_row_in_fullres_rounded = -1 * pxl_row_in_fullres_rounded +
                    dim(imgRaster(spe))[1] / scaleFactors(spe)
            ) |>
            group_by(capture_area) |>
            arrange(array_row) |>
            slice_head(n = 10) |>
            ungroup()
        coords <- coords_new |>
            select(-c(array_row, array_col)) |>
            rename(
                array_row = array_row_original, array_col = array_col_original
            )
        inter_spot_dist_px <- 277.1524

        coords_err <- .add_error_metrics(coords, coords_new, inter_spot_dist_px)

        #   Two error-metric columns should be added
        expect_equal(
            setdiff(colnames(coords_err), colnames(coords_new)),
            c("euclidean_error", "shared_neighbors")
        )

        #   Euclidean error is contrained to be between 0 and (less than) 1 spot
        expect_equal(
            all((coords_err$euclidean_error >= 0) & (coords_err$euclidean_error < 1)),
            TRUE
        )

        #   Shared neighbors must similarly be between 0 and 1 (it's a
        #   proportion). NAs are allowed when there are no neighbors originally
        temp <- coords_err$shared_neighbors[!is.na(coords_err$shared_neighbors)]
        expect_equal(all((temp >= 0) & (temp <= 1)), TRUE)
    }
)
