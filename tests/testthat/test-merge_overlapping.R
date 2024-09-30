test_that(
    "merge_overlapping",
    {
        if (!exists("spe")) {
            spe <- spatialLIBD::fetch_data(type = "visiumStitched_brain_spe")
        }

        #   Group colData by group and array coordinates
        grouped_coldata = colData(spe) |>
            as_tibble() |>
            group_by(group, array_row, array_col)

        #   Find the first 100 keys that overlap other spots and don't, respectively
        overlapping_keys = grouped_coldata |>
            filter(n() > 1) |>
            slice_head(n = 2) |>
            ungroup() |>
            slice_head(n = 100) |>
            pull(key)
        nonoverlapping_keys = grouped_coldata |>
            filter(n() == 1) |>
            ungroup() |>
            slice_head(n = 100) |>
            pull(key)

        #   Built a small SPE containing some overlaps and some non-overlapping spots
        small_spe = spe[, c(overlapping_keys, nonoverlapping_keys)]

        #   Merge overlapping spots
        small_spe_merged = merge_overlapping(small_spe)

        #   All array coordinates should have just one unique spot after merging
        spots_per_coord = colData(small_spe_merged) |>
            as_tibble() |>
            group_by(group, array_row, array_col) |>
            summarize(n = n()) |>
            pull(n)
        expect_equal(all(spots_per_coord == 1), TRUE)

        #   Grab a couple keys that overlap from different capture areas
        overlapping_keys = colData(small_spe) |>
            as_tibble() |>
            group_by(group, array_row, array_col) |>
            filter(n() == 2, length(unique(capture_area)) == 2) |>
            slice_head(n = 2) |>
            ungroup() |>
            slice_head(n = 2) |>
            pull(key)
        
        #   The key 
        dominant_key = overlapping_keys[
            small_spe$exclude_overlapping[
                match(overlapping_keys, small_spe$key)
            ]
        ]

        #   Expression should be added across overlapping spots for all
        #   genes
        expect_equal(
            rowSums(assays(small_spe)$counts[, overlapping_keys]),
            assays(small_spe_merged)$counts[, dominant_key]
        )
    }
)
