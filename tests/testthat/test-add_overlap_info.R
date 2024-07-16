test_that(
    "add_overlap_info",
    {
        spe <- fetch_data(type = "Visium_LS_spe")

        spe$my_metric = 0
        spe$my_metric[spe$capture_area == "V13B23-283_A1"] = 1
        spe$my_metric[spe$capture_area == "V13B23-283_C1"] = 2
        spe$my_metric[spe$capture_area == "V13B23-283_D1"] = 3

        spe_new = add_overlap_info(spe, "my_metric")

        #   V13B23-283_A1 should have excluded spots (lowest in my_metric) while
        #   V13B23-283_D1 should have none
        expect_equal(
            any(spe$exclude_overlapping[spe$capture_area == "V13B23-283_A1"]),
            TRUE
        )
        expect_equal(
            any(spe$exclude_overlapping[spe$capture_area == "V13B23-283_D1"]),
            FALSE
        )

        overlap_keys = colData(spe) |>
            filter(capture_area == "V13B23-283_A1", overlap_key != "") |>
            pull(overlap_key) |>
            paste(collapse = ",") |>
            str_split(',')
        
        #   All keys overlapping V13B23-283_A1 should exist in the object and
        #   belong to capture area V13B23-283_D1
        expect_equal(all(overlap_keys[[1]] %in% spe$key), TRUE)
        expect_equal(all(str_split_i(overlap_keys[[1]], '_', 2) == "D1"), TRUE)
    }
)
