test_that(
    "add_overlap_info",
    {
        spe <- fetch_data(type = "visiumStitched_brain_spe")

        spe$my_metric = 0
        spe$my_metric[spe$capture_area == "V13B23-283_A1"] = 1
        spe$my_metric[spe$capture_area == "V13B23-283_C1"] = 2
        spe$my_metric[spe$capture_area == "V13B23-283_D1"] = 3

        spe_new = add_overlap_info(spe, "my_metric")

        #   V13B23-283_A1 should have excluded spots (lowest in my_metric) while
        #   V13B23-283_D1 should have none
        expect_equal(
            any(spe_new$exclude_overlapping[spe_new$capture_area == "V13B23-283_A1"]),
            TRUE
        )
        expect_equal(
            any(spe_new$exclude_overlapping[spe_new$capture_area == "V13B23-283_D1"]),
            FALSE
        )

        #   'overlap_key' should consist of a comma-separated list of valid keys
        #   (or the empty string)
        one_key_regex = '[ACTG]+-1_V13B23-283_[ACD]1'
        expect_equal(
            all(
                grepl(
                    gsub('X', one_key_regex, '^(|X|X(,X)*)$'),
                    spe_new$overlap_key
                )
            ),
            TRUE
        )

        overlap_keys = colData(spe_new) |>
            as_tibble() |>
            filter(capture_area == "V13B23-283_A1", overlap_key != "") |>
            pull(overlap_key) |>
            paste(collapse = ",") |>
            strsplit(',')
        
        #   All keys overlapping V13B23-283_A1 should exist in the object and
        #   belong to capture area V13B23-283_D1 (or rarely V13B23-283_A1,
        #   because of precision issues in pixel coordinates)
        expect_equal(all(overlap_keys[[1]] %in% spe$key), TRUE)
        expect_equal(
            all(grepl('[ACTG]+-1_V13B23-283_[AD]1', overlap_keys[[1]])), TRUE
        )
    }
)
