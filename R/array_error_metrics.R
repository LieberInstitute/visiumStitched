.get_neighbors = function(i, coords) {
    this_array_row = coords[[i, 'array_row']]
    this_array_col = coords[[i, 'array_col']]
    this_capture_area = coords[[i, 'capture_area']]

    #   Grab the keys of the neighboring spots and return
    neighbor_keys = coords |>
        filter(
            capture_area == this_capture_area,
            (array_row == this_array_row & array_col == this_array_col - 2) |
            (array_row == this_array_row & array_col == this_array_col + 2) |
            (array_row == this_array_row - 1 & array_col == this_array_col - 1) |
            (array_row == this_array_row - 1 & array_col == this_array_col + 1) |
            (array_row == this_array_row + 1 & array_col == this_array_col - 1) |
            (array_row == this_array_row + 1 & array_col == this_array_col + 1)
        ) |>
        pull(key)
    return(neighbor_keys)
}

.get_shared_neighbors = function(coords_new, coords) {
    coords_new$shared_neighbors = sapply(
        seq_len(nrow(coords)),
        function(i) {
            n_before = .get_neighbors(i, coords)
            n_after = .get_neighbors(i, coords_new)
            return(mean(n_before %in% n_after))
        }
    )

    return(coords_new)
}

.add_error_metrics = function(coords, coords_new, inter_spot_dist_px) {
    coords_new = coords_new |>
        mutate(
            euclidean_error = (
                (pxl_col_in_fullres - pxl_col_in_fullres_rounded) ** 2 +
                (pxl_row_in_fullres - pxl_row_in_fullres_rounded) ** 2
            ) ** 0.5 / inter_spot_dist_px,
            capture_area = stringr::string_split_i(key, '^[ACTG]+-1_', 2),
        ) |>
        .get_shared_neighbors(coords) |>
        select(-capture_area)
    
    return(coords_new)
}
