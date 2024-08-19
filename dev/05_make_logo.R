library("spatialLIBD")

## Grab SpatialExperiment with normalized counts
spe_norm <- fetch_data(type = "visiumStitched_brain_spe")

## PRECAST k = 2 clusters
spe_norm$precast_k2_stitched <- factor(spe_norm$precast_k2_stitched)
pdf(here::here("dev", "logo.pdf"), height = 7)
vis_clus(
    spe_norm,
    clustervar = "precast_k2_stitched",
    is_stitched = TRUE,
    colors = c(
        "1" = "gold",
        "2" = "darkblue",
        "NA" = "white"
    ),
    spatial = FALSE
)
dev.off()
