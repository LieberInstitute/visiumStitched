# visiumStitched 1.1.1

NEW FEATURES

* `build_SpatialExperiment()` (and `add_array_coords()`) takes a new parameter `algorithm` whose default is now "LSAP". The older implementation is available via the "Euclidean" algorithm. The now LSAP approach is a significant improvement in how new array coordinates are assigned, preventing all duplicate and most empty mappings within a single capture area, leading to improvements in downstream applications like clustering.

# visiumStitched 0.99.0

NEW FEATURES

* Initial version of the `visiumStiched` package that provides utilities for
stitching together Visium capture areas and enables seamless downstream
spatially-aware clustering methods to identify clusters on the stitched data.
An example postmortem human brain dataset composed of 3 Visium capture areas
was generated to demonstrate the utility of `visiumStitched`. The dataset is
available from `spatiaLIBD::fetch_data()` version 1.17.8 or newer. The example
data is described at <https://research.libd.org/visiumStitched_brain/>.
