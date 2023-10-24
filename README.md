
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatialNAcUtils

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/LieberInstitute/spatialNAcUtils)](https://github.com/LieberInstitute/spatialNAcUtils/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/LieberInstitute/spatialNAcUtils)](https://github.com/LieberInstitute/spatialNAcUtils/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![R-CMD-check-bioc](https://github.com/LieberInstitute/spatialNAcUtils/actions/workflows/check-bioc.yml/badge.svg)](https://github.com/LieberInstitute/spatialNAcUtils/actions/workflows/check-bioc.yml)
[![Codecov test
coverage](https://codecov.io/gh/LieberInstitute/spatialNAcUtils/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/LieberInstitute/spatialNAcUtils?branch=devel)
<!-- badges: end -->

`spatialNAcUtils` is intended to provide re-usable functions helpful for
various analyses in the [spatialNAc
project](https://github.com/LieberInstitute/spatial_NAc). Code is
maintained here in an R package for ease of use, and to ensure dependent
`spatialNAc` code can easily refer to the latest functions here.

For details, check the [documentation
site](http://research.libd.org/spatialNAcUtils/).

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `spatialNAcUtils` from
[GitHub](https://github.com/LieberInstitute/spatialNAcUtils) using the
following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("LieberInstitute/spatialNAcUtils")
```

## Citation

Below is the citation output from using `citation('spatialNAcUtils')` in
R. Please run this yourself to check for any updates on how to cite
**spatialNAcUtils**.

``` r
print(citation("spatialNAcUtils"), bibtex = TRUE)
#> Warning in citation("spatialNAcUtils"): could not determine year for
#> 'spatialNAcUtils' from package DESCRIPTION file
#> To cite package 'spatialNAcUtils' in publications use:
#> 
#>   Eagles N (????). _spatialNAcUtils: Functions Useful for the
#>   spatial_NAc Project_. R package version 0.99.0.
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Manual{,
#>     title = {spatialNAcUtils: Functions Useful for the spatial_NAc Project},
#>     author = {Nicholas J. Eagles},
#>     note = {R package version 0.99.0},
#>   }
```

Please note that the `spatialNAcUtils` was only made possible thanks to
many other R and bioinformatics software authors, which are cited either
in the vignettes and/or the paper(s) describing this package.

## Code of Conduct

Please note that the `spatialNAcUtils` project is released with a
[Contributor Code of
Conduct](http://bioconductor.org/about/code-of-conduct/). By
contributing to this project, you agree to abide by its terms.

## Development tools

- Continuous code testing is possible thanks to [GitHub
  actions](https://www.tidyverse.org/blog/2020/04/usethis-1-6-0/)
  through *[usethis](https://CRAN.R-project.org/package=usethis)*,
  *[remotes](https://CRAN.R-project.org/package=remotes)*, and
  *[rcmdcheck](https://CRAN.R-project.org/package=rcmdcheck)* customized
  to use [Bioconductor’s docker
  containers](https://www.bioconductor.org/help/docker/) and
  *[BiocCheck](https://bioconductor.org/packages/3.17/BiocCheck)*.
- Code coverage assessment is possible thanks to
  [codecov](https://codecov.io/gh) and
  *[covr](https://CRAN.R-project.org/package=covr)*.
- The [documentation
  website](http://LieberInstitute.github.io/spatialNAcUtils) is
  automatically updated thanks to
  *[pkgdown](https://CRAN.R-project.org/package=pkgdown)*.
- The code is styled automatically thanks to
  *[styler](https://CRAN.R-project.org/package=styler)*.
- The documentation is formatted thanks to
  *[devtools](https://CRAN.R-project.org/package=devtools)* and
  *[roxygen2](https://CRAN.R-project.org/package=roxygen2)*.

For more details, check the `dev` directory.

This package was developed using
*[biocthis](https://bioconductor.org/packages/3.17/biocthis)*.
