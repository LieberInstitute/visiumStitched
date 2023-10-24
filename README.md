
<!-- README.md is generated from README.Rmd. Please edit that file -->

# spatialNAcUtils

<!-- badges: start -->

[![GitHub
issues](https://img.shields.io/github/issues/LieberInstitute/spatialNAcUtils)](https://github.com/LieberInstitute/spatialNAcUtils/issues)
[![GitHub
pulls](https://img.shields.io/github/issues-pr/LieberInstitute/spatialNAcUtils)](https://github.com/LieberInstitute/spatialNAcUtils/pulls)
[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![Bioc release
status](http://www.bioconductor.org/shields/build/release/bioc/spatialNAcUtils.svg)](https://bioconductor.org/checkResults/release/bioc-LATEST/spatialNAcUtils)
[![Bioc devel
status](http://www.bioconductor.org/shields/build/devel/bioc/spatialNAcUtils.svg)](https://bioconductor.org/checkResults/devel/bioc-LATEST/spatialNAcUtils)
[![Bioc downloads
rank](https://bioconductor.org/shields/downloads/release/spatialNAcUtils.svg)](http://bioconductor.org/packages/stats/bioc/spatialNAcUtils/)
[![Bioc
support](https://bioconductor.org/shields/posts/spatialNAcUtils.svg)](https://support.bioconductor.org/tag/spatialNAcUtils)
[![Bioc
history](https://bioconductor.org/shields/years-in-bioc/spatialNAcUtils.svg)](https://bioconductor.org/packages/release/bioc/html/spatialNAcUtils.html#since)
[![Bioc last
commit](https://bioconductor.org/shields/lastcommit/devel/bioc/spatialNAcUtils.svg)](http://bioconductor.org/checkResults/devel/bioc-LATEST/spatialNAcUtils/)
[![Bioc
dependencies](https://bioconductor.org/shields/dependencies/release/spatialNAcUtils.svg)](https://bioconductor.org/packages/release/bioc/html/spatialNAcUtils.html#since)
[![R-CMD-check-bioc](https://github.com/LieberInstitute/spatialNAcUtils/actions/workflows/R-CMD-check-bioc.yaml/badge.svg)](https://github.com/LieberInstitute/spatialNAcUtils/actions/workflows/R-CMD-check-bioc.yaml)
[![Codecov test
coverage](https://codecov.io/gh/LieberInstitute/spatialNAcUtils/branch/devel/graph/badge.svg)](https://app.codecov.io/gh/LieberInstitute/spatialNAcUtils?branch=devel)
<!-- badges: end -->

The goal of `spatialNAcUtils` is to …

## Installation instructions

Get the latest stable `R` release from
[CRAN](http://cran.r-project.org/). Then install `spatialNAcUtils` from
[Bioconductor](http://bioconductor.org/) using the following code:

``` r
if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
}

BiocManager::install("spatialNAcUtils")
```

And the development version from
[GitHub](https://github.com/LieberInstitute/spatialNAcUtils) with:

``` r
BiocManager::install("LieberInstitute/spatialNAcUtils")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library("spatialNAcUtils")
## basic example code
```

What is special about using `README.Rmd` instead of just `README.md`?
You can include R chunks like so:

``` r
summary(cars)
#>      speed           dist       
#>  Min.   : 4.0   Min.   :  2.00  
#>  1st Qu.:12.0   1st Qu.: 26.00  
#>  Median :15.0   Median : 36.00  
#>  Mean   :15.4   Mean   : 42.98  
#>  3rd Qu.:19.0   3rd Qu.: 56.00  
#>  Max.   :25.0   Max.   :120.00
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date.

You can also embed plots, for example:

<img src="man/figures/README-pressure-1.png" width="100%" />

In that case, don’t forget to commit and push the resulting figure
files, so they display on GitHub!

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
