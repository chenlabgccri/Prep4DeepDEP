% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/PrepCNA.r
\name{.PrepCNA}
\alias{.PrepCNA}
\title{Genomic binning and calculation of CNA scores}
\usage{
.PrepCNA(cna.original, filenames, exportTable = FALSE)
}
\description{
PrepCNA extracts 7460 bins (each with 10k bases) required by Prep4DeepDEP and each CNA (in log scale) is mapped onto these bins weighted by the percentage covered by CNA. The function is called by the main function of this package, Prep4DeepDEP.
}
