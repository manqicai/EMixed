#' Reference RNA-seq Data
#'
#' This dataset contains reference RNA-seq data for deconvolution.
#' It includes gene expression levels across different cell types, which are used
#' as a reference for the EMixed deconvolution process.
#'
#' @format A matrix with genes as rows and cell types as columns. Each entry represents
#' the expression level of a specific gene in a particular cell type.
#'
#' @examples
#' data(ref_RNA)
#' head(ref_RNA)
"ref_RNA"

#' Reference DNA Methylation Data
#'
#' This dataset contains reference DNA methylation (DNAm) data for deconvolution.
#' It includes beta value across different cell types, which are used
#' as a reference for the EMixed deconvolution process.
#'
#' @format A matrix with CpGs as rows and cell types as columns. Each entry represents
#' the beta value of a specific CpG site in a particular cell type.
#'
#' @examples
#' data(ref_DNAm)
#' head(ref_DNAm)
"ref_DNAm"

#' Bulk RNA-seq Data
#'
#' This dataset contains bulk RNA-seq data for deconvolution. It represents the mixed
#' gene expression levels from a mixture of different cell types in a biological sample.
#'
#' @format A matrix with genes as rows and samples as columns. Each entry represents
#' the expression level of a specific gene in a bulk RNA sample.
#'
#' @examples
#' data(bulk_RNA)
#' head(bulk_RNA)
"bulk_RNA"

#' Bulk DNA Methylation Data
#'
#' This dataset contains bulk DNA methylation (DNAm) data for deconvolution. It represents
#' the mixed CpG methylation levels from a mixture of different cell types in a biological sample.
#'
#' @format A matrix with CpGs as rows and samples as columns. Each entry represents
#' the methylation level of a specific CpG site in a bulk DNA sample.
#'
#' @examples
#' data(bulk_DNAm)
#' head(bulk_DNAm)
"bulk_DNAm"
