#' EMixed Deconvolution Wrapper for RNA and DNAm Data
#'
#' This function wraps around the EM_wrap function to perform deconvolution using both RNA and DNAm data. It
#' processes the bulk RNA and DNAm data, ensures there are no rows with all zero values, and calls the EM_wrap
#' function to perform the deconvolution for each data type (RNA and DNAm). The final estimates are averaged.
#'
#' @param bulk_RNA Bulk RNA-seq data.
#' @param ref_RNA Reference RNA-seq data.
#' @param bulk_DNAm Bulk DNA methylation data.
#' @param ref_DNAm Reference DNA methylation data.
#' @param parallel Logical flag indicating whether to enable parallel computation. If TRUE, parallel computation will be used.
#' @param ncore Number of cores to use for parallel computation (default is 2).
#'
#' @import foreach
#' @import tidyverse
#' @import sparseMatrixStats
#' @import purrr
#' @import progress
#' @import quadprog
#' @import parallel
#' @import doSNOW
#' @import dplyr
#' @import doParallel
#' @import qsmooth
#' @import preprocessCore
#'
#' @return A list containing the averaged deconvolution results for RNA and DNAm.
EMixed <- function(bulk_RNA, ref_RNA, bulk_DNAm, ref_DNAm, parallel, ncore = 2) {

  # Remove rows that contain all zeros in the RNA data (both bulk and reference)
  ref_RNA <- ref_RNA[rowSums(ref_RNA) != 0, ]
  bulk_RNA <- bulk_RNA[rowSums(bulk_RNA) != 0, ]

  # Find the overlap between the CpGs in the bulk and reference DNAm datasets
  overlap_cpg <- intersect(rownames(ref_DNAm), rownames(bulk_DNAm))

  # Find the overlap between the genes in the bulk and reference RNA datasets
  overlap_gene <- intersect(rownames(ref_RNA), rownames(bulk_RNA))

  # Subset RNA data based on overlapping genes
  ref_RNA <- ref_RNA[overlap_gene, ]
  bulk_RNA <- bulk_RNA[overlap_gene, ]

  # Subset DNAm data based on overlapping CpGs
  ref_DNAm <- ref_DNAm[overlap_cpg, ]
  bulk_DNAm <- bulk_DNAm[overlap_cpg, ]

  # Correct the reference RNA data for cell size by normalizing columns
  r <- colSums(ref_RNA)
  ref_RNA_correct <- t(apply(ref_RNA, 1, function(x) x / r))

  # Perform EM deconvolution using RNA-seq data
  res_rna <- EM_wrap("arrayref", bulk_beta = bulk_DNAm, ref_beta = as.matrix(ref_DNAm),
                     num_iterations = 500, bulk_rna = bulk_RNA, ref_rna = ref_RNA_correct, cell_size = r, omega = 0,
                     ncore = ncore, parallel = parallel)

  # Extract the estimated cell type proportions from the RNA deconvolution
  p_rna <- t(sapply(res_rna, function(i) i[["est_comb_new"]]$theta))

  # Perform EM deconvolution using DNAm data
  res_dnam <- EM_wrap("arrayref", bulk_beta = bulk_DNAm, ref_beta = as.matrix(ref_DNAm),
                     num_iterations = 500, bulk_rna = bulk_RNA, ref_rna = ref_RNA_correct, cell_size = r, omega = 1,
                     ncore = ncore,parallel = parallel)

  # Extract the estimated cell type proportions from the DNAm deconvolution
  p_dnam <- t(sapply(res_dnam, function(i) i[["est_comb_new"]]$theta))

  # Compute the mean of the RNA and DNAm estimates
  p_emixed <- (p_rna + p_dnam) / 2
  rownames(p_emixed) <- colnames(bulk_RNA)
  colnames(p_emixed) <- colnames(ref_DNAm)

  # Set row and column names for RNA and DNAm results
  rownames(p_dnam) <- colnames(bulk_RNA)
  colnames(p_dnam) <- colnames(ref_DNAm)

  rownames(p_rna) <- colnames(bulk_RNA)
  colnames(p_rna) <- colnames(ref_DNAm)

  # Return a list containing the combined, RNA, and DNAm deconvolution results
  return(list(EMixed_multi = p_emixed, EMixed_DNAm = p_dnam, EMixed_RNA = p_rna))
}
