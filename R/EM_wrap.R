#' Function to Impute Values 0 and 1 with Nearest Second Largest/Smallest Values
#'
#' This function replaces occurrences of 1 with the second-largest unique value and
#' occurrences of 0 with the second-smallest unique value in the dataset.
#'
#' @param data A vector or matrix containing the data with 0s and 1s to be replaced.
#'
#' @return Returns the modified dataset where 0s and 1s are imputed.
#' @export
#'
#' @examples
#' impute01(c(0, 1, 2, 3, 4, 5))  # Example usage
impute01 <- function(data){
  # Extract unique values and sort them
  unique_values <- unique(c(data))
  unique_values <- sort(unique_values)

  # Identify the second-largest and second-smallest values
  second_largest <- unique_values[length(unique_values) - 1]
  second_smallest <- unique_values[2]

  # Replace 1 with second-largest value and 0 with second-smallest value
  data[data == 1] <- second_largest
  data[data == 0] <- second_smallest

  # Return the modified data
  return(data)
}

#' Adjust Beta Values to Lie Between 0 and 1
#'
#' This function adjusts the input number to ensure it falls between 0 and 1. If the number
#' is less than 0, it returns 0, and if greater than 1, it returns 1.
#'
#' @param number A numeric value to adjust between 0 and 1.
#'
#' @return Returns the adjusted number constrained between 0 and 1.
#' @export
#'
#' @examples
#' addjust_beta(1.2)  # Will return 1
addjust_beta = function(number){
  if(number < 0){
    return(0)  # Return 0 if number is less than 0
  }
  if(number > 1){
    return(1)  # Return 1 if number is greater than 1
  }
  return(number)  # Otherwise, return the number as is
}

#' EM Algorithm Wrapper Function for Deconvolution of Methylation Data
#'
#' This function wraps around the core EM deconvolution logic, handling array or sequencing
#' data and processing the inputs. It supports parallel processing and allows for different
#' types of input data (arrays or sequencing data) to be processed accordingly.
#'
#' @param type Specifies the type of DNAm data being used ('array' or 'seq').
#' @param bulk_beta Bulk methylation beta values (if applicable).
#' @param ref_beta Reference DNAm beta values.
#' @param bulk_meth Bulk DNAm methylated counts data.
#' @param bulk_cov Bulk DNAM coverage data.
#' @param ref_meth Reference DNAm methylated counts data.
#' @param ref_cov Reference coverage data.
#' @param bulk_rna Bulk RNA data.
#' @param ref_rna Reference RNA data.
#' @param num_iterations Number of EM iterations to perform.
#' @param bulk_fac Bulk factor used in processing.
#' @param mapping Optional mapping parameter.
#' @param stopping Stopping criterion for EM convergence.
#' @param stopping_pi Stopping criterion for pi convergence.
#' @param cell_size Cell size of RNA deconvolution.
#' @param omega Parameter to choose between DNAm (1) part and RNA (0) part.
#' @param parallel Logical flag indicating whether to enable parallel computation. If TRUE, parallel computation will be used.
#' @param ncore Number of cores for parallel processing.
#' @return Returns a list of estimated deconvolution results.
#' @export
#'
#' @examples
#' # Example usage for sequencing data:
#' result <- EM_wrap("seq", bulk_beta = bulk_beta_matrix, ref_beta = ref_beta_matrix, parallel = TRUE)
EM_wrap <- function(type, bulk_beta = NULL, ref_beta = NULL, bulk_meth = NULL, bulk_cov = NULL,
                    ref_meth = NULL, ref_cov = NULL, bulk_rna = NULL, ref_rna = NULL,
                    num_iterations = 150,  bulk_fac, mapping = NULL,
                    stopping = 1e-4, stopping_pi = 5e-4, cell_size, omega = NULL, parallel = TRUE,
                    ncore = 90){  # Added parallel parameter

  library(qsmooth)

  # Process data based on type ('array' or 'seq')
  if(type == "array"){
    common_g <- intersect(rownames(bulk_beta), rownames(ref_meth))
    bulk_beta <- bulk_beta[common_g,]
    ref_meth <- ref_meth[common_g,]
    ref_cov <- ref_cov[common_g,]
    ref_beta <- ref_beta[common_g,]
    bulk_meth <- bulk_beta
    bulk_cov <- matrix(1, nrow = nrow(bulk_beta), ncol = ncol(bulk_beta))

  } else if(type == "seq"){
    common_g <- intersect(rownames(bulk_meth), rownames(ref_meth))
    bulk_meth <- bulk_meth[common_g,]
    bulk_cov <- bulk_cov[common_g,]
    ref_beta <- ref_beta[common_g,]
    ref_meth <- ref_meth[common_g,]
    ref_cov <- ref_cov[common_g,]
    bulk_beta <- bulk_meth / bulk_cov  # Calculate beta from meth and cov

  } else {
    # Default case for 'bulk' data type
    common_g <- intersect(rownames(bulk_beta), rownames(ref_beta))
    bulk_beta <- bulk_beta[common_g,]
    ref_beta <- ref_beta[common_g,]
    bulk_meth <- bulk_beta
    bulk_cov <- matrix(1, nrow = nrow(bulk_beta), ncol = ncol(bulk_beta))
    ref_meth <- ref_beta
    ref_cov <- matrix(1, nrow = nrow(ref_beta), ncol = ncol(ref_beta))
    rownames(ref_cov) <- rownames(ref_beta)
    colnames(ref_cov) <- colnames(ref_beta)
  }

  # library(preprocessCore)

  # Quantile normalization of bulk and reference beta values
  dat_qn <- normalize.quantiles(cbind(bulk_beta, ref_beta), keep.names = TRUE)
  bulk_beta <- dat_qn[, 1:ncol(bulk_beta)]
  ref_beta <- dat_qn[, (ncol(bulk_beta) + 1):ncol(dat_qn)]
  bulk_meth <- bulk_beta * bulk_cov
  ref_meth <- ref_beta * ref_cov


  if(parallel){  # Check if parallel computation is enabled
    # Progress bar setup for monitoring parallel computation
    pb <- progress_bar$new(
      format = "Current : :current [:bar] :elapsed | percent: :percent",
      total = ncol(bulk_meth),
      clear = FALSE,
      force = TRUE,
      width = 60
    )

    progress_letter <- rep(1:10, 10)  # Tokens reported in progress bar

    progress <- function(n){
      pb$tick(tokens = list(letter = progress_letter[n]))
    }

    opts <- list(progress = progress)

    cl <- makeCluster(ncore, type = "SOCK")
    registerDoSNOW(cl)

    # Run EM algorithm in parallel for each column of bulk_meth
    system.time(res <- foreach(i = 1:ncol(bulk_meth), .options.snow = opts, .errorhandling = 'pass') %dopar% {
      #library(DescTools)

      # Call the EM deconvolution function for each sample
      est_comb_new <- EM_deconv(Y = bulk_meth[,i,drop = FALSE], D_ng = bulk_cov[,i,drop = FALSE],
                                Y_ref = ref_cov, D_ng_ref = ref_cov, my_pi = ref_beta, bulk_rna = bulk_rna[,i,drop = FALSE],
                                ref_rna = ref_rna, num_iterations = num_iterations,
                                rescale = TRUE,
                                gspecific = FALSE, theta_g = FALSE, omega = omega,
                                cell_size = cell_size)

      return(list(est_comb_new = est_comb_new))
    })

    stopCluster(cl)  # Stop the cluster after computation

  } else {
    # Non-parallel implementation (sequential computation)
    res <- list()
    for(i in 1:ncol(bulk_meth)){
      est_comb_new <- EM_deconv(Y = bulk_meth[,i,drop = FALSE], D_ng = bulk_cov[,i,drop = FALSE],
                                Y_ref = ref_cov, D_ng_ref = ref_cov, my_pi = ref_beta, bulk_rna = bulk_rna[,i,drop = FALSE],
                                ref_rna = ref_rna, num_iterations = num_iterations,
                                rescale = TRUE,
                                gspecific = FALSE, theta_g = FALSE, omega = omega,
                                cell_size = cell_size)

      res[[i]] <- list(est_comb_new = est_comb_new)
    }
  }

  return(res)  # Return the deconvolution results
}
