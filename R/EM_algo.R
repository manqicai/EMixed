# Initialize parameters function for imputation
#' Impute 0 and 1 Values with Second Smallest/Largest Unique Values
#'
#' This function replaces 0s in the data with the second smallest unique value,
#' and replaces 1s with the second largest unique value.
#'
#' @param data A vector or matrix containing data with 0s and 1s to impute.
#'
#' @return A vector or matrix with imputed values.
#' @export
#'
#' @examples
#' impute01(c(0, 1, 2, 3, 4, 5))  # Example usage
impute01 <- function(data){
  unique_values <- unique(c(data))  # Extract and sort unique values from the data
  unique_values <- sort(unique_values)
  second_largest <- unique_values[length(unique_values) - 1]  # Find second largest value
  second_smallest <- unique_values[2]  # Find second smallest value

  data[data == 1] <- second_largest  # Replace 1 with second largest
  data[data == 0] <- second_smallest  # Replace 0 with second smallest
  return(data)  # Return the imputed data
}

# # Load required Rcpp and RcppArmadillo packages for C++ integration
# library(Rcpp)


# Normalize matrix rows to sum to one
#' Normalize Matrix Rows to Sum to One
#'
#' This function ensures each row of the input matrix sums to 1. It handles missing and invalid values appropriately.
#'
#' @param matrix A numeric matrix to be normalized.
#'
#' @return The matrix with rows normalized to sum to 1.
#' @export
sum_to_one <- function(matrix){

  matrix[matrix  < 0] <- 0  # Set negative values to 0
  matrix <- apply(matrix, 2, function(x){
    x[is.na(x)] = 0  # Set NA values to 0
    return(x)
  })

  if(any(rowSums(matrix) == 0)){
    matrix[rowSums(matrix) == 0,] <- 1/ncol(matrix)  # For rows that sum to 0, assign uniform distribution
  }
  matrix <- matrix / rowSums(matrix)  # Normalize rows to sum to 1
  return(matrix)  # Return normalized matrix
}

#' EM Deconvolution for DNAm and RNA Data
#'
#' This function performs deconvolution using the Expectation-Maximization (EM) algorithm.
#' It updates the parameters based on input DNAm and RNA data.
#'
#' @param Y Bulk DNAm methylated count data.
#' @param D_ng Bulk coverage data.
#' @param Y_ref Reference DNAm methylated count data.
#' @param D_ng_ref Reference coverage data.
#' @param my_pi Initial Pi matrix.
#' @param omega Weighting factor between observed and reference data.
#' @param num_iterations Maximum number of iterations for the EM algorithm.
#' @param true_theta True cell-type proportions (if available, for validation).
#' @param verbose Whether to print progress updates.
#' @param gspecific Whether to rescale D_ng and Y gene-specifically.
#' @param theta_g Whether to use the same D_ng for all genes during theta updates.
#' @param bulk_rna Bulk RNA data.
#' @param ref_rna Reference RNA data.
#' @param rescale Whether to rescale data.
#' @param stopping Stopping criterion for theta updates.
#' @param cell_size A normalization factor for cell sizes.
#'
#' @return A list containing the updated theta and performance metrics.
#' @export
#'
#' @examples
#' result <- EM_deconv(Y, D_ng, Y_ref, D_ng_ref, my_pi, omega=0.5, num_iterations=150, true_theta=true_theta_matrix)
EM_deconv <- function(Y, D_ng, Y_ref, D_ng_ref, my_pi, bulk_rna, ref_rna, omega, num_iterations = 150,
                      true_theta = NULL, verbose = TRUE, rescale = FALSE, gspecific = FALSE, theta_g = FALSE,
                      stopping = 1e-4, cell_size){

  # Rescale Y and D_ng if specified
  if(rescale){
    if(gspecific){
      Y_prime <- apply(Y, 2, function(i){
        rowMeans(D_ng_ref) * i
      })
      D_ng_prime <- apply(D_ng, 2, function(i){
        rowMeans(D_ng_ref) * i
      })
    } else {
      Y_prime <- apply(Y, 2, function(i){
        mean(D_ng_ref) * i
      })
      D_ng_prime <- apply(D_ng, 2, function(i){
        mean(D_ng_ref) * i
      })
    }
  } else {
    Y_prime <- Y
    D_ng_prime <- D_ng
  }

  # Initialize variables
  K <- ncol(my_pi)  # Number of cell types
  N <- ncol(Y)  # Number of samples
  G <- nrow(Y)  # Number of genes
  I <- nrow(bulk_rna)  # Number of genes in RNA data
  psi_rna <- array(0, dim = c(N, I, K))  # Conditional probabilities for RNA data

  theta <- matrix(1/K, nrow = N, ncol = K)  # Initialize theta matrix
  theta_dna_l <- theta_rna_l <- list()  # Store theta updates for each iteration
  colnames(theta) <- colnames(my_pi)
  rownames(theta) <- colnames(Y)
  psi_1 <- array(0, dim = c(N, G, K))  # Conditional probabilities for DNAm (psi_1)
  psi_0 <- array(0, dim = c(N, G, K))  # Conditional probabilities for DNAm (psi_0)


  # Main EM loop
  diff_theta <- NA
  a <- matrix(1/K, nrow = N, ncol = K)  # Parameter for theta updates
  b <- matrix(1/K, nrow = N, ncol = K)  # Parameter for RNA theta updates

  for (iteration in 1:num_iterations) {
    theta_old <- theta
    # E-step: Calculate conditional probabilities based on current theta
    for (n in 1:N) {
      psi_1[n,,] <- t(theta[n, ] * t(my_pi)) / rowSums(t(theta[n, ] * t(my_pi)))
      psi_0[n,,] <- t(theta[n, ] * t((1-my_pi))) / rowSums(t(theta[n, ] * t((1-my_pi))))
      psi_rna[n,,] <- t(theta[n, ] * t(ref_rna)) / rowSums(t(theta[n, ] * t(ref_rna)))
    }

    for (i in 1:I) {
      for (n in 1:N) {
        psi_rna[n,i,] <- psi_rna[n,i,] * bulk_rna[i,n]  # Adjust psi_rna using bulk RNA data
      }
    }

    # M-step: Update theta using DNAm and RNA data
    for (n in 1:N) {
      for (k in 1:K) {
        if(theta_g){
          numerator <- sum(Y_prime[,n ] * psi_1[n, ,k]) + sum((D_ng_prime[,n]-Y_prime[,n]) * psi_0[n, ,k])
          denominator <- sum(Y_prime[,n ] * psi_1[n, ,] ) + sum((D_ng_prime[,n] - Y_prime[,n ]) * psi_0[n,, ] )
          numerator_rna <- sum(psi_rna[n, ,k])
          denominator_rna <- sum(psi_rna[n, ,] )
          theta[n, k] <- (omega * numerator + (1 - omega) * numerator_rna) / (omega * denominator + (1 - omega) * denominator_rna)
        } else {
          # No gene-specific rescaling
          numerator_rna <- sum(psi_rna[n, ,k])
          denominator_rna <- sum(psi_rna[n, ,] )
          numerator <- sum(Y[,n ] * psi_1[n, ,k]) + sum((D_ng[,n]-Y[,n]) * psi_0[n, ,k])
          denominator <- sum(Y[,n ] * psi_1[n, ,] ) + sum((D_ng[,n] - Y[,n ]) * psi_0[n,, ] )
          theta[n, k] <- (omega * numerator + (1 - omega) * numerator_rna) / (omega * denominator + (1 - omega) * denominator_rna)
        }
      }
    }

    # Check for convergence
    if(norm(theta_old - theta, type = "F") < stopping){break}
  }

  # Return final theta
  return(list(theta = theta))

}
