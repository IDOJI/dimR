algorithm_2 <- function(Xs, r, epsilon = 1e-6, max_iter = 100) {
  n <- nrow(Xs[[1]])
  U <- qr.Q(qr(matrix(rnorm(n * r), n, r)))  # Random orthonormal initial U
  Ys <- lapply(Xs, function(X) t(U) %*% X %*% U)  # Initial Ys using U

  for (iter in 1:max_iter) {
    # Compute the matrix for SVD
    SVD_matrix <- Reduce("+", lapply(seq_along(Xs), function(i) {
      Ys[[i]] %*% t(U) %*% Xs[[i]]
    }))

    # Perform the SVD on the matrix
    svd_result <- svd(SVD_matrix)
    U_new <- svd_result$u %*% t(svd_result$v)

    # Compute Y^{k+1}_t
    Ys_new <- lapply(Xs, function(X) t(U_new) %*% X %*% U_new)

    # Check for convergence using the objective function f(U)
    f_U_new <- sum(sapply(Ys_new, function(Y) sum(diag(Y))^2))
    f_U <- sum(sapply(Ys, function(Y) sum(diag(Y))^2))

    if (abs(f_U_new - f_U) / f_U < epsilon) {
      break  # Convergence achieved
    }

    # Update U and Ys for the next iteration
    U <- U_new
    Ys <- Ys_new
  }

  return(list(U = U, Ys = Ys))
}

# Example usage:
# Assuming Xs is your list of covariance matrices and r is the target rank
# Replace the following with your actual data
# Xs = list(matrix(rnorm(100), 10, 10), matrix(rnorm(100), 10, 10), ...)
# r = 2
# result <- algorithm_2(Xs, r)
