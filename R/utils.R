#' Objective Function for Cox Model Optimization
#'
#' Computes the objective value, gradient, and Hessian for optimizing a Cox
#' proportional hazards model with spline-based baseline hazard estimation.
#'
#' @param par Numeric vector of parameters to be optimized, including spline
#'        coefficients for the baseline hazard (`alpha`) and regression
#'        coefficients (`beta`).
#' @param time Numeric vector of observed survival times.
#' @param status Numeric vector of event indicators.
#' @param x Numeric matrix of covariates (features) for the model.
#' @param boundary Numeric vector specifying the boundary knots for the splines.
#' @param theta Numeric vector of prior parameters for regularization.
#' @param hessian_prev Numeric matrix of prior Hessian for regularization.
#'
#' @return A list containing:
#'   \item{value}{Scalar value of the objective function.}
#'   \item{gradient}{Numeric vector of the gradient of the objective function.}
#'   \item{hessian}{Numeric matrix of the Hessian of the objective function.}
#'
#' @details
#' The objective function combines two components:
#' - A regularization term based on prior parameters (`theta`) and prior Hessian
#'   (`hessian`).
#' - A likelihood-based term for the Cox proportional hazards model, where the
#'   baseline hazard is modeled using Bernstein polynomials.
#'
#' The function calculates the loss, gradient, and Hessian for both components
#' and combines them to return the overall objective value, gradient, and
#' Hessian. The baseline hazard is integrated using helper functions
#' (`int_basehaz`) to compute cumulative hazard and its derivatives.
objective <- function(par, time, status, x, boundary, theta, hessian_prev) {
  n_parameters <- length(par)
  n_features <- ncol(x)
  n_basis <- n_parameters - n_features

  n_parameters_pre <- length(theta)
  n_basis_pre <- n_parameters_pre - n_features

  # Transform the parameters and Hessian matrix from n_basis_pre to n_basis
  prox <- prox_forward(n_basis, n_basis_pre)
  prox <- rbind(
    cbind(prox, matrix(0, nrow(prox), n_features)),
    cbind(matrix(0, n_features, ncol(prox)), diag(n_features))
  )
  theta <- qr.solve(prox, theta)
  hess1 <- t(prox) %*% hessian_prev %*% prox

  loss1 <- as.numeric((par - theta) %*% hess1 %*% (par - theta) / 2)
  grad1 <- as.vector(hess1 %*% (par - theta))

  alpha <- par[seq_len(n_basis)]
  beta <- par[(n_basis + 1):n_parameters]
  b <- if (n_basis == 1) {
    matrix(1, nrow = length(time), ncol = 1)
  } else {
    splines2::bpoly(
      time,
      degree = n_basis - 1, intercept = TRUE,
      Boundary.knots = boundary
    )
  }
  eta <- x %*% beta

  # Vectorized computation of cumulative baseline hazard and derivatives
  int_result <- int_basehaz(as.vector(time), alpha, n_basis, boundary)
  cbh <- matrix(int_result$cbh, ncol = 1)
  dcbh <- int_result$dcbh
  d2cbh <- int_result$d2cbh

  dalpha <- t(dcbh) %*% exp(eta) - t(b) %*% status
  dbeta <- t(x) %*% (cbh * exp(eta) - status)
  d2alpha <- matrix(colSums(sweep(d2cbh, 1, exp(eta), "*")), n_basis)
  d2beta <- t(x) %*% (as.vector(cbh * exp(eta)) * x)
  dalph_dabeta <- t(dcbh) %*% (as.vector(exp(eta)) * x)

  loss2 <- sum(cbh * exp(eta) - status * (b %*% alpha + eta))
  grad2 <- c(dalpha, dbeta)
  hess2 <- cbind(
    rbind(d2alpha, t(dalph_dabeta)),
    rbind(dalph_dabeta, d2beta)
  )

  res <- loss1 + loss2
  attr(res, "gradient") <- grad1 + grad2
  attr(res, "hessian") <- hess1 + hess2
  res
}

#' Compute Integrated Baseline Hazard and Its Derivatives
#'
#' Computes the integrated baseline hazard and its derivatives for multiple
#' time points using Gaussian quadrature and Bernstein polynomial basis.
#'
#' @param times Numeric vector of time points (upper limits of integration).
#' @param alpha Numeric vector of baseline hazard coefficients.
#' @param n_basis Integer number of Bernstein polynomial basis functions.
#' @param boundary Numeric vector of length 2 for boundary knots.
#' @param n_nodes Integer number of quadrature nodes (default 10).
#'
#' @return A list containing:
#'   \item{cbh}{Numeric vector of cumulative baseline hazards.}
#'   \item{dcbh}{Matrix of gradients (n_times x n_basis).}
#'   \item{d2cbh}{Matrix of vectorized Hessians (n_times x n_basis^2).}
#'
#' @keywords internal
#' @importFrom statmod gauss.quad.prob
int_basehaz <- function(times, alpha, n_basis, boundary, n_nodes = 5) {
  n_times <- length(times)

  # Get quadrature nodes for unit interval [0, 1]
  quad_unit <- gauss.quad.prob(n_nodes, "uniform", l = 0, u = 1)
  nodes_unit <- quad_unit$nodes
  weights_unit <- quad_unit$weights

  # Handle trivial case
  if (n_basis == 1) {
    bh_val <- exp(alpha[1])
    cbh <- times * bh_val
    dcbh <- matrix(cbh, ncol = 1)
    d2cbh <- matrix(cbh, ncol = 1)
    return(list(cbh = cbh, dcbh = dcbh, d2cbh = d2cbh))
  }

  # Create all quadrature nodes: nodes_all[i, j] = times[i] * nodes_unit[j]
  nodes_all <- outer(times, nodes_unit)

  # Flatten and compute basis functions in one call
  nodes_flat <- as.vector(t(nodes_all))
  b_flat <- splines2::bpoly(
    nodes_flat,
    degree = n_basis - 1, intercept = TRUE,
    Boundary.knots = boundary
  )

  # Compute baseline hazard and weighted values
  bh_flat <- exp(as.vector(b_flat %*% alpha))
  w_bh_flat <- bh_flat * rep(weights_unit, n_times)

  # Reshape weighted hazard to n_times x n_nodes
  w_bh_mat <- matrix(w_bh_flat, nrow = n_times, ncol = n_nodes, byrow = TRUE)

  # cbh: sum of weighted hazard
  cbh <- times * rowSums(w_bh_mat)

  # Group index for rowsum
  group_idx <- rep(seq_len(n_times), each = n_nodes)

  # Gradient: use rowsum
  b_weighted <- b_flat * w_bh_flat
  dcbh <- rowsum(b_weighted, group_idx) * times

  # Hessian: vectorized outer product
  b_outer <- b_flat[, rep(seq_len(n_basis), n_basis)] *
    b_flat[, rep(seq_len(n_basis), each = n_basis)]
  b_outer_weighted <- b_outer * w_bh_flat
  d2cbh <- rowsum(b_outer_weighted, group_idx) * times

  list(cbh = cbh, dcbh = dcbh, d2cbh = d2cbh)
}

#' Proximal Operator (Forward)
#'
#' Constructs a transformation matrix to map the coefficients of Bernstein
#' polynomial basis functions from one dimension to another, ensuring
#' consistency with the properties of Bernstein polynomial bases. Note that
#' the original dimension must be less than or equal to the target dimension.
#'
#'
#' @param p Integer. The number of basis functions in the original dimension.
#' @param q Integer. The number of basis functions in the target dimension.
#'
#' @return A matrix representing the transformation from the original dimension
#'         to the target dimension.
#'
#' @details This function iteratively constructs the transformation matrix by
#'          expanding the basis functions from \code{p} to \code{q}, while
#'          preserving the mathematical properties of Bernstein polynomials.
#'
#' @seealso \code{\link{prox_reverse}}
prox_forward <- function(p, q) {
  if (
    !is.numeric(p) || !is.numeric(q) || p != as.integer(p) || q != as.integer(q)
  ) {
    stop("p and q must be integers")
  }
  if (p <= 0 || q <= 0) stop("p and q must be positive integers")
  if (p > q) stop("p must be less than or equal to q")
  if (p == q) {
    return(diag(p))
  }
  prox_forward <- diag(p)
  for (i in p:(q - 1)) {
    if (i == 1) {
      prox <- matrix(1, 2, 1)
    } else {
      prox <- matrix(0, i + 1, i)
      diag(prox) <- (i - 0:(i - 1)) / i
      diag(prox[-1, ]) <- (1:i) / i
    }
    prox_forward <- prox %*% prox_forward
  }
  prox_forward
}

#' Proximal Operator (Reverse)
#'
#' Constructs a transformation matrix to map the coefficients of Bernstein
#' polynomial basis functions from one dimension to another, ensuring
#' consistency with the properties of Bernstein polynomial bases. Note that
#' the original dimension must be greater than or equal to the target dimension.
#'
#' @param q Integer. The number of basis functions in the original dimension.
#' @param p Integer. The number of basis functions in the target dimension.
#'
#' @return A matrix representing the transformation from the original dimension
#'         to the target dimension.
#'
#' @details This function iteratively constructs the transformation matrix by
#'          expanding the basis functions from \code{p} to \code{q}, while
#'          preserving the mathematical properties of Bernstein polynomials.
#'
#' @seealso \code{\link{prox_forward}}
#'
#' @importFrom MASS ginv
prox_reverse <- function(q, p) {
  if (
    !is.numeric(p) || !is.numeric(q) || p != as.integer(p) || q != as.integer(q)
  ) {
    stop("p and q must be integers")
  }
  if (p <= 0 || q <= 0) stop("p and q must be positive integers")
  if (p > q) stop("p must be less than or equal to q")
  if (p == q) {
    return(diag(p))
  }
  prox_forward <- prox_forward(p, q)
  ginv(prox_forward)
}
