# This file contains methods for anisotropic SPDEs (precision matrix and simulation)

#' @title Anisotropic field methods
#' @description
#' `r lifecycle::badge("experimental")`
#' Methods for anisotropic SPDEs and GMRFs.
#' @name fm_aniso
NULL

#' @importFrom lifecycle deprecated
#' @import methods
#' @import Matrix
#' @import fmesher
#' @importFrom stats optim
#'
#' @describeIn fm_aniso
#' Construct the (sparse) precision matrix for the basis weights of anisotropic
#' SPDE models.
#'
#' @param x A mesh object, e.g. from `fm_mesh_2d()`.
#' @param aniso List `[kappa,vec]` where `kappa` controls the (inverse) correlation range
#' and (the half angle version of) `vec` controls the main directions of the anisotropy
#' @param log_sigma The logarithm of the marginal variance sigma (is a constant)
#' @export
#' @examples
#' mesh <- fm_mesh_2d_inla(boundary = fm_extensions(cbind(2, 1), convex = 1, 2))
#' kappa <- 1
#' v <- c(1, 2)
#' aniso <- list(rep(kappa, mesh$n), matrix(v, mesh$n, 2))
#' Q <- fm_aniso_precision(mesh, aniso)
fm_aniso_precision <- function(x, aniso, log_sigma = 0) {
  sigma <- exp(log_sigma)
  scaling <- 1 / (4 * pi * sigma^2) # Calculates scaling so that Q_fem * scaling has variance sigma
  fem <- fm_fem_aniso(x, aniso)
  Q <- (fem$c0 + 2 * fem$g1 + fem$g2) * scaling
  Q
}




#' @describeIn fm_aniso
#' Simulate an anisotropic field given a mesh and
#' anisotropic parameters, and optionally evaluate at given locations.
#'
#' @param loc locations to evaluate the random field, compatible with
#' `fm_evaluate(x, loc = loc, field = ...)`
#'
#' @return `fm_aniso_sample()` returns a matrix, where each column is a sampled
#' field. If `loc` is `NULL`, the `fm_dof(mesh)` basis weights are given.
#' Otherwise, the evaluated field at the `nrow(loc)` locations `loc` are given.
#' @export

fm_aniso_sample <- function(x, aniso, n = 1, loc = NULL) {
  Q <- fm_aniso_precision(x, aniso)
  x <- fm_sample(n = n, Q = Q)
  if (!is.null(loc)) {
    x <- fm_evaluate(x, loc = loc, field = x)
  }
  x
}

#' @describeIn fm_aniso
#' Simulates the basis weights u_i in \eqn{u(x) = \sum(u_j \phi_j(x))}{u(x) = sum(u_j phi_j(x))}
#'
#' @param x A 2d mesh object.
#' @param n Number of samples (default=1).
#' @param aniso List `[kappa,vec]` where `kappa` controls the (inverse) correlation range
#' and (the half angle version of) `vec` controls the main directions of the anisotropy.
#'
#' @return `fm_aniso_sample()` returns a vector whose j_th component is a sample of the
#' weight u_j of the jth basis vector.
#' @export

fm_aniso_basis_weights_sample <- function(x, aniso, n = 1, log_sigma = 0) {
  Q <- fm_aniso_precision(x, aniso, log_sigma) # Calculate the precision
  L_solve <- function(fact, b) {
    Matrix::solve(fact, Matrix::solve(fact, b, system = "P"), system = "L")
  }
  Lt_solve <- function(fact, b) {
    Matrix::solve(fact, Matrix::solve(fact, b, system = "Lt"), system = "Pt")
  }
  # Find P and L such that Q = P' L L' P
  fact <- Matrix::Cholesky(Q, perm = TRUE, LDL = FALSE)
  # To find the value at the weights we need to solve PLu = w
  u <- Lt_solve(
    fact,
    Matrix::Matrix(stats::rnorm(n * nrow(Q)), nrow(Q), n)
  )
  return(u)
}
