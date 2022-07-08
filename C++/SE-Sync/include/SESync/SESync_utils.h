/** This file provides a convenient set of utility functions for reading in a
set of pose-graph SLAM measurements and constructing the corresponding data
matrices used in the SE-Sync algorithm.
 *
 * Copyright (C) 2016 - 2022 by David M. Rosen (dmrosen@mit.edu)
 */

#pragma once

#include <string>

#include <Eigen/Sparse>

#include "SESync/RelativePoseMeasurement.h"
#include "SESync/SESync_types.h"

namespace SESync {

/** Given the name of a file containing a description of a special Euclidean
 * synchronization problem expressed in the .g2o format (i.e. using "EDGE_SE2 or
 * EDGE_SE3:QUAT" measurements), this function constructs and returns the
 * corresponding vector of RelativePoseMeasurements, and reports the total
 * number of poses in the pose-graph */
measurements_t read_g2o_file(const std::string &filename, size_t &num_poses);

/** Given a vector of relative pose measurements, this function constructs and
 * returns the Laplacian of the rotational weight graph L(W^rho) */
SparseMatrix
construct_rotational_weight_graph_Laplacian(const measurements_t &measurements);

/** Given a vector of relative pose measurements, this function constructs and
 * returns the Laplacian of the translational weight graph L(W^tau) */
SparseMatrix construct_translational_weight_graph_Laplacian(
    const measurements_t &measurements);

/** Given a vector of relative pose measurements, this function constructs and
 * returns the corresponding rotational connection Laplacian */
SparseMatrix
construct_rotational_connection_Laplacian(const measurements_t &measurements);

/** Given a vector of relative pose measurements, this function construct and
 * returns the associated oriented incidence matrix A */
SparseMatrix
construct_oriented_incidence_matrix(const measurements_t &measurements);

/** Given a vector of relative pose measurements, this function constructs and
 * returns the associated diagonal matrix of translational measurement
 * precisions */
DiagonalMatrix
construct_translational_precision_matrix(const measurements_t &measurements);

/** Given a vector of relative pose measurements, this function constructs and
 * returns the associated matrix of raw translational measurements */
SparseMatrix
construct_translational_data_matrix(const measurements_t &measurements);

/** Given a vector of relative pose measurements, this function constructs and
 * returns the matrix B3 defined in equation (69c) of the SE-Sync tech report */
SparseMatrix construct_B3_matrix(const measurements_t &measurements);

/** Given a vector of relative pose measurements, this function constructs and
 * returns the matrices B1 and B2 defined in equation (69) of the tech report */
void construct_B1_B2_matrices(const measurements_t &measurements,
                              SparseMatrix &B2, SparseMatrix &B3);

/** Given a vector of relative pose measurements, this function constructs and
 * returns the matrix M parameterizing the translation-explicit formulation of
 * the special Euclidean synchronization problem (Problem 2 in the SE-Sync tech
 * report) */
SparseMatrix construct_M_matrix(const measurements_t &measurements);

/** Given the measurement matrix B3 defined in equation (69c) of the tech report
 * and the problem dimension d, this function computes and returns the
 * corresponding chordal initialization for the rotational states */
Matrix chordal_initialization(size_t d, const SparseMatrix &B3);

/** Given the measurement matrices B1 and B2 and a matrix R of rotational state
 * estimates, this function computes and returns the corresponding optimal
 * translation estimates */
Matrix recover_translations(const SparseMatrix &B1, const SparseMatrix &B2,
                            const Matrix &R);

/** Given a square d x d matrix, this function returns a closest element of
 * SO(d) */
Matrix project_to_SOd(const Matrix &M);

/** Given two matrices X, Y in SO(d)^n, this function computes and returns the
 * orbit distance d_S(X,Y) between them and (optionally) the optimal
 * registration G_S in SO(d) aligning Y to X, as described in Appendix C.1 of
 * the SE-Sync tech report.
 */
Scalar dS(const Matrix &X, const Matrix &Y, Matrix *G_S = nullptr);

/** Given two matrices X, Y in O(d)^n, this function computes and returns the
 * orbit distance d_O(X,Y) between them and (optionally) the optimal
 * registration G_O in O(d) aligning Y to X, as described in Appendix C.1 of the
 * SE-Sync tech report.
 */
Scalar dO(const Matrix &X, const Matrix &Y, Matrix *G_O = nullptr);

/** This function implements the fast solution verification method (Algorithm 3)
 * described in the paper "Accelerating Certifiable Estimation with
 * Preconditioned Eigensolvers".
 *
 * Given a symmetric sparse matrix S, this function returns a Boolean value
 * indicating whether the regularized matrix M := S + eta * I is
 * positive-semidefinite.  In the event that M is *not* PSD, this function
 * additionally computes a direction of negative curvature x of S, and its
 * associated Rayleight quotient theta := x'Sx < 0, using the LOBPCG method.
 *
 * Here:
 *
 * - nx is the size of the block to use in LOBPCG
 * - num_iters is the number of iterations LOBPCG executed
 * - max_iters is the maximum number of LOBPCG iterations
 * - 'max_fill_factor' and 'drop_tol' are parameters controlling the sparsity of
 *   the incomplete symmetric indefinite factorization-based preconditioner used
 *   in conjunction with LOBPCG: each column of the inexact sparse triangular
 *   factor L is guanteed to have at most max_fill_factor * (nnz(A) / dim(A))
 *   nonzero elements, and any elements l in L_k (the kth column of L)
 *   satisfying |l| <= drop_tol * |L_k|_1 will be set to 0.
 */
bool fast_verification(const SparseMatrix &S, Scalar eta, size_t nx,
                       Scalar &theta, Vector &x, size_t &num_iters,
                       size_t max_iters = 1000, Scalar max_fill_factor = 3,
                       Scalar drop_tol = 1e-3);

} // namespace SESync
