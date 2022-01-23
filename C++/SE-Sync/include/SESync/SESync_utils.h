/** This file provides a convenient set of utility functions for reading in a
set of pose-graph SLAM measurements and constructing the corresponding data
matrices used in the SE-Sync algorithm.
 *
 * Copyright (C) 2016 - 2018 by David M. Rosen (dmrosen@mit.edu)
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
 * returns the matrix B1 defined in equation (69a) of the SE-Sync tech report */
SparseMatrix construct_B1_matrix(const measurements_t &measurements);

/** Given a vector of relative pose measurements, this function constructs and
 * returns the matrices B1 and B2 defined in equation (69) of the tech report */
void construct_B2_B3_matrices(const measurements_t &measurements,
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

} // namespace SESync
