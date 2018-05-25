/** A set of typedefs describing the types of matrices and factorizations that
 * will be used in the SE-Sync algorithm
 *
 *  Copyright (C) 2016, 2017 by David M. Rosen
 */

#pragma once

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace SESync {

/** Some useful typedefs for the SE-Sync library */
typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagonalMatrix;

/** We use row-major storage order to take advantage of fast (sparse-matrix) *
 * (dense-vector) multiplications when OpenMP is available (cf. the Eigen
 * documentation page on "Eigen and Multithreading") */
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMatrix;

/** The specific formulation of special Euclidean synchronization problem to
 * solve */
enum Formulation {
  /** Construct and solve the simplified version of the special Euclidean
  * synchronization problem obtained by analytically eliminating the
  *  translational states from the estimation (cf. Problem 4 in the SE-Sync tech
  *  report).
  */
  Simplified,

  /** Construct and solve the formulation of the special Euclidean
  * synchronization problem that explicitly estimates both rotational and
  * translational states (cf. Problem 2 in the SE-Sync tech report).
  */
  Explicit
};

/** The set of available preconditioning strategies to use in the Riemannian
 * Trust Region when solving this problem */
enum Preconditioner { None, Jacobi, IncompleteCholesky, RegularizedCholesky };

} // namespace SESync
