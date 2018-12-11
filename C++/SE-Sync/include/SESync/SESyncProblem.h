/** This class encapsulates an instance of the rank-restricted Riemannian form
 * of the semidefinite relaxation solved by SE-Sync (Problem 9 in the SE-Sync
 * tech report).  It contains all of the precomputed and cached data matrices
 * necessary to describe the problem and run the optimization algorithm, as well
 * as functions for performing geometric operations on the underlying manifold
 * (tangent space projection and retraction) and evaluating the optimization
 * objective and its gradient and Hessian operator.
 *
 * Copyright (C) 2016 - 2018 by David M. Rosen (dmrosen@mit.edu)
 */

#pragma once

/** Use external matrix factorizations/linear solves provided by SuiteSparse
 * (SPQR and Cholmod) */

#include <Eigen/CholmodSupport>
#include <Eigen/Dense>
#include <Eigen/SPQRSupport>
#include <Eigen/Sparse>

#include "SESync/RelativePoseMeasurement.h"
#include "SESync/SESync_types.h"
#include "SESync/SESync_utils.h"
#include "SESync/StiefelProduct.h"

namespace SESync {

/** The type of the sparse Cholesky factorization to use in the computation of
 * the orthogonal projection operation */
typedef Eigen::CholmodDecomposition<SparseMatrix> SparseCholeskyFactorization;

/** The type of the QR decomposition to use in the computation of the orthogonal
 * projection operation */
typedef Eigen::SPQR<SparseMatrix> SparseQRFactorization;

/** The type of the incomplete Cholesky decomposition we will use for
 * preconditioning the conjugate gradient iterations in the RTR method */
typedef Eigen::IncompleteCholesky<Scalar> IncompleteCholeskyFactorization;

class SESyncProblem {
private:
  /// PROBLEM DATA

  /** The specific formulation of the SE-Sync problem to be solved
(translation-implicit, translation-explicit, or robust) */
  Formulation form_;

  /** Number of poses */
  size_t n_ = 0;

  /** Number of measurements */
  size_t m_ = 0;

  /** Dimensional parameter d for the special Euclidean group SE(d) over which
   * this problem is defined */
  size_t d_ = 0;

  /** Relaxation rank */
  size_t r_ = 0;

  /** The oriented incidence matrix A encoding the underlying measurement
   * graph for this problem */
  SparseMatrix A_;

  /** The matrices B1, B2, and B3 defined in equation (69) of the SE-Sync tech
   * report */
  SparseMatrix B1_, B2_, B3_;

  /** The matrix M parameterizing the quadratic form appearing in the Explicit
   * form of the SE-Sync problem (Problem 2 in the SE-Sync tech report) */
  SparseMatrix M_;

  /** The rotational connection Laplacian for the special Euclidean
   * synchronization problem, cf. eq. 14 of the SE-Sync tech report.  Only
   * used in Implicit mode.*/
  SparseMatrix LGrho_;

  /** The weighted reduced oriented incidence matrix Ared Omega^(1/2) (cf.
   * eq. 39 of the SE-Sync tech report).  Only used in Implicit mode. */
  SparseMatrix Ared_SqrtOmega_;

  /** The transpose of the above matrix; we cache this for computational
   * efficiency, since it's used frequently.  Only used in Implicit mode.
   */
  SparseMatrix SqrtOmega_AredT_;

  /** The weighted translational data matrix Omega^(1/2) T (cf. eqs. 22-24
   * of the SE-Sync tech report.  Only used in Implicit mode. */
  SparseMatrix SqrtOmega_T_;

  /** The transpose of the above matrix; we cache this for computational
   * efficiency, since it's used frequently.  Only used in Implicit mode. */
  SparseMatrix TT_SqrtOmega_;

  /** An Eigen sparse linear solver that encodes the Cholesky factor L used
   * in the computation of the orthogonal projection function (cf. eq. 39 of the
   * SE-Sync tech report) */
  SparseCholeskyFactorization L_;

  /** An Eigen sparse linear solver that encodes the QR factorization used in
   * the computation of the orthogonal projection function (cf. eq. 98 of the
   * SE-Sync tech report) */

  // When using Eigen::SPQR, the destructor causes a segfault if this variable
  // isn't explicitly initialized (i.e. not just default-constructed)
  SparseQRFactorization *QR_ = nullptr;

  /** A Boolean variable determining whether to use the Cholesky or QR
   * decompositions for computing the orthogonal projection */
  ProjectionFactorization projection_factorization_;

  /** The preconditioning strategy to use when running the Riemannian
   * trust-region algorithm */
  Preconditioner preconditioner_;

  /** Diagonal Jacobi preconditioner */
  DiagonalMatrix Jacobi_precon_;

  /** Incomplete Cholesky Preconditioner */
  IncompleteCholeskyFactorization *iChol_precon_ = nullptr;

  /** Tikhonov-regularized Cholesky Preconditioner */
  SparseCholeskyFactorization reg_Chol_precon_;

  /** Upper-bound on the admissible condition number of the regularized
   * approximate Hessian matrix used for Cholesky preconditioner */
  Scalar reg_Chol_precon_max_cond_;

  /** The underlying manifold in which the generalized orientations lie in the
  rank-restricted Riemannian optimization problem (Problem 9 in the SE-Sync tech
  report).*/
  StiefelProduct SP_;

public:
  /// CONSTRUCTORS AND MUTATORS

  /** Default constructor; doesn't actually do anything */
  SESyncProblem() {}

  /** Basic constructor.  Here
   *
   * - measurements is a vector of relative pose measurements defining the
          pose-graph SLAM problem to be solved.
   * - formulation is an enum type specifying whether to solve the simplified
   *      form of the SDP relaxation (in which translational states have been
   *      eliminated) or the explicit form (in which the translational states
   *      are explicitly represented).
   * - projection_factorization is an enum type specifying the kind of matrix
   *      factorization to use when computing the action of the orthogonal
   *      projection operator Pi.  Only operative when solving the Simplified
   *      formulation of the special Euclidean synchronization problem
   *  - preconditioner is an enum type specifying the preconditioning strategy
   *      to employ
   */
  SESyncProblem(
      const measurements_t &measurements,
      const Formulation &formulation = Formulation::Simplified,
      const ProjectionFactorization &projection_factorization =
          ProjectionFactorization::Cholesky,
      const Preconditioner &preconditioner = Preconditioner::IncompleteCholesky,
      Scalar reg_chol_precon_max_cond = 1e6);

  /** Set the maximum rank of the rank-restricted semidefinite relaxation */
  void set_relaxation_rank(size_t rank);

  /// ACCESSORS

  /** Returns the specific formulation of this SE-Sync problem */
  Formulation formulation() const { return form_; }

  /** Returns the type of matrix factorization used to compute the action of the
   * orthogonal projection operator Pi when solving a Simplified instance of the
   * special Euclidean synchronization problem */
  ProjectionFactorization projection_factorization() const {
    return projection_factorization_;
  }

  /** Returns the preconditioning strategy */
  Preconditioner preconditioner() const { return preconditioner_; }

  /** Returns the maximum admissible condition number for the regularized
   * Cholesky preconditioner */
  Scalar regularized_Cholesky_preconditioner_max_condition() const {
    return reg_Chol_precon_max_cond_;
  }

  /** Returns the number of poses appearing in this problem */
  size_t num_poses() const { return n_; }

  /** Returns the number of measurements in this problem */
  size_t num_measurements() const { return m_; }

  /** Returns the dimensional parameter d for the special Euclidean group SE(d)
   * over which this problem is defined */
  size_t dimension() const { return d_; }

  /** Returns the relaxation rank r of this problem */
  size_t relaxation_rank() const { return r_; }

  /** Returns the oriented incidence matrix A of the underlying measurement
   * graph over which this problem is defined */
  const SparseMatrix &oriented_incidence_matrix() const { return A_; }

  /** Returns the StiefelProduct manifold underlying this SE-Sync problem */
  const StiefelProduct &Stiefel_product_manifold() const { return SP_; }

  /// OPTIMIZATION AND GEOMETRY

  /** Given a matrix X, this function computes and returns the orthogonal
   *projection Pi * X */
  // We inline this function in order to take advantage of Eigen's ability
  // to optimize matrix expressions as compile time
  inline Matrix Pi_product(const Matrix &X) const {
    if (projection_factorization_ == ProjectionFactorization::Cholesky)
      return X - SqrtOmega_AredT_ * L_.solve(Ared_SqrtOmega_ * X);
    else {
      Matrix PiX = X;
      for (size_t c = 0; c < X.cols(); c++) {
        // Eigen's SPQR support only supports solving with vectors(!) (i.e.
        // 1-column matrices)
        PiX.col(c) = X.col(c) - SqrtOmega_AredT_ * QR_->solve(X.col(c));
      }
      return PiX;
    }
  }

  /** This function computes and returns the product QX */
  // We inline this function in order to take advantage of Eigen's ability to
  // optimize matrix expressions as compile time
  inline Matrix Q_product(const Matrix &X) const {
    return LGrho_ * X + TT_SqrtOmega_ * Pi_product(SqrtOmega_T_ * X);
  }

  /** Given a matrix Y, this function computes and returns the matrix product
  SY, where S is the matrix that determines the quadratic form defining the
  objective  F(Y) := tr(S * Y' * Y) for the SE-Sync problem.  More precisely:
  *
  * If formulation == Implicit, this returns Q * Y, where Q is as defined in
  equation (24) of the SE-Sync tech report.
  *
  * If formulation == Explicit, this returns M * Y, where M is as defined in
  equation (18) of the SE-Sync tech report. */
  Matrix data_matrix_product(const Matrix &Y) const;

  /** Given a matrix Y, this function computes and returns F(Y), the value of
   * the objective evaluated at X */
  Scalar evaluate_objective(const Matrix &Y) const;

  /** Given a matrix Y, this function computes and returns nabla F(Y), the
   * *Euclidean* gradient of F at Y. */
  Matrix Euclidean_gradient(const Matrix &Y) const;

  /** Given a matrix Y in the domain D of the SE-Sync optimization problem and
   * the *Euclidean* gradient nabla F(Y) at Y, this function computes and
   * returns the *Riemannian* gradient grad F(Y) of F at Y */
  Matrix Riemannian_gradient(const Matrix &Y, const Matrix &nablaF_Y) const;

  /** Given a matrix Y in the domain D of the SE-Sync optimization problem, this
   * function computes and returns grad F(Y), the *Riemannian* gradient of F
   * at Y */
  Matrix Riemannian_gradient(const Matrix &Y) const;

  /** Given a matrix Y in the domain D of the SE-Sync optimization problem, the
   * *Euclidean* gradient nablaF_Y of F at Y, and a tangent vector dotY in
   * T_D(Y), the tangent space of the domain of the optimization problem at Y,
   * this function computes and returns Hess F(Y)[dotY], the action of the
   * Riemannian Hessian on dotY */
  Matrix Riemannian_Hessian_vector_product(const Matrix &Y,
                                           const Matrix &nablaF_Y,
                                           const Matrix &dotY) const;

  /** Given a matrix Y in the domain D of the SE-Sync optimization problem, and
   * a tangent vector dotY in T_D(Y), the tangent space of the domain of the
   * optimization problem at Y, this function computes and returns Hess
   * F(Y)[dotY], the action of the Riemannian Hessian on dotX */
  Matrix Riemannian_Hessian_vector_product(const Matrix &Y,
                                           const Matrix &dotY) const;

  /** Given a matrix Y in the domain D of the SE-Sync optimization problem, and
   * a tangent vector dotY in T_D(Y), this function applies the selected
   * preconditioning strategy to dotY */
  Matrix precondition(const Matrix &Y, const Matrix &dotY) const;

  /** Given a matrix Y in the domain D of the SE-Sync optimization problem and a
  tangent vector dotY in T_Y(E), the tangent space of Y considered as a generic
  matrix, this function computes and returns the orthogonal projection of dotY
  onto T_D(Y), the tangent space of the domain D at Y*/
  Matrix tangent_space_projection(const Matrix &Y, const Matrix &dotY) const;

  /** Given a matrix Y in the domain D of the SE-Sync optimization problem and a
   * tangent vector dotY in T_D(Y), this function returns the point Yplus in D
   * obtained by retracting along dotY */
  Matrix retract(const Matrix &Y, const Matrix &dotY) const;

  /** Given a point Y in the domain D of the rank-r relaxation of the SE-Sync
   * optimization problem, this function computes and returns a matrix X = [t |
   * R] comprised of translations and rotations for a set of feasible poses for
   * the original estimation problem obtained by rounding the point Y */
  Matrix round_solution(const Matrix Y) const;

  /** Given a critical point Y of the rank-r relaxation of the SE-Sync
   * optimization problem, this function computes and returns a d x dn matrix
   * comprised of dxd block elements of the associated block-diagonal Lagrange
   * multiplier matrix associated with the orthonormality constraints on the
   * generalized orientations of the poses (cf. eq. (119) in the SE-Sync tech
   * report) */
  Matrix compute_Lambda_blocks(const Matrix &Y) const;

  /** Given the d x dn block matrix containing the diagonal blocks of Lambda,
   * this function computes and returns the matrix Lambda itself */
  SparseMatrix
  compute_Lambda_from_Lambda_blocks(const Matrix &Lambda_blocks) const;

  /** Given a critical point Y of the rank-r relaxation of the SE-Sync
   * optimization problem, this function computes and returns the corresponding
   * Lagrange multiplier matrix Lambda */
  SparseMatrix compute_Lambda(const Matrix &Y) const;

  /** Given a critical point Y in the domain of the optimization problem, this
   *function computes the smallest eigenvalue lambda_min of S - Lambda and its
   *associated eigenvector v.  Returns a Boolean value indicating whether the
   *Lanczos method used to estimate the smallest eigenpair converged to
   *within the required tolerance. */
  bool compute_S_minus_Lambda_min_eig(
      const Matrix &Y, Scalar &min_eigenvalue, Vector &min_eigenvector,
      size_t &num_iterations, size_t max_iterations = 10000,
      Scalar min_eigenvalue_nonnegativity_tolerance = 1e-5,
      size_t num_Lanczos_vectors = 20) const;

  /** Computes and returns the chordal initialization for the rank-restricted
   * semidefinite relaxation */
  Matrix chordal_initialization() const;

  /** Randomly samples a point in the domain for the rank-restricted
   * semidefinite relaxation */
  Matrix random_sample() const;

  ~SESyncProblem() {
    if (QR_)
      delete QR_;

    if (iChol_precon_)
      delete iChol_precon_;
  }

  /// MINIMUM EIGENVALUE COMPUTATIONS

  /** This is a lightweight struct used in conjunction with Spectra to compute
   *the minimum eigenvalue and eigenvector of S - Lambda(X); it has a single
   *nontrivial function, perform_op(x,y), that computes and returns the product
   *y = (S - Lambda + sigma*I) x */
  struct SMinusLambdaProdFunctor {
    const SESyncProblem *problem_;

    // Diagonal blocks of the matrix Lambda
    Matrix Lambda_blocks_;

    // Number of rows and columns of the matrix B - Lambda
    int rows_;
    int cols_;

    // Dimensional parameter d of the special Euclidean group SE(d) over which
    // this synchronization problem is defined
    int dim_;
    Scalar sigma_;

    // Constructor
    SMinusLambdaProdFunctor(const SESyncProblem *prob, const Matrix &Y,
                            Scalar sigma = 0);

    int rows() const { return rows_; }
    int cols() const { return cols_; }

    // Matrix-vector multiplication operation
    void perform_op(Scalar *x, Scalar *y) const;
  };
};
} // namespace SESync
