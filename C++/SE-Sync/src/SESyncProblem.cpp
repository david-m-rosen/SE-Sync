#include "Spectra/MatOp/SparseSymMatProd.h"
#include "Spectra/SymEigsSolver.h" // Spectra's symmetric eigensolver

#include "SESync/SESyncProblem.h"
#include "SESync/SESync_utils.h"

#include <random>

namespace SESync {

SESyncProblem::SESyncProblem(
    const measurements_t &measurements, const Formulation &formulation,
    const ProjectionFactorization &projection_factorization,
    const Preconditioner &precon, Scalar reg_chol_precon_max_cond)
    : form_(formulation), projection_factorization_(projection_factorization),
      preconditioner_(precon),
      reg_Chol_precon_max_cond_(reg_chol_precon_max_cond) {

  // Construct oriented incidence matrix for the underlying pose graph
  A_ = construct_oriented_incidence_matrix(measurements);

  // Construct B matrices (used to compute chordal initializations)
  construct_B_matrices(measurements, B1_, B2_, B3_);

  /// Set dimensions of the problem
  n_ = A_.rows();
  m_ = A_.cols();
  d_ = (!measurements.empty() ? measurements[0].t.size() : 0);
  r_ = d_;

  /// Set dimensions of the product of Stiefel manifolds in which the
  /// (generalized) rotational states lie
  SP_.set_k(d_);
  SP_.set_n(n_);
  SP_.set_p(r_);

  if (form_ == Formulation::Simplified) {
    /// Construct data matrices required for the implicit formulation of the
    /// SE-Sync problem

    // Construct rotational connection Laplacian
    LGrho_ = construct_rotational_connection_Laplacian(measurements);

    // Construct square root of the (diagonal) matrix of translational
    // measurement precisions
    DiagonalMatrix SqrtOmega =
        construct_translational_precision_matrix(measurements)
            .diagonal()
            .cwiseSqrt()
            .asDiagonal();

    // Construct Ared * SqrtOmega
    Ared_SqrtOmega_ = A_.topRows(n_ - 1) * SqrtOmega;

    // We cache the transpose of the above matrix as well to avoid having to
    // dynamically recompute this as an intermediate step each time the
    // transpose operator is applied
    SqrtOmega_AredT_ = Ared_SqrtOmega_.transpose();

    // Construct translational data matrix T
    SparseMatrix T = construct_translational_data_matrix(measurements);

    SqrtOmega_T_ = SqrtOmega * T;
    // Likewise, we also cache this transpose
    TT_SqrtOmega_ = SqrtOmega_T_.transpose();

    /// Construct matrices necessary to compute orthogonal projection onto the
    /// kernel of the weighted reduced oriented incidence matrix Ared_SqrtOmega
    if (projection_factorization_ == ProjectionFactorization::Cholesky) {
      // Compute and cache the Cholesky factor L of Ared * Omega * Ared^T
      L_.compute(Ared_SqrtOmega_ * SqrtOmega_AredT_);
    } else {
      // Compute the QR decomposition of Omega^(1/2) * Ared^T (cf. eq. (98) of
      // the tech report).Note that Eigen's sparse QR factorization can only be
      // called on matrices stored in compressed format
      SqrtOmega_AredT_.makeCompressed();

      QR_ = new SparseQRFactorization();
      QR_->compute(SqrtOmega_AredT_);
    }

    /** Compute and cache preconditioning matrices, if required */
    if (preconditioner_ == Preconditioner::Jacobi) {
      Vector diag = LGrho_.diagonal();
      Jacobi_precon_ = diag.cwiseInverse().asDiagonal();
    } else if (preconditioner_ == Preconditioner::IncompleteCholesky)
      iChol_precon_ = new IncompleteCholeskyFactorization(LGrho_);
    else if (preconditioner_ == Preconditioner::RegularizedCholesky) {
      // Compute maximum eigenvalue of LGrho

      // NB: Spectra's built-in SparseSymProduct matrix assumes that input
      // matrices are stored in COLUMN-MAJOR order
      Eigen::SparseMatrix<Scalar, Eigen::ColMajor> LGrho_col_major(LGrho_);

      Spectra::SparseSymMatProd<Scalar> op(LGrho_col_major);
      Spectra::SymEigsSolver<Scalar, Spectra::LARGEST_MAGN,
                             Spectra::SparseSymMatProd<Scalar>>
          max_eig_solver(&op, 1, 3);
      max_eig_solver.init();

      int max_iterations = 10000;
      Scalar tol = 1e-4; // We only require a relatively loose estimate here ...
      int nconv = max_eig_solver.compute(max_iterations, tol);

      Scalar lambda_max = max_eig_solver.eigenvalues()(0);
      reg_Chol_precon_.compute(
          LGrho_ +
          SparseMatrix(Vector::Constant(LGrho_.rows(),
                                        lambda_max / reg_Chol_precon_max_cond_)
                           .asDiagonal()));
    }

  } else {
    // form == Explicit
    M_ = construct_quadratic_form_data_matrix(measurements);

    /** Compute and cache preconditioning matrices, if required */
    if (preconditioner_ == Preconditioner::Jacobi) {
      Vector diag = M_.diagonal();
      Jacobi_precon_ = diag.cwiseInverse().asDiagonal();
    } else if (preconditioner_ == Preconditioner::IncompleteCholesky)
      iChol_precon_ = new IncompleteCholeskyFactorization(M_);
    else if (preconditioner_ == Preconditioner::RegularizedCholesky) {
      // Compute maximum eigenvalue of M

      // NB: Spectra's built-in SparseSymProduct matrix assumes that input
      // matrices are stored in COLUMN-MAJOR order
      Eigen::SparseMatrix<Scalar, Eigen::ColMajor> M_col_major(M_);

      Spectra::SparseSymMatProd<Scalar> op(M_col_major);
      Spectra::SymEigsSolver<Scalar, Spectra::LARGEST_MAGN,
                             Spectra::SparseSymMatProd<Scalar>>
          max_eig_solver(&op, 1, 3);
      max_eig_solver.init();

      int max_iterations = 10000;
      Scalar tol = 1e-4; // We only require a relatively loose estimate here ...
      int nconv = max_eig_solver.compute(max_iterations, tol);

      Scalar lambda_max = max_eig_solver.eigenvalues()(0);
      reg_Chol_precon_.compute(
          M_ +
          SparseMatrix(Vector::Constant(M_.rows(),
                                        lambda_max / reg_Chol_precon_max_cond_)
                           .asDiagonal()));
    }
  }
}

void SESyncProblem::set_relaxation_rank(size_t rank) {
  r_ = rank;
  SP_.set_p(r_);
}

Matrix SESyncProblem::data_matrix_product(const Matrix &Y) const {
  if (form_ == Formulation::Simplified)
    return Q_product(Y);
  else
    return M_ * Y;
}

Scalar SESyncProblem::evaluate_objective(const Matrix &Y) const {
  if (form_ == Formulation::Simplified)
    return (Y * Q_product(Y.transpose())).trace();
  else // form == Explicit
    return (Y * M_ * Y.transpose()).trace();
}

Matrix SESyncProblem::Euclidean_gradient(const Matrix &Y) const {
  if (form_ == Formulation::Simplified)
    return 2 * data_matrix_product(Y.transpose()).transpose();
  else // form == Explicit
    return 2 * Y * M_;
}

Matrix SESyncProblem::Riemannian_gradient(const Matrix &Y,
                                          const Matrix &nablaF_Y) const {
  return tangent_space_projection(Y, nablaF_Y);
}

Matrix SESyncProblem::Riemannian_gradient(const Matrix &Y) const {
  return tangent_space_projection(Y, Euclidean_gradient(Y));
}

Matrix SESyncProblem::Riemannian_Hessian_vector_product(
    const Matrix &Y, const Matrix &nablaF_Y, const Matrix &dotY) const {
  if (form_ == Formulation::Simplified)
    return SP_.Proj(Y, 2 * Q_product(dotY.transpose()).transpose() -
                           SP_.SymBlockDiagProduct(dotY, Y, nablaF_Y));
  else {
    // Euclidean Hessian-vector product
    Matrix H_dotY = 2 * dotY * M_;

    H_dotY.block(0, n_, r_, d_ * n_) = SP_.Proj(
        Y.block(0, n_, r_, d_ * n_),
        H_dotY.block(0, n_, r_, d_ * n_) -
            SP_.SymBlockDiagProduct(dotY.block(0, n_, r_, d_ * n_),
                                    Y.block(0, n_, r_, d_ * n_),
                                    nablaF_Y.block(0, n_, r_, d_ * n_)));
    return H_dotY;
  }
}

Matrix
SESyncProblem::Riemannian_Hessian_vector_product(const Matrix &Y,
                                                 const Matrix &dotY) const {
  return Riemannian_Hessian_vector_product(Y, Euclidean_gradient(Y), dotY);
}

Matrix SESyncProblem::precondition(const Matrix &Y, const Matrix &dotY) const {
  if (preconditioner_ == Preconditioner::None)
    return dotY;
  else if (preconditioner_ == Preconditioner::Jacobi)
    return tangent_space_projection(Y, Jacobi_precon_ * dotY);
  else if (preconditioner_ == Preconditioner::IncompleteCholesky)
    return tangent_space_projection(
        Y, iChol_precon_->solve(dotY.transpose()).transpose());
  else // preconditioner == RegularizedCholesky
    return tangent_space_projection(
        Y, reg_Chol_precon_.solve(dotY.transpose()).transpose());
}

Matrix SESyncProblem::tangent_space_projection(const Matrix &Y,
                                               const Matrix &dotY) const {
  if (form_ == Formulation::Simplified)
    return SP_.Proj(Y, dotY);
  else {
    // form == Explicit
    Matrix P(dotY.rows(), dotY.cols());

    // Projection of translational states is the identity
    P.block(0, 0, r_, n_) = dotY.block(0, 0, r_, n_);

    // Projection of generalized rotational states comes from the product of
    // Stiefel manifolds
    P.block(0, n_, r_, d_ * n_) =
        SP_.Proj(Y.block(0, n_, r_, d_ * n_), dotY.block(0, n_, r_, d_ * n_));

    return P;
  }
}

Matrix SESyncProblem::retract(const Matrix &Y, const Matrix &dotY) const {
  if (form_ == Formulation::Simplified)
    return SP_.retract(Y, dotY);
  else // form == Explicit
  {
    Matrix Yplus = Y;
    Yplus.block(0, 0, r_, n_) += dotY.block(0, 0, r_, n_);

    Yplus.block(0, n_, r_, d_ * n_) = SP_.retract(
        Y.block(0, n_, r_, d_ * n_), dotY.block(0, n_, r_, d_ * n_));

    return Yplus;
  }
}

Matrix SESyncProblem::round_solution(const Matrix Y) const {

  // First, compute a thin SVD of Y
  Eigen::JacobiSVD<Matrix> svd(Y, Eigen::ComputeThinV);

  Vector sigmas = svd.singularValues();
  // Construct a diagonal matrix comprised of the first d singular values
  DiagonalMatrix Sigma_d(d_);
  DiagonalMatrix::DiagonalVectorType &diagonal = Sigma_d.diagonal();
  for (size_t i = 0; i < d_; ++i)
    diagonal(i) = sigmas(i);

  // First, construct a rank-d truncated singular value decomposition for Y
  Matrix R = Sigma_d * svd.matrixV().leftCols(d_).transpose();

  Vector determinants(n_);

  // Compute the offset at which the rotation matrix blocks begin
  size_t rot_offset = (form_ == Formulation::Simplified ? 0 : n_);

  size_t ng0 = 0; // This will count the number of blocks whose
  // determinants have positive sign
  for (size_t i = 0; i < n_; ++i) {
    // Compute the determinant of the ith dxd block of R
    determinants(i) = R.block(0, rot_offset + i * d_, d_, d_).determinant();
    if (determinants(i) > 0)
      ++ng0;
  }

  if (ng0 < n_ / 2) {
    // Less than half of the total number of blocks have the correct sign, so
    // reverse their orientations

    // Get a reflection matrix that we can use to reverse the signs of those
    // blocks of R that have the wrong determinant
    Matrix reflector = Matrix::Identity(d_, d_);
    reflector(d_ - 1, d_ - 1) = -1;

    R = reflector * R;
  }

// Finally, project each dxd rotation block to SO(d)
#pragma omp parallel for
  for (size_t i = 0; i < n_; ++i)
    R.block(0, rot_offset + i * d_, d_, d_) =
        project_to_SOd(R.block(0, rot_offset + i * d_, d_, d_));

  if (form_ == Formulation::Explicit)
    return R;
  else // form == Explicit
  {
    // In this case, we also need to recover the corresponding translations
    Matrix X(d_, (d_ + 1) * n_);

    // Set rotational states
    X.block(0, n_, d_, d_ * n_) = R;

    // Recover translational states
    X.block(0, 0, d_, n_) = recover_translations(B1_, B2_, R);

    return X;
  }
}

Matrix SESyncProblem::compute_Lambda_blocks(const Matrix &Y) const {
  // Compute S * Y', where S is the data matrix defining the quadratic form
  // for
  // the specific version of the SE-Sync problem we're solving
  Matrix SYt = data_matrix_product(Y.transpose());

  // Preallocate storage for diagonal blocks of Lambda
  Matrix Lambda_blocks(d_, n_ * d_);

  // Index of the row/column at which the rotational blocks begin in matrix X
  size_t offset = (form_ == Formulation::Simplified ? 0 : n_);

#pragma omp parallel for
  for (size_t i = 0; i < n_; ++i) {
    Matrix P = SYt.block(offset + i * d_, 0, d_, Y.rows()) *
               Y.block(0, offset + i * d_, Y.rows(), d_);
    Lambda_blocks.block(0, i * d_, d_, d_) = .5 * (P + P.transpose());
  }
  return Lambda_blocks;
}

SparseMatrix SESyncProblem::compute_Lambda_from_Lambda_blocks(
    const Matrix &Lambda_blocks) const {

  size_t offset = (form_ == Formulation::Simplified ? 0 : n_);

  std::vector<Eigen::Triplet<Scalar>> elements;
  elements.reserve(d_ * d_ * n_);

  for (size_t i = 0; i < n_; ++i)
    for (size_t r = 0; r < d_; ++r)
      for (size_t c = 0; c < d_; ++c)
        elements.emplace_back(offset + i * d_ + r, offset + i * d_ + c,
                              Lambda_blocks(r, i * d_ + c));

  SparseMatrix Lambda(offset + d_ * n_, offset + d_ * n_);
  Lambda.setFromTriplets(elements.begin(), elements.end());
  return Lambda;
}

SparseMatrix SESyncProblem::compute_Lambda(const Matrix &Y) const {
  // First, compute the diagonal blocks of Lambda
  Matrix Lambda_blocks = compute_Lambda_blocks(Y);

  return compute_Lambda_from_Lambda_blocks(Lambda_blocks);
}

bool SESyncProblem::compute_S_minus_Lambda_min_eig(
    const Matrix &Y, Scalar &min_eigenvalue, Vector &min_eigenvector,
    size_t &num_iterations, size_t max_iterations,
    Scalar min_eigenvalue_nonnegativity_tolerance,
    size_t num_Lanczos_vectors) const {
  // First, compute the largest-magnitude eigenvalue of this matrix
  SMinusLambdaProdFunctor lm_op(this, Y);
  Spectra::SymEigsSolver<Scalar, Spectra::SELECT_EIGENVALUE::LARGEST_MAGN,
                         SMinusLambdaProdFunctor>
      largest_magnitude_eigensolver(&lm_op, 1,
                                    std::min(num_Lanczos_vectors, n_ * d_));
  largest_magnitude_eigensolver.init();

  int num_converged = largest_magnitude_eigensolver.compute(
      max_iterations, 1e-4, Spectra::SELECT_EIGENVALUE::LARGEST_MAGN);

  // Check convergence and bail out if necessary
  if (num_converged != 1)
    return false;

  Scalar lambda_lm = largest_magnitude_eigensolver.eigenvalues()(0);

  if (lambda_lm < 0) {
    // The largest-magnitude eigenvalue is negative, and therefore also the
    // minimum eigenvalue, so just return this solution
    min_eigenvalue = lambda_lm;
    min_eigenvector = largest_magnitude_eigensolver.eigenvectors(1);
    min_eigenvector.normalize(); // Ensure that this is a unit vector
    return true;
  }

  // The largest-magnitude eigenvalue is positive, and is therefore the
  // maximum  eigenvalue.  Therefore, after shifting the spectrum of S - Lambda
  // by -2*lambda_lm (by forming S - Lambda - 2*lambda_max*I), the  shifted
  // spectrum will lie in the interval [lambda_min(A) - 2*  lambda_max(A),
  // -lambda_max*A]; in particular, the largest-magnitude eigenvalue of  S -
  // Lambda - 2*lambda_max*I is lambda_min - 2*lambda_max, with  corresponding
  // eigenvector v_min; furthermore, the condition number sigma of S - Lambda
  // -2*lambda_max is then upper-bounded by 2 :-).

  SMinusLambdaProdFunctor min_shifted_op(this, Y, -2 * lambda_lm);

  Spectra::SymEigsSolver<Scalar, Spectra::SELECT_EIGENVALUE::LARGEST_MAGN,
                         SMinusLambdaProdFunctor>
      min_eigensolver(&min_shifted_op, 1,
                      std::min(num_Lanczos_vectors, n_ * d_));

  // If Y is a critical point of F, then Y^T is also in the null space of S -
  // Lambda(Y) (cf. Lemma 6 of the tech report), and therefore its rows are
  // eigenvectors corresponding to the eigenvalue 0.  In the case  that the
  // relaxation is exact, this is the *minimum* eigenvalue, and therefore the
  // rows of Y are exactly the eigenvectors that we're looking for.  On the
  // other hand, if the relaxation is *not* exact, then S - Lambda(Y) has at
  // least one strictly negative eigenvalue, and the rows of Y are *unstable
  // fixed points* for the Lanczos iterations.  Thus, we will take a slightly
  // "fuzzed" version of the first row of Y as an initialization for the Lanczos
  // iterations; this allows for rapid convergence in the case that the
  // relaxation is exact (since are starting close to a solution), while
  // simultaneously allowing the iterations to escape from this fixed point in
  // the case that the relaxation is not exact.
  Vector v0 = Y.row(0).transpose();
  Vector perturbation(v0.size());
  perturbation.setRandom();
  perturbation.normalize();
  Vector xinit = v0 + (.03 * v0.norm()) * perturbation; // Perturb v0 by ~3%

  // Use this to initialize the eigensolver
  min_eigensolver.init(xinit.data());

  // Now determine the relative precision required in the Lanczos method in
  // order to be able to estimate the smallest eigenvalue within an *absolute*
  // tolerance of 'min_eigenvalue_nonnegativity_tolerance'
  num_converged = min_eigensolver.compute(
      max_iterations, min_eigenvalue_nonnegativity_tolerance / lambda_lm,
      Spectra::SELECT_EIGENVALUE::LARGEST_MAGN);

  if (num_converged != 1)
    return false;

  min_eigenvector = min_eigensolver.eigenvectors(1);
  min_eigenvector.normalize(); // Ensure that this is a unit vector
  min_eigenvalue = min_eigensolver.eigenvalues()(0) + 2 * lambda_lm;
  num_iterations = min_eigensolver.num_iterations();
  return true;
}

Matrix SESyncProblem::chordal_initialization() const {
  Matrix Y;
  if (form_ == Formulation::Simplified) {
    Y = Matrix::Zero(r_, n_ * d_);
    Y.topRows(d_) = SESync::chordal_initialization(d_, B3_);
  } else // form == explicit
  {
    Y = Matrix::Zero(r_, n_ * (d_ + 1));

    // Compute rotations using chordal initialization
    Y.block(0, n_, d_, n_ * d_) = SESync::chordal_initialization(d_, B3_);

    // Recover corresponding translations
    Y.block(0, 0, d_, n_) =
        recover_translations(B1_, B2_, Y.block(0, n_, d_, n_ * d_));
  }

  return Y;
}

Matrix SESyncProblem::random_sample() const {
  Matrix Y;
  if (form_ == Formulation::Simplified)
    // Randomly sample a point on the Stiefel manifold
    Y = SP_.random_sample();
  else // form == Explicit
  {
    Y = Matrix::Zero(r_, n_ * (d_ + 1));

    // Randomly sample a set of elements on the Stiefel product manifold
    Y.block(0, n_, r_, n_ * d_) = SP_.random_sample();

    // Randomly sample a set of coordinates for the initial positions from the
    // standard normal distribution
    std::default_random_engine generator;
    std::normal_distribution<Scalar> g;

    for (size_t i = 0; i < r_; ++i)
      for (size_t j = 0; j < n_; ++j)
        Y(i, j) = g(generator);
  }

  return Y;
}

/// MINIMUM EIGENVALUE COMPUTATION STRUCT

SESyncProblem::SMinusLambdaProdFunctor::SMinusLambdaProdFunctor(
    const SESyncProblem *prob, const Matrix &Y, Scalar sigma)
    : problem_(prob), dim_(prob->dimension()), sigma_(sigma) {

  if (problem_->formulation() == Formulation::Simplified) {
    rows_ = problem_->dimension() * problem_->num_poses();
    cols_ = problem_->dimension() * problem_->num_poses();
  } else // mode == Explicit
  {
    rows_ = (problem_->dimension() + 1) * problem_->num_poses();
    cols_ = (problem_->dimension() + 1) * problem_->num_poses();
  }

  // Compute and cache this on construction
  Lambda_blocks_ = problem_->compute_Lambda_blocks(Y);
}

void SESyncProblem::SMinusLambdaProdFunctor::perform_op(Scalar *x,
                                                        Scalar *y) const {
  Eigen::Map<Vector> X(x, cols_);
  Eigen::Map<Vector> Y(y, rows_);

  Y = problem_->data_matrix_product(X);

  // Offset corresponding to the first index in X and Y associated with
  // rotational blocks
  size_t offset = (problem_->formulation() == Formulation::Simplified
                       ? 0
                       : problem_->num_poses());

#pragma omp parallel for
  for (size_t i = 0; i < problem_->num_poses(); ++i)
    Y.segment(offset + i * dim_, dim_) -=
        Lambda_blocks_.block(0, i * dim_, dim_, dim_) *
        X.segment(offset + i * dim_, dim_);

  if (sigma_ != 0)
    Y += sigma_ * X;
}
} // namespace SESync
