#include "SESync/SESyncProblem.h"
#include "SESync/SESync_utils.h"

#include "Optimization/LinearAlgebra/LOBPCG.h"

#include <random>

namespace SESync {

SESyncProblem::SESyncProblem(
    const measurements_t &measurements, const Formulation &formulation,
    const ProjectionFactorization &projection_factorization,
    const Preconditioner &precon, Scalar reg_chol_precon_max_cond)
    : form_(formulation), projection_factorization_(projection_factorization),
      preconditioner_(precon),
      reg_Chol_precon_max_cond_(reg_chol_precon_max_cond) {

  /// Construct oriented incidence matrix for the underlying pose graph
  A_ = construct_oriented_incidence_matrix(measurements);

  /// SET PROBLEM DIMENSIONS

  /// Set dimensions of the problem
  n_ = A_.rows();
  m_ = A_.cols();
  d_ = (!measurements.empty() ? measurements[0].R.rows() : 0);
  r_ = d_;

  /// Set dimensions of the product of Stiefel manifolds in which the
  /// (generalized) rotational states lie
  SP_.set_k(d_);
  SP_.set_n(n_);
  SP_.set_p(r_);

  /// Construct B matrices

  // Matrix B3 is required by all methods to construct chordal initializations
  B3_ = construct_B3_matrix(measurements);

  if (form_ != Formulation::SOSync) {
    // When solving the Simplified or Explicit forms of the problem, we also
    // require the matrices B1 and B2 to calculate chordal initializations
    // and/or recover the optimal assignment t(R) of the translational states
    // corresponding to the estimate for the robot orientations
    construct_B1_B2_matrices(measurements, B1_, B2_);

    // The full data matrix M is also needed for either form of
    // SE-synchronization
    M_ = construct_M_matrix(measurements);
  }

  /// Construct any additional auxiliary data matrices that are required
  /// (depending upon the specific problem formulation selected)

  if (form_ == Formulation::Simplified || form_ == Formulation::SOSync) {
    /// Construct data matrices required for the implicit or rotation-only
    /// formulations of the problem

    // Construct rotational connection Laplacian
    LGrho_ = construct_rotational_connection_Laplacian(measurements);

    if (form_ == Formulation::Simplified) {

      /// Construct the auxiliary data (matrices and cached factorizations)
      /// needed to compute products with the objective matrix Q

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
      /// kernel of the weighted reduced oriented incidence matrix
      /// Ared_SqrtOmega
      if (projection_factorization_ == ProjectionFactorization::Cholesky) {
        // Compute and cache the Cholesky factor L of Ared * Omega * Ared^T
        L_.compute(Ared_SqrtOmega_ * SqrtOmega_AredT_);
      } else {
        // Compute the QR decomposition of Omega^(1/2) * Ared^T (cf. eq. (98) of
        // the tech report).Note that Eigen's sparse QR factorization can only
        // be called on matrices stored in compressed format
        SqrtOmega_AredT_.makeCompressed();

        QR_ = new SparseQRFactorization();
        QR_->compute(SqrtOmega_AredT_);
      }
    } // if (form_ == Formulation::Simplified)
  }   // Auxiliary data matrix construction

  /// PRECONDITIONER CONSTRUCTION

  if (preconditioner_ == Preconditioner::Jacobi) {

    // We build a Jacobi (diagonal scaling) preconditioner by inverting the
    // diagonal of the data matrix D, depending upon the selected problem
    // formulation

    // We use the data matrix M for the translation-explicit case, and the
    // rotational connection Laplacian LGrho otherwise
    const SparseMatrix &D = (form_ == Formulation::Explicit ? M_ : LGrho_);

    Jacobi_precon_ = D.diagonal().cwiseInverse().asDiagonal();
  } else if (preconditioner_ == Preconditioner::RegularizedCholesky) {
    /// We will construct and cache a Cholesky factorization of the regularized
    /// data matrix P := D + lambda_reg * I, where the data matrix D depends
    /// upon the selected problem formulation formulation

    // We build the preconditioner from LGrho for SO-synchronization, and from
    // M for SE-synchronization
    const SparseMatrix &D = (form_ == Formulation::SOSync ? LGrho_ : M_);

    /// Next, we must estimate the spectral norm of D in order to determine
    /// the value of the regularization constant lambda_reg necessary to
    /// guarantee that the upper bound for the desired condition number of the
    /// preconditioner P is achieved

    // Here we use the fact that D >= 0, so that
    // ||D||_2 = lambda_max(D) = - lambda_min(-D)

    Optimization::LinearAlgebra::SymmetricLinearOperator<Matrix> neg_D_op =
        [&D](const Matrix &X) -> Matrix { return -(D * X); };

    // Estimate the algebraically-smallest eigenvalue of -D using LOBPCG

    size_t num_iters;
    size_t nc;
    Vector theta;
    Matrix X;
    std::tie(theta, X) = Optimization::LinearAlgebra::LOBPCG<Vector, Matrix>(
        neg_D_op,
        std::optional<
            Optimization::LinearAlgebra::SymmetricLinearOperator<Matrix>>(
            std::nullopt),
        std::optional<
            Optimization::LinearAlgebra::SymmetricLinearOperator<Matrix>>(
            std::nullopt),
        D.rows(), 4, 1, 100, num_iters, nc, 1e-2);

    // Extract estimated norm of M
    Scalar Dnorm = -theta(0);

    // Compute the required value of the regularization parameter lambda_reg
    Scalar lambda_reg = Dnorm / (reg_Chol_precon_max_cond_ - 1);

    /// Construct and factor the regularized data matrix P := D + lambda_reg * I

    // Construct regularized data matrix Mbar
    SparseMatrix P =
        D + SparseMatrix(Vector::Constant(D.rows(), lambda_reg).asDiagonal());

    // Compute and cache Cholesky factorization of Mbar
    reg_Chol_precon_.compute(P);
  } // Preconditioner construction
}

void SESyncProblem::set_relaxation_rank(size_t rank) {
  r_ = rank;
  SP_.set_p(r_);
}

Matrix SESyncProblem::data_matrix_product(const Matrix &Y) const {
  if (form_ == Formulation::Simplified)
    return Q_product(Y);
  else if (form_ == Formulation::Explicit)
    return M_ * Y;
  else // form_ == Formulation::SOSync
    return LGrho_ * Y;
}

Scalar SESyncProblem::evaluate_objective(const Matrix &Y) const {
  return (Y * data_matrix_product(Y.transpose())).trace();
}

Matrix SESyncProblem::Euclidean_gradient(const Matrix &Y) const {
  return 2 * data_matrix_product(Y.transpose()).transpose();
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
  if (form_ == Formulation::Simplified || form_ == Formulation::SOSync)
    return SP_.Proj(Y, 2 * data_matrix_product(dotY.transpose()).transpose() -
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
    return tangent_space_projection(Y, dotY * Jacobi_precon_);
  else {
    // preconditioner == RegularizedCholesky
    if (form_ != Formulation::Simplified) {
      return tangent_space_projection(
          Y, reg_Chol_precon_.solve(dotY.transpose()).transpose());
    } else {
      // When preconditioning the Simplified form of the problem (whose
      // objective matrix S is the generalized Schur complement of the data
      // matrix M with respect to the translational states), we make use of the
      // fact that for a block matrix of the form:
      //
      // M = [A  B]
      //     [B' C]
      //
      // and Schur complement S := C - B' * A^-1 * B, the product PYdot := S^-1
      // * Ydot is given by the second block of the solution Z := (X,PYdot) to
      // the following linear system:
      //
      //  M  [X]   =   [0]
      //   [PYdot] = [Ydot]

      // Allocate right-hand side
      Matrix rhs = Matrix::Zero(M_.rows(), r_);

      // Set second block to Ydot
      rhs.bottomRows(d_ * n_) = dotY.transpose();

      // Solve linear system
      Matrix Z = reg_Chol_precon_.solve(rhs);

      // Extract PYdot from Z and return
      return tangent_space_projection(Y, Z.bottomRows(d_ * n_).transpose());
    } // formulation == Simplified
  }   // preconditioner == RegularizedCholesky
}

Matrix SESyncProblem::tangent_space_projection(const Matrix &Y,
                                               const Matrix &dotY) const {
  if (form_ == Formulation::Simplified || form_ == Formulation::SOSync)
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
  if (form_ == Formulation::Simplified || form_ == Formulation::SOSync)
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
  size_t rot_offset =
      ((form_ == Formulation::Simplified || form_ == Formulation::SOSync) ? 0
                                                                          : n_);

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

  if ((form_ == Formulation::Explicit) || (form_ == Formulation::SOSync)) {
    // In this case, either the matrix R already includes the translation
    // estimates (Explicit), or we are solving the SO-Synchronization version of
    // the problem (so they are not necessary)
    return R;
  } else {
    // form_ == Simplified:  In this case, we also need to recover the
    // optimal translations corresponding to the estimated rotational states
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
  // for the specific version of the SE-Sync problem we're solving
  Matrix SYt = data_matrix_product(Y.transpose());

  // Preallocate storage for diagonal blocks of Lambda
  Matrix Lambda_blocks(d_, n_ * d_);

  // Index of the row/column at which the rotational blocks begin in matrix X
  size_t offset =
      ((form_ == Formulation::Simplified || form_ == Formulation::SOSync) ? 0
                                                                          : n_);

#pragma omp parallel for
  for (size_t i = 0; i < n_; ++i) {
    Matrix P = SYt.block(offset + i * d_, 0, d_, Y.rows()) *
               Y.block(0, offset + i * d_, Y.rows(), d_);
    Lambda_blocks.block(0, i * d_, d_, d_) = .5 * (P + P.transpose());
  }
  return Lambda_blocks;
}

SparseMatrix
SESyncProblem::compute_Lambda_from_Lambda_blocks(const Matrix &Lambda_blocks,
                                                 size_t offset) const {

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

SparseMatrix SESyncProblem::compute_Lambda_from_Lambda_blocks(
    const Matrix &Lambda_blocks) const {

  size_t offset = (form_ == Formulation::Explicit) ? n_ : 0;

  return compute_Lambda_from_Lambda_blocks(Lambda_blocks, offset);
}

SparseMatrix SESyncProblem::compute_Lambda(const Matrix &Y) const {
  // First, compute the diagonal blocks of Lambda
  Matrix Lambda_blocks = compute_Lambda_blocks(Y);

  return compute_Lambda_from_Lambda_blocks(Lambda_blocks);
}

bool SESyncProblem::verify_solution(const Matrix &Y, Scalar eta, size_t nx,
                                    Scalar &theta, Vector &x, size_t &num_iters,
                                    size_t max_LOBPCG_iters,
                                    Scalar max_fill_factor,
                                    Scalar drop_tol) const {

  /// Construct certificate matrix S

  Matrix Lambda_blocks;
  SparseMatrix S;

  if (form_ == Formulation::SOSync) {
    S = LGrho_ - compute_Lambda(Y);
  } else {
    // We compute the certificate matrix corresponding to the *full* (i.e.
    // translation-explicit) form of the problem
    Lambda_blocks = compute_Lambda_blocks(Y);
    S = M_ - compute_Lambda_from_Lambda_blocks(Lambda_blocks, n_);
  }

  /// Test positive-semidefiniteness of certificate matrix S using fast
  /// verification method
  bool PSD = fast_verification(S, eta, nx, theta, x, num_iters,
                               max_LOBPCG_iters, max_fill_factor, drop_tol);

  if (!PSD && (form_ == Formulation::Simplified)) {
    // Extract the (trailing) portion of the tangent vector corresponding to the
    // rotational states
    Vector v = x.tail(n_ * d_).normalized();
    x = v;

    // Compute x's Rayleight quotient with the simplified certificate matrix
    SparseMatrix Lambda = compute_Lambda_from_Lambda_blocks(Lambda_blocks);
    Vector Sx = data_matrix_product(x) - Lambda * x;
    theta = x.dot(Sx);
  }

  return PSD;
}

Matrix SESyncProblem::chordal_initialization() const {
  Matrix Y;
  if ((form_ == Formulation::Simplified) || (form_ == Formulation::SOSync)) {
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
  if ((form_ == Formulation::Simplified) || (form_ == Formulation::SOSync))
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

} // namespace SESync
