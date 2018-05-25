#include "MatOp/SparseSymMatProd.h"
#include "SymEigsSolver.h" // Spectra's symmetric eigensolver

#include "SESync/SESyncProblem.h"
#include "SESync/SESync_utils.h"

#include <random>

namespace SESync {

SESyncProblem::SESyncProblem(const measurements_t &measurements,
                             const Formulation &formulation, bool Cholesky,
                             const Preconditioner &precon,
                             double reg_chol_precon_max_cond)
    : form(formulation), use_Cholesky(Cholesky), preconditioner(precon),
      RegCholPrecon_max_cond(reg_chol_precon_max_cond) {

  // Construct oriented incidence matrix for the underlying pose graph
  A = construct_oriented_incidence_matrix(measurements);

  // Construct B matrices (used to compute chordal initializations)
  construct_B_matrices(measurements, B1, B2, B3);

  /// Set dimensions of the problem
  n = A.rows();
  m = A.cols();
  d = (!measurements.empty() ? measurements[0].t.size() : 0);
  r = d;

  /// Set dimensions of the product of Stiefel manifolds in which the
  /// (generalized) rotational states lie
  SP.set_k(d);
  SP.set_n(n);
  SP.set_p(r);

  if (form == Simplified) {
    /// Construct data matrices required for the implicit formulation of the
    /// SE-Sync problem

    // Construct rotational connection Laplacian
    LGrho = construct_rotational_connection_Laplacian(measurements);

    // Construct square root of the (diagonal) matrix of translational
    // measurement precisions
    DiagonalMatrix SqrtOmega =
        construct_translational_precision_matrix(measurements)
            .diagonal()
            .cwiseSqrt()
            .asDiagonal();

    // Construct Ared * SqrtOmega
    Ared_SqrtOmega = A.topRows(n - 1) * SqrtOmega;

    // We cache the transpose of the above matrix as well to avoid having to
    // dynamically recompute this as an intermediate step each time the
    // transpose operator is applied
    SqrtOmega_AredT = Ared_SqrtOmega.transpose();

    // Construct translational data matrix T
    SparseMatrix T = construct_translational_data_matrix(measurements);

    SqrtOmega_T = SqrtOmega * T;
    // Likewise, we also cache this transpose
    TT_SqrtOmega = SqrtOmega_T.transpose();

    /// Construct matrices necessary to compute orthogonal projection onto the
    /// kernel of the weighted reduced oriented incidence matrix Ared_SqrtOmega
    if (use_Cholesky) {
      // Compute and cache the Cholesky factor L of Ared * Omega * Ared^T
      L.compute(Ared_SqrtOmega * SqrtOmega_AredT);
    } else {
      // Compute the QR decomposition of Omega^(1/2) * Ared^T (cf. eq. (98) of
      // the tech report).Note that Eigen's sparse QR factorization can only be
      // called on matrices stored in compressed format
      SqrtOmega_AredT.makeCompressed();

      QR = new SparseQRFactorization();
      QR->compute(SqrtOmega_AredT);
    }

    /** Compute and cache preconditioning matrices, if required */
    if (preconditioner == Jacobi) {
      Eigen::VectorXd diag = LGrho.diagonal();
      JacobiPrecon = diag.cwiseInverse().asDiagonal();
    } else if (preconditioner == IncompleteCholesky)
      iCholPrecon = new IncompleteCholeskyFactorization(LGrho);
    else if (preconditioner == RegularizedCholesky) {
      // Compute maximum eigenvalue of LGrho

      // NB: Spectra's built-in SparseSymProduct matrix assumes that input
      // matrices are stored in COLUMN-MAJOR order
      Eigen::SparseMatrix<double, Eigen::ColMajor> LGrho_col_major(LGrho);

      Spectra::SparseSymMatProd<double> op(LGrho_col_major);
      Spectra::SymEigsSolver<double, Spectra::LARGEST_MAGN,
                             Spectra::SparseSymMatProd<double>>
          max_eig_solver(&op, 1, 3);
      max_eig_solver.init();

      int max_iterations = 10000;
      double tol = 1e-4; // We only require a relatively loose estimate here ...
      int nconv = max_eig_solver.compute(max_iterations, tol);

      double lambda_max = max_eig_solver.eigenvalues()(0);
      RegCholPrecon.compute(
          LGrho +
          SparseMatrix(Vector::Constant(LGrho.rows(),
                                        lambda_max / RegCholPrecon_max_cond)
                           .asDiagonal()));
    }

  } else {
    // form == Explicit
    M = construct_quadratic_form_data_matrix(measurements);

    /** Compute and cache preconditioning matrices, if required */
    if (preconditioner == Jacobi) {
      Eigen::VectorXd diag = M.diagonal();
      JacobiPrecon = diag.cwiseInverse().asDiagonal();
    } else if (preconditioner == IncompleteCholesky)
      iCholPrecon = new IncompleteCholeskyFactorization(M);
    else if (preconditioner == RegularizedCholesky) {
      // Compute maximum eigenvalue of M

      // NB: Spectra's built-in SparseSymProduct matrix assumes that input
      // matrices are stored in COLUMN-MAJOR order
      Eigen::SparseMatrix<double, Eigen::ColMajor> M_col_major(M);

      Spectra::SparseSymMatProd<double> op(M_col_major);
      Spectra::SymEigsSolver<double, Spectra::LARGEST_MAGN,
                             Spectra::SparseSymMatProd<double>>
          max_eig_solver(&op, 1, 3);
      max_eig_solver.init();

      int max_iterations = 10000;
      double tol = 1e-4; // We only require a relatively loose estimate here ...
      int nconv = max_eig_solver.compute(max_iterations, tol);

      double lambda_max = max_eig_solver.eigenvalues()(0);
      RegCholPrecon.compute(
          M + SparseMatrix(Vector::Constant(M.rows(),
                                            lambda_max / RegCholPrecon_max_cond)
                               .asDiagonal()));
    }
  }
}

void SESyncProblem::set_relaxation_rank(unsigned int rank) {
  r = rank;
  SP.set_p(r);
}

Matrix SESyncProblem::data_matrix_product(const Matrix &Y) const {
  if (form == Simplified)
    return Q_product(Y);
  else
    return M * Y;
}

double SESyncProblem::evaluate_objective(const Matrix &Y) const {
  if (form == Simplified)
    return (Y * Q_product(Y.transpose())).trace();
  else // form == Explicit
    return (Y * M * Y.transpose()).trace();
}

Matrix SESyncProblem::Euclidean_gradient(const Matrix &Y) const {
  if (form == Simplified)
    return 2 * data_matrix_product(Y.transpose()).transpose();
  else // form == Explicit
    return 2 * Y * M;
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
  if (form == Simplified)
    return SP.Proj(Y, 2 * Q_product(dotY.transpose()).transpose() -
                          SP.SymBlockDiagProduct(dotY, Y, nablaF_Y));
  else {
    // Euclidean Hessian-vector product
    Matrix H_dotY = 2 * dotY * M;

    H_dotY.block(0, n, r, d * n) =
        SP.Proj(Y.block(0, n, r, d * n),
                H_dotY.block(0, n, r, d * n) -
                    SP.SymBlockDiagProduct(dotY.block(0, n, r, d * n),
                                           Y.block(0, n, r, d * n),
                                           nablaF_Y.block(0, n, r, d * n)));
    return H_dotY;
  }
}

Matrix
SESyncProblem::Riemannian_Hessian_vector_product(const Matrix &Y,
                                                 const Matrix &dotY) const {
  return Riemannian_Hessian_vector_product(Y, Euclidean_gradient(Y), dotY);
}

Matrix SESyncProblem::precondition(const Matrix &Y, const Matrix &dotY) const {
  if (preconditioner == None)
    return dotY;
  else if (preconditioner == Jacobi)
    return tangent_space_projection(Y, JacobiPrecon * dotY);
  else if (preconditioner == IncompleteCholesky)
    return tangent_space_projection(
        Y, iCholPrecon->solve(dotY.transpose()).transpose());
  else // preconditioner == RegularizedCholesky
    return tangent_space_projection(
        Y, RegCholPrecon.solve(dotY.transpose()).transpose());
}

Matrix SESyncProblem::tangent_space_projection(const Matrix &Y,
                                               const Matrix &dotY) const {
  if (form == Simplified)
    return SP.Proj(Y, dotY);
  else {
    // form == Explicit
    Matrix P(dotY.rows(), dotY.cols());

    // Projection of translational states is the identity
    P.block(0, 0, r, n) = dotY.block(0, 0, r, n);

    // Projection of generalized rotational states comes from the product of
    // Stiefel manifolds
    P.block(0, n, r, d * n) =
        SP.Proj(Y.block(0, n, r, d * n), dotY.block(0, n, r, d * n));

    return P;
  }
}

Matrix SESyncProblem::retract(const Matrix &Y, const Matrix &dotY) const {
  if (form == Simplified)
    return SP.retract(Y, dotY);
  else // form == Explicit
  {
    Matrix Yplus = Y;
    Yplus.block(0, 0, r, n) += dotY.block(0, 0, r, n);

    Yplus.block(0, n, r, d * n) =
        SP.retract(Y.block(0, n, r, d * n), dotY.block(0, n, r, d * n));

    return Yplus;
  }
}

Matrix SESyncProblem::round_solution(const Matrix Y) const {

  // First, compute a thin SVD of Y
  Eigen::JacobiSVD<Matrix> svd(Y, Eigen::ComputeThinV);

  Eigen::VectorXd sigmas = svd.singularValues();
  // Construct a diagonal matrix comprised of the first d singular values
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> Sigma_d(d);
  Eigen::DiagonalMatrix<double, Eigen::Dynamic>::DiagonalVectorType &diagonal =
      Sigma_d.diagonal();
  for (unsigned int i = 0; i < d; i++)
    diagonal(i) = sigmas(i);

  // First, construct a rank-d truncated singular value decomposition for Y
  Eigen::MatrixXd R = Sigma_d * svd.matrixV().leftCols(d).transpose();

  Eigen::VectorXd determinants(n);

  // Compute the offset at which the rotation matrix blocks begin
  unsigned int rot_offset = (form == Simplified ? 0 : n);

  unsigned int ng0 = 0; // This will count the number of blocks whose
  // determinants have positive sign
  for (unsigned int i = 0; i < n; i++) {
    // Compute the determinant of the ith dxd block of R
    determinants(i) = R.block(0, rot_offset + i * d, d, d).determinant();
    if (determinants(i) > 0)
      ng0++;
  }

  if (ng0 < n / 2) {
    // Less than half of the total number of blocks have the correct sign, so
    // reverse their orientations

    // Get a reflection matrix that we can use to reverse the signs of those
    // blocks of R that have the wrong determinant
    Eigen::MatrixXd reflector = Eigen::MatrixXd::Identity(d, d);
    reflector(d - 1, d - 1) = -1;

    R = reflector * R;
  }

// Finally, project each dxd rotation block to SO(d)
#pragma omp parallel for
  for (unsigned int i = 0; i < n; i++)
    R.block(0, rot_offset + i * d, d, d) =
        project_to_SOd(R.block(0, rot_offset + i * d, d, d));

  if (form == Explicit)
    return R;
  else // form == Explicit
  {
    // In this case, we also need to recover the corresponding translations
    Matrix X(d, (d + 1) * n);

    // Set rotational states
    X.block(0, n, d, d * n) = R;

    // Recover translational states
    X.block(0, 0, d, n) = recover_translations(B1, B2, R);

    return X;
  }
}

Matrix SESyncProblem::compute_Lambda_blocks(const Matrix &Y) const {
  // Compute S * Y', where S is the data matrix defining the quadratic form
  // for
  // the specific version of the SE-Sync problem we're solving
  Matrix SYt = data_matrix_product(Y.transpose());

  // Preallocate storage for diagonal blocks of Lambda
  Matrix Lambda_blocks(d, n * d);

  // Index of the row/column at which the rotational blocks begin in matrix X
  unsigned int offset = (form == Simplified ? 0 : n);

#pragma omp parallel for
  for (unsigned int i = 0; i < n; i++) {
    Matrix P = SYt.block(offset + i * d, 0, d, Y.rows()) *
               Y.block(0, offset + i * d, Y.rows(), d);
    Lambda_blocks.block(0, i * d, d, d) = .5 * (P + P.transpose());
  }
  return Lambda_blocks;
}

SparseMatrix SESyncProblem::compute_Lambda_from_Lambda_blocks(
    const Matrix &Lambda_blocks) const {

  unsigned int offset = (form == Simplified ? 0 : n);

  std::vector<Eigen::Triplet<SparseMatrix::Scalar>> elements;
  elements.reserve(d * d * n);

  for (unsigned int i = 0; i < n; i++)
    for (unsigned int r = 0; r < d; r++)
      for (unsigned int c = 0; c < d; c++)
        elements.emplace_back(offset + i * d + r, offset + i * d + c,
                              Lambda_blocks(r, i * d + c));

  SparseMatrix Lambda(offset + d * n, offset + d * n);
  Lambda.setFromTriplets(elements.begin(), elements.end());
  return Lambda;
}

SparseMatrix SESyncProblem::compute_Lambda(const Matrix &Y) const {
  // First, compute the diagonal blocks of Lambda
  Matrix Lambda_blocks = compute_Lambda_blocks(Y);

  return compute_Lambda_from_Lambda_blocks(Lambda_blocks);
}

bool SESyncProblem::compute_S_minus_Lambda_min_eig(
    const Matrix &Y, double &min_eigenvalue, Eigen::VectorXd &min_eigenvector,
    unsigned int max_iterations, double min_eigenvalue_nonnegativity_tolerance,
    unsigned int num_Lanczos_vectors) const {
  // First, compute the largest-magnitude eigenvalue of this matrix
  SMinusLambdaProdFunctor lm_op(this, Y);
  Spectra::SymEigsSolver<double, Spectra::SELECT_EIGENVALUE::LARGEST_MAGN,
                         SMinusLambdaProdFunctor>
      largest_magnitude_eigensolver(&lm_op, 1,
                                    std::min(num_Lanczos_vectors, n * d));
  largest_magnitude_eigensolver.init();

  int num_converged = largest_magnitude_eigensolver.compute(
      max_iterations, 1e-4, Spectra::SELECT_EIGENVALUE::LARGEST_MAGN);

  // Check convergence and bail out if necessary
  if (num_converged != 1)
    return false;

  double lambda_lm = largest_magnitude_eigensolver.eigenvalues()(0);

  if (lambda_lm < 0) {
    // The largest-magnitude eigenvalue is negative, and therefore also the
    // minimum eigenvalue, so just return this solution
    min_eigenvalue = lambda_lm;
    min_eigenvector = largest_magnitude_eigensolver.eigenvectors(1);
    min_eigenvector.normalize(); // Ensure that this is a unit vector
    return true;
  }

  // The largest-magnitude eigenvalue is positive, and is therefore the
  // maximum
  // eigenvalue.  Therefore, after shifting the spectrum of S - Lambda by -
  // 2*lambda_lm (by forming S - Lambda - 2*lambda_max*I), the  shifted
  // spectrum
  // will lie in the interval [lambda_min(A) - 2*  lambda_max(A),
  // -lambda_max*A]; in particular, the largest-magnitude eigenvalue of  S -
  // Lambda - 2*lambda_max*I is lambda_min - 2*lambda_max, with  corresponding
  // eigenvector v_min; furthermore, the condition number sigma of S - Lambda
  // -
  // 2*lambda_max is then upper-bounded by 2 :-).

  SMinusLambdaProdFunctor min_shifted_op(this, Y, -2 * lambda_lm);

  Spectra::SymEigsSolver<double, Spectra::SELECT_EIGENVALUE::LARGEST_MAGN,
                         SMinusLambdaProdFunctor>
      min_eigensolver(&min_shifted_op, 1, std::min(num_Lanczos_vectors, n * d));

  // If Y is a critical point of F, then Y^T is also in the null space
  // of S - Lambda(Y) (cf. Lemma 6 of the tech report), and therefore its
  // rows are eigenvectors corresponding to the eigenvalue 0.  In the case
  // that
  // the relaxation is exact, this is the *minimum* eigenvalue, and therefore
  // the rows of Y are exactly the eigenvectors that we're looking for.  On
  // the other hand, if the relaxation is *not* exact, then S - Lambda(Y)
  // has at least one strictly negative eigenvalue, and the rows of Y are
  // *unstable fixed points* for the Lanczos iterations.  Thus, we will take a
  // slightly "fuzzed" version of the first row of Y as an initialization
  // for the Lanczos iterations; this allows for rapid convergence in the case
  // that the relaxation is exact (since are starting close to a solution),
  // while simultaneously allowing the iterations to escape from this fixed
  // point in the case that the relaxation is not exact.
  Eigen::VectorXd v0 = Y.row(0).transpose();
  Eigen::VectorXd perturbation(v0.size());
  perturbation.setRandom();
  perturbation.normalize();
  Eigen::VectorXd xinit =
      v0 + (.03 * v0.norm()) * perturbation; // Perturb v0 by ~3%

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
  return true;
}

Matrix SESyncProblem::chordal_initialization() const {
  Matrix Y;
  if (form == Simplified) {
    Y = Matrix::Zero(r, n * d);
    Y.topRows(d) = SESync::chordal_initialization(d, B3);
  } else // form == explicit
  {
    Y = Matrix::Zero(r, n * (d + 1));

    // Compute rotations using chordal initialization
    Y.block(0, n, d, n * d) = SESync::chordal_initialization(d, B3);

    // Recover corresponding translations
    Y.block(0, 0, d, n) = recover_translations(B1, B2, Y.block(0, n, d, n * d));
  }

  return Y;
}

Matrix SESyncProblem::random_sample() const {
  Matrix Y;
  if (form == Simplified)
    // Randomly sample a point on the Stiefel manifold
    Y = SP.random_sample();
  else // form == Explicit
  {
    Y = Matrix::Zero(r, n * (d + 1));

    // Randomly sample a set of elements on the Stiefel product manifold
    Y.block(0, n, r, n * d) = SP.random_sample();

    // Randomly sample a set of coordinates for the initial positions from the
    // standard normal distribution
    std::default_random_engine generator;
    std::normal_distribution<double> g;

    for (unsigned int i = 0; i < r; i++)
      for (unsigned int j = 0; j < n; j++)
        Y(i, j) = g(generator);
  }

  return Y;
}

/// MINIMUM EIGENVALUE COMPUTATION STRUCT

SESyncProblem::SMinusLambdaProdFunctor::SMinusLambdaProdFunctor(
    const SESyncProblem *prob, const Matrix &Y, double sigma)
    : _problem(prob), _dim(prob->dimension()), _sigma(sigma) {

  if (_problem->formulation() == Simplified) {
    _rows = _problem->dimension() * _problem->num_poses();
    _cols = _problem->dimension() * _problem->num_poses();
  } else // mode == Explicit
  {
    _rows = (_problem->dimension() + 1) * _problem->num_poses();
    _cols = (_problem->dimension() + 1) * _problem->num_poses();
  }

  // Compute and cache this on construction
  _Lambda_blocks = _problem->compute_Lambda_blocks(Y);
}

void SESyncProblem::SMinusLambdaProdFunctor::perform_op(double *x,
                                                        double *y) const {
  Eigen::Map<Eigen::VectorXd> X(x, _cols);
  Eigen::Map<Eigen::VectorXd> Y(y, _rows);

  Y = _problem->data_matrix_product(X);

  // Offset corresponding to the first index in X and Y associated with
  // rotational blocks
  unsigned int offset =
      (_problem->formulation() == Simplified ? 0 : _problem->num_poses());

#pragma omp parallel for
  for (unsigned int i = 0; i < _problem->num_poses(); i++)
    Y.segment(offset + i * _dim, _dim) -=
        _Lambda_blocks.block(0, i * _dim, _dim, _dim) *
        X.segment(offset + i * _dim, _dim);

  if (_sigma != 0)
    Y += _sigma * X;
}
} // namespace SESync
