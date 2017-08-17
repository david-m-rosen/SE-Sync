#include <SymEigsSolver.h> // Spectra's symmetric eigensolver

#include "SESyncProblem.h"
#include "SESync_utils.h"

namespace SESync {

void SESyncProblem::set_problem_data(
    const SparseMatrix &rotational_connection_Laplacian,
    const SparseMatrix &oriented_incidence_matrix,
    const SparseMatrix &translational_data_matrix,
    const DiagonalMatrix &translational_precisions_matrix, bool Cholesky,
    const Preconditioner &precon) {
  // Initialize subclass data members

  // Dimensions of the problem
  n = oriented_incidence_matrix.rows();
  m = oriented_incidence_matrix.cols();
  d = rotational_connection_Laplacian.rows() / n;

  // Matrices

  LGrho = rotational_connection_Laplacian;

  DiagonalMatrix SqrtOmega(
      translational_precisions_matrix.diagonal().cwiseSqrt());

  Ared_SqrtOmega = oriented_incidence_matrix.topRows(n - 1) * SqrtOmega;
  SqrtOmega_AredT = Ared_SqrtOmega.transpose();

  SqrtOmega_T = SqrtOmega * translational_data_matrix;
  TT_SqrtOmega = SqrtOmega_T.transpose();

  use_Cholesky = Cholesky;

  preconditioner = precon;

  if (use_Cholesky) {
    // Compute and cache the Cholesky factor L of Ared * Omega * Ared^T
    L.compute(Ared_SqrtOmega * SqrtOmega_AredT);
  } else {
    // Compute the QR decomposition of Omega^(1/2) * Ared^T (cf. eq. (98) of the
    // tech report)

    // Note that Eigen's sparse QR factorization can only be called on matrices
    // stored in compressed format
    SqrtOmega_AredT.makeCompressed();

    QR = new SparseQRFactorization();
    QR->compute(SqrtOmega_AredT);
  }

  /** Compute and cache preconditioning matrices, if required */
  if (preconditioner == Jacobi) {
    Eigen::VectorXd diag = LGrho.diagonal();
    JacobiPreconditioner = diag.cwiseInverse().asDiagonal();
  } else if (preconditioner == Preconditioner::IncompleteCholesky)
    iChol = new IncompleteCholeskyFactorization(LGrho);

  // General configuration
  SetUseGrad(true);
  SetUseHess(true);
}

void SESyncProblem::set_relaxation_rank(unsigned int rank) {
  // Record this value
  r = rank;

  // Construct the prototype Stiefel manifold
  if (Stdr)
    delete Stdr;
  Stdr = new ROPTLIB::Stiefel(r, d);

  // Use the Euclidean metric, QF retraction, vector transport by projection and
  // extrinsic representation
  Stdr->ChooseStieParamsSet3();
  // Stdr->ChooseStieParamsSet4();

  // Now construct the product manifold for this problem
  if (domain)
    delete domain;
  domain = new ROPTLIB::ProductManifold(1, Stdr, n);
  domain->SetIsIntrApproach(false);

  // Call the super function that sets the domain for this problem
  SetDomain(domain);

  // Allocate working space
  if (Y)
    delete Y;
  Y = new Matrix(r, n * d);
}

Matrix SESyncProblem::Pi_product(const Matrix &X) const {
  if (use_Cholesky)
    return X - SqrtOmega_AredT * L.solve(Ared_SqrtOmega * X);
  else {
    Eigen::MatrixXd PiX = X;
    for (unsigned int c = 0; c < X.cols(); c++) {
      // Eigen's SPQR support only supports solving with vectors(!) (i.e.
      // 1-column matrices)
      PiX.col(c) = X.col(c) - SqrtOmega_AredT * QR->solve(X.col(c));
    }
    return PiX;
  }
}

Matrix SESyncProblem::Q_product(const Matrix &X) const {
  return LGrho * X + TT_SqrtOmega * Pi_product(SqrtOmega_T * X);
}

Matrix SESyncProblem::compute_Lambda_blocks(const Matrix &Ystar) const {
  Eigen::MatrixXd QYstarT = Q_product(Ystar.transpose());
  Eigen::MatrixXd Lambda_blocks(d, n * d);

  // Preallocation of working space for computing the block elements of Lambda
  Eigen::MatrixXd B(d, d);

  for (unsigned int i = 0; i < n; i++) {
    B = QYstarT.block(i * d, 0, d, Ystar.rows()) *
        Ystar.block(0, i * d, Ystar.rows(), d);
    Lambda_blocks.block(0, i * d, d, d) = .5 * (B + B.transpose());
  }
  return Lambda_blocks;
}

bool SESyncProblem::compute_Q_minus_Lambda_min_eig(
    const Matrix &Ystar, double &min_eigenvalue,
    Eigen::VectorXd &min_eigenvector, int max_iterations, double precision,
    unsigned int num_Lanczos_vectors) const {
  // First, compute the largest-magnitude eigenvalue of this matrix
  QMinusLambdaProdFunctor lm_op(this, Ystar);
  Spectra::SymEigsSolver<double, Spectra::SELECT_EIGENVALUE::LARGEST_MAGN,
                         QMinusLambdaProdFunctor>
      largest_magnitude_eigensolver(&lm_op, 1,
                                    std::min(num_Lanczos_vectors, n * d));
  largest_magnitude_eigensolver.init();

  int num_converged = largest_magnitude_eigensolver.compute(
      max_iterations, precision, Spectra::SELECT_EIGENVALUE::LARGEST_MAGN);

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

  // The largest-magnitude eigenvalue is positive, and is therefore the maximum
  // eigenvalue.  Therefore, after shifting the spectrum of Q - Lambda by -
  // 2*lambda_lm (by forming Q - Lambda - 2*lambda_max*I), the  shifted spectrum
  // will line in the interval [lambda_min(A) - 2*  lambda_max(A),
  // -lambda_max*A]; in particular, the largest-magnitude eigenvalue of  Q -
  // Lambda - 2*lambda_max*I is lambda_min - 2*lambda_max, with  corresponding
  // eigenvector v_min; furthermore, the condition number sigma of Q - Lambda -
  // 2*lambda_max is then upper-bounded by 2 :-).

  QMinusLambdaProdFunctor min_shifted_op(this, Ystar, -2 * lambda_lm);

  Spectra::SymEigsSolver<double, Spectra::SELECT_EIGENVALUE::LARGEST_MAGN,
                         QMinusLambdaProdFunctor>
      min_eigensolver(&min_shifted_op, 1, std::min(num_Lanczos_vectors, n * d));

  // If Ystar is a critical point of F, then Ystar^T is also in the null space
  // of Q - Lambda(Ystar) (cf. Lemma 6 of the tech report), and therefore its
  // rows are eigenvectors corresponding to the eigenvalue 0.  In the case that
  // the relaxation is exact, this is the *minimum* eigenvalue, and therefore
  // the rows of Ystar are exactly the eigenvectors that we're looking for.  On
  // the other hand, if the relaxation is *not* exact, then Q - Lambda(Ystar)
  // has at least one strictly negative eigenvalue, and the rows of Ystar are
  // *unstable fixed points* for the Lanczos iterations.  Thus, we will take a
  // slightly "fuzzed" version of the first row of Ystar as an initialization
  // for the Lanczos iterations; this allows for rapid convergence in the case
  // that the relaxation is exact (since are starting close to a solution),
  // while simultaneously allowing the iterations to escape from this fixed
  // point in the case that the relaxation is not exact.
  Eigen::VectorXd v0 = Ystar.row(0).transpose();
  Eigen::VectorXd perturbation(v0.size());
  perturbation.setRandom();
  perturbation.normalize();
  Eigen::VectorXd xinit =
      v0 + (.03 * v0.norm()) * perturbation; // Perturb v0 by ~3%

  // Use this to initialize the eigensolver
  min_eigensolver.init(xinit.data());
  num_converged = min_eigensolver.compute(
      max_iterations, precision, Spectra::SELECT_EIGENVALUE::LARGEST_MAGN);

  if (num_converged != 1)
    return false;

  min_eigenvector = min_eigensolver.eigenvectors(1);
  min_eigenvector.normalize(); // Ensure that this is a unit vector
  min_eigenvalue = min_eigensolver.eigenvalues()(0) + 2 * lambda_lm;
  return true;
}

double SESyncProblem::f(ROPTLIB::Variable *x) const {
  ROPTLIB::ProductElement *X = static_cast<ROPTLIB::ProductElement *>(x);
  StiefelProd2Mat(*X, *Y);
  return ((*Y) * Q_product(Y->transpose())).trace();
}

void SESyncProblem::EucGrad(ROPTLIB::Variable *x, ROPTLIB::Vector *g) const {
  ROPTLIB::ProductElement *X = static_cast<ROPTLIB::ProductElement *>(x);
  ROPTLIB::ProductElement *G = static_cast<ROPTLIB::ProductElement *>(g);
  StiefelProd2Mat(*X, *Y);
  Mat2StiefelProd(2 * Q_product(Y->transpose()).transpose(), *G);
}

void SESyncProblem::EucHessianEta(ROPTLIB::Variable *x, ROPTLIB::Vector *v,
                                  ROPTLIB::Vector *Hv) const {
  ROPTLIB::ProductElement *X = static_cast<ROPTLIB::ProductElement *>(x);
  ROPTLIB::ProductElement *V = static_cast<ROPTLIB::ProductElement *>(v);
  ROPTLIB::ProductElement *HV = static_cast<ROPTLIB::ProductElement *>(Hv);
  StiefelProd2Mat(*V, *Y);
  Mat2StiefelProd(2 * Q_product(Y->transpose()).transpose(), *HV);
}
}
