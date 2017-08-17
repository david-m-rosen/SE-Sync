/** This class encapsulates an instance of the special Euclidean synchronization
* problem
*
* dmrosen 18 May 2017
*/

#ifndef _SESYNCPROBLEM_H_
#define _SESYNCPROBLEM_H_

/** Use external matrix factorizations/linear solves provided by SuiteSparse
 * (SPQR and Cholmod) */

#include <Eigen/CholmodSupport>
#include <Eigen/Dense>
#include <Eigen/SPQRSupport>
#include <Eigen/Sparse>

#include "RelativePoseMeasurement.h"
#include "SESync_types.h"
#include "SESync_utils.h"

#include "Problem.h"
#include "ProductManifold.h"
#include "Stiefel.h"

namespace SESync {

// Forward declaration to avoid having circular includes
class SESyncRTRNewton;

/** The type of the sparse Cholesky factorization to use in the computation of
 * the orthogonal projection operation */
typedef Eigen::CholmodDecomposition<SparseMatrix> SparseCholeskyFactorization;

/** The type of the QR decomposition to use in the computation of the orthogonal
 * projection operation */
typedef Eigen::SPQR<SparseMatrix> SparseQRFactorization;

/** the type of the incomplete Cholesky decomposition we will use for
 * preconditioning the conjugate gradient iterations in the RTR method */
typedef Eigen::IncompleteCholesky<double> IncompleteCholeskyFactorization;

class SESyncProblem : public ROPTLIB::Problem {
private:
  /** Number of poses */
  unsigned int n = 0;

  /** Number of measurements */
  unsigned int m = 0;

  /** Dimensionality of the Euclidean space */
  unsigned int d = 0;

  /** The rotational connection Laplacian for the special Euclidean
*synchronization problem, cf. eq. 14 of the tech report */
  SparseMatrix LGrho;

  /** The weighted translational data matrix Omega^(1/2) T cf. e.g. eqs. 22-24
*of the tech report */
  SparseMatrix SqrtOmega_T;

  /** The transpose of the above matrix; we cache this for computational
 * efficiency, since it's used frequently */
  SparseMatrix TT_SqrtOmega;

  /** The weighted reduced oriented incidence matrix Ared Omega^(1/2), cf. e.g.
*eq. 39 of the tech report */
  SparseMatrix Ared_SqrtOmega;

  /** The transpose of the above matrix; we cache this for computational
 * efficiency, since it's used frequently */
  SparseMatrix SqrtOmega_AredT;

  /** An Eigen sparse linear solver that encodes the Cholesky factor L used
*in the computation of the orthogonal projection function, cf. eq. (39) of the
*tech report. */
  SparseCholeskyFactorization L;

  /** An Eigen sparse linear solver that encodes the QR factorization used in
 * the computation of the orthogonal projection function, cf. eq. (98) of the
 * tech report */

  // When using Eigen::SPQR, the destructor causes a segfault if this variable
  // isn't explicitly initialized (i.e. not just default-constructed)
  SparseQRFactorization *QR = nullptr;

  /** A Boolean variable determining whether to use the Cholesky or QR
 * decompositions for computing the orthogonal projection */
  bool use_Cholesky;

  /** The preconditioning strategy to use when running the Riemannian
   * trust-region algorithm */
  Preconditioner preconditioner;

  // Diagonal Jacobi preconditioner
  DiagonalMatrix JacobiPreconditioner;

  // Incomplete Cholesky Preconditioner
  IncompleteCholeskyFactorization *iChol = nullptr;

  /** The rank of the rank-restricted relaxation */
  unsigned int r = 0;

  /** The prototype Stiefel manifold for the rank-restricted relaxation */
  ROPTLIB::Stiefel *Stdr = nullptr;

  /** The product manifold that is the domain of this problem */
  ROPTLIB::ProductManifold *domain = nullptr;

  /** Preallocated working space for ProductElement <=> Matrix conversion */
  Matrix *Y = nullptr;

public:
  /** Default constructor; doesn't actually do anything */
  SESyncProblem() {}

  /** Constructor using a vector of relative pose measurements */
  SESyncProblem(
      const std::vector<SESync::RelativePoseMeasurement> &measurements) {
    set_problem_data(construct_rotational_connection_Laplacian(measurements),
                     construct_oriented_incidence_matrix(measurements),
                     construct_translational_data_matrix(measurements),
                     construct_translational_precision_matrix(measurements));
  }

  /** This function constructs the special Euclidean synchronization problem
*from the passed data matrices */
  SESyncProblem(const SparseMatrix &rotational_connection_Laplacian,
                const SparseMatrix &oriented_incidence_matrix,
                const SparseMatrix &translational_data_matrix,
                const DiagonalMatrix &translational_precisions_matrix,
                bool Cholesky = true,
                const Preconditioner &precon = IncompleteCholesky) {
    // This is just a passthrough to the initialization function
    set_problem_data(rotational_connection_Laplacian, oriented_incidence_matrix,
                     translational_data_matrix, translational_precisions_matrix,
                     Cholesky, precon);
  }

  /** This function initializes the special Euclidean synchronization problem
*using the passed data matrices */
  void set_problem_data(const SparseMatrix &rotational_connection_Laplacian,
                        const SparseMatrix &oriented_incidence_matrix,
                        const SparseMatrix &translational_data_matrix,
                        const DiagonalMatrix &translational_precisions_matrix,
                        bool Cholesky = true,
                        const Preconditioner &precon = IncompleteCholesky);

  /** Set the maximum rank of the rank-restricted semidefinite relaxation */
  void set_relaxation_rank(unsigned int rank);

  /** Given a matrix X, this function computes and returns the orthogonal
*projection Pi * X */
  Matrix Pi_product(const Matrix &X) const;

  /** This function computes and returns the product QX */
  Matrix Q_product(const Matrix &X) const;

  /** Given Y* in St(r,d)^n, this function computes and returns a d x nd matrix
*comprised of the dxd block elements of the associated block-diagonal
*Lagrange multiplier matrix Lambda(Y) */
  Matrix compute_Lambda_blocks(const Matrix &Ystar) const;

  /** Given a critical point Y* in St(r,d)^n of the optimization problem, this
*function computes the smallest eigenvalue lambda_min of Q - Lambda and its
*associated eigenvector v.  Returns a Boolean value indicating whether the
*Lanczos method used to estimate the smallest eigenpair convergence to
*within the required tolerance. */
  bool compute_Q_minus_Lambda_min_eig(
      const Matrix &Ystar, double &min_eigenvalue,
      Eigen::VectorXd &min_eigenvector, int max_iterations = 10000,
      double precision = 1e-6, unsigned int num_Lanczos_vectors = 20) const;

  /// ACCESSORS

  unsigned int num_poses() const { return n; }

  unsigned int num_measurements() const { return m; }

  unsigned int dimension() const { return d; }

  /** Get the rank of the rank-restricted semidefinite relaxation */
  unsigned int relaxation_rank() const { return r; }

  /** Get a const pointer to the Stiefel product manifold */
  const ROPTLIB::ProductManifold *get_Stiefel_product_manifold() const {
    return domain;
  }

  /// OVERRIDDEN PURE VIRTUAL BASE CLASS (ROPTLIB::PROBLEM) FUNCTIONS

  /** Evaluates the problem objective */
  double f(ROPTLIB::Variable *x) const;

  /** Evaluates the Euclidean gradient of the function */
  void EucGrad(ROPTLIB::Variable *x, ROPTLIB::Vector *g) const;

  /** Evaluates the action of the Euclidean Hessian of the function */
  void EucHessianEta(ROPTLIB::Variable *x, ROPTLIB::Vector *v,
                     ROPTLIB::Vector *Hv) const;

  ~SESyncProblem() {
    if (Stdr)
      delete Stdr;
    if (domain)
      delete domain;
    if (Y)
      delete Y;

    if (QR)
      delete QR;

    if (iChol)
      delete iChol;
  }

  /** This is a lightweight struct used in conjunction with Spectra to compute
*the minimum eigenvalue and eigenvector of Q - Lambda; it has a single
*nontrivial function, perform_op(x,y), that computes and returns the product
*y = (Q - Lambda + sigma*I) x */
  struct QMinusLambdaProdFunctor {
    const SESyncProblem *_problem;
    Matrix _Lambda_blocks;
    int _rows;
    int _cols;
    int _dim;
    double _sigma;

    QMinusLambdaProdFunctor(const SESyncProblem *prob, const Matrix &Ystar,
                            double sigma = 0)
        : _problem(prob), _rows(prob->dimension() * prob->num_poses()),
          _cols(prob->dimension() * prob->num_poses()), _dim(prob->dimension()),
          _sigma(sigma) {
      // Compute and cache this on construction
      _Lambda_blocks = prob->compute_Lambda_blocks(Ystar);
    }
    int rows() const { return _rows; }
    int cols() const { return _cols; }

    void perform_op(double *x, double *y) const {
      Eigen::Map<Eigen::VectorXd> X(x, _cols);
      Eigen::Map<Eigen::VectorXd> Y(y, _rows);

      Y = _problem->Q_product(X);

      for (unsigned int i = 0; i < _problem->num_poses(); i++)
        Y.segment(i * _dim, _dim) -=
            _Lambda_blocks.block(0, i * _dim, _dim, _dim) *
            X.segment(i * _dim, _dim);

      if (_sigma != 0)
        Y += _sigma * X;
    }
  };

  // Declare the specialization of the RTR Newton method to be a friend so we
  // can directly access private/protected datas
  friend SESyncRTRNewton;
};
}
#endif // _SESYNCPROBLEM_H_
