/** This file provides a convenient functional interface to the SESync algorithm
*
* dmrosen 22 May 2017
*/

#ifndef _SESYNC_H_
#define _SESYNC_H_

#include <vector>

#include <Eigen/Dense>

#include "SESync_types.h"

#include "RelativePoseMeasurement.h"

namespace SESync {

/** This struct contains the various parameters that control the SESync
 * algorithm */
struct SESyncOpts {
  /// OPTIMIZATION STOPPING CRITERIA

  /** Stopping tolerance for the norm of the Riemannian gradient */
  double grad_norm_tol;

  /** Stopping criterion based upon the relative decrease in function value */
  double rel_func_decrease_tol;

  /** Maximum permitted number of (outer) iterations of the RTR algorithm */
  unsigned int max_RTR_iterations;

  /** Maximum number of inner (truncated conjugate-gradient) iterations to
 * perform per out iteration */
  unsigned int max_tCG_iterations;

  /// SE-SYNC PARAMETERS

  /** The initial level of the Riemannian Staircase */
  unsigned int r0;

  /** The maximum level of the Riemannian Staircase to explore */
  unsigned int rmax;

  /** Relative tolerance to require for the computation of the minimum
 * eigenvalue of Q - Lambda(Y*) using Lanczos */
  double eig_comp_tol;

  /** The maximum number of Lanczos iterations to admit for eigenvalue
 * computations */
  unsigned int max_eig_iterations;

  /** A numerical tolerance for acceptance of the minimum eigenvalue of Q -
 * Lambda(Y*)
 * as numerically nonnegative; this should be a small magnitude value e.g. 10^-4
 */
  double min_eig_num_tol;

  /** The number of working vectors to use in the minimum eigenvalue computation
(using the implicitly-restarted Arnoldi algorithm); must be in the range [1,
(#poses) * (#problem dimension) - 1] */
  unsigned int num_Lanczos_vectors;

  /** If no initial iterate Y0 is supplied, this boolean determines the
 * initialization strategy employed by SE-Sync: 'true' -> chordal, 'false' ->
 * random sampling */
  bool use_chordal_initialization;

  /** Whether to use the Cholesky or QR factorization when computing the
   * orthogonal projection */
  bool use_Cholesky;

  /** The preconditioning strategy to use in the Riemannian trust-region
   * algorithm*/
  Preconditioner precon;

  /** Whether to print output as the algorithm runs */
  bool verbose;

  /** Default constructor; all arguments are optional, with 'reasonable' default
 * values */
  SESyncOpts(double gradient_norm_tolerance = 1e-2,
             double relative_function_decrease_tolerance = 1e-6,
             unsigned int max_RTR_iters = 100,
             unsigned int max_tCG_iters = 2000, unsigned int init_RS_level = 5,
             unsigned int max_RS_level = 7,
             double eigenvalue_computation_tolerance = 1e-10,
             unsigned int max_eigenvalue_iterations = 10000,
             double eigenvalue_nonnegativity_numerical_tolerance = 1e-5,
             unsigned int num_vectors = 20, bool chordal_initialization = true,
             bool Cholesky = true,
             const Preconditioner &preconditioner = IncompleteCholesky,
             bool verbose_output = false)
      : grad_norm_tol(gradient_norm_tolerance),
        rel_func_decrease_tol(relative_function_decrease_tolerance),
        max_RTR_iterations(max_RTR_iters), max_tCG_iterations(max_tCG_iters),
        r0(init_RS_level), rmax(max_RS_level),
        eig_comp_tol(eigenvalue_computation_tolerance),
        max_eig_iterations(max_eigenvalue_iterations),
        min_eig_num_tol(eigenvalue_nonnegativity_numerical_tolerance),
        num_Lanczos_vectors(num_vectors),
        use_chordal_initialization(chordal_initialization),
        use_Cholesky(Cholesky), precon(preconditioner),
        verbose(verbose_output) {}
};

/** These enumerations describe the termination status of the SE-Sync algorithm
 */
enum SESyncStatus {
  /** The algorithm converged to a certified global optimum */
  GLOBAL_OPT,

  /** The algorithm converged to a saddle point, but the backtracking line
     search was unable to escape it */
  SADDLE_POINT,

  /** The algorithm converged to a first-order critical point, but the
     minimum-eigenvalue computation did not converge to sufficient precision to
     enable its characterization */
  EIG_IMPRECISION,

  /** The algorithm exhausted the maximum number of iterations of the Riemannian
     Staircase before finding an optimal solution */
  RS_ITER_LIMIT
};

/** This struct contains the output of the SESync algorithm */
struct SESyncResult {
  /** The optimal value of the SDP relaxation */
  double SDPval;

  /** A minimizer Yopt \in St(r,d)^n of the rank-restricted relaxation */
  Matrix Yopt;

  /** The norm of the Riemannian gradient at Yopt */
  double gradnorm;

  /** The minimum eigenvalue of the matrix Q - Lambda(Yopt) */
  double lambda_min;

  /** The corresponding eigenvector of the minimum eigenvalue */
  Vector v_min;

  /** The optimal value of the reprojected solution */
  double Fxhat;

  /** The reprojected rotational state estimates Rhat \in SO(d)^n */
  Matrix Rhat;

  /** The elapsed computation time used to compute the initialization for the
   * Riemannian Staircase */
  double initialization_time;

  /** The number of iterations performed by the Riemannian Trust Region
   * optimization method at each level of the Riemannian Staircase */
  std::vector<unsigned int> RTR_iterations;

  /** The sequence of function values obtained by the Riemannian Staircase as
   * the algorithm runs */
  std::vector<double> function_values;

  /** The sequence of norms of the gradients obtained by the Riemannian
   * Staircase as the algorithm runs */
  std::vector<double> gradient_norm_values;

  /** The number of Hessian-vector multiplications performed by the Riemannian
   * Trust Region algorithm at each level of the Riemannian Staircase */
  std::vector<unsigned int> Hessian_multiplications;

  /** The elapsed computation time since the start of the Riemannian Staircase
   * algorithm at which the returned function values were obtained */
  std::vector<double> elapsed_optimization_times;

  /** The optimal translational estimates corresponding to the rotational
   * estimate Rhat */
  Matrix that;

  SESyncStatus status;
};

/** Given a vector of relative pose measurements specifying a special Euclidean
 * synchronization problem, performs synchronization using the SESync algorithm
 */
SESyncResult SESync(const std::vector<RelativePoseMeasurement> &measurements,
                    const SESyncOpts &options = SESyncOpts(),
                    const Matrix &Y0 = Matrix());
} // namespace SESync

#endif // _SESYNC_H_
