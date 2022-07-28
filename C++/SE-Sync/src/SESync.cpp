#include <functional>

#include "SESync/SESync.h"
#include "SESync/SESyncProblem.h"
#include "SESync/SESync_types.h"
#include "SESync/SESync_utils.h"

#include "Optimization/Riemannian/TNT.h"

#include <algorithm>

namespace SESync {

SESyncResult SESync(SESyncProblem &problem, const SESyncOpts &options,
                    const Matrix &Y0) {

  /// INPUT SANITATION

  if (options.r0 < problem.dimension())
    throw std::invalid_argument("Initial relaxation rank must be at least "
                                "the dimension of the estimation problem.");

  if (options.rmax < options.r0)
    throw std::invalid_argument("Maximum relaxation rank must be greater than "
                                "or equal to initial relaxation rank.");

  if (options.max_computation_time <= 0)
    throw std::invalid_argument(
        "Maximum computation time must be a positive value");

  if (options.min_eig_num_tol <= 0)
    throw std::invalid_argument("Numerical tolerance for minimum eigenvalue "
                                "nonnegativity must be a positive value");

  if (options.LOBPCG_block_size < 1)
    throw std::invalid_argument("LOBPCG block size must be a positive integer");

  if (options.LOBPCG_max_fill_factor <= 0)
    throw std::invalid_argument("Maximum fill factor for LOBPCG preconditioner "
                                "must be a positive value");

  if (options.LOBPCG_drop_tol <= 0 || options.LOBPCG_drop_tol > 1)
    throw std::invalid_argument("Drop tolerance for LOBPCG preconditioner must "
                                "be a positive value in the range (0, 1]");

  if (options.LOBPCG_max_iterations <= 0)
    throw std::invalid_argument(
        "Maximum number of LOBPCG iterations must be a positive value");

  /// ALGORITHM DATA

  // The current iterate in the Riemannian Staircase
  Matrix Y;

  // A cache variable to store the *Euclidean* gradient at the current iterate Y
  Matrix NablaF_Y;

  // The output results struct that we will return
  SESyncResult sesync_result;
  sesync_result.status = MaxRank;

  /// OPTION PARSING AND OUTPUT TO USER

  if (options.verbose) {
    std::cout << "========= SE-Sync ==========" << std::endl << std::endl;

    std::cout << "ALGORITHM SETTINGS:" << std::endl << std::endl;
    std::cout << "SE-Sync settings:" << std::endl;
    std::cout << " SE-Sync problem formulation: ";
    if (problem.formulation() == Formulation::Simplified)
      std::cout << "Simplified";
    else if (problem.formulation() == Formulation::Explicit)
      std::cout << "Explicit";
    else // formulation == SOSync
      std::cout << "SO-Sync";
    std::cout << std::endl;
    std::cout << " Initial level of Riemannian staircase: " << options.r0
              << std::endl;
    std::cout << " Maximum level of Riemannian staircase: " << options.rmax
              << std::endl;
    std::cout << " Tolerance for accepting an eigenvalue as numerically "
                 "nonnegative in optimality verification: "
              << options.min_eig_num_tol << std::endl;
    std::cout << " LOBPCG block size: " << options.LOBPCG_block_size
              << std::endl;
    std::cout << " LOBPCG preconditioner maximum fill factor: "
              << options.LOBPCG_max_fill_factor << std::endl;
    std::cout << " LOBPCG preconditioner drop tolerance: "
              << options.LOBPCG_drop_tol << std::endl;

    std::cout << " Maximum number of LOBPCG iterations for escape direction "
                 "computation: "
              << options.LOBPCG_max_iterations << std::endl;

    if (problem.formulation() == Formulation::Simplified) {
      std::cout << " Using "
                << (problem.projection_factorization() ==
                            ProjectionFactorization::Cholesky
                        ? "Cholesky"
                        : "QR")
                << " decomposition to compute orthogonal projections"
                << std::endl;
    }
    std::cout << " Initialization method: "
              << (options.initialization == Initialization::Chordal ? "chordal"
                                                                    : "random")
              << std::endl;
    if (options.log_iterates)
      std::cout << " Logging entire sequence of Riemannian Staircase iterates"
                << std::endl;
#if defined(_OPENMP)
    std::cout << " Running SE-Sync with " << options.num_threads << " threads"
              << std::endl;
#endif
    std::cout << std::endl;

    std::cout << "Riemannian trust-region settings:" << std::endl;
    std::cout << " Stopping tolerance for norm of Riemannian gradient: "
              << options.grad_norm_tol << std::endl;
    std::cout << " Stopping tolerance for norm of preconditioned Riemannian "
                 "gradient: "
              << options.preconditioned_grad_norm_tol << std::endl;
    std::cout << " Stopping tolerance for relative function decrease: "
              << options.rel_func_decrease_tol << std::endl;
    std::cout << " Stopping tolerance for the norm of an accepted update step: "
              << options.stepsize_tol << std::endl;
    std::cout << " Maximum number of trust-region iterations: "
              << options.max_iterations << std::endl;
    std::cout << " Maximum number of truncated conjugate gradient iterations "
                 "per outer iteration: "
              << options.max_tCG_iterations << std::endl;
    std::cout << " STPCG fractional gradient tolerance (kappa): "
              << options.STPCG_kappa << std::endl;
    std::cout << " STPCG target q-superlinear convergence rate (1 + theta): "
              << (1 + options.STPCG_theta) << std::endl;
    std::cout
        << " Preconditioning the truncated conjugate gradient method using ";
    if (problem.preconditioner() == Preconditioner::None)
      std::cout << "the identity preconditioner";
    else if (problem.preconditioner() == Preconditioner::Jacobi)
      std::cout << "Jacobi preconditioner";
    else if (problem.preconditioner() == Preconditioner::RegularizedCholesky)
      std::cout << "regularized Cholesky preconditioner with maximum condition "
                   "number "
                << problem.regularized_Cholesky_preconditioner_max_condition();

    std::cout << std::endl << std::endl;
  } // if (options.verbose)

  /// ALGORITHM START
  auto SESync_start_time = Stopwatch::tick();

// Set number of threads
#if defined(_OPENMP)
  omp_set_num_threads(options.num_threads);
#endif

  /// SET UP OPTIMIZATION

  /// Function handles required by the TNT optimization algorithm

  // Objective
  Optimization::Objective<Matrix, Scalar, Matrix> F =
      [&problem](const Matrix &Y, const Matrix &NablaF_Y) {
        return problem.evaluate_objective(Y);
      };

  // Local quadratic model constructor
  Optimization::Riemannian::QuadraticModel<Matrix, Matrix, Matrix> QM =
      [&problem](const Matrix &Y, Matrix &grad,
                 Optimization::Riemannian::LinearOperator<Matrix, Matrix,
                                                          Matrix> &HessOp,
                 Matrix &NablaF_Y) {
        // Compute and cache Euclidean gradient at the current iterate
        NablaF_Y = problem.Euclidean_gradient(Y);

        // Compute Riemannian gradient from Euclidean gradient
        grad = problem.Riemannian_gradient(Y, NablaF_Y);

        // Define linear operator for computing Riemannian Hessian-vector
        // products (cf. eq. (44) in the SE-Sync tech report)
        HessOp = [&problem](const Matrix &Y, const Matrix &Ydot,
                            const Matrix &NablaF_Y) {
          return problem.Riemannian_Hessian_vector_product(Y, NablaF_Y, Ydot);
        };
      };

  // Riemannian metric

  // We consider a realization of the product of Stiefel manifolds as an
  // embedded submanifold of R^{r x dn}; consequently, the induced Riemannian
  // metric is simply the usual Euclidean inner product
  Optimization::Riemannian::RiemannianMetric<Matrix, Matrix, Scalar, Matrix>
      metric = [&problem](const Matrix &Y, const Matrix &V1, const Matrix &V2,
                          const Matrix &NablaF_Y) {
        return (V1 * V2.transpose()).trace();
      };

  // Retraction operator
  Optimization::Riemannian::Retraction<Matrix, Matrix, Matrix> retraction =
      [&problem](const Matrix &Y, const Matrix &Ydot, const Matrix &NablaF_Y) {
        return problem.retract(Y, Ydot);
      };

  // Preconditioning operator (optional)
  std::optional<
      Optimization::Riemannian::LinearOperator<Matrix, Matrix, Matrix>>
      precon;
  if (options.preconditioner == Preconditioner::None)
    precon = std::nullopt;
  else {
    Optimization::Riemannian::LinearOperator<Matrix, Matrix, Matrix> precon_op =
        [&problem](const Matrix &Y, const Matrix &Ydot,
                   const Matrix &NablaF_Y) {
          return problem.precondition(Y, Ydot);
        };
    precon = precon_op;
  }

  /// INITIALIZATION
  if (options.verbose)
    std::cout << "INITIALIZATION:" << std::endl;

  problem.set_relaxation_rank(options.r0);

  if (Y0.size() != 0) {
    if (options.verbose)
      std::cout << " Using user-supplied initial iterate Y0" << std::endl;

    Y = Y0;
  } else {
    if (options.initialization == Initialization::Chordal) {
      if (options.verbose)
        std::cout << " Computing chordal initialization ... ";

      auto chordal_init_start_time = Stopwatch::tick();
      Y = problem.chordal_initialization();
      double chordal_init_elapsed_time =
          Stopwatch::tock(chordal_init_start_time);
      if (options.verbose)
        std::cout << "elapsed computation time: " << chordal_init_elapsed_time
                  << " seconds" << std::endl;

    } else {
      if (options.verbose)
        std::cout << " Sampling a random initialization ... " << std::endl;
      Y = problem.random_sample();
    }
  }

  sesync_result.initialization_time = Stopwatch::tock(SESync_start_time);
  if (options.verbose)
    std::cout << " SE-Sync initialization finished; elapsed time: "
              << sesync_result.initialization_time << " seconds" << std::endl
              << std::endl;

  if (options.verbose) {
    // Compute and display the initial objective value
    std::cout << "Initial objective value: " << problem.evaluate_objective(Y)
              << std::endl;
  }

  /// RIEMANNIAN STAIRCASE

  // Configure optimization parameters
  Optimization::Riemannian::TNTParams<Scalar> params;
  params.gradient_tolerance = options.grad_norm_tol;
  params.preconditioned_gradient_tolerance =
      options.preconditioned_grad_norm_tol;
  params.relative_decrease_tolerance = options.rel_func_decrease_tol;
  params.stepsize_tolerance = options.stepsize_tol;
  params.max_iterations = options.max_iterations;
  params.max_TPCG_iterations = options.max_tCG_iterations;
  params.kappa_fgr = options.STPCG_kappa;
  params.theta = options.STPCG_theta;
  params.log_iterates = options.log_iterates;
  params.verbose = options.verbose;

  auto riemannian_staircase_start_time = Stopwatch::tick();

  for (size_t r = options.r0; r <= options.rmax; r++) {
    // The elapsed time from the start of the Riemannian Staircase algorithm
    // until the start of this iteration of RTR
    double RTR_iteration_start_time =
        Stopwatch::tock(riemannian_staircase_start_time);

    /// Test temporal stopping condition

    if (RTR_iteration_start_time >= options.max_computation_time) {
      sesync_result.status = ElapsedTime;
      break;
    }

    // Set  maximum permitted computation time for this level of the
    // Riemannian Staircase
    params.max_computation_time =
        options.max_computation_time - RTR_iteration_start_time;

    if (options.verbose)
      std::cout << std::endl
                << std::endl
                << "====== RIEMANNIAN STAIRCASE (level r = " << r
                << ") ======" << std::endl
                << std::endl;

    /// Run optimization!
    Optimization::Riemannian::TNTResult<Matrix, Scalar> tnt_result =
        Optimization::Riemannian::TNT<Matrix, Matrix, Scalar, Matrix>(
            F, QM, metric, retraction, Y, NablaF_Y, precon, params,
            options.user_function);

    // Extract the results
    sesync_result.Yopt = tnt_result.x;
    sesync_result.SDPval = tnt_result.f;
    sesync_result.gradnorm =
        problem.Riemannian_gradient(sesync_result.Yopt).norm();

    // Record sequence of function values
    sesync_result.function_values.push_back(tnt_result.objective_values);

    // Record sequence of gradient norms
    sesync_result.gradient_norms.push_back(tnt_result.gradient_norms);

    // Record sequence of preconditioned gradient norms
    sesync_result.preconditioned_gradient_norms.push_back(
        tnt_result.preconditioned_gradient_norms);

    // Record sequence of (# Hessian-vector products)
    sesync_result.Hessian_vector_products.push_back(
        tnt_result.inner_iterations);

    // Record sequence of update step norms
    sesync_result.update_step_norms.push_back(tnt_result.update_step_norms);

    // Record sequence of update step M-norms
    sesync_result.update_step_M_norms.push_back(tnt_result.update_step_M_norms);

    // Record sequence of gain ratios for the update steps
    sesync_result.gain_ratios.push_back(tnt_result.gain_ratios);

    // Record sequence of elapsed optimization times
    sesync_result.elapsed_optimization_times.push_back(tnt_result.time);

    // Record sequence of pose estimates, if requested
    if (options.log_iterates)
      sesync_result.iterates.push_back(tnt_result.iterates);

    /// Check TNT termination status
    if (tnt_result.status == Optimization::Riemannian::TNTStatus::ElapsedTime) {
      sesync_result.status = SESyncStatus::ElapsedTime;
      break;
    }

    if (options.verbose) {
      // Display some output to the user
      std::cout << std::endl
                << "Found first-order critical point with value F(Y) = "
                << sesync_result.SDPval
                << "!  Elapsed computation time: " << tnt_result.elapsed_time
                << " seconds" << std::endl
                << std::endl;
      std::cout << "Checking second order optimality ... " << std::endl;
    }

    /// Check second-order optimality

    size_t num_lobpcg_iters;
    auto verification_start_time = Stopwatch::tick();

    Vector v;     // Escape direction
    Scalar theta; // Curvature of certificate matrix along escape direction

    bool global_opt = problem.verify_solution(
        sesync_result.Yopt, options.min_eig_num_tol, options.LOBPCG_block_size,
        theta, v, num_lobpcg_iters, options.LOBPCG_max_iterations,
        options.LOBPCG_max_fill_factor, options.LOBPCG_drop_tol);
    double verification_elapsed_time = Stopwatch::tock(verification_start_time);

    // Check eigenvalue convergence
    if (!global_opt && theta >= -options.min_eig_num_tol / 2) {
      if (options.verbose)
        std::cout
            << "WARNING! ESCAPE DIRECTION COMPUTATION DID NOT CONVERGE TO "
               "DESIRED PRECISION!"
            << std::endl;
      sesync_result.status = EigImprecision;
      break;
    }

    // Record results of eigenvalue computation
    sesync_result.escape_direction_curvatures.push_back(theta);
    sesync_result.LOBPCG_iters.push_back(num_lobpcg_iters);
    sesync_result.verification_times.push_back(verification_elapsed_time);

    if (global_opt) {
      // results.Yopt is a second-order critical point (global optimum)!
      if (options.verbose)
        std::cout
            << "Found second-order critical point! Elapsed computation time: "
            << verification_elapsed_time << " seconds." << std::endl;
      sesync_result.status = GlobalOpt;
      break;
    } // global optimality
    else {

      /// ESCAPE FROM SADDLE!
      if (options.verbose) {
        std::cout << "Saddle point detected! Curvature along escape direction: "
                  << theta << ".  Elapsed computation time: "
                  << verification_elapsed_time << " seconds ("
                  << num_lobpcg_iters << " LOBPCG iterations)." << std::endl;
      }

      // Augment the rank of the rank-restricted semidefinite relaxation in
      // preparation for ascending to the next level of the Riemannian
      // Staircase
      problem.set_relaxation_rank(r + 1);

      Matrix Yplus;
      bool escape_success = escape_saddle(
          problem, sesync_result.Yopt, theta, v, options.grad_norm_tol,
          options.preconditioned_grad_norm_tol, Yplus);

      if (escape_success) {
        // Update initialization point for next level in the Staircase
        Y = Yplus;
      } else {
        if (options.verbose)
          std::cout
              << "WARNING!  BACKTRACKING LINE SEARCH FAILED TO ESCAPE FROM "
                 "SADDLE POINT!  (Try decreasing the preconditioned "
                 "gradient norm tolerance)"
              << std::endl;
        sesync_result.status = SaddlePoint;
        break;
      }
    } // saddle point
  }   // Riemannian Staircase

  /// POST-PROCESSING

  if (options.verbose) {
    std::cout << std::endl
              << std::endl
              << "===== END RIEMANNIAN STAIRCASE =====" << std::endl
              << std::endl;

    switch (sesync_result.status) {
    case GlobalOpt:
      std::cout << "Found global optimum!" << std::endl;
      break;
    case EigImprecision:
      std::cout << "WARNING: Escape direction computation did not achieve "
                   "sufficient accuracy; solution may not be globally optimal!"
                << std::endl;
      break;
    case SaddlePoint:
      std::cout << "WARNING: Line search was unable to escape saddle point!  "
                   "Solution is not globally optimal!"
                << std::endl;
      break;
    case MaxRank:
      std::cout << "WARNING: Riemannian Staircase reached the maximum "
                   "permitted level before finding global optimum!"
                << std::endl;
      break;
    case ElapsedTime:
      std::cout << "WARNING: Algorithm exhausted the allotted computation "
                   "time before finding global optimum!"
                << std::endl;
      break;
    }
  } // if (options.verbose)

  if (options.verbose)
    std::cout << std::endl << "Rounding solution ... ";

  // Round solution
  auto rounding_start_time = Stopwatch::tick();
  // Recover the complete pose matrix X = [t | R]
  sesync_result.xhat = problem.round_solution(sesync_result.Yopt);
  double rounding_elapsed_time = Stopwatch::tock(rounding_start_time);

  if (options.verbose)
    std::cout << "elapsed computation time: " << rounding_elapsed_time
              << " seconds" << std::endl
              << std::endl;

  sesync_result.total_computation_time = Stopwatch::tock(SESync_start_time);

  /// Compute some additional interesting bits of data

  // Evaluate objective function at ROUNDED solution.
  // Note that since xhat contains the *complete* set of pose estimates, we must
  // extract only the *rotational* elements of xhat if the SE synchronization
  // problem was solved using the simplified formulation
  sesync_result.Fxhat =
      (problem.formulation() == Formulation::Simplified
           ? problem.evaluate_objective(sesync_result.xhat.block(
                 0, problem.num_states(), problem.dimension(),
                 problem.dimension() * problem.num_states()))
           : problem.evaluate_objective(sesync_result.xhat));

  // Compute the primal optimal SDP solution Lambda and its objective value
  Matrix Lambda_blocks = problem.compute_Lambda_blocks(sesync_result.Yopt);

  sesync_result.trLambda = 0;
  for (size_t i = 0; i < problem.num_states(); i++)
    sesync_result.trLambda +=
        Lambda_blocks
            .block(0, i * problem.dimension(), problem.dimension(),
                   problem.dimension())
            .trace();

  sesync_result.Lambda =
      problem.compute_Lambda_from_Lambda_blocks(Lambda_blocks);

  // Get the duality gap for the primal-dual pair (Y'*Y, Lambda) of SDP
  // estimates

  sesync_result.duality_gap = sesync_result.SDPval - sesync_result.trLambda;

  // Get an upper bound on the (global) suboptimality of the recovered (rounded)
  // pose estimates
  sesync_result.suboptimality_bound =
      sesync_result.Fxhat - sesync_result.trLambda;

  /// FINAL OUTPUT

  if (options.verbose) {
    std::cout << "SDP RESULTS:" << std::endl;
    std::cout << "Value of dual SDP solution F(Y): " << sesync_result.SDPval
              << std::endl;
    std::cout << "Norm of Riemannian gradient grad F(Y): "
              << sesync_result.gradnorm << std::endl;
    std::cout << "Value of primal SDP solution tr(Lambda): "
              << sesync_result.trLambda << std::endl;
    std::cout << "SDP duality gap: " << sesync_result.duality_gap << std::endl
              << std::endl;
    std::cout << "SE-SYNCHRONIZATION RESULTS:" << std::endl;
    std::cout << "Value of rounded pose estimates F(x): " << sesync_result.Fxhat
              << std::endl;
    std::cout << "Suboptimality bound F(x) - tr(Lambda) of recovered pose "
                 "estimate: "
              << sesync_result.suboptimality_bound << std::endl
              << std::endl;
    std::cout << "Total elapsed computation time: "
              << sesync_result.total_computation_time << " seconds" << std::endl
              << std::endl;

    std::cout << "===== END SE-SYNC =====" << std::endl << std::endl;
  } // if (options.verbose)
  return sesync_result;
}

SESyncResult SESync(const measurements_t &measurements,
                    const SESyncOpts &options, const Matrix &Y0) {
  if (options.verbose)
    std::cout << "Constructing SE-Sync problem instance ... ";

  auto problem_construction_start_time = Stopwatch::tick();
  SESyncProblem problem(
      measurements, options.formulation, options.projection_factorization,
      options.preconditioner, options.reg_Cholesky_precon_max_condition_number);
  double problem_construction_elapsed_time =
      Stopwatch::tock(problem_construction_start_time);
  if (options.verbose)
    std::cout << "elapsed computation time: "
              << problem_construction_elapsed_time << " seconds" << std::endl
              << std::endl;

  return SESync(problem, options, Y0);
}

bool escape_saddle(const SESyncProblem &problem, const Matrix &Y, Scalar theta,
                   const Vector &v, Scalar gradient_tolerance,
                   Scalar preconditioned_gradient_tolerance, Matrix &Yplus) {

  /** v is an eigenvector corresponding to a negative eigenvalue of Q - Lambda,
   * so the KKT conditions for the semidefinite relaxation are not satisfied;
   * this implies that Y is a saddle point of the rank-restricted semidefinite
   * optimization.  Fortunately, v_min can be used to compute a descent
   * direction from this saddle point, as described in Theorem 3.9 of the paper
   * "A Riemannian Low-Rank Method for Optimization over Semidefinite  Matrices
   * with Block-Diagonal Constraints". Define the vector Ydot := e_{r+1} * v';
   * this is a tangent vector to the domain of the SDP and provides a direction
   * of negative curvature */

  // Function value at current iterate (saddle point)
  Scalar FY = problem.evaluate_objective(Y);

  // Relaxation rank at the NEXT level of the Riemannian Staircase, i.e. we
  // require that r = Y.rows() + 1
  size_t r = problem.relaxation_rank();

  // Construct the corresponding representation of the saddle point Y in the
  // next level of the Riemannian Staircase by adding a row of 0's
  Matrix Y_augmented = Matrix::Zero(r, Y.cols());
  Y_augmented.topRows(r - 1) = Y;

  Matrix Ydot = Matrix::Zero(r, Y.cols());
  Ydot.bottomRows<1>() = v.transpose();

  // Set the initial step length to the greater of 10 times the distance needed
  // to arrive at a trial point whose gradient is large enough to avoid
  // triggering the gradient norm tolerance stopping condition (according to the
  // local second-order model), or at least 2^4 times the minimum admissible
  // steplength,
  Scalar alpha_min = 1e-6; // Minimum stepsize
  Scalar alpha =
      std::max(16 * alpha_min, 10 * gradient_tolerance / fabs(theta));

  // Vectors of trial stepsizes and corresponding function values
  std::vector<double> alphas;
  std::vector<double> fvals;

  /// Backtracking line search
  Matrix Ytest;
  while (alpha >= alpha_min) {

    // Retract along the given tangent vector using the given stepsize
    Ytest = problem.retract(Y_augmented, alpha * Ydot);

    // Ensure that the trial point Ytest has a lower function value than
    // the current iterate Y, and that the gradient at Ytest is
    // sufficiently large that we will not automatically trigger the
    // gradient tolerance stopping criterion at the next iteration
    Scalar FYtest = problem.evaluate_objective(Ytest);
    Matrix grad_FYtest = problem.Riemannian_gradient(Ytest);
    Scalar grad_FYtest_norm = grad_FYtest.norm();
    Scalar preconditioned_grad_FYtest_norm =
        problem.precondition(Ytest, grad_FYtest).norm();

    // Record trial stepsize and function value
    alphas.push_back(alpha);
    fvals.push_back(FYtest);

    if ((FYtest < FY) && (grad_FYtest_norm > gradient_tolerance) &&
        (preconditioned_grad_FYtest_norm > preconditioned_gradient_tolerance)) {
      // Accept this trial point and return success
      Yplus = Ytest;
      return true;
    }
    alpha /= 2;
  }

  // If control reaches here, we failed to find a trial point that satisfied
  // *both* the function decrease *and* gradient bounds.  In order to make
  // forward progress, we will fall back to accepting the trial point that
  // simply minimized the objective value, provided that it strictly *decreased*
  // the objective from the current (saddle) point

  // Find minimum function value from among the trial points
  auto fmin_iter = std::min_element(fvals.begin(), fvals.end());
  auto min_idx = std::distance(fvals.begin(), fmin_iter);

  double f_min = fvals[min_idx];
  double a_min = alphas[min_idx];

  if (f_min < FY) {
    // If this trial point strictly decreased the objective value, accept it and
    // return success
    Yplus = problem.retract(Y_augmented, a_min * Ydot);
    return true;
  } else {
    // NO trial point decreased the objective value: we were unable to escape
    // the saddle point!
    return false;
  }
}
} // namespace SESync
