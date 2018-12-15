#include <functional>

#include "SESync/SESync.h"
#include "SESync/SESyncProblem.h"
#include "SESync/SESync_types.h"
#include "SESync/SESync_utils.h"

#include "Optimization/Smooth/TNT.h"

namespace SESync {

SESyncResult SESync(SESyncProblem &problem, const SESyncOpts &options,
                    const Matrix &Y0) {

  /// ALGORITHM DATA

  // The current iterate in the Riemannian Staircase
  Matrix Y;

  // A cache variable to store the *Euclidean* gradient at the current iterate Y
  Matrix NablaF_Y;

  // The output results struct that we will return
  SESyncResult SESyncResults;
  SESyncResults.status = RS_ITER_LIMIT;

  /// OPTION PARSING AND OUTPUT TO USER

  if (options.verbose) {
    std::cout << "========= SE-Sync ==========" << std::endl << std::endl;

    std::cout << "ALGORITHM SETTINGS:" << std::endl << std::endl;
    std::cout << "SE-Sync settings:" << std::endl;
    std::cout << " SE-Sync problem formulation: ";
    if (problem.formulation() == Formulation::Simplified)
      std::cout << "Simplified";
    else // Explicit)
      std::cout << "Explicit";
    std::cout << std::endl;
    std::cout << " Initial level of Riemannian staircase: " << options.r0
              << std::endl;
    std::cout << " Maximum level of Riemannian staircase: " << options.rmax
              << std::endl;
    std::cout << " Number of Lanczos vectors to use in minimum eigenvalue "
                 "computation: "
              << options.num_Lanczos_vectors << std::endl;
    std::cout << " Maximum number of iterations for eigenvalue computation: "
              << options.max_eig_iterations << std::endl;
    std::cout << " Tolerance for accepting an eigenvalue as numerically "
                 "nonnegative in optimality verification: "
              << options.min_eig_num_tol << std::endl;
    if (options.formulation == Formulation::Simplified) {
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
    else if (problem.preconditioner() == Preconditioner::IncompleteCholesky)
      std::cout << "incomplete Cholesky preconditioner";
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
  Optimization::Smooth::QuadraticModel<Matrix, Matrix, Matrix> QM =
      [&problem](
          const Matrix &Y, Matrix &grad,
          Optimization::Smooth::LinearOperator<Matrix, Matrix, Matrix> &HessOp,
          Matrix &NablaF_Y) {
        // Compute and cache Euclidean gradient at the current iterate
        NablaF_Y = problem.Euclidean_gradient(Y);

        // Compute Riemannian gradient from Euclidean gradient
        grad = problem.Riemannian_gradient(Y, NablaF_Y);

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
  Optimization::Smooth::RiemannianMetric<Matrix, Matrix, Scalar, Matrix>
      metric = [&problem](const Matrix &Y, const Matrix &V1, const Matrix &V2,
                          const Matrix &NablaF_Y) {
        return (V1 * V2.transpose()).trace();
      };

  // Retraction operator
  Optimization::Smooth::Retraction<Matrix, Matrix, Matrix> retraction =
      [&problem](const Matrix &Y, const Matrix &Ydot, const Matrix &NablaF_Y) {
        return problem.retract(Y, Ydot);
      };

  // Preconditioning operator (optional)
  std::experimental::optional<
      Optimization::Smooth::LinearOperator<Matrix, Matrix, Matrix>>
      precon;
  if (options.preconditioner == Preconditioner::None)
    precon = std::experimental::nullopt;
  else {
    Optimization::Smooth::LinearOperator<Matrix, Matrix, Matrix> precon_op =
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

  SESyncResults.initialization_time = Stopwatch::tock(SESync_start_time);
  if (options.verbose)
    std::cout << " SE-Sync initialization finished; elapsed time: "
              << SESyncResults.initialization_time << " seconds" << std::endl
              << std::endl;

  if (options.verbose) {
    // Compute and display the initial objective value
    std::cout << "Initial objective value: " << problem.evaluate_objective(Y)
              << std::endl;
  }

  /// RIEMANNIAN STAIRCASE

  // Configure optimization parameters
  Optimization::Smooth::TNTParams<Scalar> params;
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
      SESyncResults.status = ELAPSED_TIME;
      break;
    }

    // Set  maximum permitted computation time for this level of the Riemannian
    // Staircase
    params.max_computation_time =
        options.max_computation_time - RTR_iteration_start_time;

    if (options.verbose)
      std::cout << std::endl
                << std::endl
                << "====== RIEMANNIAN STAIRCASE (level r = " << r
                << ") ======" << std::endl
                << std::endl;

    /// Run optimization!
    Optimization::Smooth::TNTResult<Matrix, Scalar> TNTResults =
        Optimization::Smooth::TNT<Matrix, Matrix, Scalar, Matrix>(
            F, QM, metric, retraction, Y, NablaF_Y, precon, params,
            options.user_function);

    // Extract the results
    SESyncResults.Yopt = TNTResults.x;
    SESyncResults.SDPval = TNTResults.f;
    SESyncResults.gradnorm =
        problem.Riemannian_gradient(SESyncResults.Yopt).norm();

    // Record sequence of function values
    SESyncResults.function_values.push_back(TNTResults.objective_values);

    // Record sequence of gradient norm values
    SESyncResults.gradient_norms.push_back(TNTResults.gradient_norms);

    // Record sequence of (# Hessian-vector products)
    SESyncResults.Hessian_vector_products.push_back(
        TNTResults.inner_iterations);

    // Record sequence of elapsed optimization times
    SESyncResults.elapsed_optimization_times.push_back(TNTResults.time);

    // Record sequence of pose estimates, if requested
    if (options.log_iterates)
      SESyncResults.iterates.push_back(TNTResults.iterates);

    if (options.verbose) {
      // Display some output to the user
      std::cout << std::endl
                << "Found first-order critical point with value F(Y) = "
                << SESyncResults.SDPval
                << "!  Elapsed computation time: " << TNTResults.elapsed_time
                << " seconds" << std::endl
                << std::endl;
      std::cout << "Checking second order optimality ... " << std::endl;
    }

    /// Check second-order optimality

    // Compute the minimum eigenvalue lambda and corresponding eigenvector
    // of Q - Lambda
    size_t num_min_eig_iterations;
    auto eig_start_time = Stopwatch::tick();
    bool eigenvalue_convergence = problem.compute_S_minus_Lambda_min_eig(
        SESyncResults.Yopt, SESyncResults.lambda_min, SESyncResults.v_min,
        num_min_eig_iterations, options.max_eig_iterations,
        options.min_eig_num_tol, options.num_Lanczos_vectors);
    double eig_elapsed_time = Stopwatch::tock(eig_start_time);

    // Check eigenvalue convergence
    if (!eigenvalue_convergence) {
      if (options.verbose)
        std::cout << "WARNING!  EIGENVALUE COMPUTATION DID NOT CONVERGE TO "
                     "DESIRED PRECISION!"
                  << std::endl;
      SESyncResults.status = EIG_IMPRECISION;
      break;
    }

    // Record results of eigenvalue computation
    SESyncResults.minimum_eigenvalues.push_back(SESyncResults.lambda_min);
    SESyncResults.minimum_eigenvalue_matrix_ops.push_back(
        num_min_eig_iterations);
    SESyncResults.minimum_eigenvalue_computation_times.push_back(
        eig_elapsed_time);

    // Test nonnegativity of minimum eigenvalue
    if (SESyncResults.lambda_min > -options.min_eig_num_tol) {
      // results.Yopt is a second-order critical point (global optimum)!
      if (options.verbose)
        std::cout << "Found second-order critical point! Minimum eigenvalue: "
                  << SESyncResults.lambda_min
                  << ".  Elapsed computation time: " << eig_elapsed_time
                  << " seconds (" << num_min_eig_iterations
                  << " matrix-vector multiplications)." << std::endl;
      SESyncResults.status = GLOBAL_OPT;
      break;
    } // global optimality
    else {

      /// ESCAPE FROM SADDLE!
      if (options.verbose) {
        std::cout << "Saddle point detected! Minimum eigenvalue: "
                  << SESyncResults.lambda_min
                  << ".  Elapsed computation time: " << eig_elapsed_time
                  << " seconds (" << num_min_eig_iterations
                  << " matrix-vector multiplications)." << std::endl;

        std::cout << "Computing escape direction ... " << std::endl;
      }

      // Augment the rank of the rank-restricted semidefinite relaxation in
      // preparation for ascending to the next level of the Riemannian Staircase
      problem.set_relaxation_rank(r + 1);

      Matrix Yplus;
      bool escape_success =
          escape_saddle(problem, SESyncResults.Yopt, SESyncResults.lambda_min,
                        SESyncResults.v_min, options.grad_norm_tol,
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
        SESyncResults.status = SADDLE_POINT;
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

    switch (SESyncResults.status) {
    case GLOBAL_OPT:
      std::cout << "Found global optimum!" << std::endl;
      break;
    case EIG_IMPRECISION:
      std::cout << "WARNING: Minimum eigenvalue computation did not achieve "
                   "sufficient accuracy; solution may not be globally optimal!"
                << std::endl;
      break;
    case SADDLE_POINT:
      std::cout << "WARNING: Line search was unable to escape saddle point!  "
                   "Solution is not globally optimal!"
                << std::endl;
      break;
    case RS_ITER_LIMIT:
      std::cout << "WARNING: Riemannian Staircase reached the maximum "
                   "permitted level before finding global optimum!"
                << std::endl;
      break;
    case ELAPSED_TIME:
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
  SESyncResults.xhat = problem.round_solution(SESyncResults.Yopt);
  double rounding_elapsed_time = Stopwatch::tock(rounding_start_time);

  if (options.verbose)
    std::cout << "elapsed computation time: " << rounding_elapsed_time
              << " seconds" << std::endl
              << std::endl;

  SESyncResults.total_computation_time = Stopwatch::tock(SESync_start_time);

  /// Compute some additional interesting bits of data

  // Evaluate objective function at ROUNDED solution.
  // Note that since xhat contains the *complete* set of pose estimates, we must
  // extract only the *rotational* elements of xhat if the SE synchronization
  // problem was solved using the simplified formulation
  SESyncResults.Fxhat =
      (options.formulation == Formulation::Simplified
           ? problem.evaluate_objective(SESyncResults.xhat.block(
                 0, problem.num_poses(), problem.dimension(),
                 problem.dimension() * problem.num_poses()))
           : problem.evaluate_objective(SESyncResults.xhat));

  // Compute the primal optimal SDP solution Lambda and its objective value
  Matrix Lambda_blocks = problem.compute_Lambda_blocks(SESyncResults.Yopt);

  SESyncResults.trace_Lambda = 0;
  for (size_t i = 0; i < problem.num_poses(); i++)
    SESyncResults.trace_Lambda +=
        Lambda_blocks
            .block(0, i * problem.dimension(), problem.dimension(),
                   problem.dimension())
            .trace();

  SESyncResults.Lambda =
      problem.compute_Lambda_from_Lambda_blocks(Lambda_blocks);

  // Get the duality gap for the primal-dual pair (Y'*Y, Lambda) of SDP
  // estimates

  SESyncResults.SDP_duality_gap =
      SESyncResults.SDPval - SESyncResults.trace_Lambda;

  // Get an upper bound on the (global) suboptimality of the recovered (rounded)
  // pose estimates
  SESyncResults.suboptimality_upper_bound =
      SESyncResults.Fxhat - SESyncResults.trace_Lambda;

  /// FINAL OUTPUT

  if (options.verbose) {
    std::cout << "SDP RESULTS:" << std::endl;
    std::cout << "Value of dual SDP solution F(Y): " << SESyncResults.SDPval
              << std::endl;
    std::cout << "Norm of Riemannian gradient grad F(Y): "
              << SESyncResults.gradnorm << std::endl;
    std::cout << "Value of primal SDP solution tr(Lambda): "
              << SESyncResults.trace_Lambda << std::endl;
    std::cout << "Minimum eigenvalue of certificate matrix S - Lambda: "
              << SESyncResults.lambda_min << std::endl;
    std::cout << "SDP duality gap: " << SESyncResults.SDP_duality_gap
              << std::endl
              << std::endl;
    std::cout << "SE-SYNCHRONIZATION RESULTS:" << std::endl;
    std::cout << "Value of rounded pose estimates F(x): " << SESyncResults.Fxhat
              << std::endl;
    std::cout
        << "Suboptimality bound F(x) - tr(Lambda) of recovered pose estimate: "
        << SESyncResults.suboptimality_upper_bound << std::endl
        << std::endl;
    std::cout << "Total elapsed computation time: "
              << SESyncResults.total_computation_time << " seconds" << std::endl
              << std::endl;

    std::cout << "===== END SE-SYNC =====" << std::endl << std::endl;
  } // if (options.verbose)
  return SESyncResults;
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

bool escape_saddle(const SESyncProblem &problem, const Matrix &Y,
                   Scalar lambda_min, const Vector &v_min,
                   Scalar gradient_tolerance,
                   Scalar preconditioned_gradient_tolerance, Matrix &Yplus) {

  /** v_min is an eigenvector corresponding to a negative eigenvalue of Q -
   * Lambda, so the KKT conditions for the semidefinite relaxation are not
   * satisfied; this implies that Y is a saddle point of the rank-restricted
   * semidefinite  optimization.  Fortunately, v_min can be used to compute a
   * descent  direction from this saddle point, as described in Theorem 3.9
   * of the paper "A Riemannian Low-Rank Method for Optimization over
   * Semidefinite  Matrices with Block-Diagonal Constraints". Define the vector
   * Ydot := e_{r+1} * v'; this is a tangent vector to the domain of the SDP
   * and provides a direction of negative curvature */

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
  Ydot.bottomRows<1>() = v_min.transpose();

  // Set the initial step length to 100 times the distance needed to
  // arrive at a trial point whose gradient is large enough to avoid
  // triggering the gradient norm tolerance stopping condition,
  // according to the local second-order model
  Scalar alpha = 2 * 1000 * gradient_tolerance / fabs(lambda_min);
  Scalar alpha_min = 1e-32; // Minimum stepsize

  // Initialize line search
  bool escape_success = false;

  Matrix Ytest;
  do {
    alpha /= 2;

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

    if ((FYtest < FY) && (grad_FYtest_norm > gradient_tolerance) &&
        (preconditioned_grad_FYtest_norm > preconditioned_gradient_tolerance))
      escape_success = true;
  } while (!escape_success && (alpha > alpha_min));
  if (escape_success) {
    // Update initialization point for next level in the Staircase
    Yplus = Ytest;
    return true;
  } else {
    // If control reaches here, we exited the loop without finding a suitable
    // iterate, i.e. we failed to escape the saddle point
    return false;
  }
}
} // namespace SESync
