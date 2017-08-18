#include <chrono>
#include <functional>

#include "SESync.h"
#include "SESyncProblem.h"
#include "SESync_types.h"
#include "SESync_utils.h"

#include "SESyncRTRNewton.h"

namespace SESync {

SESyncResult SESync(const std::vector<RelativePoseMeasurement> &measurements,
                    const SESyncOpts &options, const Eigen::MatrixXd &Y0) {

  /// ALGORITHM DATA
  SparseMatrix LGrho;   // Rotational connection Laplacian
  SparseMatrix A;       // Oriented incidence matrix
  DiagonalMatrix Omega; // Translational measurement precision matrix
  SparseMatrix T;       // Translational data matrix
  Matrix Y;             // The current iterate in the Riemannian Staircase

  SparseMatrix B1, B2, B3; // The measurement matrices B1, B2, B3 defined in
                           // equations (69) of the tech report

  SESyncResult results;
  results.status = RS_ITER_LIMIT;

  if (options.verbose) {
    std::cout << "========= SE-Sync ==========" << std::endl << std::endl;

    std::cout << "ALGORITHM SETTINGS:" << std::endl << std::endl;
    std::cout << "SE-Sync settings:" << std::endl;
    std::cout << " Initial level of Riemannian staircase: " << options.r0
              << std::endl;
    std::cout << " Maximum level of Riemannian staircase: " << options.rmax
              << std::endl;
    std::cout << " Relative tolerance for minimum eigenvalue computation in "
                 "optimality verification: "
              << options.eig_comp_tol << std::endl;
    std::cout << " Number of Lanczos vectors to use in minimum eigenvalue "
                 "computation: "
              << options.num_Lanczos_vectors << std::endl;
    std::cout << " Maximum number of iterations for eigenvalue computation: "
              << options.max_eig_iterations << std::endl;
    std::cout << " Tolerance for accepting an eigenvalue as numerically "
                 "nonnegative in optimality verification: "
              << options.min_eig_num_tol << std::endl;
    std::cout << " Using " << (options.use_Cholesky ? "Cholseky" : "QR")
              << " decomposition to compute orthogonal projections"
              << std::endl;
    std::cout << " Initialization method: "
              << (options.use_chordal_initialization ? "chordal" : "random")
              << std::endl
              << std::endl;

    std::cout << "ROPTLIB settings:" << std::endl;
    std::cout << " Stopping tolerance for norm of Riemannian gradient: "
              << options.grad_norm_tol << std::endl;
    std::cout << " Stopping tolerance for relative function decrease: "
              << options.rel_func_decrease_tol << std::endl;
    std::cout << " Maximum number of trust-region iterations: "
              << options.max_RTR_iterations << std::endl;
    std::cout << " Maximum number of truncated conjugate gradient iterations "
                 "per outer iteration: "
              << options.max_tCG_iterations << std::endl;
    std::cout
        << " Preconditioning the truncated conjugate gradient method using ";
    if (options.precon == None)
      std::cout << "the identity preconditioner";
    else if (options.precon == Jacobi)
      std::cout << "Jacobi preconditioning";
    else
      std::cout << "incomplete Cholesky preconditioning";

    std::cout << std::endl << std::endl;
  }

  /// ALGORITHM START
  auto SESync_start_time = std::chrono::high_resolution_clock::now();

  /// INITIALIZATION
  if (options.verbose) {
    std::cout << "INITIALIZATION:" << std::endl;
    std::cout << " Constructing auxiliary data matrices ..." << std::endl;
  }

  // Construct rotational connection Laplacian
  if (options.verbose)
    std::cout << " Constructing rotational connection Laplacian L(G^rho) ... ";
  auto LGrho_start_time = std::chrono::high_resolution_clock::now();
  LGrho = construct_rotational_connection_Laplacian(measurements);
  auto LGrho_counter =
      std::chrono::high_resolution_clock::now() - LGrho_start_time;
  double LGrho_elapsed_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(LGrho_counter)
          .count() /
      1000.0;
  if (options.verbose)
    std::cout << "elapsed computation time: " << LGrho_elapsed_time
              << " seconds" << std::endl;

  // Construct oriented incidence matrix
  if (options.verbose)
    std::cout << " Constructing oriented incidence matrix A ... ";
  auto A_start_time = std::chrono::high_resolution_clock::now();
  A = construct_oriented_incidence_matrix(measurements);
  auto A_counter = std::chrono::high_resolution_clock::now() - A_start_time;
  double A_elapsed_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(A_counter).count() /
      1000.0;
  if (options.verbose)
    std::cout << "elapsed computation time: " << A_elapsed_time << " seconds"
              << std::endl;

  // Construct translational measurement precision matrix
  if (options.verbose)
    std::cout << " Constructing translational measurement precision matrix "
                 "Omega ... ";
  auto Omega_start_time = std::chrono::high_resolution_clock::now();
  Omega = construct_translational_precision_matrix(measurements);
  auto Omega_counter =
      std::chrono::high_resolution_clock::now() - Omega_start_time;
  double Omega_elapsed_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(Omega_counter)
          .count() /
      1000.0;
  if (options.verbose)
    std::cout << "elapsed computation time: " << Omega_elapsed_time
              << " seconds" << std::endl;

  // Construct translational data matrix
  if (options.verbose)
    std::cout << " Constructing translational data matrix T ... ";
  auto T_start_time = std::chrono::high_resolution_clock::now();
  T = construct_translational_data_matrix(measurements);
  auto T_counter = std::chrono::high_resolution_clock::now() - T_start_time;
  double T_elapsed_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(T_counter).count() /
      1000.0;
  if (options.verbose)
    std::cout << "elapsed computation time: " << T_elapsed_time << " seconds"
              << std::endl;

  // Construct measurement matrices B1, B2, B3
  if (options.verbose)
    std::cout << " Constructing measurement matrices B1, B2, B3 ... ";
  auto B_start_time = std::chrono::high_resolution_clock::now();
  construct_B_matrices(measurements, B1, B2, B3);
  auto B_counter = std::chrono::high_resolution_clock::now() - B_start_time;
  double B_elapsed_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(B_counter).count() /
      1000.0;
  if (options.verbose)
    std::cout << "elapsed computation time: " << B_elapsed_time << " seconds"
              << std::endl
              << std::endl;

  /// Construct SESync problem instance

  if (options.verbose)
    std::cout << "Constructing SE-Sync problem instance ... ";

  auto problem_construction_start_time =
      std::chrono::high_resolution_clock::now();
  SESyncProblem problem(LGrho, A, T, Omega, options.use_Cholesky,
                        options.precon);
  auto problem_construction_counter =
      std::chrono::high_resolution_clock::now() -
      problem_construction_start_time;
  auto problem_construction_elapsed_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(
          problem_construction_counter)
          .count() /
      1000.0;
  if (options.verbose)
    std::cout << "elapsed computation time: "
              << problem_construction_elapsed_time << " seconds" << std::endl
              << std::endl;

  // Set initial relaxation rank
  problem.set_relaxation_rank(options.r0);

  /// Construct initial iterate

  if (Y0.size() != 0) {
    if (options.verbose)
      std::cout << "Using user-supplied initial iterate Y0" << std::endl;

    // The user supplied an initial iterate, so check that it has the correct
    // size, and if so, use this
    assert((Y0.rows() == options.r0) &&
           (Y0.cols() == problem.dimension() * problem.num_poses()));

    Y = Y0;
  } else {
    if (options.use_chordal_initialization) {
      if (options.verbose)
        std::cout << "Computing chordal initialization ... ";

      auto chordal_init_start_time = std::chrono::high_resolution_clock::now();
      Matrix Rinit = chordal_initialization(problem.dimension(), B3);
      auto chordal_init_counter =
          std::chrono::high_resolution_clock::now() - chordal_init_start_time;
      double chordal_init_elapsed_time =
          std::chrono::duration_cast<std::chrono::milliseconds>(
              chordal_init_counter)
              .count() /
          1000.0;
      if (options.verbose)
        std::cout << "elapsed computation time: " << chordal_init_elapsed_time
                  << " seconds" << std::endl;

      Y.resize(options.r0, problem.dimension() * problem.num_poses());
      Y.setZero();
      Y.topRows(problem.dimension()) = Rinit;
    } else {
      if (options.verbose)
        std::cout << "Sampling a random point on St(" << options.r0 << ","
                  << problem.dimension() << ")^" << problem.num_poses()
                  << std::endl;

      ROPTLIB::StieVariable St(options.r0, problem.dimension());
      ROPTLIB::ProductElement Yinit(1, &St, problem.num_poses());
      Yinit.RandInManifold();

      Y.resize(options.r0, problem.dimension() * problem.num_poses());
      StiefelProd2Mat(Yinit, Y);
    }
  }

  auto initialization_counter =
      std::chrono::high_resolution_clock::now() - SESync_start_time;
  double initialization_elapsed_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(
          initialization_counter)
          .count() /
      1000.0;

  results.initialization_time = initialization_elapsed_time;

  if (options.verbose)
    std::cout << "SE-Sync initialization finished; elapsed time: "
              << initialization_elapsed_time << " seconds" << std::endl
              << std::endl;

  if (options.verbose) {
    // Compute and display the initial objective value
    double F0 = (Y * problem.Q_product(Y.transpose())).trace();
    std::cout << "Initial objective value: " << F0;
  }

  /// RIEMANNIAN STAIRCASE
  auto riemannian_staircase_start_time =
      std::chrono::high_resolution_clock::now();

  for (unsigned int r = options.r0; r <= options.rmax; r++) {

    // The elapsed time from the start of the Riemannian Staircase algorithm
    // until the start of this iteration of RTR
    double RTR_iteration_start_time =
        std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::high_resolution_clock::now() -
            riemannian_staircase_start_time)
            .count() /
        1000.0;

    if (options.verbose)
      std::cout << std::endl
                << std::endl
                << "====== RIEMANNIAN STAIRCASE (level r = " << r
                << ") ======" << std::endl
                << std::endl;

    // Update rank of relaxation
    problem.set_relaxation_rank(r);

    // Allocate storage for a new iterate
    ROPTLIB::StieVariable St(r, problem.dimension());
    ROPTLIB::ProductElement Yinit_ropt(1, &St, problem.num_poses());
    Yinit_ropt.NewMemoryOnWrite();

    // Initialize the value
    Mat2StiefelProd(Y, Yinit_ropt);

    /// Set up RTR solver!
    SESyncRTRNewton RTR(&problem, &Yinit_ropt, options.grad_norm_tol,
                        options.rel_func_decrease_tol);

    RTR.Max_Iteration = options.max_RTR_iterations;
    RTR.Max_Inner_Iter = options.max_tCG_iterations;
    RTR.maximum_Delta = 1e4;
    RTR.Debug = (options.verbose ? ROPTLIB::DEBUGINFO::ITERRESULT
                                 : ROPTLIB::DEBUGINFO::NOOUTPUT);

    /// RUN RTR!
    RTR.Run();

    // Extract the results
    const ROPTLIB::ProductElement *Yopt_ropt =
        static_cast<const ROPTLIB::ProductElement *>(RTR.GetXopt());

    results.Yopt.resize(r, problem.num_poses() * problem.dimension());
    StiefelProd2Mat(*Yopt_ropt, results.Yopt);
    results.SDPval = RTR.Getfinalfun();
    results.gradnorm = RTR.Getnormgf();

    // Record some interesting info about the solving process

    // Number of iterations in the RTR algorithm
    results.RTR_iterations.push_back(RTR.GetIter());

    // Number of Hessian-vector multiplications
    results.Hessian_multiplications.push_back(RTR.GetnH());

    // Sequence of function values
    results.function_values.insert(results.function_values.end(),
                                   RTR.GetfunSeries(),
                                   RTR.GetfunSeries() + RTR.GetlengthSeries());

    // Sequence of gradient norm values
    results.gradient_norm_values.insert(
        results.gradient_norm_values.end(), RTR.GetgradSeries(),
        RTR.GetgradSeries() + RTR.GetlengthSeries());

    // Elapsed time since the start of the Riemannian Staircase at which these
    // values were obtained
    std::vector<double> RTR_iteration_function_times(
        RTR.GettimeSeries(), RTR.GettimeSeries() + RTR.GetlengthSeries());
    for (unsigned int i = 0; i < RTR_iteration_function_times.size(); i++)
      RTR_iteration_function_times[i] += RTR_iteration_start_time;
    results.elapsed_optimization_times.insert(
        results.elapsed_optimization_times.end(),
        RTR_iteration_function_times.begin(),
        RTR_iteration_function_times.end());

    if (options.verbose) {
      // Display some output to the user
      std::cout << std::endl
                << "Found first-order critical point with value F(Y) = "
                << results.SDPval
                << "!  Elapsed computation time: " << RTR.GetComTime()
                << " seconds" << std::endl
                << std::endl;
      std::cout << "Checking second order optimality ... " << std::endl;
    }

    // Compute the minimum eigenvalue lambda and corresponding eigenvector of Q
    // - Lambda
    auto eig_start_time = std::chrono::high_resolution_clock::now();
    bool eigenvalue_convergence = problem.compute_Q_minus_Lambda_min_eig(
        results.Yopt, results.lambda_min, results.v_min,
        options.max_eig_iterations, options.eig_comp_tol,
        options.num_Lanczos_vectors);
    auto eig_counter =
        std::chrono::high_resolution_clock::now() - eig_start_time;
    double eig_elapsed_time =
        std::chrono::duration_cast<std::chrono::milliseconds>(eig_counter)
            .count() /
        1000.0;

    /// Tests for eigenvalue precision
    if (!eigenvalue_convergence) {
      std::cout << "WARNING!  EIGENVALUE COMPUTATION DID NOT CONVERGE TO "
                   "DESIRED PRECISION!"
                << std::endl;
      results.status = EIG_IMPRECISION;
      break;
    }

    // Test nonnegativity of minimum eigenvalue
    if (results.lambda_min > -options.min_eig_num_tol) {
      // results.Yopt is a second-order critical point!
      if (options.verbose)
        std::cout << "Found second-order critical point! (minimum eigenvalue = "
                  << results.lambda_min
                  << "). Elapsed computation time: " << eig_elapsed_time
                  << " seconds" << std::endl;
      results.status = GLOBAL_OPT;
      break;
    } else {

      /// ESCAPE FROM SADDLE!
      ///
      /// We will perform a backtracking line search along the direction of the
      /// negative eigenvector to escape from the saddle point

      // Set the initial step length to 100 times the distance needed to arrive
      // at a trial point whose gradient is large enough to avoid triggering the
      // RTR gradient norm tolerance stopping condition, according to the local
      // second-order model
      double alpha = 2 * 100 * options.grad_norm_tol / fabs(results.lambda_min);
      double alpha_min = 1e-6; // Minimum stepsize

      if (options.verbose) {
        std::cout << "Saddle point detected (minimum eigenvalue = "
                  << results.lambda_min
                  << "). Elapsed computation time: " << eig_elapsed_time
                  << " seconds" << std::endl;

        std::cout << "Computing escape direction ... " << std::endl;
      }

      /**  lambda_min is a negative eigenvalue of Q - Lambda, so the KKT
 * conditions for the semidefinite relaxation are not satisfied; this
 * implies that Y is a saddle point of the rank-restricted
 * semidefinite optimization.  Fortunately, the eigenvector v_min
 * corresponding to lambda_min can be used to provide a descent
 * direction from this saddle point, as described in Theorem 3.9 ofthe
 * paper "A Riemannian Low-Rank Method for Optimization overSemidefinite
 * Matrices with Block-Diagonal Constraints". Define the vector Ydot :=
 * e_{r+1} * v'; this is tangent to the manifold St(d, r+1)^n at Y,
 * and provides a direction ofnegative curvature */

      // Construct a new matrix Y by augmenting Yopt with a new row of zeros
      // at the bottom
      Y = Eigen::MatrixXd::Zero(r + 1,
                                problem.dimension() * problem.num_poses());
      Y.topRows(r) = results.Yopt;
      Eigen::MatrixXd Ydot = Eigen::MatrixXd::Zero(
          r + 1, problem.num_poses() * problem.dimension());
      Ydot.bottomRows<1>() = results.v_min.transpose();

      // Update the rank of the relaxation
      problem.set_relaxation_rank(r + 1);

      // Allocate new product manifold and tangent variables of the
      // appropriate size
      ROPTLIB::StieVariable St(r + 1, problem.dimension());
      ROPTLIB::StieVector StVec(r + 1, problem.dimension());

      ROPTLIB::ProductElement Xprod(1, &St, problem.num_poses());
      Mat2StiefelProd(Y, Xprod);

      ROPTLIB::ProductElement EtaProd(1, &StVec, problem.num_poses());
      ROPTLIB::ProductElement Ytest(1, &St, problem.num_poses());

      // Initialize line search
      bool escape_success = false;
      do {
        alpha /= 2;

        // Retract along the given tangent vector using the given stepsize
        Mat2StiefelProd(alpha * Ydot, EtaProd);
        problem.GetDomain()->Retraction(&Xprod, &EtaProd, &Ytest);
        Eigen::MatrixXd YtestMat(r + 1,
                                 problem.num_poses() * problem.dimension());
        StiefelProd2Mat(Ytest, YtestMat);

        // Ensure that the trial point Ytest has a lower function value than the
        // current iterate Xprod, and that the gradient at Ytest is sufficiently
        // negative that we will not automatically trigger the gradient
        // tolerance stopping criterion at the next iteration
        double FYtest = problem.f(&Ytest);
        double FYtest_gradnorm =
            2 * problem.Q_product(YtestMat.transpose()).norm();

        if ((FYtest < results.SDPval) &&
            (FYtest_gradnorm > options.grad_norm_tol))
          escape_success = true;
      } while (!escape_success && (alpha > alpha_min));
      if (escape_success) {
        // Update initialization point for next level in the Staircase
        StiefelProd2Mat(Ytest, Y);
      } else {
        std::cout << "WARNING!  BACKTRACKING LINE SEARCH FAILED TO ESCAPE FROM "
                     "SADDLE POINT!"
                  << std::endl;
        results.status = SADDLE_POINT;
        break;
      }
    } // if (saddle point)
  }   // Riemannian Staircase

  /// POST-PROCESSING

  if (options.verbose) {
    std::cout << std::endl
              << std::endl
              << "===== END RIEMANNIAN STAIRCASE =====" << std::endl
              << std::endl;

    switch (results.status) {
    case GLOBAL_OPT:
      std::cout << "Found global optimum!" << std::endl;
      break;
    case EIG_IMPRECISION:
      std::cout << "WARNING: Minimum eigenvalue computation did not achieve "
                   "sufficient accuracy; solution may not be globally optimal!"
                << std::endl;
      break;
    case SADDLE_POINT:
      std::cout << "WARNING: Line-search was unable to escape saddle point!"
                << std::endl;
      break;
    case RS_ITER_LIMIT:
      std::cout << "WARNING:  Riemannian Staircase reached the maximum level "
                   "before finding global optimum!"
                << std::endl;
      break;
    }
  }

  if (options.verbose)
    std::cout << std::endl << "Rounding solution ... ";

  auto rounding_start_time = std::chrono::high_resolution_clock::now();
  results.Rhat = round_solution(results.Yopt, problem.dimension());
  auto rounding_counter =
      std::chrono::high_resolution_clock::now() - rounding_start_time;
  double rounding_elapsed_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(rounding_counter)
          .count() /
      1000.0;

  if (options.verbose)
    std::cout << "elapsed computation time: " << rounding_elapsed_time
              << " seconds" << std::endl;

  results.Fxhat =
      (results.Rhat * problem.Q_product(results.Rhat.transpose())).trace();

  if (options.verbose)
    std::cout << "Recovering translational state estimates ... ";

  auto translations_recovery_start_time =
      std::chrono::high_resolution_clock::now();
  results.that = recover_translations(
      B1, B2, results.Rhat.block(0, 0, problem.dimension(), problem.dimension())
                      .inverse() *
                  results.Rhat);
  auto translations_recovery_counter =
      std::chrono::high_resolution_clock::now() -
      translations_recovery_start_time;
  double translations_recovery_elapsed_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(
          translations_recovery_counter)
          .count() /
      1000.0;

  if (options.verbose)
    std::cout << "elapsed computation time: "
              << translations_recovery_elapsed_time << " seconds" << std::endl
              << std::endl;

  auto SESync_counter =
      std::chrono::high_resolution_clock::now() - SESync_start_time;
  double SESync_elapsed_time =
      std::chrono::duration_cast<std::chrono::milliseconds>(SESync_counter)
          .count() /
      1000.0;

  if (options.verbose) {
    std::cout << "Value of SDP solution F(Y): " << results.SDPval << std::endl;
    std::cout << "Norm of Riemannian gradient grad F(Y): " << results.gradnorm
              << std::endl;
    std::cout << "Minimum eigenvalue of Q - Lambda(Y): " << results.lambda_min
              << std::endl;
    std::cout << "Value of rounded pose estimate Rhat: " << results.Fxhat
              << std::endl;
    std::cout << "Suboptimality bound of recovered pose estimate: "
              << results.Fxhat - results.SDPval << std::endl;
    std::cout << "Total elapsed computation time: " << SESync_elapsed_time
              << std::endl
              << std::endl;

    std::cout << "===== END SE-SYNC =====" << std::endl << std::endl;
  }
  return results;
}
}
