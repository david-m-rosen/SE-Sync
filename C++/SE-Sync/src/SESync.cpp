#include <chrono>
#include <functional>

#include "SESync.h"
#include "SESyncProblem.h"
#include "SESync_types.h"
#include "SESync_utils.h"

#include "RTRNewton.h"
#include "SolversLS.h"

namespace SESync {

/** Helper function; trigger the stopping criterion based upon the gradient norm
 * tolerance stored in the SESyncProblem object*/
bool gradient_norm_stopping_criterion(ROPTLIB::Variable *x, ROPTLIB::Vector *gf,
                                      double f, double ngf, double ngf0,
                                      const ROPTLIB::Problem *problem,
                                      const ROPTLIB::Solvers *solver) {
  const SESyncProblem *SESync_problem_ptr =
      static_cast<const SESyncProblem *>(problem);
  return ngf < SESync_problem_ptr->RTR_gradient_norm_tolerance;
}

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
              << options.tolgradnorm << std::endl;
    std::cout << " Stopping tolerance for relative function decrease: "
              << options.rel_func_decrease_tol << std::endl;
    std::cout << " Maximum number of trust-region iterations: "
              << options.max_RTR_iterations << std::endl;
    std::cout << " Maximum number of truncated conjugate gradient iterations "
                 "per outer iteration: "
              << options.max_tCG_iterations << std::endl
              << std::endl;
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
  SESyncProblem problem(LGrho, A, T, Omega, options.use_Cholesky);
  problem.RTR_gradient_norm_tolerance = options.tolgradnorm;
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
      // Matrix Rinit = chordal_initialization(LGrho,
      // problem.dimension(),options.max_eig_iterations, 100 *
      // options.eig_comp_tol);
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

  if (options.verbose) {
    // Compute and display the initial objective value
    double F0 = (Y * problem.Q_product(Y.transpose())).trace();
    std::cout << "Initial objective value: " << F0;
  }

  /// RIEMANNIAN STAIRCASE

  for (unsigned int r = options.r0; r <= options.rmax; r++) {

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
    ROPTLIB::RTRNewton RTR(&problem, &Yinit_ropt);

    // Set stopping criteria
    RTR.Stop_Criterion = ROPTLIB::StopCrit::FUN_REL;
    RTR.Tolerance = options.rel_func_decrease_tol;
    // Note that this custom stopping criterion is called before, and IN
    // ADDITION TO, the relative function decrease tolerance; thus, by setting
    // both, we enforce stopping base both upon gradient norm AND relative
    // function decrease

    RTR.StopPtr = &gradient_norm_stopping_criterion;
    RTR.Max_Iteration = options.max_RTR_iterations;
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
      double alpha = 2 * 100 * options.tolgradnorm / fabs(results.lambda_min);
      double alpha_min = 1e-6; // Minimum stepsize

      //            // First, double-check that the eigenvalue computation
      //            converged with sufficient accuracy ...
      //            Eigen::VectorXd QminusLambda_vmin(results.v_min.size());
      //            SESyncProblem::QMinusLambdaProdFunctor
      //            QminusLambda(&problem, results.Yopt);
      //            QminusLambda.perform_op(results.v_min.data(),
      //            QminusLambda_vmin.data());
      //            double rayleigh_quotient =
      //            results.v_min.dot(QminusLambda_vmin);

      //            // If the rayleigh quotient has the wrong sign, or if the
      //            Rayleigh quotient is more than an order of magnitude off of
      //            the estimate, the eigenvector is likely not accurate enough
      //            to enable us to escape the saddle
      //            if (rayleigh_quotient > 0 || fabs(rayleigh_quotient /
      //            results.lambda_min) < .1) {
      //                results.status = EIG_IMPRECISION;
      //                if (options.verbose) {
      //                    std::cout << "WARNING!  Eigencomputation is not
      //                    sufficiently accurate to enable escape from saddle
      //                    point!" << std::endl;
      //                    std::cout << "Estimated eigenvalue: " <<
      //                    results.lambda_min << std::endl;
      //                    std::cout << "Rayleigh quotient: " <<
      //                    rayleigh_quotient << std::endl
      //                              << std::endl;
      //                }
      //                break;
      //            }

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
            (FYtest_gradnorm > options.tolgradnorm))
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
