#include "SESync/RelativePoseMeasurement.h"
#include "SESync/SESync.h"
#include "SESync/SESyncProblem.h"
#include "SESync/SESync_types.h"
#include "SESync/SESync_utils.h"

#include <tuple>

#include <pybind11/eigen.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(PySESync, m) {

  m.doc() = "A library for certifiably correct synchronization over the "
            "special Euclidean group.";

  /// Bindings for SE-Sync enum classes

  // Problem formulation
  py::enum_<SESync::Formulation>(
      m, "Formulation",
      "The specific formulation of the synchronization problem to solve")
      .value("Simplified", SESync::Formulation::Simplified,
             "Translational states have been analytically eliminated")
      .value("Explicit", SESync::Formulation::Explicit,
             "Translational states are explicitly included in the optimization")
      .value("SOSync", SESync::Formulation::SOSync,
             "Rotation synchronization (rotation averaging) only: ignore all "
             "translational data");

  // Matrix factorization to use when computing orthogonal projections
  py::enum_<SESync::ProjectionFactorization>(
      m, "ProjectionFactorization",
      "The type of cached matrix factorization to use when computing "
      "orthogonal projections")
      .value("Cholesky", SESync::ProjectionFactorization::Cholesky)
      .value("QR", SESync::ProjectionFactorization::QR);

  // Preconditioner type
  py::enum_<SESync::Preconditioner>(m, "Preconditioner")
      .value("None", SESync::Preconditioner::None)
      .value("Jacobi", SESync::Preconditioner::Jacobi)
      .value("RegularizedCholesky",
             SESync::Preconditioner::RegularizedCholesky);

  // Initialization method
  py::enum_<SESync::Initialization>(
      m, "Initialization",
      "The initialization method to use constructing an initial estimate, if "
      "none is provided")
      .value("Chordal", SESync::Initialization::Chordal)
      .value("Random", SESync::Initialization::Random);

  // SE-Sync algorithm termination status
  py::enum_<SESync::SESyncStatus>(
      m, "SESyncStatus", "Termination status flag for the SE-Sync algorithm")
      .value("GlobalOpt", SESync::SESyncStatus::GlobalOpt,
             "The returned estimate is a certified globally optimal solution")
      .value("SaddlePoint", SESync::SESyncStatus::SaddlePoint,
             "SE-Sync converged to a first-order saddle point, but could not "
             "escape via backtracking line-search ")
      .value("EigImprecision", SESync::SESyncStatus::EigImprecision,
             "The algorithm converged to a first-order critical point, but the "
             "minimum-eigenvalue computation did not converge with sufficient "
             "precision to determine its global optimality")
      .value("MaxRank", SESync::SESyncStatus::MaxRank,
             "The algorithm exhausted the maximum number of Riemannian "
             "Staircase iterations before finding an optimal solution")
      .value("ElapsedTime", SESync::SESyncStatus::ElapsedTime,
             "The algorithm exhausted the alloted computation time before "
             "finding an optimal solution");

  /// Bindings for the RelativePoseMeasurement struct

  py::class_<SESync::RelativePoseMeasurement>(m, "RelativePoseMeasurement")
      .def(py::init<>(),
           "Default constructor: produces an ininitialized measurement")
      .def(py::init<size_t, size_t, SESync::Matrix, SESync::Vector,
                    SESync::Scalar, SESync::Scalar>())
      .def_readwrite("i", &SESync::RelativePoseMeasurement::i,
                     "Index of first pose (0-based)")
      .def_readwrite("j", &SESync::RelativePoseMeasurement::j,
                     "Index of second pose (0-based)")
      .def_readwrite("R", &SESync::RelativePoseMeasurement::R,
                     "Rotational measurement")
      .def_readwrite("t", &SESync::RelativePoseMeasurement::t,
                     "Translational measurement")
      .def_readwrite("kappa", &SESync::RelativePoseMeasurement::kappa,
                     "Rotational measurement precision")
      .def_readwrite("tau", &SESync::RelativePoseMeasurement::tau,
                     "Translational measurement precision");

  /// Bindings for the SESyncOpts struct

  py::class_<SESync::SESyncOpts>(
      m, "SESyncOpts", "Options structure for the main SE-Sync algorithm")
      .def(py::init<>())
      .def_readwrite(
          "grad_norm_tol", &SESync::SESyncOpts::grad_norm_tol,
          "Stopping tolerance for the norm of the Riemannian gradient")
      .def_readwrite("preconditioned_grad_norm_tol",
                     &SESync::SESyncOpts::preconditioned_grad_norm_tol,
                     "Stopping tolerance for the norm of the preconditioned "
                     "Riemannian gradient")
      .def_readwrite("rel_func_decrease_tol",
                     &SESync::SESyncOpts::rel_func_decrease_tol,
                     "Stopping tolerance based upon the relative decrease in "
                     "function value between accepted iterations")
      .def_readwrite(
          "stepsize_tol", &SESync::SESyncOpts::stepsize_tol,
          "Stopping criterion based upon the norm of an accepted update step")
      .def_readwrite("max_time", &SESync::SESyncOpts::max_computation_time)

      .def_readwrite(
          "max_iterations", &SESync::SESyncOpts::max_iterations,
          "Maximum permitted number of (outer) iterations of the Riemannian "
          "Trust-Region method when solving the local optimization problem at "
          "each level of the Riemannian Staircase")
      .def_readwrite("max_tCG_iterations",
                     &SESync::SESyncOpts::max_tCG_iterations,
                     "Maximum number of inner (truncated conjugate-gradient) "
                     "iterations to perform per outer iteration")

      .def_readwrite("STPCG_kappa", &SESync::SESyncOpts::STPCG_kappa)
      .def_readwrite("STPCG_theta", &SESync::SESyncOpts::STPCG_theta)

      .def_readwrite(
          "formulation", &SESync::SESyncOpts::formulation,
          "The specific form of the synchronization problem to solve")
      .def_readwrite("r0", &SESync::SESyncOpts::r0,
                     "Initial level of the Riemannian Staircase")
      .def_readwrite("rmax", &SESync::SESyncOpts::rmax,
                     "Maximum level of the Riemannian Staircase to explore")
      .def_readwrite("min_eig_num_tol", &SESync::SESyncOpts::min_eig_num_tol,
                     "Numerical tolerance for accepting the minimum eigenvalue "
                     "of the certificate matrix as nonnegative; this should be "
                     "a small positive constant.")
      .def_readwrite("LOBPCG_block_size",
                     &SESync::SESyncOpts::LOBPCG_block_size,
                     "Block size to use in LOBPCG when computing a minimum "
                     "eigenpair of the certificate matrix")
      .def_readwrite("LOBPCG_max_fill_factor",
                     &SESync::SESyncOpts::LOBPCG_max_fill_factor)
      .def_readwrite("LOBPCG_drop_tol", &SESync::SESyncOpts::LOBPCG_drop_tol)
      .def_readwrite("LOBPCG_max_iterations",
                     &SESync::SESyncOpts::LOBPCG_max_iterations,
                     "Maximum number of LOBPCG iterations to permit for the "
                     "minimum-eigenpair computation")

      .def_readwrite("projection_factorization",
                     &SESync::SESyncOpts::projection_factorization,
                     "Type of cached matrix factorization to use for computing "
                     "orthogonal projections")
      .def_readwrite("preconditioner", &SESync::SESyncOpts::preconditioner,
                     "The preconditioning strategy to use in the Riemannian "
                     "trust-region algorithm")
      .def_readwrite(
          "reg_Chol_precon_max_cond",
          &SESync::SESyncOpts::reg_Cholesky_precon_max_condition_number)

      .def_readwrite("initialization", &SESync::SESyncOpts::initialization,
                     "Initialization method to use for calculating an initial "
                     "iterate Y0, if none was provided ")

      .def_readwrite("verbose", &SESync::SESyncOpts::verbose,
                     "Boolean value indicating whether to print output as the "
                     "algorithm runs")
      .def_readwrite(
          "log_iterates", &SESync::SESyncOpts::log_iterates,
          "If this value is true, SE-Sync will log and return the entire "
          "sequence of iterates generated by the Riemannian Staircase")
      .def_readwrite("num_threads", &SESync::SESyncOpts::num_threads,
                     "Number of threads to use for parallel parallelization");

  /// Bindings for the SESyncResult struct

  py::class_<SESync::SESyncResult>(m, "SESyncResult")
      .def(py::init<>())
      .def_readwrite("Yopt", &SESync::SESyncResult::Yopt,
                     "Estimate of a low-rank factor Y for the global minimizer "
                     "Z = YY' of the SDP relaxation")
      .def_readwrite("SDPval", &SESync::SESyncResult::SDPval,
                     "Value of the SDP objective F(YY')")
      .def_readwrite("gradnorm", &SESync::SESyncResult::gradnorm,
                     "Norm of the Riemannian gradient at Yopt")
      .def_readwrite(
          "Lambda", &SESync::SESyncResult::Lambda,
          "The Lagrange multiplier matrix Lambda corresponding to Yopt")
      .def_readwrite("trLambda", &SESync::SESyncResult::trLambda,
                     "The trace of Lambda; this is the value of the (primal) "
                     "SDP relaxation ")
      .def_readwrite("duality_gap", &SESync::SESyncResult::duality_gap,
                     "Duality gap between the estimates for the primal and "
                     "dual SDP solutions")
      .def_readwrite("Fxhat", &SESync::SESyncResult::Fxhat,
                     "The objective value of the rounded solution xhat")
      .def_readwrite("xhat", &SESync::SESyncResult::xhat,
                     "The rounded solution xhat in SE(d)^n")
      .def_readwrite("suboptimality_bound",
                     &SESync::SESyncResult::suboptimality_bound,
                     "Upper bound on the global suboptimality of the returned "
                     "(rounded) estimates")
      .def_readwrite("total_computation_time",
                     &SESync::SESyncResult::total_computation_time,
                     "Total elapsed computation time for the SE-Sync algorithm")
      .def_readwrite("initialization_time",
                     &SESync::SESyncResult::initialization_time,
                     "Elapsed time needed to compute an initial estimate for "
                     "the Riemannian Staircase")
      .def_readwrite(
          "function_values", &SESync::SESyncResult::function_values,
          "A vector containing the sequence of function values obtained during "
          "optimization at each level of the Riemannian Staircase")
      .def_readwrite("gradient_norms", &SESync::SESyncResult::gradient_norms,
                     "A vector containing the sequence of Riemannian gradient "
                     "norms obtained during optimization at each level of the "
                     "Riemannian Staircase")
      .def_readwrite("preconditioned_gradient_norms",
                     &SESync::SESyncResult::preconditioned_gradient_norms,
                     "A vector containing the sequence of preconditioned "
                     "Riemannian gradient norms obtained during optimization "
                     "at each level of the Riemannian Staircase")
      .def_readwrite(
          "Hessian_vector_products",
          &SESync::SESyncResult::Hessian_vector_products,
          "A vector containing the sequence of "
          "(# Hessian-vector products required) at each level of the "
          "Riemannian Staircase")
      .def_readwrite("update_step_norms",
                     &SESync::SESyncResult::update_step_norms,
                     "A vector containing the sequence of norms of the update "
                     "steps computed during the optimization at each level of "
                     "the Riemannian Staircase")
      .def_readwrite(
          "update_step_M_norms", &SESync::SESyncResult::update_step_M_norms,
          "A vector containing the sequence of M-norms | h_k |_{M_k} of the "
          "update steps h_k computed during the optimization at each level of "
          "the Riemannian Staircase, where M_k is the preconditioner "
          "constructed in the kth iteration")
      .def_readwrite("gain_ratios", &SESync::SESyncResult::gain_ratios,
                     "A vector containing the sequence of gain ratios for the "
                     "update steps computed during optimization at each level "
                     "of the Riemannian Staircase")
      .def_readwrite(
          "elapsed_optimization_times",
          &SESync::SESyncResult::elapsed_optimization_times,
          "A vector containing the sequence of elapsed times in the "
          "optimization at each level of the Riemannian Staircase at which the "
          "corresponding function values and gradients were obtained")
      .def_readwrite(
          "min_eigs", &SESync::SESyncResult::escape_direction_curvatures,
          "A vector containing the sequence of minimum eigenvalues of the "
          "certificate matrix constructed at the critical point recovered from "
          "optimization at each level of the Riemannian Staircase")
      .def_readwrite("LOBPCG_iters", &SESync::SESyncResult::LOBPCG_iters)
      .def_readwrite(
          "verification_times", &SESync::SESyncResult::verification_times,
          "A vector containing the elapsed time of the minimum eigenvalue "
          "computation at each level of the Riemannian Staircase")
      .def_readwrite("iterates", &SESync::SESyncResult::iterates,
                     "If log_iterates = true, this will contain the sequence "
                     "of iterates generated by the TNT method at each level of "
                     "the Riemannian Staircase")
      .def_readwrite("status", &SESync::SESyncResult::status,
                     "Termination status of the SE-Sync algorithm");

  /// Bindings for the SESync_utils functions

  // NB:  Here we are actually binding an anonymous lambda function that
  // simply wraps the C++ read_g2o_file.  We do this so that we can modify
  // the signature of the Python binding for this function: the binding will
  // accept a single string providing the path of the file to load, and
  // return a tuple consisting of (1) the list of relative pose
  // measurements, and (2) the number of poses in the pose graph

  m.def(
      "read_g2o_file",

      [](const std::string &filename)
          -> std::pair<SESync::measurements_t, size_t> {
        size_t num_poses;
        SESync::measurements_t measurements =
            SESync::read_g2o_file(filename, num_poses);

        return std::pair<SESync::measurements_t, size_t>(measurements,
                                                         num_poses);
      },
      "Given the name of a file containing a description of a special "
      "Euclidean synchronization problem expressed in the .g2o format (i.e. "
      "using 'EDGE_SE2' or 'EDGE_SE3:QUAT' measurements), this function "
      "constructs and returns a pair consisting of (1) the corresponding "
      "vector of "
      "RelativePoseMeasurements and (2) the total number of poses in the "
      "pose-graph");

  m.def(
      "construct_rotational_weight_graph_Laplacian",
      &SESync::construct_rotational_weight_graph_Laplacian,
      "Given a vector of relative pose measurements, this function constructs "
      "and returns the Laplacian of the rotational weight graph L(W^rho)");

  m.def(
      "construct_translational_weight_graph_Laplacian",
      &SESync::construct_translational_weight_graph_Laplacian,
      "Given a vector of relative pose measurements, this function constructs "
      "and returns the Laplacian of the translational weight graph L(W^tau)");

  m.def("construct_rotational_connection_Laplacian",
        &SESync::construct_rotational_connection_Laplacian,
        "Given a vector of relative pose measurements, this function "
        "constructs "
        "and  returns the corresponding rotational connection Laplacian ");
  m.def("construct_oriented_incidence_matrix",
        &SESync::construct_oriented_incidence_matrix,
        "Given a list of relative pose measurements, this function constructs "
        "and returns the associated oriented incidence matrix A");
  m.def("construct_translational_precision_matrix",
        &SESync::construct_translational_precision_matrix,
        "Given a list of relative pose measurements, this function constructs "
        "and returns the associated diagonal matrix of translational "
        "measurement precisions");
  m.def("construct_translational_data_matrix",
        &SESync::construct_translational_data_matrix,
        "Given a list of relative pose measurements, this function constructs "
        "and returns the associated matrix of raw translational measurements");
  m.def("construct_M_matrix", &SESync::construct_M_matrix,
        "Given a list of relative pose measurements, this function "
        "constructs and returns the matrix M parameterizing the "
        "translation-explicit formulation ofthe special Euclidean "
        "synchronization problem");
  m.def("project_to_SOd", &SESync::project_to_SOd,
        "Given a square d x d matrix, this function returns a closest element "
        "of SO(d)");

  m.def("constuct_B3_matrix", &SESync::construct_B3_matrix,
        "Given a list of relative pose measurements, this function "
        "constructs and returns the matrix B3 defined in equation (69c) of the "
        "SE-Sync tech report");

  m.def(
      "construct_B1_B2_matrices",
      [](const SESync::measurements_t &measurements)
          -> std::pair<SESync::SparseMatrix, SESync::SparseMatrix> {
        SESync::SparseMatrix B1, B2;

        SESync::construct_B1_B2_matrices(measurements, B1, B2);

        return std::make_pair(B1, B2);
      },
      "Given a list of relative pose measurements, this function "
      "constructs and returns the matrices B1 and B2 defined in equation (69) "
      "of the SE-Sync tech report");
  m.def(
      "chordal_initialization", &SESync::chordal_initialization,
      "Given the measurement matrix B3 defined in equation (69c) of the tech "
      "report and the problem dimension d, this function computes and returns "
      "the corresponding chordal initialization for the rotational states");
  m.def("recover_translations", &SESync::recover_translations,
        "Given the measurement matrices B1 and B2 and a matrix R of rotational "
        "state estimates, this function computes and returns the "
        "corresponding optimal translation estimates");

  m.def(
      "dS",
      [](const SESync::Matrix &X,
         const SESync::Matrix &Y) -> std::pair<SESync::Scalar, SESync::Matrix> {
        SESync::Matrix G_S;

        SESync::Scalar dS = SESync::dS(X, Y, &G_S);

        return std::make_pair(dS, G_S);
      },
      " Given two matrices X, Y in SO(d)^n, this function computes and returns "
      "a pair consisting of (1) the orbit distance d_S(X,Y) between them and "
      "(2) the optimal registration G_S in SO(d) aligning Y to X");

  m.def(
      "dO",
      [](const SESync::Matrix &X,
         const SESync::Matrix &Y) -> std::pair<SESync::Scalar, SESync::Matrix> {
        SESync::Matrix G_O;

        SESync::Scalar dO = SESync::dS(X, Y, &G_O);

        return std::make_pair(dO, G_O);
      },
      " Given two matrices X, Y in O(d)^n, this function computes and returns "
      "a pair consisting of (1) the orbit distance d_O(X,Y) between them and "
      "(2) the optimal registration G_O in O(d) aligning Y to X");

  m.def(
      "fast_verification",
      [](const SESync::SparseMatrix &S, SESync::Scalar eta, size_t nx,
         size_t max_iters, SESync::Scalar max_fill_factor,
         SESync::Scalar drop_tol)
          -> std::tuple<bool, SESync::Scalar, SESync::Vector, size_t> {
        SESync::Scalar theta;
        SESync::Vector x;
        size_t num_iters;

        bool PSD =
            SESync::fast_verification(S, eta, nx, theta, x, num_iters,
                                      max_iters, max_fill_factor, drop_tol);

        return std::make_tuple(PSD, theta, x, num_iters);
      },
      "Given a symmetric sparse matrix S and a numerical tolerance eta, this "
      "function returns a tuple consisting of a boolean value indicating "
      "whether S + eta * I is positive-semidefinite, and if it is not, an "
      "estimate (theta, x) for the minimum eigenpair of S, and the number of "
      "iterations required by LOBPCG to estimate this minimum eigenpair ");

  /// Bindings for SESyncProblem class
  py::class_<SESync::SESyncProblem>(
      m, "SESyncProblem",
      "This class encapsulates an instance of the rank-restricted Riemannian "
      "form of the semidefinite relaxation of a special Euclidean.  It "
      "contains all of the precomputed and cached data matrices necessary to "
      "describe the problem and run the optimization algorithm, as well as "
      "functions for performing geometric operations on the underlying "
      "manifold (tangent space projection and retraction) and evaluating the "
      "optimization objective and its gradient and Hessian operator. ")
      .def(py::init<>(), "Default constructor: produces an empyt "
                         "(uninitialized) problem instance")
      .def(py::init<SESync::measurements_t, SESync::Formulation,
                    SESync::ProjectionFactorization, SESync::Preconditioner,
                    SESync::Scalar>(),
           py::arg("measurements"),
           py::arg("formulation") = SESync::Formulation::Simplified,
           py::arg("projection_factorization") =
               SESync::ProjectionFactorization::Cholesky,
           py::arg("preconditioner") =
               SESync::Preconditioner::RegularizedCholesky,
           py::arg("reg_chol_precon_max_cond") = 1e6, "Basic constructor.")
      .def("set_relaxation_rank", &SESync::SESyncProblem::set_relaxation_rank,
           "Set maximum rank of the rank-restricted semidefinite relaxation.")
      .def("formulation", &SESync::SESyncProblem::formulation,
           "Get the specific formulation of this problem")
      .def(
          "projection_factorization",
          &SESync::SESyncProblem::projection_factorization,
          "Returns the type of matrix factorization used to compute the action "
          "of the orthogonal projection operator Pi when solving a Simplified "
          "instance of the special Euclidean synchronization problem")
      .def("preconditioner", &SESync::SESyncProblem::preconditioner,
           "Get the preconditioning strategy")
      .def("num_states", &SESync::SESyncProblem::num_states,
           "Get the number of states (poses or rotations) appearing in this "
           "problem")
      .def("num_measurements", &SESync::SESyncProblem::num_measurements,
           "Get the number of measurements appearing in this problem")
      .def("dimension", &SESync::SESyncProblem::dimension,
           "Get the dimension d of the group SE(d) or SO(d) over which this "
           "problem is defined")
      .def("relaxation_rank", &SESync::SESyncProblem::relaxation_rank,
           "Get the current relaxation rank r of this problem")
      .def("oriented_incidence_matrix",
           &SESync::SESyncProblem::oriented_incidence_matrix,
           "Returns the oriented incidence matrix A of the underlying "
           "measurement graph over which the problem is defined")
      .def("Pi_product", &SESync::SESyncProblem::Pi_product,
           "Given a matrix X, this function computes and returns the "
           "orthogonal projection Pi*X")
      .def("Q_product", &SESync::SESyncProblem::Q_product,
           "Given a matrix X, computes the product Q*X")
      .def("data_matrix_product", &SESync::SESyncProblem::data_matrix_product,
           "Given a matrix Y, this function computes and returns the matrix "
           "product SY, where S is the symmetric matrix parameterizing the "
           "quadratic objective")
      .def("evaluate_objective", &SESync::SESyncProblem::evaluate_objective,
           "Evaluate the objective of the rank-restricted relaxation")
      .def("Euclidean_gradient", &SESync::SESyncProblem::Euclidean_gradient,
           "Evaluate the Euclidean gradient of the objective")
      .def("Riemannian_gradient",
           static_cast<SESync::Matrix (SESync::SESyncProblem::*)(
               const SESync::Matrix &) const>(
               &SESync::SESyncProblem::Riemannian_gradient),
           "Evaluate the Riemannian gradient of the objective")
      .def("Riemannian_Hessian_vector_product",
           static_cast<SESync::Matrix (SESync::SESyncProblem::*)(
               const SESync::Matrix &, const SESync::Matrix &) const>(
               &SESync::SESyncProblem::Riemannian_Hessian_vector_product),
           "Given a matrix Y in the domain D of the relaxation and a tangent "
           "vector dotY in T_Y(D), this function computes and returns Hess "
           "F(Y)[dotY], the action of the Riemannian Hessian on dotY")
      .def("precondition", &SESync::SESyncProblem::precondition,
           "Given a matrix Y in the domain D of the relaxation and a tangent "
           "vector dotY in T_D(Y), this function applies the selected "
           "preconditioning strategy to dotY")
      .def("tangent_space_projection",
           &SESync::SESyncProblem::tangent_space_projection,
           "Given a matrix Y in the domain D of the relaxation and a tangent "
           "vector dotY of the ambient Eucliean space at Y, this function "
           "computes and returns the orthogonal projection of dotY onto T_D(Y)")
      .def("retract", &SESync::SESyncProblem::retract,
           "Given a matrix Y in the domain of the relaxation and a tangent "
           "vector dotY in T_Y(D), this function computes the retraction of Y "
           "along dotY")
      .def("round_solution", &SESync::SESyncProblem::round_solution,
           "Given a point Y in the domain D of the rank-r relaxation, this "
           "function computes and returns a matrix X = [t|R] composed of "
           "translations and rotations for a set of feasible poses for the "
           "original estimation problem obtained by rounding the point Y")
      .def("compute_Lambda", &SESync::SESyncProblem::compute_Lambda,
           "Given a critical point Y of the rank-r relaxation, this function "
           "computes and returns the corresponding Lagrange multiplier matrix "
           "Lambda")
      .def("chordal_initialization",
           &SESync::SESyncProblem::chordal_initialization,
           "This function computes and returns a chordal initialization for "
           "the rank-restricted semidefinite relaxation")
      .def("random_sample", &SESync::SESyncProblem::random_sample,
           "Randomly sample a point in the domain of the rank-restricted "
           "semidefinite relaxation");

  /// Bindings for the main SESync driver
  m.def(
      "SESync",
      [](const SESync::measurements_t &measurements,
         const SESync::SESyncOpts &options,
         const SESync::Matrix &Y0) -> SESync::SESyncResult {
        // Redirect emitted output from (C++) stdout to (Python) sys.stdout
        py::scoped_ostream_redirect stream(
            std::cout, py::module::import("sys").attr("stdout"));
        return SESync::SESync(measurements, options, Y0);
      },
      py::arg("measurements"), py::arg("options") = SESync::SESyncOpts(),
      py::arg("Y0") = SESync::Matrix(),
      "Main SE-Sync function:  Given a list of relative pose measurements "
      "specifying a special Euclidean synchronization problem, this "
      "function "
      "computes and returns an estimated solution using the SE-Sync "
      "algorithm ");

  m.def(
      "SESync",
      [](SESync::SESyncProblem &problem, const SESync::SESyncOpts &options,
         const SESync::Matrix &Y0) -> SESync::SESyncResult {
        // Redirect emitted output from (C++) stdout to (Python) sys.stdout
        py::scoped_ostream_redirect stream(
            std::cout, py::module::import("sys").attr("stdout"));
        return SESync::SESync(problem, options, Y0);
      },
      py::arg("problem"), py::arg("options") = SESync::SESyncOpts(),
      py::arg("Y0") = SESync::Matrix(),
      "Main SE-Sync function:  Given an SESyncProblem instance, this "
      "function computes and returns an estimated solution using the SE-Sync "
      "algorithm ");
}
