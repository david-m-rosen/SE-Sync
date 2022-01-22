#include "SESync/RelativePoseMeasurement.h"
#include "SESync/SESync.h"
#include "SESync/SESync_types.h"
#include "SESync/SESync_utils.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

PYBIND11_MODULE(sesync, m) {

  /// Bindings for SE-Sync enum classes

  // Problem formulation
  py::enum_<SESync::Formulation>(m, "Formulation")
      .value("Simplified", SESync::Formulation::Simplified)
      .value("Explicit", SESync::Formulation::Explicit);

  // Matrix factorization to use when computing orthogonal projections
  py::enum_<SESync::ProjectionFactorization>(m, "ProjectionFactorization")
      .value("Cholesky", SESync::ProjectionFactorization::Cholesky)
      .value("QR", SESync::ProjectionFactorization::QR);

  // Preconditioner type
  py::enum_<SESync::Preconditioner>(m, "Preconditioner")
      .value("None", SESync::Preconditioner::None)
      .value("Jacobi", SESync::Preconditioner::Jacobi)
      .value("IncompleteCholesky", SESync::Preconditioner::IncompleteCholesky)
      .value("RegularizedCholesky",
             SESync::Preconditioner::RegularizedCholesky);

  // Initialization method
  py::enum_<SESync::Initialization>(m, "Initialization")
      .value("Chordal", SESync::Initialization::Chordal)
      .value("Random", SESync::Initialization::Random);

  // SE-Sync algorithm termination status
  py::enum_<SESync::SESyncStatus>(m, "SESyncStatus")
      .value("GlobalOpt", SESync::SESyncStatus::GlobalOpt)
      .value("SaddlePoint", SESync::SESyncStatus::SaddlePoint)
      .value("EigImprecision", SESync::SESyncStatus::EigImprecision)
      .value("MaxRank", SESync::SESyncStatus::MaxRank);

  /// Bindings for the RelativePoseMeasurement struct

  py::class_<SESync::RelativePoseMeasurement>(m, "RelativePoseMeasurement")
      .def(py::init<>())
      .def_readwrite("i", &SESync::RelativePoseMeasurement::i)
      .def_readwrite("j", &SESync::RelativePoseMeasurement::j)
      .def_readwrite("R", &SESync::RelativePoseMeasurement::R)
      .def_readwrite("t", &SESync::RelativePoseMeasurement::t)
      .def_readwrite("kappa", &SESync::RelativePoseMeasurement::kappa)
      .def_readwrite("tau", &SESync::RelativePoseMeasurement::tau);

  /// Bindings for the SESyncOpts struct

  py::class_<SESync::SESyncOpts>(m, "SESyncOpts")
      .def(py::init<>())
      .def_readwrite("grad_norm_tol", &SESync::SESyncOpts::grad_norm_tol)
      .def_readwrite("preconditioned_grad_norm_tol",
                     &SESync::SESyncOpts::preconditioned_grad_norm_tol)
      .def_readwrite("rel_func_decrease_tol",
                     &SESync::SESyncOpts::rel_func_decrease_tol)
      .def_readwrite("stepsize_tol", &SESync::SESyncOpts::stepsize_tol)
      .def_readwrite("max_time", &SESync::SESyncOpts::max_computation_time)

      .def_readwrite("max_iterations", &SESync::SESyncOpts::max_iterations)
      .def_readwrite("max_tCG_iterations",
                     &SESync::SESyncOpts::max_tCG_iterations)

      .def_readwrite("STPCG_kappa", &SESync::SESyncOpts::STPCG_kappa)
      .def_readwrite("STPCG_theta", &SESync::SESyncOpts::STPCG_theta)
      .def_readwrite("grad_norm_tol", &SESync::SESyncOpts::grad_norm_tol)

      .def_readwrite("formulation", &SESync::SESyncOpts::formulation)
      .def_readwrite("r0", &SESync::SESyncOpts::r0)
      .def_readwrite("rmax", &SESync::SESyncOpts::rmax)
      .def_readwrite("max_eig_iters", &SESync::SESyncOpts::max_eig_iterations)
      .def_readwrite("min_eig_num_tol", &SESync::SESyncOpts::min_eig_num_tol)
      .def_readwrite("num_Lanczos_vectors",
                     &SESync::SESyncOpts::num_Lanczos_vectors)

      .def_readwrite("projection_factorization",
                     &SESync::SESyncOpts::projection_factorization)
      .def_readwrite(
          "reg_Chol_precon_max_cond",
          &SESync::SESyncOpts::reg_Cholesky_precon_max_condition_number)

      .def_readwrite("initialization", &SESync::SESyncOpts::initialization)

      .def_readwrite("verbose", &SESync::SESyncOpts::verbose)
      .def_readwrite("log_iterates", &SESync::SESyncOpts::log_iterates)
      .def_readwrite("num_threads", &SESync::SESyncOpts::num_threads);

  /// Bindings for the SESyncResult struct

  py::class_<SESync::SESyncResult>(m, "SESyncResult")
      .def(py::init<>())
      .def_readwrite("Yopt", &SESync::SESyncResult::Yopt)
      .def_readwrite("SDPval", &SESync::SESyncResult::SDPval)
      .def_readwrite("gradnorm", &SESync::SESyncResult::gradnorm)
      .def_readwrite("Lambda", &SESync::SESyncResult::Lambda)
      .def_readwrite("trLambda", &SESync::SESyncResult::trLambda)
      .def_readwrite("duality_gap", &SESync::SESyncResult::duality_gap)
      .def_readwrite("lambda_min", &SESync::SESyncResult::lambda_min)
      .def_readwrite("vmin", &SESync::SESyncResult::v_min)
      .def_readwrite("Fxhat", &SESync::SESyncResult::Fxhat)
      .def_readwrite("xhat", &SESync::SESyncResult::xhat)
      .def_readwrite("suboptimality_bound",
                     &SESync::SESyncResult::suboptimality_bound)
      .def_readwrite("total_computation_time",
                     &SESync::SESyncResult::total_computation_time)
      .def_readwrite("initialization_time",
                     &SESync::SESyncResult::initialization_time)
      .def_readwrite("function_values", &SESync::SESyncResult::function_values)
      .def_readwrite("gradient_norms", &SESync::SESyncResult::gradient_norms)
      .def_readwrite("Hessian_vector_products",
                     &SESync::SESyncResult::Hessian_vector_products)
      .def_readwrite("elapsed_optimization_times",
                     &SESync::SESyncResult::elapsed_optimization_times)
      .def_readwrite("min_eigs", &SESync::SESyncResult::minimum_eigenvalues)
      .def_readwrite("min_eig_mat_ops", &SESync::SESyncResult::min_eig_mat_ops)
      .def_readwrite("min_eig_comp_times",
                     &SESync::SESyncResult::min_eig_comp_times)
      .def_readwrite("iterates", &SESync::SESyncResult::iterates)
      .def_readwrite("status", &SESync::SESyncResult::status);

  /// Bindings for (some) of the SESync_utils functions

  // NB:  Here we are actually binding an anonymous lambda function that
  // simply wraps the C++ read_g2o_file.  We do this so that we can modify
  // the signature of the Python binding for this function: the binding will
  // accept a single string providing the path of the file to load, and
  // return a tuple consisting of (1) the list of relative pose
  // measurements, and (2) the number of poses in the pose graph

  m.def("read_g2o_file",

        [](const std::string &filename)
            -> std::pair<SESync::measurements_t, size_t> {
          size_t num_poses;
          SESync::measurements_t measurements =
              SESync::read_g2o_file(filename, num_poses);

          return std::pair<SESync::measurements_t, size_t>(measurements,
                                                           num_poses);
        });

  m.def("construct_oriented_incidence_matrix",
        &SESync::construct_oriented_incidence_matrix);

  m.def("project_to_SOd", &SESync::project_to_SOd);

  m.def("orbit_distance_dS", &SESync::orbit_distance_dS);

  m.def("orbit_distance_dO", &SESync::orbit_distance_dO);

  /// Bindings for the main SESync driver
  m.def(
      "SESync",
      [](const SESync::measurements_t &measurements,
         const SESync::SESyncOpts &options,
         const SESync::Matrix &Y0) -> SESync::SESyncResult {
        return SESync::SESync(measurements, options, Y0);
      },
      py::arg("measurements"), py::arg("options") = SESync::SESyncOpts(),
      py::arg("Y0") = SESync::Matrix());
}
