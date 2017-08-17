#include "SESyncRTRNewton.h"

#include "SESync_utils.h"

namespace SESync {

bool SESyncRTRNewton::IsStopped() {

  // Stopping criterion based upon elapsed computation time
  if (static_cast<double>(ROPTLIB::getTickCount() - starttime) / CLK_PS >
      TimeBound)
    return true;

  // Stopping criterion based upon gradient norm
  if (ngf < grad_norm_tol)
    return true;

  // Stopping criterion based upon relative norm tolerance
  double relative_decrease = (f1 - f2) / (fabs(f1) + 1e-6);
  if (!is_initial_iteration && relative_decrease < rel_dec_tol)
    return true;

  // We are executing this function, so the next time it is called cannot be its
  // first execution
  is_initial_iteration = false;

  return false;
}

void SESyncRTRNewton::PreConditioner(ROPTLIB::Variable *x, ROPTLIB::Vector *eta,
                                     ROPTLIB::Vector *result) {
  if (problem->preconditioner == None) {
    // Identity preconditioning; this is simply a passthrough
    eta->CopyTo(result);
  } else {
    // Convert the passed tangent vector to an Eigen matrix
    Matrix preconditioned_eta_mat;
    StiefelProd2Mat(*static_cast<ROPTLIB::ProductElement *>(eta), eta_mat);

    if (problem->preconditioner == Jacobi) {
      // Apply Jacobi preconditioning
      preconditioned_eta_mat = eta_mat * problem->JacobiPreconditioner;
    } else {
      // Incomplete Cholesky preconditioning
      preconditioned_eta_mat =
          problem->iChol->solve(eta_mat.transpose()).transpose();
    }

    // Convert back to a ProductElement
    Mat2StiefelProd(preconditioned_eta_mat, *preconditioned_eta);

    // Reproject this onto the tangent space at x
    Prob->GetDomain()->Projection(x, preconditioned_eta, result);
  } // else (problem->preconditioner == None)
}

} // namespace SESync
