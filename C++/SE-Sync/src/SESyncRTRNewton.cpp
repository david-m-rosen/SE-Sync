#include "SESyncRTRNewton.h"

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
}
