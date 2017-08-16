/** This file provides a specialized subclass of ROPTLIB::RTRNewton in order to
* enable overriding of the base class's stopping criterion and preconditioning
* functions
*
* dmrosen 16 August 2017
*
*/

#ifndef _SESYNCRTRNEWTON_H_
#define _SESYNCRTRNEWTON_H_

#include "RTRNewton.h"
#include "SolversLS.h"

namespace SESync {

class SESyncRTRNewton : public ROPTLIB::RTRNewton {

protected:
  /// Stopping criteria for the RTR Newton method

  /** Stopping criterion for relative decrease in objective value*/
  double rel_dec_tol;

  /** Stopping criterion for the norm of the Riemannian gradient */
  double grad_norm_tol;

  /** This flag simply keeps track of whether or not the isStopped() function is
   * being called prior to accepting a step; we use this to prevent the relative
   * decrease stopping criterion from immediately triggering */
  bool is_initial_iteration = true;

public:
  typedef ROPTLIB::RTRNewton Base; // Typedef for base class

  SESyncRTRNewton(const ROPTLIB::Problem *prob, const ROPTLIB::Variable *x0,
                  double gradient_norm_tolerance,
                  double relative_decrease_tolerance)
      : Base(prob, x0), grad_norm_tol(gradient_norm_tolerance),
        rel_dec_tol(relative_decrease_tolerance) {}

  /** Overridden virtual function from ROPTLIB::Solvers base class;
   * reimplemented here to enable stopping based upon gradient tolerance OR
   * relative decrease */
  bool IsStopped();

  ~SESyncRTRNewton() {} // Nothing to do here
};

} // namespace SESync

#endif // _SESYNCRTRNEWTON_H_
