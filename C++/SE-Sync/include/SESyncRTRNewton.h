#ifndef _SESYNCRTRNEWTON_H_
#define _SESYNCRTRNEWTON_H_

#include "RTRNewton.h"
#include "SolversLS.h"

namespace SESync {

class SESyncRTRNewton : public ROPTLIB::RTRNewton {

public:
  typedef ROPTLIB::RTRNewton Base; // Typedef for base class

  SESyncRTRNewton(const ROPTLIB::Problem *prob, const ROPTLIB::Variable *x0)
      : Base(prob, x0) {}

  ~SESyncRTRNewton() {} // Nothing to do here
};

} // namespace SESync

#endif // _SESYNCRTRNEWTON_H_
