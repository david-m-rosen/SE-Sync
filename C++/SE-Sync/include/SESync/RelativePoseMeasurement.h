/** A lightweight struct encapsulating a single relative pose-graph SLAM
 * measurement sampled from the generative model given in equation (10) of the
 * SE-Sync tech report.
 *
 * Copyright (C) 2016, 2017 by David M. Rosen
 */

#pragma once

#include <iostream>

#include <Eigen/Dense>

namespace SESync {

/** A simple struct that contains the elements of a relative pose measurement */
struct RelativePoseMeasurement {

  /** 0-based index of first pose */
  size_t i;

  /** 0-based index of second pose */
  size_t j;

  /** Rotational measurement */
  Eigen::MatrixXd R;

  /** Translational measurement */
  Eigen::VectorXd t;

  /** Rotational measurement precision */
  double kappa;

  /** Translational measurement precision */
  double tau;

  /** Simple default constructor; does nothing */
  RelativePoseMeasurement() {}

  /** Basic constructor */
  RelativePoseMeasurement(size_t first_pose, size_t second_pose,
                          const Eigen::MatrixXd &relative_rotation,
                          const Eigen::VectorXd &relative_translation,
                          double rotational_precision,
                          double translational_precision)
      : i(first_pose), j(second_pose), R(relative_rotation),
        t(relative_translation), kappa(rotational_precision),
        tau(translational_precision) {}

  /** A utility function for streaming Nodes to cout */
  inline friend std::ostream &
  operator<<(std::ostream &os, const RelativePoseMeasurement &measurement) {
    os << "i: " << measurement.i << std::endl;
    os << "j: " << measurement.j << std::endl;
    os << "R: " << std::endl << measurement.R << std::endl;
    os << "t: " << std::endl << measurement.t << std::endl;
    os << "Kappa: " << measurement.kappa << std::endl;
    os << "Tau: " << measurement.tau << std::endl;

    return os;
  }
};
}
