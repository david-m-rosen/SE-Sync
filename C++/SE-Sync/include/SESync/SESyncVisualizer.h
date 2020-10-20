/**
 * @file SESyncVisualizer.h
 * @brief Simple 3D visualization tool based on Pangolin.
 * @author Tonio Teran, David Rosen, {teran,dmrosen}@mit.edu
 * Copyright (C) 2016 - 2020 by David M. Rosen (dmrosen@mit.edu)
 */

#pragma once

#include <pangolin/pangolin.h>
#include <pangolin/scene/axis.h>
#include <pangolin/scene/scenehandler.h>

#include <Eigen/Dense>
#include <tuple>
#include <vector>

#include "SESync/SESync.h"

namespace SESync {

/// Trajectory consisiting of vector of Eigen-aligned 4x4 SE(3) matrices.
typedef std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>>
    Trajectory3;

/**
 * @class SESyncVisualizer
 * @brief Class for wrapping OpenGL and Pangoling to visualize in 3D.
 */
class SESyncVisualizer {
public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  /**
   * @brief Default constructor.
   * @param result  Structure with SE-Sync results.
   */
  explicit SESyncVisualizer(const size_t num_poses,
                            const measurements_t &measurements,
                            const SESyncOpts &options);

  /**
   * @brief Default destructor.
   */
  ~SESyncVisualizer();

  /**
   * @brief Main visualization for Simple3dWorld that does all the drawing.
   */
  void RenderWorld();

private:
  /**
   * @brief Renders the trajectory as a sequence of triads.
   * @param[in] trajectory Eigen-aligned vector of 3D poses.
   */
  void DrawTrajectory(const Trajectory3 &trajectory,
                      const double axesLength = 0.2) const;

  /**
   * @brief Parse an Xhat solution matrix into individual poses.
   * @param[in] xhat  Solution matrix with [translations|Rotations].
   * @return Vector with parsed poses.
   */
  Trajectory3 ParseXhatToVector(const Matrix &xhat) const;

  /**
   * @brief Anchor solution at the origin.
   * @param[in] xhat  Solution matrix with [translations|Rotations].
   * @return Solution matrix with first pose at the origin.
   */
  Matrix AnchorSolution(const Matrix &xhat) const;

  std::vector<Matrix> iterates_;       ///< Rounded iterates for visualization.
  std::vector<Trajectory3> solutions_; ///< Parsed solutions.

  Trajectory3 est_; ///< Current state estimate trajectory.
  Trajectory3 tgt_; ///< Current trajectory estimate for the target.

  float w_ = 1200.0f; ///< Width of the screen [px].
  float h_ = 800.0f;  ///< Heigh of the screen [px].
  float f_ = 300.0f;  ///< Focal distance of the visualization camera [px].

  measurements_t measurements_;            ///< Relative pose measurements data.
  size_t num_poses_;                       ///< Total number of poses.
  size_t dim_;                             ///< Dimension of the problem.
  SESyncOpts options_;                     ///< Initial description of problem.
  SESyncResult result_;                    ///< Bundle of magic is all here.
  std::shared_ptr<SESyncProblem> problem_; ///< Problem instance.
};

} // namespace SESync
