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
/// 3D Pose with axes length (1st double) and line width (2nd double) for viz.
typedef std::tuple<Eigen::Matrix4d, double, double> VizPose;

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
  explicit SESyncVisualizer(const measurements_t &measurements,
                            const SESyncOpts &options);

  /**
   * @brief Default destructor.
   */
  ~SESyncVisualizer();

  /**
   * @brief Main visualization for Simple3dWorld that does all the drawing.
   */
  void RenderWorld();

  /**
   * @brief Add a visualization pose element.
   * @param[in] vpose   Visualization tuple with pose, axes length, and width.
   */
  void AddVizPose(const VizPose &vpose) { vposes_.push_back(vpose); }

  /**
   * @brief Add a visualization pose element.
   * @param[in] pose     3D pose of triad to visualize.
   * @param[in] length   Length of the pose axes.
   * @param[in] width    Width of the pose axes..
   */
  void AddVizPose(const Eigen::Matrix4d &pose, const double length,
                  const double width) {
    AddVizPose(std::make_tuple(pose, length, width));
  }

private:
  /**
   * @brief Renders the trajectory as a sequence of triads.
   * @param[in] trajectory Eigen-aligned vector of 3D poses.
   */
  void DrawTrajectory(const Trajectory3 &trajectory,
                      const double axesLength = 0.2) const;

  // Manually-modifiable variables.
  std::vector<VizPose> vposes_;  ///< Manually added poses to visualize.
  std::vector<Matrix> iterates_; ///< Rounded iterates for visualization.

  Trajectory3 est_; ///< Current state estimate trajectory.
  Trajectory3 tgt_; ///< Current trajectory estimate for the target.

  float w_ = 1200.0f; ///< Width of the screen [px].
  float h_ = 800.0f;  ///< Heigh of the screen [px].
  float f_ = 300.0f;  ///< Focal distance of the visualization camera [px].

  measurements_t measurements_;            ///< Relative pose measurements data.
  SESyncOpts options_;                     ///< Initial description of problem.
  SESyncResult result_;                    ///< Bundle of magic is all here.
  std::shared_ptr<SESyncProblem> problem_; ///< Problem instance.
};

} // namespace SESync
