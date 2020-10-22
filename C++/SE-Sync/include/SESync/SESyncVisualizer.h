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
#include <unordered_map>
#include <vector>

#include "SESync/SESync.h"

namespace SESync {

/// Trajectory consisiting of vector of Eigen-aligned 4x4 SE(3) matrices.
typedef std::vector<Eigen::Matrix4d, Eigen::aligned_allocator<Eigen::Matrix4d>>
    Trajectory3;

struct VisualizationOpts {
  std::string img_name = "SE-Sync-Iter";   ///< Name for output screenshots.
  std::string img_dir = "SE-Sync-Images";  ///< Output directory name.
  double delay = 0.5;                      ///< Delay between iterates [s].
};

/**
 * @class SESyncVisualizer
 * @brief Class for wrapping OpenGL and Pangoling to visualize in SESync in 3D.
 */
class SESyncVisualizer {
 public:
  /**
   * @brief Default constructor.
   * @param num_poses     Total number of poses in the problem.
   * @param measurements  Vector with relative pose measurements.
   * @param options       Structure with SE-Sync initialization options.
   */
  explicit SESyncVisualizer(
      const size_t num_poses, const measurements_t &measurements,
      const SESyncOpts &options,
      const VisualizationOpts &vopts = VisualizationOpts());

  /**
   * @brief Default destructor.
   */
  ~SESyncVisualizer();

  /**
   * @brief Main visualization for enjoying the synchronization process.
   */
  void RenderSynchronization();

 private:
  /**
   * @brief Renders the solution as points and lines.
   * @param[in] trajectory  Eigen-aligned vector of 3D poses.
   * @param[in] lcs         Vector of 3D position pairs for loop-closing lines.
   * @param[in] marker      Whether to draw point markers or not.
   * @param[in] loops       Whether to draw loop closing lines or not.
   */
  void DrawIterate(const Trajectory3 &trajectory,
                   const std::vector<Eigen::Vector3d> &lcs, const bool marker,
                   const bool loops) const;

  /**
   * @brief Parse an Xhat solution matrix into individual poses.
   * @param[in] xhat  Solution matrix with [translations|Rotations].
   * @return Vector with parsed poses.
   */
  Trajectory3 ParseXhatToVector(const Matrix &xhat) const;

  /**
   * @brief Rotates a solution such that the first pose had identity R.
   * @param[in] xhat  Solution matrix with [translations|Rotations].
   * @return Solution matrix with R0 = I.
   */
  Matrix RotateSolution(const Matrix &xhat) const;

  /**
   * @brief Draws text to screen with iterate number and staircase level.
   * @param[in] iter   Number of the current iterate.
   * @param[in] bkgnd  If true, black background enabled, otherwise white.
   */
  void DrawInfoText(const size_t iter, const bool bkgnd) const;

  /**
   * @brief Creates a padded string with specified number of digits.
   * @param[in] iter  Number to be padded with zeros.
   */
  std::string GetScreenshotName(const size_t iter,
                                const size_t digits = 3) const;

  std::vector<Matrix> iterates_;        ///< Rounded iterates for visualization.
  std::vector<Trajectory3> solutions_;  ///< Parsed solutions, raw.
  std::vector<Trajectory3> solnspind_;  ///< Parsed solutions, pinned and rot'd.
  std::vector<std::vector<Eigen::Vector3d>> lcs_;  ///< Loop closures, natural.
  std::vector<std::vector<Eigen::Vector3d>> lcspind_;  ///< Loop closures, rot.
  std::unordered_map<size_t, size_t> staircase_;  ///< iter -> staircase lvl.

  float w_ = 1200.0f;  ///< Width of the screen [px].
  float h_ = 800.0f;   ///< Heigh of the screen [px].
  float f_ = 300.0f;   ///< Focal distance of the visualization camera [px].

  measurements_t measurements_;  ///< Relative pose measurements data.
  size_t num_poses_;             ///< Total number of poses.
  size_t num_iters_;             ///< Total number of iterates.
  size_t dim_;                   ///< Dimension of the problem.
  VisualizationOpts vopts_;      ///< Visualization related parameters.
  SESyncOpts options_;           ///< Initial description of problem.
  SESyncResult result_;          ///< Bundle of magic is all here.
  std::shared_ptr<SESyncProblem> problem_;  ///< Problem instance.

 public:
  EIGEN_MAKE_ALIGNED_OPERATOR_NEW
};

}  // namespace SESync
