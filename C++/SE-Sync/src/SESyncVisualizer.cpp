/**
 * @file SESyncVisualizer.h
 * @brief Simple 3D visualization tool based on Pangolin.
 * @author Tonio Teran, David Rosen, {teran,dmrosen}@mit.edu
 * Copyright (C) 2016 - 2020 by David M. Rosen (dmrosen@mit.edu)
 */

#include "SESync/SESyncVisualizer.h"

#include <Eigen/StdVector>

namespace SESync {

/* *************************************************************************  */
SESyncVisualizer::SESyncVisualizer(const size_t num_poses,
                                   const measurements_t &measurements,
                                   const SESyncOpts &options)
    : num_poses_(num_poses), measurements_(measurements), options_(options) {
  // Build an SE-Sync problem.
  problem_ = std::make_shared<SESyncProblem>(
      measurements_, options_.formulation, options_.projection_factorization,
      options_.preconditioner,
      options_.reg_Cholesky_precon_max_condition_number);

  // Solve.
  result_ = SESync(*problem_, options_);

  // Round all iterates using approriate rank.
  std::cout << "Rounding iterates..." << std::endl;
  std::vector<std::vector<Matrix>> iterates = result_.iterates;
  size_t rank_counter = options_.r0;
  for (size_t l = 0; l < iterates.size(); l++) {
    problem_->set_relaxation_rank(rank_counter);
    for (const Matrix &Y : iterates[l]) {
      iterates_.push_back(problem_->round_solution(Y));
    }
    rank_counter++;
  }

  // Get problem's dimension from optimal solution.
  dim_ = result_.xhat.rows();

  // Parse all solutions into vectors of matrices.
  std::cout << "Parsing rounded solutions..." << std::endl;
  for (const Matrix &xhat : iterates_) {
    solutions_.push_back(ParseXhatToVector(xhat));
    solnspind_.push_back(ParseXhatToVector(RotateSolution(xhat)));
    lcs_.push_back(std::vector<Eigen::Vector3d>());     // Empty.
    lcspind_.push_back(std::vector<Eigen::Vector3d>()); // Empty.
  }

  // Detect loop closures for drawing.
  std::cout << "Analyzing loop closures..." << std::endl;
  size_t M = measurements_.size();
  for (size_t m = 0; m < M; m++) {
    int i = measurements_[m].i;
    int j = measurements_[m].j;
    if (std::abs(j - i) != 1) {
      for (size_t v = 0; v < solutions_.size(); v++) {
        // "Natural" solutions in the parameter space.
        lcs_[v].push_back(solutions_[v][i].block(0, 3, 3, 1)); // i-th position.
        lcs_[v].push_back(solutions_[v][j].block(0, 3, 3, 1)); // j-th position.
        // Pinned and rotated solutions.
        lcspind_[v].push_back(solnspind_[v][i].block(0, 3, 3, 1));
        lcspind_[v].push_back(solnspind_[v][j].block(0, 3, 3, 1));
      }
    }
  }
}

/* *************************************************************************  */
SESyncVisualizer::~SESyncVisualizer() {}

/* *************************************************************************  */
void SESyncVisualizer::RenderSynchronization() {
  std::cout << "Starting the visualization thread." << std::endl;
  pangolin::CreateWindowAndBind("SE-Sync Super Official Viewer", w_, h_);
  glEnable(GL_DEPTH_TEST);

  // Define Projection and initial ModelView matrix.
  pangolin::OpenGlMatrix proj =
      pangolin::ProjectionMatrix(w_, h_, f_, f_, w_ / 2.0, h_ / 2.0, 0.2, 1000);
  pangolin::OpenGlRenderState s_cam(
      proj,
      pangolin::ModelViewLookAt(1.0, 1.0, 1.0, 0.0, 0.0, 0.0, pangolin::AxisZ));

  // Create the camera for the view.
  pangolin::View &d_cam = pangolin::CreateDisplay()
                              .SetBounds(0.0, 1.0, 0.0, 1.0)
                              .SetHandler(new pangolin::Handler3D(s_cam));

  // Real-time toggles using key presses.
  bool show_z0 = true; // For the XY-plane grid.
  pangolin::RegisterKeyPressCallback('z', [&]() { show_z0 = !show_z0; });
  bool restart = false; // For restarting the visualization.
  pangolin::RegisterKeyPressCallback('r', [&]() { restart = !restart; });
  bool pin = false; // For pinning and rotating iterates to the origin.
  pangolin::RegisterKeyPressCallback('p', [&]() { pin = !pin; });

  // Book-keeping variables for advancing between iterates.
  auto clock = Stopwatch::tick();
  double time = 0; // [s].
  size_t counter = 0, soln_idx = 0, num_iters = solutions_.size();
  double dt = 0.5; // The desired time between iterate visualization [s].

  while (!pangolin::ShouldQuit()) {
    // Clear screen and activate view to render into
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    d_cam.Activate(s_cam);
    glClearColor(1.0f, 1.0f, 1.0f, 0.0f); // Background color.

    // Show either the "natural" solution or the rotated and pinned one.
    if (pin) { // Rotated and pinned.
      DrawIterate(solnspind_[soln_idx], lcspind_[soln_idx]);
    } else { // In parameter space.
      DrawIterate(solutions_[soln_idx], lcs_[soln_idx]);
    }

    s_cam.Apply();
    glColor3f(1.0, 1.0, 1.0);
    if (show_z0) // Draw XY plane, if desired (toggle with 'z' key).
      pangolin::glDraw_z0(10.0, 10);

    time = Stopwatch::tock(clock);
    counter++;
    // Advance to next iterate at the desired time (`dt`).
    if (fmod(time, dt) < 1e-4 && counter > 100) {
      counter = 0;
      soln_idx++;
      if (soln_idx == num_iters)
        soln_idx = num_iters - 1;
    }

    if (restart) { // Restart the iterates (toggle with 'r' key).
      restart = false;
      soln_idx = 0;
    }

    pangolin::FinishFrame(); // Swap frames and process events.
  }
}

/* ************************************************************************** */
void SESyncVisualizer::DrawIterate(
    const Trajectory3 &trajectory,
    const std::vector<Eigen::Vector3d> &lcs) const {
  std::vector<Eigen::Vector3d> positions;

  // Get all positions.
  for (const Eigen::Matrix4d &p : trajectory) {
    positions.push_back(p.block<3, 1>(0, 3));
    pangolin::glDrawAxis(p, 0.1); // TEMP.
  }

  // Draw a line connecting all poses.
  glColor4f(0.2, 0.2, 1.0, 0.8);
  glLineWidth(1.0);
  pangolin::glDrawLineStrip(positions); // Odometry lines.
  glPointSize(2.0);
  glColor4f(0.6, 0.6, 1.0, 0.3);
  pangolin::glDrawPoints(positions); // Position points.

  // Draw loop closures.
  glColor4f(0.7, 0.7, 1.0, 0.3);
  glLineWidth(1.0);
  pangolin::glDrawLines(lcs); // Loop-closing lines.
  glLineWidth(1.0);
}

/* ************************************************************************** */
Trajectory3 SESyncVisualizer::ParseXhatToVector(const Matrix &xhat) const {
  Trajectory3 poses;
  poses.reserve(num_poses_);

  for (size_t i = 0; i < num_poses_; i++) {
    Eigen::Matrix4d p = Eigen::Matrix4d::Identity();
    p.block(0, 0, dim_, dim_) =
        xhat.block(0, num_poses_ + i * dim_, dim_, dim_);
    p.block(0, 3, dim_, 1) = xhat.block(0, i, dim_, 1);
    poses.push_back(p);
  }

  return poses;
}

/* ************************************************************************** */
Matrix SESyncVisualizer::AnchorSolution(const Matrix &xhat) const {
  // Rotate.
  Matrix xhat_anchored =
      xhat.block(0, num_poses_, dim_, dim_).transpose() * xhat;

  // Pin to origin.
  Matrix anchored_t =
      xhat.block(0, 0, dim_, num_poses_).colwise() - xhat.col(0);
  xhat_anchored.block(0, 0, dim_, num_poses_) = anchored_t;

  return xhat_anchored;
}

/* ************************************************************************** */
Matrix SESyncVisualizer::RotateSolution(const Matrix &xhat) const {
  return xhat.block(0, num_poses_, dim_, dim_).transpose() * xhat;
}

} // namespace SESync
