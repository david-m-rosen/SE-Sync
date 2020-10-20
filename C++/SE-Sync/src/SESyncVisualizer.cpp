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

  // Solve (Empty initialization).
  result_ = SESync(*problem_, options_, /*Y0 = */ Matrix());

  // Round all iterates using approriate rank.
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
  for (const Matrix &xhat : iterates_) {
    solutions_.push_back(ParseXhatToVector(AnchorSolution(xhat)));
  }
}

/* *************************************************************************  */
SESyncVisualizer::~SESyncVisualizer() {}

/* *************************************************************************  */
void SESyncVisualizer::RenderWorld() {
  std::cout << "Starting the visualization thread." << std::endl;
  pangolin::CreateWindowAndBind("SE-Sync Super Official Viewer", w_, h_);
  glEnable(GL_DEPTH_TEST);

  // Define Projection and initial ModelView matrix.
  pangolin::OpenGlMatrix proj =
      pangolin::ProjectionMatrix(w_, h_, f_, f_, w_ / 2.0, h_ / 2.0, 0.2, 1000);
  pangolin::OpenGlRenderState s_cam(
      proj,
      pangolin::ModelViewLookAt(1.0, 1.0, 1.0, 0.0, 0.0, 0.0, pangolin::AxisZ));

  // Create the individual cameras for each view.
  pangolin::View &d_cam = pangolin::CreateDisplay()
                              .SetBounds(0.0, 1.0, 0.0, 1.0)
                              .SetHandler(new pangolin::Handler3D(s_cam));
  // Create the image viewer for either mono or stereo.
  pangolin::View &left_cam =
      pangolin::CreateDisplay().SetBounds(0.05, 0.3, 0.05, 0.5);
  pangolin::View &right_cam =
      pangolin::CreateDisplay().SetBounds(0.05, 0.3, 0.5, 0.95);

  // Real-time toggles using key presses.
  bool show_z0 = true;
  pangolin::RegisterKeyPressCallback('z', [&]() { show_z0 = !show_z0; });

  bool restart = false;
  pangolin::RegisterKeyPressCallback('r', [&]() { restart = !restart; });

  // Manage the size of the points.
  glPointSize(3.5); // Default is 1.
  // Useful identity.
  Eigen::Matrix4d I_4x4 = Eigen::Matrix4d::Identity();

  auto clock = Stopwatch::tick();
  double time = 0, dt = 0.5;
  size_t counter = 0, pose_idx = 0, num_iters = solutions_.size();

  while (!pangolin::ShouldQuit()) {
    // Clear screen and activate view to render into
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    // ------------
    // -- 2D plots.

    // -----------
    // -- 3D view.
    d_cam.Activate(s_cam);
    // Background color.
    glClearColor(0.9f, 0.9f, 0.9f, 0.0f);
    // Default line width.
    glLineWidth(1.0);

    DrawTrajectory(solutions_[pose_idx]);

    s_cam.Apply();
    glColor3f(1.0, 1.0, 1.0);
    if (show_z0)
      pangolin::glDraw_z0(1.0, 2);
    glLineWidth(7.5);
    if (show_z0)
      pangolin::glDrawAxis(I_4x4, 0.11);
    glLineWidth(1.0);

    time = Stopwatch::tock(clock);
    counter++;
    // std::cout << time << std::endl;
    // std::cout << fmod(time, 3.0) << std::endl;
    if (fmod(time, dt) < 1e-4 && counter > 100) {
      std::cout << "SHOULD TRIGGER HERE!" << std::endl;
      std::cout << "counter: " << counter << std::endl;
      counter = 0;
      pose_idx++;
      if (pose_idx == num_iters)
        pose_idx = num_iters - 1;
      std::cout << "pose idx: " << pose_idx << std::endl;
      std::cout << "num poses: " << num_iters << std::endl;
    }

    if (restart) {
      restart = false;
      pose_idx = 0;
    }

    // Swap frames and Process Events
    pangolin::FinishFrame();
  }
}

/* ************************************************************************** */
void SESyncVisualizer::DrawTrajectory(const Trajectory3 &trajectory,
                                      const double axesLength) const {
  std::vector<Eigen::Vector3d> positions;

  // Draw all poses.
  for (const Eigen::Matrix4d &p : trajectory) {
    positions.push_back(p.block<3, 1>(0, 3));
  }
  glLineWidth(3.0);
  if (trajectory.size())
    pangolin::glDrawAxis(trajectory.back(), axesLength);

  // Draw a line connecting all poses.
  glColor4f(0.7, 0.7, 0.7, 0.1);
  glLineWidth(3.0);
  pangolin::glDrawLineStrip(positions);
  glLineWidth(1.0);
  glColor3f(1.0, 1.0, 1.0);
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

} // namespace SESync
