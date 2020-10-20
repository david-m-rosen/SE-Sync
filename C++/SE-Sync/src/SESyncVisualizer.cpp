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
SESyncVisualizer::SESyncVisualizer(const measurements_t &measurements,
                                   const SESyncOpts &options)
    : measurements_(measurements), options_(options) {
  // Build an SE-Sync problem.
  problem_ = std::make_shared<SESyncProblem>(
      measurements_, options_.formulation, options_.projection_factorization,
      options_.preconditioner,
      options_.reg_Cholesky_precon_max_condition_number);

  // Solve (Empty initialization).
  result_ = SESync(*problem_, options_, /*Y0 = */ Matrix());

  // TEMP
  std::vector<std::vector<Matrix>> iterates = result_.iterates;
  std::cout << "Checking out the iterates." << std::endl;
  std::cout << "Staircase levels: " << iterates.size() << std::endl;
  size_t rank_counter = options_.r0;
  for (size_t l = 0; l < iterates.size(); l++) {
    problem_->set_relaxation_rank(rank_counter);
    for (const Matrix &Y : iterates[l]) {
      iterates_.push_back(problem_->round_solution(Y));
      std::cout << "xhat: " << iterates_.back().rows() << "x"
                << iterates_.back().cols() << std::endl;
    }
    std::cout << " - level " << rank_counter << ": " << iterates[l].size()
              << std::endl;
    rank_counter++;
  }
  Matrix Yopt = result_.Yopt;
  std::cout << "Yopt dimensions: " << Yopt.rows() << "x" << Yopt.cols()
            << std::endl;
  std::cout << "Rank counter: " << rank_counter << std::endl;
  std::cout << "recovered iterates: " << iterates_.size() << std::endl;
  // TEMP

  // Round all iterates using approriate rank.
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

  bool show_gtsam = true;
  pangolin::RegisterKeyPressCallback('g', [&]() { show_gtsam = !show_gtsam; });

  bool show_manual = true;
  pangolin::RegisterKeyPressCallback('m',
                                     [&]() { show_manual = !show_manual; });

  // Manage the size of the points.
  glPointSize(3.5); // Default is 1.
  // Useful identity.
  Eigen::Matrix4d I_4x4 = Eigen::Matrix4d::Identity();

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

    if (show_manual) {
      for (const auto &vp : vposes_) {
        glLineWidth(std::get<2>(vp));
        pangolin::glDrawAxis(std::get<0>(vp), std::get<1>(vp));
        glLineWidth(1.0);
      }
    }

    s_cam.Apply();
    glColor3f(1.0, 1.0, 1.0);
    if (show_z0)
      pangolin::glDraw_z0(1.0, 2);
    glLineWidth(7.5);
    if (show_z0)
      pangolin::glDrawAxis(I_4x4, 0.11);
    glLineWidth(1.0);

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

} // namespace SESync
