/**
 * @file SESyncVisualizer.h
 * @brief Simple 3D visualization tool based on Pangolin.
 * @author Tonio Teran, David Rosen, {teran,dmrosen}@mit.edu
 * Copyright (C) 2016 - 2020 by David M. Rosen (dmrosen@mit.edu)
 */

#include "SESync/SESyncVisualizer.h"

#include <Eigen/StdVector>
#include <cstdlib>

namespace SESync {

/* *************************************************************************  */
SESyncVisualizer::SESyncVisualizer(const size_t num_poses,
                                   const measurements_t &measurements,
                                   const SESyncOpts &options,
                                   const VisualizationOpts &vopts)
    : num_poses_(num_poses),
      measurements_(measurements),
      options_(options),
      vopts_(vopts) {
  // Build an SE-Sync problem.
  options_.log_iterates = true;  // Ensure this is ON at all times.
  problem_ = std::make_shared<SESyncProblem>(
      measurements_, options_.formulation, options_.projection_factorization,
      options_.preconditioner,
      options_.reg_Cholesky_precon_max_condition_number);

  // Solve.
  result_ = SESync(*problem_, options_);

  // Round all iterates using approriate rank.
  std::cout << "Rounding iterates..." << std::endl;
  std::vector<std::vector<Matrix>> iterates = result_.iterates;
  size_t rank_counter = options_.r0, iter_counter = 0;
  for (size_t l = 0; l < iterates.size(); l++) {
    problem_->set_relaxation_rank(rank_counter);
    for (const Matrix &Y : iterates[l]) {
      iterates_.push_back(problem_->round_solution(Y));
      staircase_.insert({iter_counter++, rank_counter});
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
    lcs_.push_back(std::vector<Eigen::Vector3d>());      // Empty.
    lcspind_.push_back(std::vector<Eigen::Vector3d>());  // Empty.
  }
  // Cache the total number of iterates.
  num_iters_ = solutions_.size();

  // Detect loop closures for drawing.
  std::cout << "Analyzing loop closures..." << std::endl;
  size_t M = measurements_.size();
  for (size_t m = 0; m < M; m++) {
    int i = measurements_[m].i;
    int j = measurements_[m].j;
    if (std::abs(j - i) != 1) {
      for (size_t v = 0; v < solutions_.size(); v++) {
        // "Natural" solutions in the parameter space.
        lcs_[v].push_back(
            solutions_[v][i].block(0, 3, 3, 1));  // i-th position.
        lcs_[v].push_back(
            solutions_[v][j].block(0, 3, 3, 1));  // j-th position.
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
  bool restart = false;  // For restarting the visualization.
  pangolin::RegisterKeyPressCallback('r', [&]() { restart = !restart; });
  bool pin = true;  // For pinning and rotating iterates to the origin.
  pangolin::RegisterKeyPressCallback('p', [&]() { pin = !pin; });
  bool save = false;  // For saving a full round of screenshots.
  pangolin::RegisterKeyPressCallback('s', [&]() { save = !save; });
  bool markers = false;  // For toggling the point markers on and off.
  pangolin::RegisterKeyPressCallback('m', [&]() { markers = !markers; });
  bool bkgnd = false;  // For toggling between dark and light backgrounds.
  pangolin::RegisterKeyPressCallback('b', [&]() { bkgnd = !bkgnd; });
  bool text = true;  // For rendering text information.
  pangolin::RegisterKeyPressCallback('t', [&]() { text = !text; });
  bool loops = true;  // For toggling loop closing lines.
  pangolin::RegisterKeyPressCallback('l', [&]() { loops = !loops; });

  bool take_screenshot = false;  // Auxiliary variable for saving frames.
  std::string mkdir_string = "mkdir -p " + vopts_.img_dir;

  // Book-keeping variables for advancing between iterates.
  auto clock = Stopwatch::tick();
  double time = 0;      // To keep track of elapsed time [s].
  size_t soln_idx = 0;  // Current iterate index [-].
  bool iters_complete = false;

  while (!pangolin::ShouldQuit()) {
    // Clear screen and activate view to render into
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    d_cam.Activate(s_cam);
    if (bkgnd) {
      glClearColor(0.0f, 0.0f, 0.0f, 1.0f);  // Black background.
    } else {
      glClearColor(1.0f, 1.0f, 1.0f, 1.0f);  // White background.
    }

    // Show either the "natural" solution or the rotated and pinned one.
    if (pin) {  // Rotated and pinned.
      DrawIterate(solnspind_[soln_idx], lcspind_[soln_idx], markers, loops);
    } else {  // In parameter space.
      DrawIterate(solutions_[soln_idx], lcs_[soln_idx], markers, loops);
    }

    // Show iterate number and staircase level.
    if (text) DrawInfoText(soln_idx, bkgnd);

    // Advance to next iterate at the desired time.
    time = Stopwatch::tock(clock);
    if (time > (soln_idx + 1) * vopts_.delay && !iters_complete) {
      if (take_screenshot) d_cam.SaveOnRender(GetScreenshotName(soln_idx));
      soln_idx++;
      if (soln_idx == num_iters_) {
        soln_idx = num_iters_ - 1;
        iters_complete = true;
        take_screenshot = false;
      }
    }

    if (save) {                // On save, restart the loop and save on render.
      save = false;            // Debounce.
      restart = true;          // Trigger the restart below.
      take_screenshot = true;  // Save img at each iterate.
      const int err = system(mkdir_string.c_str());  // Create img dir.
    }

    if (restart) {                // Restart the iterates (toggle with 'r' key).
      restart = false;            // Debounce.
      iters_complete = false;     // Restart iterations.
      clock = Stopwatch::tick();  // Restart clock.
      soln_idx = 0;               // Go back to beginning.
    }

    pangolin::FinishFrame();  // Swap frames and process events.
  }
}

/* ************************************************************************** */
void SESyncVisualizer::DrawIterate(const Trajectory3 &trajectory,
                                   const std::vector<Eigen::Vector3d> &lcs,
                                   const bool marker, const bool loops) const {
  std::vector<Eigen::Vector3d> positions;

  // Get all positions.
  for (const Eigen::Matrix4d &p : trajectory) {
    positions.push_back(p.block<3, 1>(0, 3));
  }

  // Draw a line connecting all poses.
  glColor4f(0.2, 0.2, 1.0, 0.8);
  glLineWidth(1.0);
  pangolin::glDrawLineStrip(positions);  // Odometry lines.
  glPointSize(2.0);
  glColor4f(0.6, 0.6, 1.0, 0.3);
  if (marker) pangolin::glDrawPoints(positions);  // Position points.

  // Draw loop closures.
  glColor4f(0.7, 0.7, 1.0, 0.3);
  glLineWidth(1.0);
  if (loops) pangolin::glDrawLines(lcs);  // Loop-closing lines.
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
Matrix SESyncVisualizer::RotateSolution(const Matrix &xhat) const {
  return xhat.block(0, num_poses_, dim_, dim_).transpose() * xhat;
}

/* ************************************************************************** */
void SESyncVisualizer::DrawInfoText(const size_t iter, const bool bkgnd) const {
  // Save previous value.
  GLboolean gl_blend_enabled;
  glGetBooleanv(GL_BLEND, &gl_blend_enabled);

  // Ensure that blending is enabled for rendering text.
  glEnable(GL_BLEND);
  glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

  if (bkgnd) {
    glColor3f(1.0f, 1.0f, 1.0f);  // White text on black background.
  } else {
    glColor3f(0.0f, 0.0f, 0.0f);  // Black text on white background.
  }

  pangolin::GlFont::I()
      .Text("Iterate: %d/%d", iter + 1, num_iters_)
      .DrawWindow(10, 22);
  pangolin::GlFont::I()
      .Text("Staircase level: %d", staircase_.at(iter))
      .DrawWindow(10, 10);

  // Restore previous value.
  if (!gl_blend_enabled) glDisable(GL_BLEND);
}

/* ************************************************************************** */
std::string SESyncVisualizer::GetScreenshotName(const size_t iter,
                                                const size_t digits) const {
  std::string numstr = std::to_string(iter);
  std::string padstr = std::string(digits - numstr.length(), '0') + numstr;
  return vopts_.img_dir + std::string("/") + vopts_.img_name + padstr;
}

}  // namespace SESync
