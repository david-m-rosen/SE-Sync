
#include <Eigen/QR>
#include <Eigen/SVD>

#include "SESync/StiefelProduct.h"
namespace SESync {

Matrix StiefelProduct::project(const Matrix &A) const {

  // We use a generalization of the well-known SVD-based projection for the
  // orthogonal and special orthogonal groups; see for example Proposition 7
  // in the paper "Projection-Like Retractions on Matrix Manifolds" by Absil
  // and Malick.

  Matrix P(p_, k_ * n_);

#pragma omp parallel for
  for (size_t i = 0; i < n_; ++i) {
    // Compute the (thin) SVD of the ith block of A
    Eigen::JacobiSVD<Matrix> SVD(A.block(0, i * k_, p_, k_),
                                 Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Set the ith block of P to the SVD-based projection of the ith block of A
    P.block(0, i * k_, p_, k_) = SVD.matrixU() * SVD.matrixV().transpose();
  }
  return P;
}

Matrix StiefelProduct::SymBlockDiagProduct(const Matrix &A, const Matrix &B,
                                           const Matrix &C) const {
  // Preallocate result matrix
  Matrix R(p_, k_ * n_);

#pragma omp parallel for
  for (size_t i = 0; i < n_; ++i) {
    // Compute block product Bi' * Ci
    Matrix P =
        B.block(0, i * k_, p_, k_).transpose() * C.block(0, i * k_, p_, k_);
    // Symmetrize this block
    Matrix S = .5 * (P + P.transpose());
    // Compute Ai * S and set corresponding block of R
    R.block(0, i * k_, p_, k_) = A.block(0, i * k_, p_, k_) * S;
  }
  return R;
}

Matrix StiefelProduct::retract(const Matrix &Y, const Matrix &V) const {

  // We use projection-based retraction, as described in "Projection-Like
  // Retractions on Matrix Manifolds" by Absil and Malick

  return project(Y + V);
}

Matrix StiefelProduct::random_sample(
    const std::default_random_engine::result_type &seed) const {
  // Generate a matrix of the appropriate dimension by sampling its elements
  // from the standard Gaussian
  std::default_random_engine generator(seed);
  std::normal_distribution<Scalar> g;

  Matrix R(p_, k_ * n_);
  for (size_t r = 0; r < p_; ++r)
    for (size_t c = 0; c < k_ * n_; ++c)
      R(r, c) = g(generator);
  return project(R);
}
} // namespace SESync
