
#include <Eigen/QR>
#include <Eigen/SVD>

#include <random> // For sampling random points on the manifold

#include "SESync/StiefelProduct.h"

Eigen::MatrixXd StiefelProduct::project(const Eigen::MatrixXd &A) const {

  // We use a generalization of the well-known SVD-based projection for the
  // orthogonal and special orthogonal groups; see for example Proposition 7
  // in the paper "Projection-Like Retractions on Matrix Manifolds" by Absil
  // and Malick.

  Eigen::MatrixXd P(_p, _k * _n);

#pragma omp parallel for
  for (unsigned int i = 0; i < _n; i++) {
    // Compute the (thin) SVD of the ith block of A
    Eigen::JacobiSVD<Eigen::MatrixXd> SVD(
        A.block(0, i * _k, _p, _k), Eigen::ComputeThinU | Eigen::ComputeThinV);

    // Set the ith block of P to the SVD-based projection of the ith block of A
    P.block(0, i * _k, _p, _k) = SVD.matrixU() * SVD.matrixV().transpose();
  }
  return P;
}

Eigen::MatrixXd
StiefelProduct::SymBlockDiagProduct(const Eigen::MatrixXd &A,
                                    const Eigen::MatrixXd &B,
                                    const Eigen::MatrixXd &C) const {
  // Preallocate result matrix
  Eigen::MatrixXd R(_p, _k * _n);

#pragma omp parallel for
  for (unsigned int i = 0; i < _n; i++) {
    // Compute block product Bi' * Ci
    Eigen::MatrixXd P =
        B.block(0, i * _k, _p, _k).transpose() * C.block(0, i * _k, _p, _k);
    // Symmetrize this block
    Eigen::MatrixXd S = .5 * (P + P.transpose());
    // Compute Ai * S and set corresponding block of R
    R.block(0, i * _k, _p, _k) = A.block(0, i * _k, _p, _k) * S;
  }
  return R;
}

Eigen::MatrixXd StiefelProduct::retract(const Eigen::MatrixXd &Y,
                                        const Eigen::MatrixXd &V) const {

  // We use projection-based retraction, as described in "Projection-Like
  // Retractions on Matrix Manifolds" by Absil and Malick

  return project(Y + V);
}

Eigen::MatrixXd StiefelProduct::random_sample() const {
  // Generate a matrix of the appropriate dimension by sampling its elements
  // from the standard Gaussian
  std::default_random_engine generator;
  std::normal_distribution<double> g;

  Eigen::MatrixXd R(_p, _k * _n);
  for (unsigned int r = 0; r < _p; r++)
    for (unsigned int c = 0; c < _k * _n; c++)
      R(r, c) = g(generator);
  return project(R);
}
