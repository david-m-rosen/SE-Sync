/** This lightweight class models the geometry of M = St(k, p)^n, the n-fold
 * product of the Stiefel manifold St(k,p) (the manifold of orthonormal k-frames
 * in p-dimensional space).  Elements of this manifold (and its tangent spaces)
 * are represented as p x kn matrices of type 'MatrixType'
 *
 * Copyright (C) 2017 by David M. Rosen
 */

#pragma once

#include <Eigen/Dense>

class StiefelProduct {

private:
  // Number of vectors in each orthonormal k-frame
  unsigned int _k;

  // Dimension of ambient Euclidean space containing the frames
  unsigned int _p;

  // Number of copies of St(k,p) in the product
  unsigned int _n;

public:
  /// CONSTRUCTORS AND MUTATORS

  // Default constructor -- sets all dimensions to 0
  StiefelProduct() {}

  StiefelProduct(unsigned int k, unsigned int p, unsigned int n)
      : _k(k), _p(p), _n(n) {}

  void set_k(unsigned int k) { _k = k; }
  void set_p(unsigned int p) { _p = p; }
  void set_n(unsigned int n) { _n = n; }

  /// ACCESSORS
  unsigned int get_k() const { return _k; }
  unsigned int get_p() const { return _p; }
  unsigned int get_n() const { return _n; }

  /// GEOMETRY

  /** Given a generic matrix A in R^{p x kn}, this function computes the
   * projection of A onto R (closest point in the Frobenius norm sense).  */
  Eigen::MatrixXd project(const Eigen::MatrixXd &A) const;

  /** Helper function -- this computes and returns the product
   *
   *  P = A * SymBlockDiag(B^T * C)
   *
   * where A, B, and C are p x kn matrices (cf. eq. (5) in the SE-Sync tech
   * report).
   */
  Eigen::MatrixXd SymBlockDiagProduct(const Eigen::MatrixXd &A,
                                      const Eigen::MatrixXd &B,
                                      const Eigen::MatrixXd &C) const;

  /** Given an element Y in M and a matrix V in T_X(R^{p x kn}) (that is, a (p
   * x kn)-dimensional matrix V considered as an element of the tangent space to
   * the *entire* ambient Euclidean space at X), this function computes and
   * returns the projection of V onto T_X(M), the tangent space of M at X (cf.
   * eq. (42) in the SE-Sync tech report).*/
  Eigen::MatrixXd Proj(const Eigen::MatrixXd &Y,
                       const Eigen::MatrixXd &V) const {
    return V - SymBlockDiagProduct(Y, Y, V);
  }

  /** Given an element Y in M and a tangent vector V in T_Y(M), this function
   * computes the retraction along V at Y using the QR-based retraction
   * specified in eq. (4.8) of Absil et al.'s  "Optimization Algorithms on
   * Matrix Manifolds").
   */
  Eigen::MatrixXd retract(const Eigen::MatrixXd &Y,
                          const Eigen::MatrixXd &V) const;

  /** Sample a random point on M */
  Eigen::MatrixXd random_sample() const;
};
