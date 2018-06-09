/** This lightweight class models the geometry of M = St(k, p)^n, the n-fold
 * product of the Stiefel manifold St(k,p) (the manifold of orthonormal k-frames
 * in p-dimensional space).  Elements of this manifold (and its tangent spaces)
 * are represented as p x kn matrices of type 'MatrixType'
 *
 * Copyright (C) 2017 by David M. Rosen
 */

#pragma once

#include <Eigen/Dense>

#include "SESync/SESync_types.h"

namespace SESync {

class StiefelProduct {

private:
  // Number of vectors in each orthonormal k-frame
  unsigned int k_;

  // Dimension of ambient Euclidean space containing the frames
  unsigned int p_;

  // Number of copies of St(k,p) in the product
  unsigned int n_;

public:
  /// CONSTRUCTORS AND MUTATORS

  // Default constructor -- sets all dimensions to 0
  StiefelProduct() {}

  StiefelProduct(unsigned int k, unsigned int p, unsigned int n)
      : k_(k), p_(p), n_(n) {}

  void set_k(unsigned int k) { k_ = k; }
  void set_p(unsigned int p) { p_ = p; }
  void set_n(unsigned int n) { n_ = n; }

  /// ACCESSORS
  unsigned int get_k() const { return k_; }
  unsigned int get_p() const { return p_; }
  unsigned int get_n() const { return n_; }

  /// GEOMETRY

  /** Given a generic matrix A in R^{p x kn}, this function computes the
   * projection of A onto R (closest point in the Frobenius norm sense).  */
  Matrix project(const Matrix &A) const;

  /** Helper function -- this computes and returns the product
   *
   *  P = A * SymBlockDiag(B^T * C)
   *
   * where A, B, and C are p x kn matrices (cf. eq. (5) in the SE-Sync tech
   * report).
   */
  Matrix SymBlockDiagProduct(const Matrix &A, const Matrix &B,
                             const Matrix &C) const;

  /** Given an element Y in M and a matrix V in T_X(R^{p x kn}) (that is, a (p
   * x kn)-dimensional matrix V considered as an element of the tangent space to
   * the *entire* ambient Euclidean space at X), this function computes and
   * returns the projection of V onto T_X(M), the tangent space of M at X (cf.
   * eq. (42) in the SE-Sync tech report).*/
  Matrix Proj(const Matrix &Y, const Matrix &V) const {
    return V - SymBlockDiagProduct(Y, Y, V);
  }

  /** Given an element Y in M and a tangent vector V in T_Y(M), this function
   * computes the retraction along V at Y using the QR-based retraction
   * specified in eq. (4.8) of Absil et al.'s  "Optimization Algorithms on
   * Matrix Manifolds").
   */
  Matrix retract(const Matrix &Y, const Matrix &V) const;

  /** Sample a random point on M */
  Matrix random_sample() const;
};

} // namespace SESync
