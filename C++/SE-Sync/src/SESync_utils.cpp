#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#include <Eigen/Geometry>
#include <Eigen/SPQRSupport>

#include "SESync/SESync_utils.h"

namespace SESync {

measurements_t read_g2o_file(const std::string &filename, size_t &num_poses) {

  // Preallocate output vector
  measurements_t measurements;

  // A single measurement, whose values we will fill in
  SESync::RelativePoseMeasurement measurement;

  // A string used to contain the contents of a single line
  std::string line;

  // A string used to extract tokens from each line one-by-one
  std::string token;

  // Preallocate various useful quantities
  Scalar dx, dy, dz, dtheta, dqx, dqy, dqz, dqw, I11, I12, I13, I14, I15, I16,
      I22, I23, I24, I25, I26, I33, I34, I35, I36, I44, I45, I46, I55, I56, I66;

  size_t i, j;

  // Open the file for reading
  std::ifstream infile(filename);

  num_poses = 0;

  while (std::getline(infile, line)) {
    // Construct a stream from the string
    std::stringstream strstrm(line);

    // Extract the first token from the string
    strstrm >> token;

    if (token == "EDGE_SE2") {
      // This is a 2D pose measurement

      /** The g2o format specifies a 2D relative pose measurement in the
       * following form:
       *
       * EDGE_SE2 id1 id2 dx dy dtheta, I11, I12, I13, I22, I23, I33
       *
       */

      // Extract formatted output
      strstrm >> i >> j >> dx >> dy >> dtheta >> I11 >> I12 >> I13 >> I22 >>
          I23 >> I33;

      // Fill in elements of this measurement

      // Pose ids
      measurement.i = i;
      measurement.j = j;

      // Raw measurements
      measurement.t = Eigen::Matrix<Scalar, 2, 1>(dx, dy);
      measurement.R = Eigen::Rotation2D<Scalar>(dtheta).toRotationMatrix();

      Eigen::Matrix<Scalar, 2, 2> TranCov;
      TranCov << I11, I12, I12, I22;
      measurement.tau = 2 / TranCov.inverse().trace();

      measurement.kappa = I33;

    } else if (token == "EDGE_SE3:QUAT") {
      // This is a 3D pose measurement

      /** The g2o format specifies a 3D relative pose measurement in the
       * following form:
       *
       * EDGE_SE3:QUAT id1, id2, dx, dy, dz, dqx, dqy, dqz, dqw
       *
       * I11 I12 I13 I14 I15 I16
       *     I22 I23 I24 I25 I26
       *         I33 I34 I35 I36
       *             I44 I45 I46
       *                 I55 I56
       *                     I66
       */

      // Extract formatted output
      strstrm >> i >> j >> dx >> dy >> dz >> dqx >> dqy >> dqz >> dqw >> I11 >>
          I12 >> I13 >> I14 >> I15 >> I16 >> I22 >> I23 >> I24 >> I25 >> I26 >>
          I33 >> I34 >> I35 >> I36 >> I44 >> I45 >> I46 >> I55 >> I56 >> I66;

      // Fill in elements of the measurement

      // Pose ids
      measurement.i = i;
      measurement.j = j;

      // Raw measurements
      measurement.t = Eigen::Matrix<Scalar, 3, 1>(dx, dy, dz);
      measurement.R =
          Eigen::Quaternion<Scalar>(dqw, dqx, dqy, dqz).toRotationMatrix();

      // Compute precisions

      // Compute and store the optimal (information-divergence-minimizing) value
      // of the parameter tau
      Eigen::Matrix<Scalar, 3, 3> TranCov;
      TranCov << I11, I12, I13, I12, I22, I23, I13, I23, I33;
      measurement.tau = 3 / TranCov.inverse().trace();

      // Compute and store the optimal (information-divergence-minimizing value
      // of the parameter kappa

      Eigen::Matrix<Scalar, 3, 3> RotCov;
      RotCov << I44, I45, I46, I45, I55, I56, I46, I56, I66;
      measurement.kappa = 3 / (2 * RotCov.inverse().trace());

    } else if ((token == "VERTEX_SE2") || (token == "VERTEX_SE3:QUAT")) {
      // This is just initialization information, so do nothing
      continue;
    } else {
      std::cout << "Error: unrecognized type: " << token << "!" << std::endl;
      assert(false);
    }

    // Update maximum value of poses found so far
    size_t max_pair = std::max<size_t>(measurement.i, measurement.j);

    num_poses = ((max_pair > num_poses) ? max_pair : num_poses);
    measurements.push_back(measurement);
  } // while

  infile.close();

  num_poses++; // Account for the use of zero-based indexing

  return measurements;
}

SparseMatrix
construct_rotational_connection_Laplacian(const measurements_t &measurements) {

  size_t num_poses = 0; // We will use this to keep track of the largest pose
  // index encountered, which in turn provides the number
  // of poses

  size_t d = (!measurements.empty() ? measurements[0].t.size() : 0);

  // Each measurement contributes 2*d elements along the diagonal of the
  // connection Laplacian, and 2*d^2 elements on a pair of symmetric
  // off-diagonal blocks

  size_t measurement_stride = 2 * (d + d * d);

  std::vector<Eigen::Triplet<Scalar>> triplets;
  triplets.reserve(measurement_stride * measurements.size());

  size_t i, j, max_pair;
  for (const SESync::RelativePoseMeasurement &measurement : measurements) {
    i = measurement.i;
    j = measurement.j;

    // Elements of ith block-diagonal
    for (size_t k = 0; k < d; k++)
      triplets.emplace_back(d * i + k, d * i + k, measurement.kappa);

    // Elements of jth block-diagonal
    for (size_t k = 0; k < d; k++)
      triplets.emplace_back(d * j + k, d * j + k, measurement.kappa);

    // Elements of ij block
    for (size_t r = 0; r < d; r++)
      for (size_t c = 0; c < d; c++)
        triplets.emplace_back(i * d + r, j * d + c,
                              -measurement.kappa * measurement.R(r, c));

    // Elements of ji block
    for (size_t r = 0; r < d; r++)
      for (size_t c = 0; c < d; c++)
        triplets.emplace_back(j * d + r, i * d + c,
                              -measurement.kappa * measurement.R(c, r));

    // Update num_poses
    max_pair = std::max<size_t>(i, j);

    if (max_pair > num_poses)
      num_poses = max_pair;
  }

  num_poses++; // Account for 0-based indexing

  // Construct and return a sparse matrix from these triplets
  SparseMatrix LGrho(d * num_poses, d * num_poses);
  LGrho.setFromTriplets(triplets.begin(), triplets.end());

  return LGrho;
}

SparseMatrix
construct_oriented_incidence_matrix(const measurements_t &measurements) {
  std::vector<Eigen::Triplet<Scalar>> triplets;
  triplets.reserve(2 * measurements.size());

  size_t num_poses = 0;
  size_t max_pair;
  for (size_t m = 0; m < measurements.size(); m++) {
    triplets.emplace_back(measurements[m].i, m, -1);
    triplets.emplace_back(measurements[m].j, m, 1);

    max_pair = std::max<size_t>(measurements[m].i, measurements[m].j);
    if (max_pair > num_poses)
      num_poses = max_pair;
  }
  num_poses++; // Account for zero-based indexing

  SparseMatrix A(num_poses, measurements.size());
  A.setFromTriplets(triplets.begin(), triplets.end());
  return A;
}

DiagonalMatrix
construct_translational_precision_matrix(const measurements_t &measurements) {

  // Allocate output matrix
  DiagonalMatrix Omega(measurements.size());

  DiagonalMatrix::DiagonalVectorType &diagonal = Omega.diagonal();

  for (size_t m = 0; m < measurements.size(); m++)
    diagonal[m] = measurements[m].tau;

  return Omega;
}

SparseMatrix
construct_translational_data_matrix(const measurements_t &measurements) {

  size_t num_poses = 0;

  size_t d = (!measurements.empty() ? measurements[0].t.size() : 0);

  std::vector<Eigen::Triplet<Scalar>> triplets;
  triplets.reserve(d * measurements.size());

  size_t max_pair;
  for (size_t m = 0; m < measurements.size(); m++) {
    for (size_t k = 0; k < d; k++)
      triplets.emplace_back(m, d * measurements[m].i + k,
                            -measurements[m].t(k));

    max_pair = std::max<size_t>(measurements[m].i, measurements[m].j);
    if (max_pair > num_poses)
      num_poses = max_pair;
  }
  num_poses++; // Account for zero-based indexing

  SparseMatrix T(measurements.size(), d * num_poses);
  T.setFromTriplets(triplets.begin(), triplets.end());

  return T;
}

void construct_B_matrices(const measurements_t &measurements, SparseMatrix &B1,
                          SparseMatrix &B2, SparseMatrix &B3) {
  // Clear input matrices
  B1.setZero();
  B2.setZero();
  B3.setZero();

  size_t num_poses = 0;
  size_t d = (!measurements.empty() ? measurements[0].t.size() : 0);

  std::vector<Eigen::Triplet<Scalar>> triplets;

  // Useful quantities to cache
  size_t d2 = d * d;
  size_t d3 = d * d * d;

  size_t i, j; // Indices for the tail and head of the given measurement
  Scalar sqrttau;
  size_t max_pair;

  /// Construct the matrix B1 from equation (69a) in the tech report
  triplets.reserve(2 * d * measurements.size());

  for (size_t e = 0; e < measurements.size(); e++) {
    i = measurements[e].i;
    j = measurements[e].j;
    sqrttau = sqrt(measurements[e].tau);

    // Block corresponding to the tail of the measurement
    for (size_t l = 0; l < d; l++) {
      triplets.emplace_back(e * d + l, i * d + l,
                            -sqrttau); // Diagonal element corresponding to tail
      triplets.emplace_back(e * d + l, j * d + l,
                            sqrttau); // Diagonal element corresponding to head
    }

    // Keep track of the number of poses we've seen
    max_pair = std::max<size_t>(i, j);
    if (max_pair > num_poses)
      num_poses = max_pair;
  }
  num_poses++; // Account for zero-based indexing

  B1.resize(d * measurements.size(), d * num_poses);
  B1.setFromTriplets(triplets.begin(), triplets.end());

  /// Construct matrix B2 from equation (69b) in the tech report
  triplets.clear();
  triplets.reserve(d2 * measurements.size());

  for (size_t e = 0; e < measurements.size(); e++) {
    i = measurements[e].i;
    sqrttau = sqrt(measurements[e].tau);
    for (size_t k = 0; k < d; k++)
      for (size_t r = 0; r < d; r++)
        triplets.emplace_back(d * e + r, d2 * i + d * k + r,
                              -sqrttau * measurements[e].t(k));
  }

  B2.resize(d * measurements.size(), d2 * num_poses);
  B2.setFromTriplets(triplets.begin(), triplets.end());

  /// Construct matrix B3 from equation (69c) in the tech report
  triplets.clear();
  triplets.reserve((d3 + d2) * measurements.size());

  for (size_t e = 0; e < measurements.size(); e++) {
    Scalar sqrtkappa = std::sqrt(measurements[e].kappa);
    const Matrix &R = measurements[e].R;

    for (size_t r = 0; r < d; r++)
      for (size_t c = 0; c < d; c++) {
        i = measurements[e].i; // Tail of measurement
        j = measurements[e].j; // Head of measurement

        // Representation of the -sqrt(kappa) * Rt(i,j) \otimes I_d block
        for (size_t l = 0; l < d; l++)
          triplets.emplace_back(e * d2 + d * r + l, i * d2 + d * c + l,
                                -sqrtkappa * R(c, r));
      }

    for (size_t l = 0; l < d2; l++)
      triplets.emplace_back(e * d2 + l, j * d2 + l, sqrtkappa);
  }

  B3.resize(d2 * measurements.size(), d2 * num_poses);
  B3.setFromTriplets(triplets.begin(), triplets.end());
}

SparseMatrix
construct_quadratic_form_data_matrix(const measurements_t &measurements) {

  size_t num_poses = 0;
  size_t d = (!measurements.empty() ? measurements[0].t.size() : 0);

  std::vector<Eigen::Triplet<Scalar>> triplets;

  /// Useful quantities to cache
  size_t d2 = d * d;

  // Number of nonzero elements contributed to L(W^tau) by each measurement
  size_t LWtau_nnz_per_measurement = 4;

  // Number of nonzero elements contributed to V by each measurement
  size_t V_nnz_per_measurement = 2 * d;

  // Number of nonzero elements contributed to L(G^rho) by each measurement
  size_t LGrho_nnz_per_measurement = 2 * d + 2 * d2;

  // Number of nonzero elements contributed to Sigma by each measurement
  size_t Sigma_nnz_per_measurement = d2;

  // Number of nonzero elements contributed to the entire matrix M by each
  // measurement
  size_t num_nnz_per_measurement =
      LWtau_nnz_per_measurement + 2 * V_nnz_per_measurement +
      LGrho_nnz_per_measurement + Sigma_nnz_per_measurement;

  /// Working space
  size_t i, j; // Indices for the tail and head of the given measurement
  size_t max_pair;

  triplets.reserve(num_nnz_per_measurement * measurements.size());

  // Scan through the set of measurements to determine the total number of poses
  // in this problem
  for (const SESync::RelativePoseMeasurement &measurement : measurements) {
    max_pair = std::max<size_t>(measurement.i, measurement.j);
    if (max_pair > num_poses)
      num_poses = max_pair;
  }
  num_poses++; // Account for zero-based indexing

  // Now scan through the measurements again, using knowledge of the total
  // number of poses to compute offsets as appropriate

  for (const SESync::RelativePoseMeasurement &measurement : measurements) {

    i = measurement.i; // Tail of measurement
    j = measurement.j; // Head of measurement

    // Add elements for L(W^tau)
    triplets.emplace_back(i, i, measurement.tau);
    triplets.emplace_back(j, j, measurement.tau);
    triplets.emplace_back(i, j, -measurement.tau);
    triplets.emplace_back(j, i, -measurement.tau);

    // Add elements for V (upper-right block)
    for (size_t k = 0; k < d; k++)
      triplets.emplace_back(i, num_poses + i * d + k,
                            measurement.tau * measurement.t(k));
    for (size_t k = 0; k < d; k++)
      triplets.emplace_back(j, num_poses + i * d + k,
                            -measurement.tau * measurement.t(k));

    // Add elements for V' (lower-left block)
    for (size_t k = 0; k < d; k++)
      triplets.emplace_back(num_poses + i * d + k, i,
                            measurement.tau * measurement.t(k));
    for (size_t k = 0; k < d; k++)
      triplets.emplace_back(num_poses + i * d + k, j,
                            -measurement.tau * measurement.t(k));

    // Add elements for L(G^rho)
    // Elements of ith block-diagonal
    for (size_t k = 0; k < d; k++)
      triplets.emplace_back(num_poses + d * i + k, num_poses + d * i + k,
                            measurement.kappa);

    // Elements of jth block-diagonal
    for (size_t k = 0; k < d; k++)
      triplets.emplace_back(num_poses + d * j + k, num_poses + d * j + k,
                            measurement.kappa);

    // Elements of ij block
    for (size_t r = 0; r < d; r++)
      for (size_t c = 0; c < d; c++)
        triplets.emplace_back(num_poses + i * d + r, num_poses + j * d + c,
                              -measurement.kappa * measurement.R(r, c));

    // Elements of ji block
    for (size_t r = 0; r < d; r++)
      for (size_t c = 0; c < d; c++)
        triplets.emplace_back(num_poses + j * d + r, num_poses + i * d + c,
                              -measurement.kappa * measurement.R(c, r));

    // Add elements for Sigma
    for (size_t r = 0; r < d; r++)
      for (size_t c = 0; c < d; c++)
        triplets.emplace_back(num_poses + i * d + r, num_poses + i * d + c,
                              measurement.tau * measurement.t(r) *
                                  measurement.t(c));
  }

  SparseMatrix M((d + 1) * num_poses, (d + 1) * num_poses);
  M.setFromTriplets(triplets.begin(), triplets.end());

  return M;
}

Matrix chordal_initialization(size_t d, const SparseMatrix &B3) {
  size_t d2 = d * d;
  size_t num_poses = B3.cols() / d2;

  /// We want to find a minimizer of
  /// || B3 * r ||
  ///
  /// For the purposes of initialization, we can simply fix the first pose to
  /// the origin; this corresponds to fixing the first d^2 elements of r to
  /// vec(I_d), and slicing off the first d^2 columns of B3 to form
  ///
  /// min || B3red * rred + c ||, where
  ///
  /// c = B3(1:d^2) * vec(I_3)

  SparseMatrix B3red = B3.rightCols((num_poses - 1) * d2);
  // Must be in compressed format to use Eigen::SparseQR!
  B3red.makeCompressed();

  // Vectorization of I_d
  Matrix Id = Matrix::Identity(d, d);
  Eigen::Map<Vector> Id_vec(Id.data(), d2);

  Vector cR = B3.leftCols(d2) * Id_vec;

  Vector rvec;
  Eigen::SPQR<SparseMatrix> QR(B3red);
  rvec = -QR.solve(cR);

  Matrix Rchordal(d, d * num_poses);
  Rchordal.leftCols(d) = Id;
  Rchordal.rightCols((num_poses - 1) * d) =
      Eigen::Map<Matrix>(rvec.data(), d, (num_poses - 1) * d);

  for (size_t i = 1; i < num_poses; i++)
    Rchordal.block(0, i * d, d, d) =
        project_to_SOd(Rchordal.block(0, i * d, d, d));
  return Rchordal;
}

Matrix recover_translations(const SparseMatrix &B1, const SparseMatrix &B2,
                            const Matrix &R) {
  size_t d = R.rows();
  size_t n = R.cols() / d;

  /// We want to find a minimizer of
  /// || B1 * t + B2 * vec(R) ||
  ///
  /// For the purposes of initialization, we can simply fix the first pose to
  /// the origin; this corresponds to fixing the first d elements of t to 0,
  /// and
  /// slicing off the first d columns of B1 to form
  ///
  /// min || B1red * tred + c) ||, where
  ///
  /// c = B2 * vec(R)

  // Vectorization of R matrix
  Eigen::Map<Vector> rvec(const_cast<Scalar *>(R.data()), d * d * n);

  // Form the matrix comprised of the right (n-1) block columns of B1
  SparseMatrix B1red = B1.rightCols(d * (n - 1));

  Vector c = B2 * rvec;

  // Solve
  Eigen::SPQR<SparseMatrix> QR(B1red);
  Vector tred = -QR.solve(c);

  // Reshape this result into a d x (n-1) matrix
  Eigen::Map<Matrix> tred_mat(tred.data(), d, n - 1);

  // Allocate output matrix
  Matrix t = Matrix::Zero(d, n);

  // Set rightmost n-1 columns
  t.rightCols(n - 1) = tred_mat;

  return t;
}

Matrix project_to_SOd(const Matrix &M) {
  // Compute the SVD of M
  Eigen::JacobiSVD<Matrix> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);

  Scalar detU = svd.matrixU().determinant();
  Scalar detV = svd.matrixV().determinant();

  if (detU * detV > 0) {
    return svd.matrixU() * svd.matrixV().transpose();
  } else {
    Matrix Uprime = svd.matrixU();
    Uprime.col(Uprime.cols() - 1) *= -1;
    return Uprime * svd.matrixV().transpose();
  }
}

Scalar orbit_distance_dS(const Matrix &X, const Matrix &Y, Matrix *G_S) {
  size_t d = X.rows();
  size_t n = X.cols() / d;

  // Compute orbit distance and optimal registration G_S according to Theorem 5
  // in the SE-Sync tech report
  Matrix XYt = X * Y.transpose();

  // Compute singular value decomposition of XY^T
  Eigen::JacobiSVD<Matrix> svd(XYt, Eigen::ComputeFullU | Eigen::ComputeFullV);

  // Compute UVt
  Matrix UVt = svd.matrixU() * svd.matrixV().transpose();

  // Construct diagonal matrix Xi = diag(1, ..., 1, det(UV^T))
  Vector Xi_diag = Vector::Constant(d, 1.0);
  // Note that since U and V are orthogonal matrices, then det(UV^T) = +/- 1
  Xi_diag(d - 1) = std::copysign(1.0, UVt.determinant());

  // Compute orbit distance dS
  Scalar dS = sqrt(fabs(
      2 * d * n - 2 * (Xi_diag.array() * svd.singularValues().array()).sum()));

  if (G_S) {
    // Compute optimal registration G_S registering Y to X
    *G_S = svd.matrixU() * Xi_diag.asDiagonal() * svd.matrixV().transpose();
  }
  return dS;
}

Scalar orbit_distance_dO(const Matrix &X, const Matrix &Y, Matrix *G_O) {
  size_t d = X.rows();
  size_t n = X.cols() / d;

  // Compute orbit distance and optimal registration G_O according to Theorem 5
  // in the SE-Sync tech report
  Matrix XYt = X * Y.transpose();

  // Compute singular value decomposition of XY^T
  Eigen::JacobiSVD<Matrix> svd(XYt, Eigen::ComputeFullU | Eigen::ComputeFullV);

  // Compute orbit distance dO
  Scalar dO = sqrt(fabs(2 * d * n - 2 * svd.singularValues().sum()));

  if (G_O) {
    // Compute optimal registration G_O registering Y to X
    *G_O = svd.matrixU() * svd.matrixV().transpose();
  }
  return dO;
}
} // namespace SESync
