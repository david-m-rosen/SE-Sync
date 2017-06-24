#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>

#include <Eigen/Geometry>
#include <Eigen/SPQRSupport>

#include <MatOp/SparseSymMatProd.h> // Spectra's built-in class for handling Eigen's sparse symmetric matrix type
#include <SymEigsSolver.h> // Spectra's symmetric eigensolver

#include "SESync_utils.h"
#include "StieVariable.h"

namespace SESync {

std::vector<SESync::RelativePoseMeasurement>
read_g2o_file(const std::string &filename, size_t &num_poses) {

  // Preallocate output vector
  std::vector<SESync::RelativePoseMeasurement> measurements;

  // A single measurement, whose values we will fill in
  SESync::RelativePoseMeasurement measurement;

  // A string used to contain the contents of a single line
  std::string line;

  // A string used to extract tokens from each line one-by-one
  std::string token;

  // Preallocate various useful quantities
  double dx, dy, dz, dtheta, dqx, dqy, dqz, dqw, I11, I12, I13, I14, I15, I16,
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
      measurement.t = Eigen::Vector2d(dx, dy);
      measurement.R = Eigen::Rotation2Dd(dtheta).toRotationMatrix();

      Eigen::Matrix2d TranCov;
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
      measurement.t = Eigen::Vector3d(dx, dy, dz);
      measurement.R = Eigen::Quaterniond(dqw, dqx, dqy, dqz).toRotationMatrix();

      // Compute precisions

      // Compute and store the optimal (information-divergence-minimizing) value
      // of the parameter tau
      Eigen::Matrix3d TranCov;
      TranCov << I11, I12, I13, I12, I22, I23, I13, I23, I33;
      measurement.tau = 3 / TranCov.inverse().trace();

      // Compute and store the optimal (information-divergence-minimizing value
      // of the parameter kappa

      Eigen::Matrix3d RotCov;
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
    size_t max_pair = std::max<double>(measurement.i, measurement.j);

    num_poses = ((max_pair > num_poses) ? max_pair : num_poses);
    measurements.push_back(measurement);
  } // while

  infile.close();

  num_poses++; // Account for the use of zero-based indexing

  return measurements;
}

SparseMatrix construct_rotational_connection_Laplacian(
    const std::vector<SESync::RelativePoseMeasurement> &measurements) {

  size_t num_poses = 0; // We will use this to keep track of the largest pose
  // index encountered, which in turn provides the number
  // of poses

  size_t d = (!measurements.empty() ? measurements[0].t.size() : 0);

  // Each measurement contributes 2*d elements along the diagonal of the
  // connection Laplacian, and 2*d^2 elements on a pair of symmetric
  // off-diagonal blocks

  size_t measurement_stride = 2 * (d + d * d);

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(measurement_stride * measurements.size());

  size_t i, j, max_pair;
  for (const SESync::RelativePoseMeasurement &measurement : measurements) {
    i = measurement.i;
    j = measurement.j;

    // Elements of ith block-diagonal
    for (unsigned int k = 0; k < d; k++)
      triplets.emplace_back(d * i + k, d * i + k, measurement.kappa);

    // Elements of jth block-diagonal
    for (unsigned int k = 0; k < d; k++)
      triplets.emplace_back(d * j + k, d * j + k, measurement.kappa);

    // Elements of ij block
    for (unsigned r = 0; r < d; r++)
      for (unsigned int c = 0; c < d; c++)
        triplets.emplace_back(i * d + r, j * d + c,
                              -measurement.kappa * measurement.R(r, c));

    // Elements of ji block
    for (unsigned int r = 0; r < d; r++)
      for (unsigned int c = 0; c < d; c++)
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
  std::vector<Eigen::Triplet<double>> triplets;
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

  std::vector<Eigen::Triplet<double>> triplets;
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

void construct_B_matrices(
    const std::vector<RelativePoseMeasurement> &measurements, SparseMatrix &B1,
    SparseMatrix &B2, SparseMatrix &B3) {
  // Clear input matrices
  B1.setZero();
  B2.setZero();
  B3.setZero();

  size_t num_poses = 0;
  size_t d = (!measurements.empty() ? measurements[0].t.size() : 0);

  std::vector<Eigen::Triplet<double>> triplets;

  // Useful quantities to cache
  unsigned int d2 = d * d;
  unsigned int d3 = d * d * d;

  unsigned int i, j; // Indices for the tail and head of the given measurement
  double sqrttau, sqrtkappa;
  size_t max_pair; // Used for keeping track of the maximum pose that we've
                   // encountered so far

  /// Construct the matrix B1 from equation (69a) in the tech report
  triplets.reserve(2 * d * measurements.size());

  for (unsigned int e = 0; e < measurements.size(); e++) {
    i = measurements[e].i;
    j = measurements[e].j;
    sqrttau = sqrt(measurements[e].tau);

    // Block corresponding to the tail of the measurement
    for (unsigned int l = 0; l < d; l++) {
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

  for (unsigned int e = 0; e < measurements.size(); e++) {
    i = measurements[e].i;
    for (unsigned int k = 0; k < d; k++)
      for (unsigned int r = 0; r < d; r++)
        triplets.emplace_back(d * e + r, d2 * i + d * k + r,
                              -sqrttau * measurements[e].t(k));
  }

  B2.resize(d * measurements.size(), d2 * num_poses);
  B2.setFromTriplets(triplets.begin(), triplets.end());

  /// Construct matrix B3 from equation (69c) in the tech report
  triplets.clear();
  triplets.reserve((d3 + d2) * measurements.size());

  for (unsigned int e = 0; e < measurements.size(); e++) {
    double sqrtkappa = sqrt(measurements[e].kappa);
    const Eigen::MatrixXd &R = measurements[e].R;

    for (unsigned int r = 0; r < d; r++)
      for (unsigned int c = 0; c < d; c++) {
        i = measurements[e].i; // Tail of measurement
        j = measurements[e].j; // Head of measurement

        // Representation of the -sqrt(kappa) * Rt(i,j) \otimes I_d block
        for (unsigned int l = 0; l < d; l++)
          triplets.emplace_back(e * d2 + d * r + l, i * d2 + d * c + l,
                                -sqrtkappa * R(c, r));
      }

    for (unsigned l = 0; l < d2; l++)
      triplets.emplace_back(e * d2 + l, j * d2 + l, sqrtkappa);
  }

  B3.resize(d2 * measurements.size(), d2 * num_poses);
  B3.setFromTriplets(triplets.begin(), triplets.end());
}

Matrix
chordal_initialization_eig(const SparseMatrix &rotational_connection_Laplacian,
                           unsigned int d, unsigned int max_iterations,
                           double precision) {
  Spectra::SparseSymMatProd<double> op(rotational_connection_Laplacian);
  Spectra::SymEigsSolver<double, Spectra::SELECT_EIGENVALUE::SMALLEST_ALGE,
                         Spectra::SparseSymMatProd<double>>
  eigensolver(&op, d, std::min<unsigned int>(
                          10 * d, rotational_connection_Laplacian.rows()));
  eigensolver.init();
  eigensolver.compute(max_iterations, precision,
                      Spectra::SELECT_EIGENVALUE::SMALLEST_ALGE);

  // Reproject the blocks of this matrix onto SO(d)
  return round_solution(eigensolver.eigenvectors().transpose(), d);
}

Matrix chordal_initialization(unsigned int d, const SparseMatrix &B3) {
  unsigned int d2 = d * d;
  unsigned int num_poses = B3.cols() / d2;

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
  B3red
      .makeCompressed(); // Must be in compressed format to use Eigen::SparseQR!

  // Vectorization of I_d
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(d, d);
  Eigen::Map<Eigen::VectorXd> Id_vec(Id.data(), d2);

  Eigen::VectorXd cR = B3.leftCols(d2) * Id_vec;

  Eigen::VectorXd rvec;
  Eigen::SPQR<SparseMatrix> QR(B3red);
  rvec = -QR.solve(cR);

  Eigen::MatrixXd Rchordal(d, d * num_poses);
  Rchordal.leftCols(d) = Id;
  Rchordal.rightCols((num_poses - 1) * d) =
      Eigen::Map<Eigen::MatrixXd>(rvec.data(), d, (num_poses - 1) * d);

  for (unsigned int i = 1; i < num_poses; i++)
    Rchordal.block(0, i * d, d, d) =
        project_to_SOd(Rchordal.block(0, i * d, d, d));
  return Rchordal;
}

Matrix recover_translations(const SparseMatrix &B1, const SparseMatrix &B2,
                            const Matrix &R) {
  unsigned int d = R.rows();
  unsigned int n = R.cols() / d;

  /// We want to find a minimizer of
  /// || B1 * t + B2 * vec(R) ||
  ///
  /// For the purposes of initialization, we can simply fix the first pose to
  /// the origin; this corresponds to fixing the first d elements of t to 0, and
  /// slicing off the first d columns of B1 to form
  ///
  /// min || B1red * tred + c) ||, where
  ///
  /// c = B2 * vec(R)

  // Vectorization of R matrix
  Eigen::Map<Eigen::VectorXd> rvec((double *)R.data(), d * d * n);

  // Form the matrix comprised of the right (n-1) block columns of B1
  SparseMatrix B1red = B1.rightCols(d * (n - 1));

  Eigen::VectorXd c = B2 * rvec;

  // Solve
  Eigen::SPQR<SparseMatrix> QR(B1red);
  Eigen::VectorXd tred = -QR.solve(c);

  // Reshape this result into a d x (n-1) matrix
  Eigen::Map<Eigen::MatrixXd> tred_mat(tred.data(), d, n - 1);

  // Allocate output matrix
  Eigen::MatrixXd t = Eigen::MatrixXd::Zero(d, n);

  // Set rightmost n-1 columns
  t.rightCols(n - 1) = tred_mat;

  return t;
}

Matrix project_to_SOd(const Matrix &M) {
  // Compute the SVD of M
  Eigen::JacobiSVD<Matrix> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);

  double detU = svd.matrixU().determinant();
  double detV = svd.matrixV().determinant();

  if (detU * detV > 0) {
    return svd.matrixU() * svd.matrixV().transpose();
  } else {
    Eigen::MatrixXd Uprime = svd.matrixU();
    Uprime.col(Uprime.cols() - 1) *= -1;
    return Uprime * svd.matrixV().transpose();
  }
}

Matrix round_solution(const Matrix &Y, unsigned int d) {
  // First, compute a thin SVD of Y
  Eigen::JacobiSVD<Matrix> svd(Y, Eigen::ComputeThinV);

  Eigen::VectorXd sigmas = svd.singularValues();
  // Construct a diagonal matrix comprised of the first d singular values
  Eigen::DiagonalMatrix<double, Eigen::Dynamic> Xid(d);
  Eigen::DiagonalMatrix<double, Eigen::Dynamic>::DiagonalVectorType &diagonal =
      Xid.diagonal();
  for (unsigned int i = 0; i < d; i++)
    diagonal(i) = sigmas(i);

  Eigen::MatrixXd R = Xid * svd.matrixV().leftCols(d).transpose();

  unsigned int n = Y.cols() / d;

  Eigen::VectorXd determinants(n);

  unsigned int ng0 = 0; // This will count the number of blocks whose
  // determinants have positive sign
  for (unsigned int i = 0; i < n; i++) {
    // Compute the determinant of the ith dxd block of R
    determinants(i) = R.block(0, i * d, d, d).determinant();
    if (determinants(i) > 0)
      ng0++;
  }

  if (ng0 < n / 2) {
    // Less than half of the total number of blocks have the correct sign, so
    // reverse their orientations

    // Get a reflection matrix that we can use to reverse the signs of those
    // blocks of R that have the wrong determinant
    Eigen::MatrixXd reflector = Eigen::MatrixXd::Identity(d, d);
    reflector(d - 1, d - 1) = -1;

    R = reflector * R;
  }

  // Finally, project each dxd block of R to SO(d)
  for (unsigned int i = 0; i < n; i++)
    R.block(0, i * d, d, d) = project_to_SOd(R.block(0, i * d, d, d));

  return R;
}

void StiefelProd2Mat(const ROPTLIB::ProductElement &product_element,
                     Eigen::MatrixXd &Y) {

  const int *sizes = product_element.GetElement(0)->Getsize();
  unsigned int r = sizes[0];
  unsigned int d = sizes[1];
  unsigned int n = product_element.GetNumofElement();

  for (unsigned int j = 0; j < n; j++) {
    // Set the ith block of the output matrix
    Y.block(0, j * d, r, d) = Eigen::Map<Eigen::MatrixXd>(
        (double *)product_element.GetElement(j)->ObtainReadData(), r, d);
  }
}

void Mat2StiefelProd(const Eigen::MatrixXd &Y,
                     ROPTLIB::ProductElement &product_element) {
  const int *sizes = product_element.GetElement(0)->Getsize();
  unsigned int r = sizes[0];
  unsigned int d = sizes[1];
  unsigned int n = product_element.GetNumofElement();
  unsigned int data_stride = r * d;

  const double *matrix_data = Y.data();

  for (unsigned int j = 0; j < n; j++) {
    double *element_data =
        product_element.GetElement(j)->ObtainWriteEntireData();
    memcpy(element_data, matrix_data + j * data_stride,
           sizeof(double) * data_stride);
  }
}
};
