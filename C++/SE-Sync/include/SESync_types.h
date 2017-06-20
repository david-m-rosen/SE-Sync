#ifndef SESYNC_TYPES_H
#define SESYNC_TYPES_H

#include <Eigen/Dense>
#include <Eigen/Sparse>

namespace SESync {

/** Some useful typedefs for the SE-Sync library */

typedef Eigen::VectorXd Vector;
typedef Eigen::MatrixXd Matrix;
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic> DiagonalMatrix;

/** We use row-major storage order to take advantage of fast (sparse-matrix) * (dense-vector) multiplications when OpenMP is available */
typedef Eigen::SparseMatrix<double, Eigen::RowMajor> SparseMatrix;

} // namespace SESync

#endif // SESYNC_TYPES_H
