#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include <Eigen/Dense>

/* class Matrix<double, Eigen::Dynamic, Eigen::Dynamic>; */
/* class Matrix<double, Eigen::Dynamic, 1>; */

/* namespace Eigen { */
/*   template<typename _Scalar, int _Rows, int _Cols> class Matrix; */
/*   template<> class Matrix<double, -1, -1>; */
/*   template<> class Matrix<double, -1, -1>; */
/*   typedef class Matrix<double, -1, -1> MatrixXd; */
/*   typedef class Matrix<double, -1, 1> VectorXd; */
/*   /\\* template<typename T, int, int, int, int, int> class Matrix; *\/ */
/*   /\* template<> class Matrix<double, -1, -1, 4,4,4>; *\/ */
/*   /\* template<> class Matrix<double, -1, -1, 4,4,4>; *\/ */
/*   /\* typedef class Matrix<double, -1, -1, 4,4,4> MatrixXd; *\/ */
/*   /\* typedef class Matrix<double, -1, 1, 4,4,4> VectorXd; *\/ */
/*   /\* class MatrixXd; *\/ */
/*   /\* class VectorXd; *\/ */
/* } */

/* typedef class Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, int, int, int> MatrixXd; */
/* typedef class Eigen::Matrix<double, Eigen::Dynamic, 1, int, int, int> VectorXd; */

class Linear_Solver {
 public:
  virtual Eigen::VectorXd solve_linear_system(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Eigen::VectorXd& x0) const = 0;
  Linear_Solver(){}
  virtual ~Linear_Solver(){}
};

#endif










