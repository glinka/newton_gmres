#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H
// #include "common_header.h"

typedef class Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXd;
typedef class Eigen::Matrix<double, Eigen::Dynamic, 1> VectorXd;

class Linear_Solver {
 public:
  virtual Eigen::VectorXd solve_linear_system(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Eigen::VectorXd& x0) = 0;
  Linear_Solver(){}
  virtual ~Linear_Solver(){}
};

#endif










