#ifndef LINEAR_SOLVER_H
#define LINEAR_SOLVER_H

#include <Eigen/Dense>

class Linear_Solver {
 public:
  virtual Eigen::VectorXd solve_linear_system(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Eigen::VectorXd& x0) const = 0;
  Linear_Solver(){}
  virtual ~Linear_Solver(){}
};

#endif










