#include "test_fns.h"

#include <Eigen/Dense>

/*
  f1 = x^2 - 1
  f2 = y^2 - 1
*/

Eigen::VectorXd F(const Eigen::VectorXd& x) {
  return x.array().square() - 1;
}

Eigen::MatrixXd DF(const Eigen::VectorXd& x) {
  return Eigen::MatrixXd((2*x).asDiagonal());
}

