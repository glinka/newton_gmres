#ifndef JACOBIAN_FN
#define JACOBIAN_FN

#include <Eigen/Dense>

class Jacobian_Fn {
 public:
  virtual Eigen::MatrixXd operator()(const Eigen::VectorXd& x) const = 0;
};

#endif
