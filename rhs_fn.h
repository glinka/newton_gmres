#ifndef RHS_FN
#define RHS_FN

#include <Eigen/Dense>

class RHS_Fn {
 public:
  virtual Eigen::VectorXd operator()(const Eigen::VectorXd& x) const = 0;
};

#endif
