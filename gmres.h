#ifndef GMRES_H
#define GMRES_H
#include "linear_solver.h"

class GMRES : public Linear_Solver {
 public:
  GMRES(const double tol, const double kmax);
  ~GMRES(){}
  Eigen::VectorXd solve_linear_system(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Eigen::VectorXd& x0);
 private: 
  const double tol_;
  const int kmax_;
};

#endif
