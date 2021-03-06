#ifndef GMRES_H
#define GMRES_H
#include "linear_solver.h"

class RHS_Fn;

class GMRES : public Linear_Solver {
 public:
  GMRES(const double tol, const int kmax);
  ~GMRES(){}
  Eigen::VectorXd solve_linear_system(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Eigen::VectorXd& x0) const;
  Eigen::VectorXd solve_linear_system(const RHS_Fn& F, const Eigen::VectorXd& x, const Eigen::VectorXd& x0, const double dx) const;
 private: 
  const double tol_;
  const int kmax_;
};

#endif
