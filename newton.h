#ifndef NEWTON_H_
#define NEWTON_H_

#include <Eigen/Dense>

class GMRES;
class Linear_Solver;
class RHS_Fn;
class Jacobian_Fn;

class Newton {
 public:
  Newton(const double tol_abs, const double tol_rel, const int itermax);
  ~Newton() {}
  Eigen::VectorXd find_zero(const RHS_Fn& F, const Jacobian_Fn& DF, const Eigen::VectorXd& x0, const Linear_Solver& ls) const;
  Eigen::VectorXd find_zero(const RHS_Fn& F, const Eigen::VectorXd& x0, const double dx, const GMRES& ls) const;
 private:
  const double tol_abs_;
  const double tol_rel_;
  const int itermax_;
};

#endif
    
    
    
