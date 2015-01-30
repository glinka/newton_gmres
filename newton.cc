#include <Eigen/Dense>
#include "newton.h"
#include "gmres.h"
#include "rhs_fn.h"
#include "jacobian_fn.h"
#include "iters_exception.h"

Newton::Newton(const double tol_abs, const double tol_rel, const int itermax): tol_abs_(tol_abs), tol_rel_(tol_rel), itermax_(itermax) {}

Eigen::VectorXd Newton::find_zero(const RHS_Fn& F, const Jacobian_Fn& DF, const Eigen::VectorXd& x0, const Linear_Solver& ls) const {
  const int n = x0.size();
  double r0 = F(x0).norm();
  double r = r0;
  Eigen::VectorXd zeros = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd x = x0;
  int iters = 0;
  while(r > r0*tol_rel_ + tol_abs_ && iters < itermax_) {
    // default to initial x = {0, 0, ..., 0}
    x -= ls.solve_linear_system(DF(x), F(x), zeros);
    r = F(x).norm();
    iters++;
  }
  // probably a better way to do this
  Iters_Exception::test_iters(iters, itermax_);
  return x;
}

Eigen::VectorXd Newton::find_zero(const RHS_Fn& F, const Eigen::VectorXd& x0, const double dx, const GMRES& ls) const {
  const int n = x0.size();
  double r0 = F(x0).norm();
  double r = r0;
  Eigen::VectorXd zeros = Eigen::VectorXd::Zero(n);
  Eigen::VectorXd x = x0;
  int iters = 0;
  while(r > r0*tol_rel_ + tol_abs_ && iters < itermax_) {
    // default to initial x = {0, 0, ..., 0}
    x += ls.solve_linear_system(F, x, zeros, dx);
    r = F(x).norm();
    iters++;
  }
  // probably a better way to do this
  Iters_Exception::test_iters(iters, itermax_);
  return x;
}
