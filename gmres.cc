#include <cmath>
// #include "la_ops.h"
#include <Eigen/Dense>
#include "gmres.h"
#include "iters_exception.h"

GMRES::GMRES(const double tol, const int kmax): tol_(tol), kmax_(kmax) {}

Eigen::VectorXd GMRES::solve_linear_system(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Eigen::VectorXd& x0) const {
  const int n = b.size();
  Eigen::VectorXd r = b - A*x0;
  /*
    V will store the basis vectors for
    the Krylov subspace as column vectors
  */
  Eigen::MatrixXd V(n, kmax_ + 1);
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(kmax_ + 1, kmax_);
  // V[0] is really x0, and not part of the Krylov subspace
  V.col(0) = r/r.norm();
  double rho = r.norm();
  const double b_norm = b.norm();
  int k = -1;
  Eigen::VectorXd g = Eigen::VectorXd::Zero(kmax_ + 1);
  Eigen::VectorXd h = Eigen::VectorXd::Zero(kmax_ + 1);
  g(0) = rho;
  double c[kmax_];
  double s[kmax_];
  // Eigen::MatrixXd G = Eigen::MatrixXd::Identity(kmax_ + 1, kmax_ + 1);
  while(rho > tol_*b_norm && k < kmax_-1) {
    k++;
    V.col(k+1) = A*V.col(k);
    for(int j = 0; j < k+1; j++) {
      h(j) = V.col(k+1).dot(V.col(j));
      V.col(k+1) -= h(j)*V.col(j);
    }
    h(k+1) = V.col(k+1).norm();
    // test for orthogonality
    V.col(k+1) /= h(k+1);
    // perform first loop separately to avoid if?
    // will not execute during first loop, as desired
    for(int i = 0; i < k; i++) {
      double temp1 = c[i]*h(i) - s[i]*h(i+1);
      double temp2 = s[i]*h(i) + c[i]*h(i+1);
      h(i) = temp1;
      h(i+1) = temp2;
      // executed every time, wasteful
      // G(k-1, k-1) = 1;
      // G(k-1, k) = 0;
      // G(k, k-1) = 0;
    }
    double nu = std::sqrt(h(k)*h(k) + h(k+1)*h(k+1));
    c[k] = h(k)/nu;
    s[k] = -h(k+1)/nu;
    h(k) = c[k]*h(k) - s[k]*h(k+1);
    h(k+1) = 0;
    // set Givens matrix
    // G(k, k) = c[k];
    // G(k, k+1) = -s[k];
    // G(k+1, k) = s[k];
    // G(k+1, k+1) = c[k];
    // actually do this by hand
    double temp1 = c[k]*g(k) - s[k]*g(k+1);
    double temp2 = s[k]*g(k) + c[k]*g(k+1);
    g(k) = temp1;
    g(k+1) = temp2;
    H.block(0, k, k+1, 1) = h.head(k+1);
    rho = std::abs(g(k+1));
  }

  Eigen::VectorXd y(k+1);
  // solve for y^k
  for(int i = k; i >= 0; i--) {
    double sum = 0;
    for(int j = k; j > i; j--) {
      sum += y(j)*H(i,j);
    }
    y(i) = (g(i) - sum)/H(i,i);
  }
  return x0 + V.leftCols(k+1)*y;
}

namespace la {

  // for use when x, w \neq 0
  Eigen::VectorXd v_next(const Eigen::VectorXd& x, const Eigen::VectorXd& w, const double h, Eigen::VectorXd (*F)(const Eigen::VectorXd&)) {
    return w.norm()*(F(x + h*x.norm()*w/(w.norm())) - F(x))/(h*x.norm());
  }

  // for use when x = 0, w \neq 0
  Eigen::VectorXd v_next_zero(const Eigen::VectorXd& x, const Eigen::VectorXd& w, const double h, Eigen::VectorXd (*F)(const Eigen::VectorXd&)) {
    return w.norm()*(F(h*w/(w.norm())) - F(x))/h;
  }

}


Eigen::VectorXd GMRES::solve_linear_system(Eigen::VectorXd (*F)(const Eigen::VectorXd&), const Eigen::VectorXd& x, const Eigen::VectorXd& x0, const double dx) const {
  Eigen::VectorXd b = -F(x);
  const int n = b.size();
  /*
    V will store the basis vectors for
    the Krylov subspace as column vectors
  */
  Eigen::MatrixXd V(n, kmax_ + 1);
  Eigen::MatrixXd H = Eigen::MatrixXd::Zero(kmax_ + 1, kmax_);
  // V[0] is really x0, and not part of the Krylov subspace
  Eigen::VectorXd (*v_next)(const Eigen::VectorXd& x, const Eigen::VectorXd& w, const double h, Eigen::VectorXd (*F)(const Eigen::VectorXd&));
  Eigen::VectorXd r;
  // assume w.norm (= V.col(k).norm) \neq 0 in further iterations
  // otherwise we should have ended earlier
  if(x.norm() == 0) {
    v_next = la::v_next_zero;
  }
  else {
    v_next = la::v_next;
  }
  // however, x0.norm could still be zero, which
  // we account for below
  if (x0.norm() == 0) {
    r = b;
  }
  else {
    r = b - v_next(x, x0, dx, F);
  }    
  V.col(0) = r/r.norm();
  double rho = r.norm();
  const double b_norm = b.norm();
  int k = -1;
  Eigen::VectorXd g = Eigen::VectorXd::Zero(kmax_ + 1);
  Eigen::VectorXd h = Eigen::VectorXd::Zero(kmax_ + 1);
  g(0) = rho;
  double c[kmax_];
  double s[kmax_];
  // Eigen::MatrixXd G = Eigen::MatrixXd::Identity(kmax_ + 1, kmax_ + 1);
  // select appropriate method for determining v_{k+1}
  while(rho > tol_*b_norm && k < kmax_-1) {
    k++;
    V.col(k+1) = v_next(x, V.col(k), dx, F);
    for(int j = 0; j < k+1; j++) {
      h(j) = V.col(k+1).dot(V.col(j));
      V.col(k+1) -= h(j)*V.col(j);
    }
    h(k+1) = V.col(k+1).norm();
    // test for orthogonality
    V.col(k+1) /= h(k+1);
    // perform first loop separately to avoid if?
    // will not execute during first loop, as desired
    for(int i = 0; i < k; i++) {
      double temp1 = c[i]*h(i) - s[i]*h(i+1);
      double temp2 = s[i]*h(i) + c[i]*h(i+1);
      h(i) = temp1;
      h(i+1) = temp2;
      // executed every time, wasteful
      // G(k-1, k-1) = 1;
      // G(k-1, k) = 0;
      // G(k, k-1) = 0;
    }
    double nu = sqrt(h(k)*h(k) + h(k+1)*h(k+1));
    c[k] = h(k)/nu;
    s[k] = -h(k+1)/nu;
    h(k) = c[k]*h(k) - s[k]*h(k+1);
    h(k+1) = 0;
    // set Givens matrix
    // G(k, k) = c[k];
    // G(k, k+1) = -s[k];
    // G(k+1, k) = s[k];
    // G(k+1, k+1) = c[k];
    double temp1 = c[k]*g(k) - s[k]*g(k+1);
    double temp2 = s[k]*g(k) + c[k]*g(k+1);
    g(k) = temp1;
    g(k+1) = temp2;
    H.block(0, k, k+1, 1) = h.head(k+1);
    rho = std::abs(g(k+1));
  }

  Eigen::VectorXd y(k+1);
  // solve for y^k
  for(int i = k; i >= 0; i--) {
    double sum = 0;
    for(int j = k; j > i; j--) {
      sum += y(j)*H(i,j);
    }
    y(i) = (g(i) - sum)/H(i,i);
  }
  return x0 + V.leftCols(k+1)*y;

}
