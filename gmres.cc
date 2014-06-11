#include <cmath>
#include "la_ops.h"
#include "gmres.h"

GMRES::GMRES(const double tol, const double kmax): tol_(tol), kmax_(kmax) {}

Eigen::VectorXd& GMRES::solve_linear_system(const Eigen::MatrixXd& A, const Eigen::VectorXd& b, const Eigen::VectorXd& x0) {
  const int n = b.get_size();
  Eigen::VectorXd r = b - A*x0;
  std::vector<int> V_max_dims(2);
  V_max_dims[0] = n;
  V_max_dims[1] = kmax_ + 1;
  std::vector<int> H_max_dims(2);
  H_max_dims[0] = kmax_ + 1;
  H_max_dims[1] = kmax_;
  /*
    V will store the basis vectors for
    the Krylov subspace as row vectors
  */
  Eigen::MatrixXd V(V_max_dims, 0);
  la::matrix H(H_max_dims, 0);
  // V[0] is really x0, and not part of the Krylov subspace
  V[0] = r/l2_norm(r);
  double rho = l2_norm(r);
  const double beta = rho;
  const double b_norm = l2_norm(b);
  int k = -1;
  la::vector g(kmax+1);
  g[0] = rho;
  double c[kmax_];
  double s[kmax_];
  la::vector h(kmax_ + 1);
  while(rho > tol_*b_norm && k < kmax-1) {
    k++;
    V[k+1] = A*V[k];
    for(int j = 0, j < k+1; j++) {
      h[j] = V[k+1]*V[k];
      V[k+1] = V[k+1] - h[j]*V[j];
    }
    h[k+1] = l2_norm(V[k+1]);
    // test for orthogonality
    V[k+1] = V[k+1]/h[k];
    // perform first loop separately to avoid if?
    // will not execute during first loop, as desired
    for(int i = 0; i < k; i++) {
      double temp1 = c[i]*h[i] - s[i+1]*h[i+1];
      double temp2 = s[i]*h[i] + c[i+1]*h[i+1];
      h[i] = temp1;
      h[i+1] = temp2;
    }
    double nu = sqrt(h[k]*h[k] + h[k+1]*h[k+1]);
    c[k] = h[k]/nu;
    s[k] = -h[k+1]/nu;
    h[k] = c[k]*h[k] - s[k]*h[k+1];
    h[k+1] = 0;
    g = calc_G(k, c[k], s[k])*g;
    for(int i = 0; i < k; i++) {
      H[i, k] = h[i];p
    }
    rho = abs(g[k+1]);
  }
  la::vector y(k+2);
  y[k+1] = 1;
  // solve for y^k
  for(int i = k; i >= 0; i--) {
    double sum = 0;
    for(int j = k; j > i; j--) {
      sum += y[j]*H[i][j];
    }
    y[i] = (g[i] - sum)/H[i][i];
  }
  la::vector x(n);
  for(int i = 0; i < n; i++) {
    double sum = 0;
    for(int j = 0; j < k+2; j++) {
      // because V is transposed
      sum += V[j][i]*y[j];
    }
    x[i] = x0[i] + sum;
  }

  // cleanup, can move to one line?
  delete r;
  delete V;
  delete H;
  delete g;
  delete h;
  delete y;

  return x;
}

#include<iostream>

int main(int argc, char **argv) {
  const int n = 3;
  la::matrix A(std::vector<int>(2, n), 0);
  la::vector b(n);
  la::vector x0(n, 1);
  for(int i = 0; i < n; i++) {
    A[i][i] = 1;
    b[i] = i;
  }
  GMRES solver = GMRES(1e-4, 3);
  la::vector x = solver.solve_linear_system(A, b, x0);
  for(int i = 0; i < n; i++) {
    std::cout << x[i] << std::endl;
  }
  return 0;
}
