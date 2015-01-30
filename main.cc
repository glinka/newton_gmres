#include <iostream>
#include "gmres.h"
#include "newton.h"
#include "rhs_fn.h"
#include "jacobian_fn.h"

class F : public RHS_Fn {
public:
  F(){}
  ~F(){}
  Eigen::VectorXd operator()(const Eigen::VectorXd& x) const {
    return x.array().square() - 1;
  }
};

class DF : public Jacobian_Fn {
public:
  DF(){}
  ~DF(){}
  Eigen::MatrixXd operator()(const Eigen::VectorXd& x) const {
    return Eigen::MatrixXd((2*x).asDiagonal());
  }
};

#include <Eigen/Dense>

int main(int argc, char* argv[]) {
  const int n = 10;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(10, 10);
  Eigen::VectorXd b = Eigen::VectorXd::Random(10);
  Eigen::VectorXd x0 = Eigen::VectorXd::Ones(10);
  Newton newton(1e-8, 1e-9, 100);
  GMRES gmres(1e-8, n);
  Eigen::VectorXd x_sln = gmres.solve_linear_system(A, b, x0);
  std::cout << "||Ax -b|| = " << (A*x_sln - b).norm() << std::endl;
  const double dx = 0.01;
  F f;
  x_sln = newton.find_zero(f, b, dx, gmres);
  std::cout << "||F(x*)|| = " << f(x_sln).norm() << std::endl;
}
