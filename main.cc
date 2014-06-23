#include <iostream>
#include "gmres.h"
#include "newton.h"
#include "test_fns.h"

#include <Eigen/Dense>

int main(int argc, char* argv[]) {
  Eigen::VectorXd x0 = Eigen::VectorXd::Ones(2);
  x0(0) = 2;
  Eigen::VectorXd gmres_x0 = Eigen::VectorXd::Zero(2);
  Newton newton(1e-8, 1e-9, 100);
  GMRES gmres(1e-8, 2);
  std::cout << gmres.solve_linear_system(DF(x0), F(x0), x0) << std::endl;
  std::cout << F(x0) << std::endl;
  std::cout << DF(x0) << std::endl;
  std::cout << newton.find_zero(F, DF, x0, gmres) << std::endl;
}
  
