#ifndef EIGEN_SOLVER_UTILS_H_
#define EIGEN_SOLVER_UTILS_H_

#include <cmath>
#include <complex>

namespace eigen_solver_utils {

  double givens_c(const double diag, const double subdiag) {
    return diag/sqrt(diag*diag + subdiag*subdiag);
  }
  double givens_s(const double diag, const double subdiag) {
    return -subdiag/sqrt(diag*diag + subdiag*subdiag);
  }
  std::complex<double> wilk_shift(const double a, const double b, const double c) {
    return c - copysign(1, (a - c)/2.0)*b*b/(abs((a - c)/2.0) + sqrt(pow((a - c)/2.0, 2) + b*b));
  }

}

#endif
