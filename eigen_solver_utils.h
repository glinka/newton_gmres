#ifndef EIGEN_SOLVER_UTILS_H_
#define EIGEN_SOLVER_UTILS_H_

#include <cmath>
#include <complex>
#include <utility>
#include <vector>
#include <algorithm>

namespace eigen_solver_utils {

  double givens_c(const double diag, const double subdiag) {
    return diag/sqrt(diag*diag + subdiag*subdiag);
  }

  double givens_s(const double diag, const double subdiag) {
    return -subdiag/sqrt(diag*diag + subdiag*subdiag);
  }

  std::complex<double> wilk_shift(const double a, const double b, const double c) {
    std::complex<double> num(copysign(1, (a - c)/2.0)*b*b, 0);
    std::complex<double> denom_i(abs((a - c)/2.0), 0);
    std::complex<double> d(c, 0);
    std::complex<double> denom_ii(pow((a - c)/2.0, 2) + b*b, 0);
    return d - num/(denom_i + std::sqrt(denom_ii));
  }

  bool pair_comp(const std::pair<double, int>& a, const std::pair<double, int>& b) {
    return (a.first < b.first);
  }

  std::vector<int> argsort(const std::vector< std::pair<double, int> >& to_sort) {
    std::vector< std::pair<double, int> > pairs = to_sort;
    std::sort(pairs.begin(), pairs.end(), pair_comp);
    std::vector<int> sorted_indices(to_sort.size());
    for(int i = 0; i < to_sort.size(); i++) {
      sorted_indices[i] = pairs[i].second;
    }
    return sorted_indices;
  }
    

}

#endif
