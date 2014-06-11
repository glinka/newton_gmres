#ifndef LA_OPS_H
#define LA_OPS_H
#include "la_array.h"

namespace la {
  /* template<typename T, unsigned int N, unsigned int M> static Array<T, N+M-2>& operator*(const Array<T, N>& A, const Array<T, M>& B); */
  template <typename T> static Array<T, 2>& operator*(const Array<T, 2>& A, const Array<T, 2>& B);
  template <typename T> static Array<T, 2>& dot(const Array<T, 2>& A, const Array<T, 2>& B);

  template <typename T> static Array<T, 1>& operator*(const Array<T, 2>& A, const Array<T, 1>& x);
  template <typename T> static Array<T, 1>& dot(const Array<T, 2>& A, const Array<T, 1>& x);

  template <typename T> static double operator*(const Array<T, 1>& x, const Array<T, 1>& y);
  template <typename T> static double dot(const Array<T, 1>& x, const Array<T, 1>& y);

  template <typename T> static double l2_norm(const Array<T, 1>&x);


  template <typename T, unsigned int N> static Array<T, N>& operator*(const Array<T, N>& x, const double c);
  template <typename T> static double operator*(const Array<T, 1>& x, const double c);
}

#include "la_ops.tpp"

#endif
