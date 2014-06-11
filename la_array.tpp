#include "dims_exception.h"

namespace la {

  // eventually need to check whether N = dims.size()

  template <typename T, unsigned int N> Array<T, N>::Array(const std::vector<int> dims): Array(dims.begin()) {
    Dims_Exception::test_dims(dims, N);
  }

  template <typename T, unsigned int N> Array<T, N>::Array(std::vector<int>::const_iterator dims): dims_(std::vector<int>(dims, dims+N)) {
    const int n = *dims;
    dims++;
    A = new Array<T, N-1>*[n];
    for(int i = 0; i < n; i++) {
      A[i] = new Array<T, N-1>(dims);
    }
  }

  template <typename T, unsigned int N> Array<T, N>::Array(std::vector<int> dims, const T init_val): Array(dims.begin(), init_val)  {
    Dims_Exception::test_dims(dims, N);
  }

  template <typename T, unsigned int N> Array<T, N>::Array(std::vector<int>::const_iterator dims, const T init_val): dims_(std::vector<int>(dims, dims+N)) {
    const int n = *dims;
    dims++;
    A = new Array<T, N-1>*[n];
    for(int i = 0; i < n; i++) {
      A[i] = new Array<T, N-1>(dims, init_val);
    }
  }

  template <typename T, unsigned int N>  Array<T, N-1>& Array<T, N>::operator[](unsigned int i) {
    return *A[i];
  }

  template <typename T, unsigned int N>  const Array<T, N-1>& Array<T, N>::operator[](unsigned int i) const {
    return *A[i];
  }

  template <typename T> Array<T, 1>::Array(std::vector<int>::const_iterator dims): size_(*dims) {
    V = new T[size_];
  }

  template <typename T> Array<T, 1>::Array(const std::vector<int> dims): size_(dims[0]) {
    Dims_Exception::test_dims(dims, 1);
    V = new T[size_];
  }

  template <typename T> Array<T, 1>::Array(std::vector<int>::const_iterator dims, const T init_val): size_(*dims) {
    V = new T[size_];
    for(int i = 0; i < size_; i++) {
      V[i] = init_val;
    }
  }

  template <typename T> Array<T, 1>::Array(const std::vector<int> dims, const T init_val): size_(dims[0]) {
    Dims_Exception::test_dims(dims, 1);
    V = new T[size_];
    for(int i = 0; i < size_; i++) {
      V[i] = init_val;
    }
  }

  template <typename T> Array<T, 1>::Array(const unsigned int size): size_(size) {
    V = new T[size_];
  }
  
  template <typename T> Array<T, 1>::Array(const unsigned int size, const T init_val): size_(size) {
    V = new T[size_];
    for(int i = 0; i < size_; i++) {
      V[i] = init_val;
    }
  }

  template <typename T>  T& Array<T, 1>::operator[](unsigned int i) {
    return V[i];
  }

  template <typename T>  const T& Array<T, 1>::operator[](unsigned int i) const {
    return V[i];
  }

}
