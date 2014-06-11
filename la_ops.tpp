#include <cmath>

namespace la {
  // template<typename T, unsigned int N, unsigned int M> static Array<T, N+M-2>& operator*(const Array<T, N>& A, const Array<T, M>& B) {
  //   Dims_Exception::test_dims(A, B);
  //   const std::vector<int> A_dims = A.get_dims();
  //   const std::vector<int> B_dims = B.get_dims();
  //   std::vector<int> new_dims(N+M-2);
  //   std::vector<int>::iterator new_dims_it = new_dims.begin();
  //   for(std::vector<int>::const_iterator val = A_dims().begin(); val != A_dims().end()-1; val++) {
  //     *new_dims_it = *val;
  //     new_dims_it++;
  //   }
  //   for(std::vector<int>::const_iterator val = B_dims().begin()+1; val != B_dims().end(); val++) {
  //     *new_dims_it = *val;
  //     new_dims_it++;
  //   }
  //   Array<T, N+M-2>* C = new Array<T, N+M-2>(new_dims);
  //   int dim_counter = N+M-2;
  //   // for(int i = A_dims
  // }
  template <typename T> static Array<T, 2>& operator*(const Array<T, 2>& A, const Array<T, 2>& B) {
    Dims_Exception::test_dims(A, B);
    const std::vector<int> A_dims = A.get_dims();
    const std::vector<int> B_dims = B.get_dims();
    std::vector<int> new_dims(2);
    new_dims[0] = A_dims.front();
    new_dims[1] = B_dims.back();
    Array<T, 2>* C = new Array<T, 2>(new_dims, 0);
    for(int i = 0; i < A_dims[0]; i++) {
      for(int j = 0; j < B_dims[0]; j++) {
	for(int k = 0; k < A_dims[1]; k++) {
	  (*C)[i][j] += A[i][k]*B[k][j];
	}
      }
    }
    return *C;
  }

  template <typename T> static Array<T, 2>& dot(const Array<T, 2>& A, const Array<T, 2>& B) {
    return A*B;
  }

  template <typename T> static Array<T, 1>& operator*(const Array<T, 2>& A, const Array<T, 1>& x) {
    Dims_Exception::test_dims(A, x);
    const std::vector<int> A_dims = A.get_dims();
    const int new_size = A_dims.front();
    Array<T, 1>* b = new Array<T, 1>(new_size, 0);
    for(int i = 0; i < A_dims[0]; i++) {
      (*b)[i] = 0;
      for(int j = 0; j < A_dims[1]; j++) {
	(*b)[i] += A[i][j]*x[j];
      }
    }
    return *b;
  }

  template <typename T> static Array<T, 1>& dot(const Array<T, 2>& A, const Array<T, 1>& x) {
    return A*x;
  }
	  
  template <typename T> static double operator*(const Array<T, 1>& x, const Array<T, 1>& y) {
    Dims_Exception::test_dims(x, y);
    const int n = x.get_dims()[0];
    double dot_prod = 0;
    for(int i = 0; i < n; i++) {
      dot_prod += x[i]*y[i];
    }
    return dot_prod;
  }

  template <typename T> static double dot(const Array<T, 1>& x, const Array<T, 1>& y) {
    return x*y;
  }

  template <typename T> static double l2_norm(const Array<T, 1>&x) {
    return sqrt(dot(x, x));
  }

  template <typename T, unsigned int N> static double operator*(const Array<T, N>& x, const double c) {
    
    const int n = x.get_dims()[0];
    double dot_prod = 0;
    for(int i = 0; i < n; i++) {
      dot_prod += x[i]*y[i];
    }
    return dot_prod;
  }

}	  

