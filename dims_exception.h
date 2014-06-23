#include <exception>
#include "la_array.h"

class Dims_Exception : public std::exception {

  virtual const char* what() const throw() {
    return "dims.size() != N: number of dimensions mismatched";
 }
  
 public:
  static void test_dims(const std::vector<int> dims, const unsigned int N) {
    if(dims.size() != N) {
      throw Dims_Exception();
    }
  }

  template <typename T, unsigned int N, unsigned int M> static void test_dims(const la::Array<T, N>& A, const la::Array<T, M>& B) {
    int a = A.get_dims().back();
    int b = B.get_dims().back();
    /* if(A.get_dims().back() != B.get_dims().front()) { */
    if(a != b) {
      throw Dims_Exception();
    }
  }
};

#include <exception>
#include "la_array.h"

class Dims_Exception : public std::exception {

  virtual const char* what() const throw() {
    return "dims.size() != N: number of dimensions mismatched";
 }
  
 public:
  static void test_dims(const std::vector<int> dims, const unsigned int N) {
    if(dims.size() != N) {
      throw Dims_Exception();
    }
  }

  template <typename T, unsigned int N, unsigned int M> static void test_dims(const la::Array<T, N>& A, const la::Array<T, M>& B) {
    int a = A.get_dims().back();
    int b = B.get_dims().back();
    /* if(A.get_dims().back() != B.get_dims().front()) { */
    if(a != b) {
      throw Dims_Exception();
    }
  }
};
