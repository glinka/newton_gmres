#ifndef ITERS_EXCEPTION_H_
#define ITERS_EXCEPTION_H_

#include <exception>

class Iters_Exception : public std::exception {

  virtual const char* what() const throw() {
    return "error: maximum iterations exceeded";
 }

 public:
  static void test_iters(const int iters, const int maxiters) {
    if(iters >= maxiters) {
      throw Iters_Exception();
    }
  }
};

#endif
