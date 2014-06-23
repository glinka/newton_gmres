#ifndef NEWTON_H_
#define NEWTON_H_

class GMRES;
class Linear_Solver;
/* class Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>; */
/* class Eigen::Matrix<double, Eigen::Dynamic, 1>; */

#include <Eigen/Dense>
/* namespace Eigen { */
/*   template<typename _Scalar, int _Rows, int _Cols> class Matrix; */
/*   template<> class Matrix<double, -1, -1>; */
/*   template<> class Matrix<double, -1, -1>; */
/*   typedef class Matrix<double, -1, -1> MatrixXd; */
/*   typedef class Matrix<double, -1, 1> VectorXd; */
/* } */

/* typedef class Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, int, int, int> MatrixXd; */
/* typedef class Eigen::Matrix<double, Eigen::Dynamic, 1, int, int, int> VectorXd; */

class Newton {
 public:
  Newton(const double tol_abs, const double tol_rel, const int itermax);
  ~Newton() {}
  Eigen::VectorXd find_zero(Eigen::VectorXd (*F)(const Eigen::VectorXd&), Eigen::MatrixXd (*DF)(const Eigen::VectorXd&), const Eigen::VectorXd& x0, const Linear_Solver& ls);
  Eigen::VectorXd find_zero(Eigen::VectorXd (*F)(const Eigen::VectorXd&), const Eigen::VectorXd& x0, const double dx, const GMRES& ls);
 private:
  const double tol_abs_;
  const double tol_rel_;
  const int itermax_;
};

#endif
    
    
    
