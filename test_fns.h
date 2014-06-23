#ifndef TEST_FNS_H_
#define TEST_FNS_H_

#include <Eigen/Dense>

/* class Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>; */
/* class Eigen::Matrix<double, Eigen::Dynamic, 1>; */

/* typedef class Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, int, int, int> MatrixXd; */
/* typedef class Eigen::Matrix<double, Eigen::Dynamic, 1, int, int, int> VectorXd; */

/* namespace Eigen { */
/*   template<typename _Scalar, int _Rows, int _Cols> class Matrix; */
/*   template<> class Matrix<double, -1, -1>; */
/*   template<> class Matrix<double, -1, -1>; */
/*   typedef class Matrix<double, -1, -1> MatrixXd; */
/*   typedef class Matrix<double, -1, 1> VectorXd; */
/*   /\* template<> class Matrix<double, -1, -1, -1, -1, -1>; *\/ */
/*   /\* template<> class Matrix<double, -1, -1, -1, -1, -1>; *\/ */
/*   /\* typedef class Matrix<double, -1, -1, -1, -1, -1> MatrixXd; *\/ */
/*   /\* typedef class Matrix<double, -1, 1, -1, -1, -1> VectorXd; *\/ */
/* } */

Eigen::VectorXd F(const Eigen::VectorXd& x);

Eigen::MatrixXd DF(const Eigen::VectorXd& x);
    
#endif
