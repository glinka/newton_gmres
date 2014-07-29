#ifndef EIGEN_SOLVER_H_
#define EIGEN_SOLVER_H_

#include <Eigen/Dense>
#include <vector>
#include <string>

namespace eigen_solver {

  int qr_tridiag(const Eigen::MatrixXd& T, Eigen::MatrixXd& Q, Eigen::MatrixXd& R);
  int arnoldi_iter(const Eigen::MatrixXd& A, const Eigen::MatrixXd& V_init, const Eigen::MatrixXd& H_init, const Eigen::VectorXd& f_init, Eigen::MatrixXd& V, Eigen::MatrixXd& H, Eigen::VectorXd& f, const int m, const double zero_tol=1e-14);
  int arnoldi_iter(const Eigen::MatrixXd& A, const Eigen::VectorXd& v, Eigen::MatrixXd& V, Eigen::MatrixXd& H, Eigen::VectorXd& f, const int m, const double zero_tol=1e-14);

  int qr_impshift_tridiag(const Eigen::MatrixXd& T, Eigen::MatrixXd& S, Eigen::VectorXd& thetas, const int maxiter, const double eigen_tol=1e-11, const double zero_tol=1e-13, const std::string shift_type="eigval")
;
  int arnoldi_method_imprestart_hermitian(const Eigen::MatrixXd& A, const Eigen::VectorXd& v, Eigen::MatrixXd& V_ritz, Eigen::VectorXd& l_ritz, const int k, const int p, const int iram_maxiter, const int qr_maxiter, const double f_tol=1e-14, const double qr_eigen_tol=1e-11, const double qr_zero_tol=1e-13);

}
    
#endif
