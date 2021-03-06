#include "eigen_solvers.h"
#include <iostream>
#include <cmath>
#include <cstdlib>

double below_diag_sum(const Eigen::MatrixXd& A) {
  double below_diag_sum = 0;
  for(int i = 0; i < A.rows(); i++) {
    for(int j = 0; j < i; j++) {
      below_diag_sum += A(i,j);
    }
  }
  return below_diag_sum;
}




void threshold(Eigen::MatrixXd& A) {
  const double zero_tol = 1e-15;
  for(int i = 0; i < A.rows(); i++) {
    for(int j = 0; j < A.cols(); j++) {
      if(abs(A(i,j)) < zero_tol) {
	A(i,j) = 0;
      }
    }
  }
}

int main(int argc, char** argv) {
  const int n = std::atoi(argv[1]);
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(n, n);
  A += A.transpose().eval();
  Eigen::VectorXd v = Eigen::VectorXd::Random(n);
  Eigen::MatrixXd V, H;
  Eigen::VectorXd f;
  const int m = std::atoi(argv[2]);
  eigen_solver::arnoldi_iter(A, v, V, H, f, m);
  Eigen::MatrixXd Q, R;
  eigen_solver::qr_tridiag(H, Q, R);
  Eigen::MatrixXd S;
  Eigen::VectorXd thetas;
  int qr_maxiter = 25;
  eigen_solver::qr_impshift_tridiag(H, S, thetas, qr_maxiter);

  std::cout << "------------------------------" << std::endl;
  std::cout << "V orthogonality test:" << std::endl;
  std::cout << "summation of off-diagonal elements in V.T*V: " << (V.transpose()*V).trace() - (V.transpose()*V).sum() << std::endl;
  std::cout << "diagonal elements in V.T*V: " << (V.transpose()*V).diagonal().transpose() << std::endl;
  std::cout << "------------------------------" << std::endl;
  std::cout << "Q orthogonality test:" << std::endl;
  std::cout << "summation of off-diagonal elements in Q.T*Q: " << (Q.transpose()*Q).trace() - (Q.transpose()*Q).sum() << std::endl;
  std::cout << "diagonal elements in Q.T*Q: " << (Q.transpose()*Q).diagonal().transpose() << std::endl;
  std::cout << "------------------------------" << std::endl;
  std::cout << "QR accuracy test:" << std::endl;
  std::cout << "residual in QR - H: " << (Q*R - H).sum() << std::endl;
  std::cout << "residual below diagonal of R: " << below_diag_sum(R) << std::endl;
  std::cout << "------------------------------" << std::endl;
  std::cout << "Imp-shift QR accuracy test:" << std::endl;
  std::cout << "residual in HV - VL: " << (H*S - S*(thetas.asDiagonal())).norm() << std::endl;
  std::cout << "------------------------------" << std::endl;

  Eigen::MatrixXd Vt1, Vt2, Ht1, Ht2;
  Eigen::VectorXd ft;
  eigen_solver::arnoldi_iter(A, v, Vt1, Ht1, ft, m);
  eigen_solver::arnoldi_iter(A, v, Vt2, Ht2, ft, 5);
  eigen_solver::arnoldi_iter(A, Eigen::MatrixXd(Vt2), Eigen::MatrixXd(Ht2), Eigen::VectorXd(ft), Vt2, Ht2, ft, m);
  std::cout << "Arnoldi iteration test:" << std::endl;
  std::cout << "residual in AV1 - AV2: " << (A*Vt1 - A*Vt2).norm() << std::endl;
  std::cout << "------------------------------" << std::endl;

  Eigen::MatrixXd V_r;
  Eigen::VectorXd l_r;
  const int k = 100, p=200, iram_maxiter=100;
  qr_maxiter=1000;
  eigen_solver::arnoldi_method_imprestart_hermitian(A, Eigen::VectorXd::Ones(n), V_r, l_r, k, p, iram_maxiter, qr_maxiter, 1e-11, 1e-11);
  std::cout << "Arnoldi method test:" << std::endl;
  std::cout << "residual in AV - VL: " << (A*V_r - V_r*(l_r.asDiagonal())).norm() << std::endl;
  std::cout << "------------------------------" << std::endl;
  

  return 1;
}
