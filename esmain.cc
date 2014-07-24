#include "eigen_solvers.h"
#include <iostream>
#include <cmath>

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
  const int n = 20;
  Eigen::MatrixXd A = Eigen::MatrixXd::Random(n, n);
  A += A.transpose().eval();
  Eigen::VectorXd v = Eigen::VectorXd::Random(n);
  Eigen::MatrixXd V, H;
  Eigen::VectorXd f;
  const int m = 15;
  eigen_solver::arnoldi_iter(A, v, V, H, f, m);
  Eigen::MatrixXd Q, R;
  eigen_solver::qr_tridiag(H, Q, R);
  Eigen::MatrixXd S;
  Eigen::VectorXd thetas;
  const int qr_maxiter = 25;
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
  std::cout << "------------------------------" << std::endl;
  std::cout << "Imp-shift QR accuracy test:" << std::endl;
  std::cout << "residual in HV - VL: " << (H*S - S*(thetas.asDiagonal())).sum() << std::endl;
  std::cout << "------------------------------" << std::endl;
  return 1;
}
