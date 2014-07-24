#include "eigen_solvers.h"
#include "eigen_solver_utils.h"
#include <iostream>

namespace eigen_solver {

  int qr_tridiag(const Eigen::MatrixXd& T, Eigen::MatrixXd& Q, Eigen::MatrixXd& R) {
    const int n = T.rows();
    if(n == 1) {
      Q = Eigen::MatrixXd(1,1);
      R = Eigen::MatrixXd(1,1);
      Q(0,0) = 1;
      R(0,0) = T(0,0);
      return 1;
    }
    Q = Eigen::MatrixXd::Identity(n, n);
    R = T;
    // compute first givens matrix
    Q(0,0) = eigen_solver_utils::givens_c(R(0,0), R(1,0));
    Q(1,1) = Q(0,0);
    Q(0,1) = eigen_solver_utils::givens_s(R(0,0), R(1,0));
    Q(1,0) = -Q(0,1);
    double coeffs[2];
    R = Q.transpose()*R;
    for(int k = 1; k < n-2; k++) {
      coeffs[0] = eigen_solver_utils::givens_c(R(k,k), R(k+1,k));
      coeffs[1] = eigen_solver_utils::givens_s(R(k,k), R(k+1,k));      
      double temp[2][k+1];
      for(int j = k; j < k+2; j++) {
	for(int i = 0; i < k+1; i++) {
	  temp[j-k][i] = Q(i, k)*coeffs[j-k];
	}
      }
      for(int j = k; j < k+2; j++) {
	for(int i = 0; i < k+1; i++) {
	  Q(i,j) = temp[j-k][i];
	}
      }
      Q(k+1,k) = -coeffs[1];
      Q(k+1,k+1) = coeffs[0];
      for(int j = k; j < k+3; j++) {
	double temp1 = R(k,j);
	double temp2 = R(k+1,j);
	R(k,j) = coeffs[0]*temp1 - coeffs[1]*temp2;
	R(k+1,j) = coeffs[1]*temp1 + coeffs[0]*temp2;
      }
    }
    // no if statements, dammit
    int k = n-2;
    coeffs[0] = eigen_solver_utils::givens_c(R(k,k), R(k+1,k));
    coeffs[1] = eigen_solver_utils::givens_s(R(k,k), R(k+1,k));      
    double temp[2][k+1];
    for(int j = k; j < k+2; j++) {
      for(int i = 0; i < k+1; i++) {
	temp[j-k][i] = Q(i, k)*coeffs[j-k];
      }
    }
    for(int j = k; j < k+2; j++) {
      for(int i = 0; i < k+1; i++) {
	Q(i,j) = temp[j-k][i];
      }
    }
    Q(k+1,k) = -coeffs[1];
    Q(k+1,k+1) = coeffs[0];
    for(int j = k; j < k+2; j++) {
      double temp1 = R(k,j);
      double temp2 = R(k+1,j);
      R(k,j) = coeffs[0]*temp1 - coeffs[1]*temp2;
      R(k+1,j) = coeffs[1]*temp1 + coeffs[0]*temp2;
    }
    return 1;
  }

  int arnoldi_iter(const Eigen::MatrixXd& A, const Eigen::MatrixXd& V_init, const Eigen::MatrixXd& H_init, const Eigen::VectorXd& f_init, Eigen::MatrixXd& V, Eigen::MatrixXd& H, Eigen::VectorXd& f, const int m, const double zero_tol) {
    const int n = A.rows();
    const int k_init = V_init.rows();

    V = Eigen::MatrixXd(n, m);
    H = Eigen::MatrixXd::Zero(m, m);
    f = f_init;
    V.block(0, 0, n, k_init) = V_init;
    H.block(0, 0, k_init, k_init) = H_init;
    
    for(int k = k_init; k < m; k++) {
      double fnorm = f.norm();
      if(fnorm < zero_tol) {
	std::cout << "arnoldi iteration failed" << std::endl;
	std::cout << "initial vector only spans " << k << " dimensional subspace" << std::endl;
	std::cout << "|| f || = " << fnorm << std::endl;
	return 0;
      }
      V.col(k) = f/fnorm;
      H(k, k-1) = fnorm;
      f = A*V.col(k);
      for(int i = 0; i < k+1; i++) {
	H(i,k) = V.col(i).dot(f);
	f = f - H(i,k)*V.col(i);
      }
    }
    return 1;
  }

  int arnoldi_iter(const Eigen::MatrixXd& A, const Eigen::VectorXd& v, Eigen::MatrixXd& V, Eigen::MatrixXd& H, Eigen::VectorXd& f, const int m, const double zero_tol) {
    const int n = A.rows();

    V = Eigen::MatrixXd(n, m);
    H = Eigen::MatrixXd::Zero(m, m);
    V.block(0, 0, n, 1) = v/v.norm();
    f = A*V.col(0);
    H(0,0) = V.col(0).dot(f);
    f = f - H(0,0)*V.col(0);
    
    for(int k = 1; k < m; k++) {
      double fnorm = f.norm();
      if(fnorm < zero_tol) {
	std::cout << "arnoldi iteration failed" << std::endl;
	std::cout << "initial vector only spans " << k << " dimensional subspace" << std::endl;
	std::cout << "|| f || = " << fnorm << std::endl;
	return 0;
      }
      V.col(k) = f/fnorm;
      H(k, k-1) = fnorm;
      f = A*V.col(k);
      for(int i = 0; i < k+1; i++) {
	H(i,k) = V.col(i).dot(f);
	f = f - H(i,k)*V.col(i);
      }
    }
    return 1;
  }

  int qr_impshift_tridiag(const Eigen::MatrixXd& T, Eigen::MatrixXd& S, Eigen::VectorXd& thetas, const int maxiter, const double eigen_tol, const double zero_tol) {
    const int n = T.rows();
    Eigen::MatrixXd Tk = T;
    
    // in case n == 1 || n == 2
    // calculate eigenpairs exactly
    if(n == 1) {
      S = Eigen::MatrixXd(1,1);
      thetas = Eigen::VectorXd(1);
      S(0,0) = 1;
      thetas(0) = T(0,0);
      return 1;
    }
    else if(n == 2) {
      double a = T(0,0);
      double b = T(0,1);
      double c = T(1,1);
      double d = a*c - b*b;
      thetas = Eigen::VectorXd(2);
      thetas(0) = (a+c + sqrt(pow(a+c, 2) - 4*d))/2.0;
      thetas(1) = (a+c - sqrt(pow(a+c, 2) - 4*d))/2.0;
      S = Eigen::MatrixXd(2,2);
      S(1,0) = 1;
      S(0,0) = (thetas(0) - b - c)*S(1,0)/(a + b - thetas(0));
      S.col(0) = S.col(0)/(S.col(0).norm());
      S(1,1) = -1;
      S(0,1) = (thetas(1) - b - c)*S(1,1)/(a + b - thetas(1));
      S.col(1) = S.col(1)/(S.col(1).norm());
      return 1;
    }
    
    Eigen::VectorXd gersh_rings(n);
    S = Eigen::MatrixXd::Identity(n, n);
    double err = 1;
    int iters = 0;
    const double IMAG_TOL = 1e-16;
    while(err > eigen_tol && iters < maxiter) {
      std::complex<double> shift = eigen_solver_utils::wilk_shift(Tk(n-2,n-2), Tk(n-2,n-1), Tk(n-1,n-1));
      
      if(abs(std::imag(shift)) > IMAG_TOL) {
	/* 
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	          Still untested
	   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	*/
	// speed may be increased by manually multiplying
	// matrices and vectors given tridiagonal structure
	// todo: just address the blocks directly through Eigen
	std::cout << "deploy le francis shift" << std::endl;
	double realshift = std::real(shift);
	double s = 2*realshift;
	double t = std::norm(shift);
	// reflection vector
	Eigen::VectorXd refl_v(3);
	refl_v(0) = Tk(0,0)*Tk(0,0) + Tk(0,1)*Tk(0,1) - s*Tk(0,0) + t;
	refl_v(1) = Tk(0,1)*(Tk(0,0) + Tk(1,1) - s);
	refl_v(2) = Tk(0,1)*Tk(2,1);
	refl_v(0) += copysign(refl_v.norm(), refl_v(0));
	// reflection matrix
	Eigen::MatrixXd refl_m = Eigen::MatrixXd::Identity(n,n);
	refl_m.block(0, 0, 3, 3) -=  2*refl_v*refl_v.transpose()/refl_v.squaredNorm();
	Tk = refl_m.transpose()*Tk*refl_m;
	S = refl_m*S;
	for(int i = 0; i < n-3; i++) {
	  refl_v = Tk.block(i+1, i, 3, 1);
	  refl_v(0) += copysign(refl_v.norm(), refl_v(0));
	  refl_m = Eigen::MatrixXd::Identity(n,n);
	  refl_m.block(i+1, i+1, 3, 3) -= 2*refl_v*refl_v.transpose()/refl_v.squaredNorm();
	  Tk = refl_m.transpose()*Tk*refl_m;
	  S = refl_m*S;
	}
	refl_v = Tk.block(n-2, n-3, 2, 1);
	refl_v(0) += copysign(refl_v.norm(), refl_v(0));
	refl_m = Eigen::MatrixXd::Identity(n,n);
	refl_m.block(n-2, n-2, 2, 2) -= 2*refl_v*refl_v.transpose()/refl_v.squaredNorm();
	Tk = refl_m.transpose()*Tk*refl_m;
	S = refl_m*S;
      }
      else {
	// wilkinson shift
	Eigen::MatrixXd Qi, Ri;
	double wilk_shift = std::real(shift);
	qr_tridiag(Tk - wilk_shift*Eigen::MatrixXd::Identity(n, n), Qi, Ri);
	S = S*Qi;
	Tk = Ri*Qi + wilk_shift*Eigen::MatrixXd::Identity(n, n);
      }
      // check for zeros
      // there is a more sophisticated method
      // of doing this in the true "francis step"
      for(int i = 0; i < n-1; i++) {
	if(abs(Tk(i,i+1)) < zero_tol) {
	  Tk(i, i+1) = 0;
	  Tk(i+1, i) = 0;
	  Eigen::MatrixXd S_top, S_bottom;
	  Eigen::VectorXd thetas_top, thetas_bottom;
	  qr_impshift_tridiag(Tk.block(0, 0, i+1, i+1), S_top, thetas_top, maxiter, eigen_tol, zero_tol);
	  qr_impshift_tridiag(Tk.block(i+1, i+1, n-i-1, n-i-1), S_bottom, thetas_bottom, maxiter, eigen_tol, zero_tol);
	  thetas = Eigen::VectorXd(n);
	  thetas.head(i+1) = thetas_top;
	  thetas.tail(n-i-1) = thetas_bottom;
	  S.block(0, 0, i+1, i+1) *= S_top;
	  S.block(i+1, i+1, n-i-1, n-i-1) *= S_bottom;
	  return 1;
	}
      }
      gersh_rings(0) = abs(Tk(0,1));
      gersh_rings(n-1) = abs(Tk(n-1,n-2));
      for(int i = 1; i < n-1; i++) {
	gersh_rings(i) = abs(Tk(i, i-1)) + abs(T(i, i+1));
      }
      err = gersh_rings.maxCoeff();
      iters++;
    }
    if(iters == maxiter) {
      std::cout << "qr failed to converge with n = " << n << std::endl;
      return 0;
    }
    thetas = Tk.diagonal();
    return 1;
  }





  int arnoldi_method_imprestart_hermitian(const Eigen::MatrixXd& A, const Eigen::VectorXd& v, Eigen::MatrixXd& V_ritz, Eigen::VectorXd& l_ritz, const int k, const int p, const double iram_maxiter, const int qr_maxiter, const double f_tol, const double qr_eigen_tol, const double qr_zero_tol){return 1;}

}
