#include "eigen_solvers.h"
#include "eigen_solver_utils.h"
#include <random>
#include <chrono>
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
    const int k_init = V_init.cols();

    V = Eigen::MatrixXd(n, m);
    H = Eigen::MatrixXd::Zero(m, m);
    f = f_init;
    V.block(0, 0, n, k_init) = V_init;
    H.block(0, 0, k_init, k_init) = H_init;
    
    for(int k = k_init; k < m; k++) {
      double fnorm = f.norm();
      if(fnorm < zero_tol) {
	// std::cout << "arnoldi iteration failed" << std::endl;
	// std::cout << "initial vector only spans " << k << " dimensional subspace" << std::endl;
	// std::cout << "|| f || = " << fnorm << std::endl;
	// std::cout << "caution, entering untested code, expect problems" << std::endl;
	std::cout << "caution, subspace has been truncated to " << k-1 << " dimensions" << std::endl;
	if(k == 1) {
	  return k-1;
	}
	else {
	  f = H(k-1, k-2)*V.col(k-1).eval();
	  V = V.block(0, 0, n, k-1).eval();
	  H = H.block(0, 0, k-1, k-1).eval();
	  return k-1;
	}
      }
      V.col(k) = f/fnorm;
      H(k, k-1) = fnorm;
      f = A*V.col(k);
      for(int i = 0; i < k+1; i++) {
	H(i,k) = V.col(i).dot(f);
	f = f - H(i,k)*V.col(i);
      }
    }
    return m;
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
	// std::cout << "arnoldi iteration failed" << std::endl;
	// std::cout << "caution, entering untested code, expect problems" << std::endl;
	// std::cout << "|| f || = " << fnorm << std::endl;
	std::cout << "caution, subspace has been truncated to " << k-1 << " dimensions" << std::endl;
	if(k == 1) {
	  return k-1;
	}
	else {
	  f = H(k-1, k-2)*V.col(k-1).eval();
	  V = V.block(0, 0, n, k-1).eval();
	  H = H.block(0, 0, k-1, k-1).eval();
	  return k-1;
	}
      }
      V.col(k) = f/fnorm;
      H(k, k-1) = fnorm;
      f = A*V.col(k);
      for(int i = 0; i < k+1; i++) {
	H(i,k) = V.col(i).dot(f);
	f = f - H(i,k)*V.col(i);
      }
    }
    return m;
  }

  int qr_impshift_tridiag(const Eigen::MatrixXd& T, Eigen::MatrixXd& S, Eigen::VectorXd& thetas, const int maxiter, const double eigen_tol, const double zero_tol, const std::string shift_type) {
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
    while(err > eigen_tol && iters < maxiter) {
      if(shift_type == "wilkinson") {
	const double IMAG_TOL = 1e-16;
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
	  S *= Qi;
	  Tk = Ri*Qi + wilk_shift*Eigen::MatrixXd::Identity(n, n);
	}
      }
      else if(shift_type == "eigval") {
	double eigval_shift = Tk(n-1, n-1);
	Eigen::MatrixXd Qi, Ri;
	qr_tridiag(Tk - eigval_shift*Eigen::MatrixXd::Identity(n, n), Qi, Ri);
	S *= Qi;
	Tk = Ri*Qi + eigval_shift*Eigen::MatrixXd::Identity(n, n);
      }

      // check for zeros
      // there is a more sophisticated method
      // of doing this in the true "francis step"
      for(int i = 0; i < n-1; i++) {
	if(std::abs(Tk(i,i+1)) < zero_tol) {
	  Tk(i, i+1) = 0;
	  Tk(i+1, i) = 0;
	  Eigen::MatrixXd S_top, S_bottom;
	  Eigen::VectorXd thetas_top, thetas_bottom;
	  qr_impshift_tridiag(Tk.block(0, 0, i+1, i+1), S_top, thetas_top, maxiter, eigen_tol, zero_tol);
	  qr_impshift_tridiag(Tk.block(i+1, i+1, n-i-1, n-i-1), S_bottom, thetas_bottom, maxiter, eigen_tol, zero_tol);
	  thetas = Eigen::VectorXd(n);
	  thetas.block(0, 0, i+1, 1) = thetas_top;
	  thetas.block(i+1, 0, n-i-1, 1) = thetas_bottom;
	  S.block(0, 0, n, i+1) *= S_top;
	  S.block(0, i+1, n, n-i-1) *= S_bottom;
	  return 1;
	}
      }
      gersh_rings(0) = std::abs(Tk(0,1));
      gersh_rings(n-1) = std::abs(Tk(n-1,n-2));
      for(int i = 1; i < n-1; i++) {
	gersh_rings(i) = std::abs(Tk(i, i-1)) + std::abs(T(i, i+1));
      }
      err = gersh_rings.maxCoeff();
      iters++;
    }
    if(iters == maxiter && err > eigen_tol) {
      // std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      // std::cout << "qr failed to converge with n = " << n << std::endl;
      // std::cout << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      thetas = Tk.diagonal();
      return 0;
    }
    thetas = Tk.diagonal();
    // std::cout << "qr finished via gersh with err = " << err << std::endl;
    return 1;
  }





  int arnoldi_method_imprestart_hermitian(const Eigen::MatrixXd& A, const Eigen::VectorXd& v, Eigen::MatrixXd& V_ritz, Eigen::VectorXd& l_ritz, const int k, const int p, const int iram_maxiter, const int qr_maxiter, const double f_tol, const double qr_eigen_tol, const double qr_zero_tol){
    const int n = A.rows();
    int m = k + p;
    Eigen::MatrixXd H;
    Eigen::VectorXd f, f_old;
    int arnoldi_dims = arnoldi_iter(A, v/v.norm(), V_ritz, H, f_old, k);

    // pretty questionable method of dealing with the case
    // when v does not span a k-dimensional subspace
    const int max_reinits = 10;
    int reinits = 0;
    int maxdims = arnoldi_dims;
    while(arnoldi_dims != k && reinits < max_reinits) {
      unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
      std::mt19937 mt(seed);
      double normalization = (double) (mt.max()+1);
      Eigen::VectorXd V_new(n);
      for(int i = 0; i < n; i++) {
	V_new(i) = mt()/normalization;
      }
      arnoldi_dims = arnoldi_iter(A, V_new/V_new.norm(), V_ritz, H, f_old, k);
      if(arnoldi_dims > maxdims) {
	maxdims = arnoldi_dims;
      }
      reinits++;
    }
    if(reinits == max_reinits) {
      std::cout << "inital arnoldi-iteration vector never spanned " << k << "-dimensional subspace" << std::endl;
      std::cout << "inital arnoldi-iteration vector spanned max of " << maxdims << " dimensions" << std::endl;      
      return 0;
    }


    double err = 1;
    int iters = 0;
    std::vector<int> eigval_sorted_indices;
    Eigen::VectorXd thetas;
    Eigen::MatrixXd V, S;
    int qr_success = 0;
    while(err > f_tol && iters < iram_maxiter) {
      arnoldi_dims = arnoldi_iter(A, V_ritz, Eigen::MatrixXd(H.block(0, 0, k, k)), f_old, V, H, f, m);


      // more questionable error handling
      if(arnoldi_dims < m) {
	if(arnoldi_dims < k) {
	  std::cout << "restarting subtracted too many dimensions" << std::endl;
	  return 0;
	}
	else {
	  m = arnoldi_dims;
	}
      }


      qr_success = qr_impshift_tridiag(H, S, thetas, qr_maxiter);
      std::vector< std::pair<double, int> > to_sort(m);
      for(int i = 0; i < m; i++) {
	to_sort[i] = std::pair<double, int>(std::abs(thetas(i)), i);
      }
      eigval_sorted_indices = eigen_solver_utils::argsort(to_sort);
      Eigen::VectorXd errors(k);
      double fnorm = f.norm();
      for(int i = 1; i < k+1; i++) {
	// e_m * S[:, sorted[m-i]]
	errors(i-1) = fnorm*std::abs(S.col(eigval_sorted_indices[m-i])(m-1));
      }
      err = errors.maxCoeff();
      if(err > f_tol) {
	Eigen::MatrixXd Q = Eigen::MatrixXd::Identity(m, m);
	for(int i = 0; i < m-k; i++) {
	  double shift = thetas(eigval_sorted_indices[i]);
	  Eigen::MatrixXd Qi, Ri;
	  qr_tridiag(H - shift*Eigen::MatrixXd::Identity(m, m), Qi, Ri);
	  H = Qi.transpose()*H*Qi;
	  Q = Q*Qi;
	}
	// if we can resize, do it
	// otherwise, perform one last, extended (100*qr_maxiter)
	// attempt to find eigvals and exit
	// (continuing to iterate as normal will have no added benefit)
	if(m - k > 0) {
	  double beta = H(k,k-1);
	  double sigma = Q(m-1, k-1);
	  f_old = beta*V.col(k) + sigma*f;
	  V_ritz = V*Q.block(0, 0, m, k);
	}
	else {
	  qr_success = qr_impshift_tridiag(H, S, thetas, 100000*qr_maxiter);
	  std::vector< std::pair<double, int> > to_sort(m);
	  for(int i = 0; i < m; i++) {
	    to_sort[i] = std::pair<double, int>(std::abs(thetas(i)), i);
	  }
	  eigval_sorted_indices = eigen_solver_utils::argsort(to_sort);
	  Eigen::VectorXd errors(k);
	  double fnorm = f.norm();
	  for(int i = 1; i < k+1; i++) {
	    // e_m * S[:, sorted[m-i]]
	    errors(i-1) = fnorm*std::abs(S.col(eigval_sorted_indices[m-i])(m-1));
	  }
	  err = errors.maxCoeff();
	  l_ritz = Eigen::VectorXd(k);
	  for(int i = 1; i < k+1; i++) {
	    V_ritz.col(i-1) = V*S.col(eigval_sorted_indices[m-i]);
	    l_ritz(i-1) = thetas(eigval_sorted_indices[m-i]);
	  }
	  if(err > f_tol) {
	    std::cout << "iram failed on last desperate attempt of converge" << std::endl;
	    return 0;
	  }
	  return 1;
	}
	iters++;
      }
    }

    l_ritz = Eigen::VectorXd(k);
    for(int i = 1; i < k+1; i++) {
      V_ritz.col(i-1) = V*S.col(eigval_sorted_indices[m-i]);
      l_ritz(i-1) = thetas(eigval_sorted_indices[m-i]);
    }
    // do some half-assed error checking
    if(iters == iram_maxiter && err > f_tol) {
      // std::cout << "!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      std::cout << "iram failed to converge" << std::endl;
      // std::cout << "!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
      return 0;
    }
    if(qr_success != 1) {
      std::cout << "----- Caution: QR failed to converge during most recent iteration -----" << std::endl;
      return 0;
    }
    return 1;
  }

}
