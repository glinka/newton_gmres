#ifndef TEST_FNS_H_
#define TEST_FNS_H_

#include <Eigen/Dense>

Eigen::VectorXd F(const Eigen::VectorXd& x);

Eigen::MatrixXd DF(const Eigen::VectorXd& x);
    
#endif
