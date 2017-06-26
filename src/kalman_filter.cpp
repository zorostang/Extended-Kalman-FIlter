#include "kalman_filter.h"
#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using Eigen::Matrix3f;
using namespace std;
KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Predict() {
  /**
  predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
    z is our sensor measurments (polar coords)
    x_ is our predicted state [x',y',vx',vy'] 
    use h(x_) eqs to map cartesian coords to polar
    y=z−h(x')
  */
  float ro = sqrt(x_[0]*x_[0]+x_[1]*x_[1]);
  if (abs(ro) < 0.001) {
    cout << ro << " ** ro **" << endl;
  }
  float phi = atan2(x_[1], x_[0]);
  if (abs(phi) > M_PI) {
    cout << phi << " ** phi **" << endl;
  }
  float rodot = (x_[0]*x_[2]+x_[1]*x_[3])/ro;
  VectorXd z_pred(3);
  z_pred << ro,
            phi,
            rodot;
  VectorXd y = z - z_pred;
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;
  
  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
