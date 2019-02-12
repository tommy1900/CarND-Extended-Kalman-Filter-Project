#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/*
 * Please note that the Eigen library does not initialize
 *   VectorXd or MatrixXd objects with zeros upon creation.
 */

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
   * TODO: predict the state
   */
   x_ = F_ * x_; // x(k+1) = F * x(k) + v -> predictor
   MatrixXd Ft = F_.transpose();
   P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
   * TODO: update the state by using Kalman Filter equations
   */
   VectorXd z_pred = H_ * x_; // measurements = H * x(k+1) + w
   VectorXd y = z - z_pred; // error = true - est
   MatrixXd Ht = H_.transpose();
   MatrixXd S = H_ * P_ * Ht + R_; // Sensor uncertainty Matrix
   MatrixXd Si = S.inverse();
   MatrixXd PHt = P_ * Ht;
   MatrixXd K = PHt * Si; // KF gain

   // New estimation
   x_ = x_ + (K * y); // x_new = x_old + correction
   long x_size = x_.size();
   MatrixXd I = MatrixXd::Identity(x_size,x_size);
   P_ = (I - K * H_) * P_; //P_new = P_old - correction
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
   * TODO: update the state by using Extended Kalman Filter equations
   */
  MatrixXd Hj(3,4);

  // recover state parameters
  double px = x_(0);
  double py = x_(1);
  double vx = x_(2);
  double vy = x_(3);

  // set up measurement parameters from kinematic
  double rho = sqrt(px*px + py*py);
  double phi = atan2(py , px);
  double rho_dot = (px*vx + py*vy) / rho;

  // since the measurement of Radar is nonlinear:
  VectorXd h = VectorXd(3);
  h << rho, phi, rho_dot;
  VectorXd y = z - h;

  // check for the rho angle (+- pi)
  while(y(1) > M_PI || y(1) < -M_PI){
    if(y(1) > M_PI) y(1) -= M_PI;
    else y(1) += M_PI;
  }

  // update with the nonlinear y:
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K =  P_ * Ht * Si;

  // New estimation
  x_ = x_ + (K * y); // x_new = x_old + correction
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size,x_size);
  P_ = (I - K * H_) * P_; //P_new = P_old - correction

}
