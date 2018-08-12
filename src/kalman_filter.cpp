#include "kalman_filter.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using namespace std;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

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

void KalmanFilter::Predict(const double delta_T) {
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

  // new estimate
  // x_ is new state matrix expressing most likely state (2d position and 2d velocity)
  // P_ is new covarience matrix expressing uncertainty of the state
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
    // x_ is predicted state matrix set by caller before this function is called
    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);

    // Convert state matrix from Cartesian coordinate to Polar coordinate
    // to match with measurement domain
    VectorXd z_pred = VectorXd(3);
    z_pred << sqrt(px*px + py*py),
              atan(py / px),
              (px*vx + py*vy) / sqrt(px*px + py*py);

    // Adjust the result of atan(py/px) to have range of -pi to pi
    // Note: atan2 should remove the necessity of this step but atan2 was not available to my environment
    if(px<0) {
        if(py>=0) {
            z_pred(1) += 3.1415926535;
        } else {
            z_pred(1) -= 3.1415926535;
        }
    }

    // Calculate y vector
    // Adjust the operation result to have range between -pi to pi
    VectorXd y = z - z_pred;

    if(y(1) >= 3.1415926535) {
        y(1) -= 2*3.1415926535;
    } else if (y(1) <= -3.1415926535) {
        y(1) += 2*3.1415926535;
    }

    // Calculate Jacobian matrix
    Tools tools;
    MatrixXd Hj = tools.CalculateJacobian(x_);
    MatrixXd Ht = Hj.transpose();

    // Use Jacobian matrix instead of transformation Matrix H comparing to Normal Kalman Filter
    // Note that Jacobian matrix makes the transformation linear
    // though Cartesian to Polar coordinate transformation is not linear
    // Non linear conversion here breaks Gaussian property that breaks
    // Kalman Filter prerequisite to work.
    MatrixXd S = Hj * P_ * Ht + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHt = P_ * Ht;
    MatrixXd K = PHt * Si;

    // new estimate
    // x_ is new state matrix expressing most likely state (2d position and 2d velocity)
    // P_ is new covarience matrix expressing uncertainty of the state
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * Hj) * P_;
}
