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
  /**
  TODO:
    * predict the state
  */
    x_ = F_ * x_;
    MatrixXd Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;

//    cout << "*** Predict ***" << endl;
//    cout << "F_ = " << endl << F_ << endl;
//    cout << "P_ = " << endl << P_ << endl;
//    cout << "Q_ = " << endl << Q_ << endl;
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
//  cout << "*** Update ***" << endl;

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
    TODO:
      * update the state by using Extended Kalman Filter equations
    */
//    cout << "*** Update EKF ***" << endl;

    Tools tools;
    MatrixXd Hj = tools.CalculateJacobian(x_);

//    cout << "Hj = " << endl << Hj << endl;

    float px = x_(0);
    float py = x_(1);
    float vx = x_(2);
    float vy = x_(3);

    VectorXd z_pred = VectorXd(3);
    z_pred << sqrt(px*px + py*py),
              atan(py / px),
              (px*vx + py*vy) / sqrt(px*px + py*py);

//    cout << "z_pred = " << endl << z_pred << endl;
//    cout << "z = " << endl << z << endl;

    if(px<0) {
        if(py>=0) {
            z_pred(1) += 3.1415926535;
        } else {
            z_pred(1) -= 3.1415926535;
        }
    }
//    cout << "phi = " << z_pred(1) << " (" << px << ", " << py << ")" << endl;

    VectorXd y = z - z_pred;

    if(y(1) >= 3.1415926535) {
        y(1) -= 2*3.1415926535;
    } else if (y(1) <= -3.1415926535) {
        y(1) += 2*3.1415926535;
    }

//    cout << "y = " << endl << y << endl;

    MatrixXd Ht = Hj.transpose();
//    cout << "Ht = " << endl << Ht << endl;

    MatrixXd S = Hj * P_ * Ht;
//    cout << "S = " << endl << S << endl;
//    cout << "R_ = " << endl << R_ << endl;
    S = S + R_;
//    cout << "S = " << endl << S << endl;

    MatrixXd Si = S.inverse();
//    cout << "Si = " << endl << Si << endl;

    MatrixXd PHt = P_ * Ht;
//    cout << "PHt = " << endl << PHt << endl;

    MatrixXd K = PHt * Si;
//    cout << "K = " << endl << K << endl;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * Hj) * P_;
}
