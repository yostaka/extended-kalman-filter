#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    /**
    TODO:
      * Initialize the state ekf_.x_ with the first measurement.
      * Create the covariance matrix.
      * Remember: you'll need to convert radar from polar to cartesian coordinates.
    */
    // first measurement
    cout << "Initialize state vector x" << endl;
    ekf_.x_ = VectorXd(4);

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      cout << "Initialize state vector x with Radar measurement" << endl;
      double ro = measurement_pack.raw_measurements_[0];
      double theta = measurement_pack.raw_measurements_[1];
      double ro_dot = measurement_pack.raw_measurements_[2];

      double px = ro * cos(theta);
      double py = ro * sin(theta);
      double vx = ro_dot * cos(theta);
      double vy = ro_dot * sin(theta);
      ekf_.x_ << px, py, vx, vy;
      cout << "x_ = " << ekf_.x_ << endl;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      cout << "Initialize state vector x with Laser measurement" << endl;
      double px = measurement_pack.raw_measurements_[0];
      double py = measurement_pack.raw_measurements_[1];
      ekf_.x_ << px, py, 0, 0;
      cout << "x_ = " << endl << ekf_.x_ << endl;

    }

    cout << "Create the covariance matrix P " << endl;
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 1000, 0, 0, 0,
            0, 1000, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;
    cout << "P_ = " << endl << ekf_.P_ << endl;

    previous_timestamp_ = measurement_pack.timestamp_;
    cout << "timestamp = " << previous_timestamp_ << endl;

//    ekf_.R_ = MatrixXd(2, 2);
//    ekf_.R_ << 0.0225, 0,
//               0, 0.0225;

    ekf_.H_ = MatrixXd(2, 4);
    ekf_.H_ << 1, 0, 0, 0,
               0, 1, 0, 0;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }


//   Skip radar for debugging purpose
//  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
//    return;
//  }
//  if (measurement_pack.sensor_type_ != MeasurementPackage::RADAR) {
//    return;
//  }


  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  long long current_timestamp = measurement_pack.timestamp_;
  double dt = (current_timestamp - previous_timestamp_) / 1000000.0;
  double dt_2 = dt * dt;
  double dt_3 = dt_2 * dt;
  double dt_4 = dt_3 * dt;

//  cout << "dt = " << dt << endl;
//  cout << "current: " << current_timestamp << endl;
//  cout << "previous:" << previous_timestamp_ << endl;

  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;

  ekf_.Q_ = MatrixXd(4, 4);
  float noise_ax = 9, noise_ay = 9;

  ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
          0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
          dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
          0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

//  cout << "Q_ = " << endl << ekf_.Q_ << endl;

  ekf_.Predict(dt);

  previous_timestamp_ = current_timestamp;

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
