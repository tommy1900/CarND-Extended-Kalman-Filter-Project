#include "FusionEKF.h"
#include <iostream>
#include <sstream>
#include <vector>
#include "measurement_package.h"
#include "Eigen/Dense"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;
using std::vector;


/**
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
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */
   H_laser_ << 1, 0, 0, 0,
               0, 1, 0, 0;

   // Hj_ do not have to initialize Hj here

   // Init the F matrix
   ekf_.F_ = MatrixXd(4,4);
   ekf_.F_ << 1, 0, 1, 0,
              0, 1, 0, 1,
              0, 0, 1, 0,
              0, 0, 0, 1;

   // Init the uncertainty
   ekf_.P_ = MatrixXd(4,4);
   ekf_.P_ << 1, 0, 0, 0,
             0, 1, 0, 0,
             0, 0, 1000, 0,
             0, 0, 0, 1000;
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
  /**
   * Initialization
   */
  if (!is_initialized_) {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // first measurement
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      // TODO: Convert radar from polar to cartesian coordinates
      //         and initialize state.
      double rho = measurement_pack.raw_measurements_[0];
      double phi = measurement_pack.raw_measurements_[1];
      double rho_dot = measurement_pack.raw_measurements_[2];
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);

      // check for 0 measurements
      if(fabs(px) < 0.0001) px = 0.0001;
      if(fabs(py) < 0.0001) py = 0.0001;

      ekf_.x_ << px,py,vx,vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      // TODO: Initialize state.
      ekf_.x_ << measurement_pack.raw_measurements_[0],
                 measurement_pack.raw_measurements_[1],
                 0,
                 0;
      // px,py,0,0 directly from laser sensor
    }
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // Saving first timestamp in seconds
    // NOTE VERY IMPORTANT: If we dont have an initialized timestamp at the begining,
    // the first x_in will mess up and we will get a huge position RMSE
    previous_timestamp_ = measurement_pack.timestamp_ ;
    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /**
   * Prediction
   */

  /**
   * TODO: Update the state transition matrix F according to the new elapsed time.
   * Time is measured in seconds.
   * TODO: Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

   float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
   previous_timestamp_ = measurement_pack.timestamp_;
   float dt_2 = dt   * dt;
   float dt_3 = dt_2 * dt;
   float dt_4 = dt_3 * dt;

  //********************************************************
  // NOTE WE STILL USING A LINEAR MODEL HERE F STAY THE SAME
  //********************************************************

  // State transition matrix update -> same as the linear one
  ekf_.F_ = MatrixXd(4, 4);
  ekf_.F_ << 1, 0, dt, 0,
             0, 1, 0, dt,
             0, 0, 1, 0,
             0, 0, 0, 1;
  // Noise coveriance matrix -> same as the linear one
  float noise_ax = 9;
  float noise_ay = 9;

  //set the process covariance matrix Q
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
             0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
             dt_3 / 2 * noise_ax, 0, dt_2*noise_ax, 0,
             0, dt_3 / 2 * noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // TODO: Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_); //update the input x_ into Jaco
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_); // re-direct to the kalman_filter.cpp

  } else {
    // TODO: Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_); // linear measurement here
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
