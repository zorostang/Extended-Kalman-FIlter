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
  //previous_timestamp_ = 0;
  counter = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2); //measurement covariance matrix - laser

  R_laser_ << 0.0225, 0,
              0, 0.0225;

  R_radar_ = MatrixXd(3, 3); //measurement covariance matrix - radar

  R_radar_ << 0.09, 0, 0,
              0, 0.0009, 0,
              0, 0, 0.09;

  H_laser_ = MatrixXd(2, 4);
  H_laser_ << 1, 0, 0, 0,
              0, 1, 0, 0;

  ekf_.F_ = MatrixXd(4,4);
  ekf_.F_ << 1, 0, 1, 0,
             0, 1, 0, 1,
             0, 0, 1, 0,
             0, 0, 0, 1;

  ekf_.P_ = MatrixXd(4,4);
  ekf_.P_ <<  1, 0, 0, 0,
              0, 1, 0, 0,
              0, 0, 1000, 0,
              0, 0, 0, 1000;

  ekf_.x_ = VectorXd(4);
  ekf_.x_ << 1,1,1,1;

  ekf_.Q_ = MatrixXd(4,4); // process covarience matrix
  Hj_ = MatrixXd(3, 4);
  noise_ax = 9;
  noise_ay = 9;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
    float x_cart;
    float y_cart;
    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {

      if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {

        // Convert radar from polar to cartesian coordinates and initialize state.
        x_cart = measurement_pack.raw_measurements_[0] * cos(measurement_pack.raw_measurements_[1]);
        y_cart = measurement_pack.raw_measurements_[0] * sin(measurement_pack.raw_measurements_[1]);
        
        if (x_cart == 0 || y_cart == 0) {
          x_cart = 0.1;
          y_cart = 0.1;
        }
        
        // Initialize state.
        ekf_.x_ << x_cart, y_cart, 0, 0;
      }
      else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
        if (measurement_pack.raw_measurements_[0] == 0 or measurement_pack.raw_measurements_[1] == 0){
           x_cart = 0.1;
           y_cart = 0.1;
        }
        // Initialize state.
        ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
      }

      // done initializing, no need to predict or update
      previous_timestamp_ = measurement_pack.timestamp_;
      is_initialized_ = true;
      return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    // compute the time elapsed between the current and previous measurements
    // dt - expressed in seconds
    double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	
    cout << dt << "** dt" << endl;
    double dt_2 = dt * dt;
    double dt_3 = dt_2 * dt;
    double dt_4 = dt_3 * dt;
    if (dt != 0)  {
    // update previous timestamp for next iteration
    previous_timestamp_ = measurement_pack.timestamp_;

    //1. Modify the F matrix so that the time is integrated
    ekf_.F_(0,2) = dt;
    ekf_.F_(1,3) = dt;
  
    //2. Set the process covariance matrix Q
    ekf_.Q_ <<  (dt_4/4)*noise_ax, 0, (dt_3/2)*noise_ax, 0,
                0, (dt_4/4)*noise_ay, 0, (dt_3/2)*noise_ay,
                (dt_3/2)*noise_ax, 0, dt_2*noise_ax, 0,
                0, (dt_3/2)*noise_ay, 0, dt_2*noise_ay;
    ekf_.Predict();
    }
  /*****************************************************************************
   *  Update
   ****************************************************************************/
  
  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates
    ekf_.R_ = R_radar_;
    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  if (counter < 6) {
    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
    counter++;
  }
  
}
