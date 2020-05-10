#include "FusionEKF.h"
#include <iostream>
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
FusionEKF::FusionEKF()
{
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  // Measurement covariance matrix - LIDAR
  R_laser_ << 0.0225, 0,
      0, 0.0225;

  // Measurement covariance matrix - RADAR
  R_radar_ << 0.09, 0, 0,
      0, 0.0009, 0,
      0, 0, 0.09;

  /**
   * TODO: Finish initializing the FusionEKF.
   * TODO: Set the process and measurement noises
   */

  // Define: (x, P, F, H, R, Q)

  // x: First measurement. 4D state vector.  Initialized in ProcessMeasurement function.

  // State covariance matrix P
  ekf_.P_ = MatrixXd(4, 4);
  ekf_.P_ << 1, 0, 0, 0,
      0, 1, 0, 0,
      0, 0, 1000, 0,
      0, 0, 0, 1000;

  // Initial transition matrix F.  Time integrated below.
  ekf_.F_ = MatrixXd(4, 4);

  // H: Measurement matrices.  RADAR defined by Jacobian matrix Hj_
  // LIDAR
  H_laser_ << 1, 0, 0, 0,
      0, 1, 0, 0;

  // R: Measurement covariance matrices for LIDAR and RADAR.  Defined above (R_laser_ and R_radar_).

  // Process covariance matrix Q.  Computed below
  ekf_.Q_ = MatrixXd(4, 4);
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  /**
   * Initialization
   */
  if (!is_initialized_)
  {
    /**
     * TODO: Initialize the state ekf_.x_ with the first measurement.
     * TODO: Create the covariance matrix.
     * You'll need to convert radar from polar to cartesian coordinates.
     */

    // First measurement. 4D state vector.
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      // TODO: Convert RADAR data from polar to cartesian coordinates
      //         and initialize state.
      cout << "EKF Measuring RADAR:" << endl;

      // Collect data
      double rho = measurement_pack.raw_measurements_[0];     // Range
      double phi = measurement_pack.raw_measurements_[1];     // Bearing
      double rho_dot = measurement_pack.raw_measurements_[2]; // Range rate (e.g., Doppler or Radial velocity)

      // Calculate positions & velocities
      double px = rho * cos(phi);
      double py = rho * sin(phi);
      double vx = rho_dot * cos(phi);
      double vy = rho_dot * sin(phi);

      // Set the state vector x
      ekf_.x_ << px, py, vx, vy;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
    {
      // TODO: Initialize state.
      cout << "EKF Measuring LIDAR:" << endl;

      // Set the state with the initial location and zero velocity
      ekf_.x_ << measurement_pack.raw_measurements_[0],
          measurement_pack.raw_measurements_[1],
          0,
          0;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

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

  // compute the time elapsed between the current and previous measurements
  // dt: expressed in seconds, dtu: expressed in microseconds
  double dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  // double dtu = (measurement_pack.timestamp_ - previous_timestamp_);
  previous_timestamp_ = measurement_pack.timestamp_;

  // Modify the F matrix so that the time is integrated
  ekf_.F_ << 1, 0, dt, 0,
      0, 1, 0, dt,
      0, 0, 1, 0,
      0, 0, 0, 1;

  // Acceleration noise components
  double noise_ax = 9.0; // Sigma squared ax
  double noise_ay = 9.0; // Sigma squared ay

  // Set the process covariance matrix Q
  ekf_.Q_ << (pow(dt, 4) / 4.0) * noise_ax, 0.0, (pow(dt, 3) / 2.0) * noise_ax, 0.0,
      0.0, (pow(dt, 4) / 4.0) * noise_ay, 0.0, (pow(dt, 3) / 2.0) * noise_ay,
      (pow(dt, 3) / 2.0) * noise_ax, 0.0, (pow(dt, 2)) * noise_ax, 0.0,
      0.0, (pow(dt, 3) / 2.0) * noise_ay, 0.0, (pow(dt, 2)) * noise_ay;

  ekf_.Predict();

  /**
   * Update
   */

  /**
   * TODO:
   * - Use the sensor type to perform the update step.
   * - Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
  {
    // TODO: Radar updates
    cout << "EKF Updating RADAR:" << endl;

    Hj_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.H_ = Hj_;
    ekf_.R_ = R_radar_;
    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  }
  else
  {
    // TODO: Laser updates
    cout << "EKF Updating LIDAR:" << endl;

    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // Print the output
  cout << "x_ = " << endl
       << ekf_.x_ << endl;
  cout << "P_ = " << endl
       << ekf_.P_ << endl;
}
