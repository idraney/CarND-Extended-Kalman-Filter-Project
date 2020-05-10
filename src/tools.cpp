#include "tools.h"
#include <iostream>

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
  /**
   * TODO: Calculate the RMSE here.
   */

  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  if (estimations.size() == 0)
  {
    std::cout << "The estimation vector size is zero. Exiting" << std::endl;

    return rmse;
  }

  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size())
  {
    std::cout << "The estimation vector size and ground truth vector size are not the same. Exiting" << std::endl;
    std::cout << "Estimation vector size:\t" << estimations.size() << std::endl;
    std::cout << "Ground truth vector size:\t" << ground_truth.size() << std::endl;

    return rmse;
  }

  // TODO: accumulate squared residuals
  for (int i = 0; i < estimations.size(); ++i)
  {
    VectorXd residual = estimations[i] - ground_truth[i];

    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // TODO: calculate the mean
  rmse = rmse / estimations.size();

  // TODO: calculate the squared root
  rmse = rmse.array().sqrt();

  std::cout << std::endl
            << "RMSE: \n"
            << rmse << std::endl
            << std::endl;

  // return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state)
{
  /**
   * TODO:
   * Calculate a Jacobian here.
   */
  MatrixXd Hj(3, 4);

  // recover state parameters
  double px = x_state(0);
  double py = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  // TODO: YOUR CODE HERE

  double px2pluspy2 = pow(px, 2) + pow(py, 2);

  // Prevent division by zero
  if (px2pluspy2 < 0.000001)
  {
    std::cout << "The square root of px and py is (very close to) zero!  Cannot compute Jacobian matrix." << std::endl;

    return Hj;
  }

  double sqrt_px2pluspy2 = sqrt(px2pluspy2);
  double vxpy = (vx * py);
  double vypx = (vy * px);

  // Compute the Jacobian matrix
  Hj << (px / sqrt_px2pluspy2), (py / sqrt_px2pluspy2), 0.0, 0.0,
      (-1.0 * py / px2pluspy2), (px / px2pluspy2), 0.0, 0.0,
      (py * (vxpy - vypx) / pow(px2pluspy2, 3.0 / 2.0)),
      (px * (vypx - vxpy) / pow(px2pluspy2, 3.0 / 2.0)),
      (px / sqrt_px2pluspy2),
      (py / sqrt_px2pluspy2);

  return Hj;
}
