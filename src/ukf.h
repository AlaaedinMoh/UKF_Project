#ifndef UKF_H
#define UKF_H

#include "Eigen/Dense"
#include "measurement_package.h"
#include <iostream>

class UKF {
  private:
  //intialize state vector and covariance matrix from radar measurments
  void InitializeFromRadar(MeasurementPackage* rdr_meas_package);
  //intialize state vector and covariance matrix from lidar measurments
  void InitializeFromLidar(MeasurementPackage* ldr_meas_package);
  //initialize noise matrix
  void InitializeNoise();
  void InitAugmentedStateVect(Eigen::VectorXd& x_aug);
  void InitAugmentedCovarMat(Eigen::MatrixXd& p_aug);
  void Init_X_sig_mat(Eigen::MatrixXd& Xsig_aug, const Eigen::VectorXd& x_aug, const Eigen::MatrixXd& sqrt_p_aug);
  void DoPredictSigmaPoints(const Eigen::MatrixXd& xSigMat, double dt);
  void PredictMeanAndCovar();
  void PredictRadarMeasurements(Eigen::VectorXd& rdr_meas_p, Eigen::MatrixXd& rdr_cov_p);
  void PredictLaserMeasurements(Eigen::VectorXd& lsr_meas_p, Eigen::MatrixXd& lsr_cov_p);
  void CalculateCorrel_Mat(Eigen::MatrixXd& corrMat);
 public:
  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * ProcessMeasurement
   * @param meas_package The latest measurement data of either radar or laser
   */
  void ProcessMeasurement(MeasurementPackage meas_package);

  /**
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   * @param delta_t Time between k and k+1 in s
   */
  void Prediction(double delta_t);

  /**
   * Updates the state and the state covariance matrix using a laser measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateLidar(MeasurementPackage meas_package);

  /**
   * Updates the state and the state covariance matrix using a radar measurement
   * @param meas_package The measurement at k+1
   */
  void UpdateRadar(MeasurementPackage meas_package);


  // initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  // if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  // if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  // state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  Eigen::VectorXd x_;

  // state covariance matrix
  Eigen::MatrixXd P_;

  // predicted sigma points matrix
  Eigen::MatrixXd Xsig_pred_;

  // time when the state is true, in us
  long long time_us_;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  // Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  // Laser measurement noise standard deviation position1 in m
  double std_laspx_;

  // Laser measurement noise standard deviation position2 in m
  double std_laspy_;

  // Radar measurement noise standard deviation radius in m
  double std_radr_;

  // Radar measurement noise standard deviation angle in rad
  double std_radphi_;

  // Radar measurement noise standard deviation radius change in m/s
  double std_radrd_ ;

  // Weights of sigma points
  Eigen::VectorXd weights_;

  // State dimension
  int n_x_;

  // Augmented state dimension
  int n_a_;

  // Sigma point spreading parameter
  double lambda_;

  Eigen::MatrixXd sig_noise_;

  int n_noise_;
};

#endif  // UKF_H