#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::cout;
using std::endl;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  //User defined paramters
  n_x_ = 5;

  n_a_ = 7;

  Xsig_pred_ = MatrixXd(n_x_, 2*n_a_ + 1);

  lambda_ = 3 - n_a_;

  n_noise_ = 2;

  sig_noise_ = MatrixXd(n_noise_, n_noise_);

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 3.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 2.0;

  /**
   * DO NOT MODIFY measurement noise values below.
   * These are provided by the sensor manufacturer.
   */

  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;

  /**
   * End DO NOT MODIFY section for measurement noise values
   */

  /**
   * TODO: Complete the initialization. See ukf.h for other member properties.
   * Hint: one or more values initialized above might be wildly off...
   */
  is_initialized_ = false;
  weights_ = VectorXd(2*n_a_ + 1);
  weights_(0) = lambda_/(lambda_ + n_a_);
  for(int i = 1; i < 2*n_a_ + 1; i++)
  {
    weights_(i)= 0.5/(lambda_ + n_a_);
  }
}

UKF::~UKF() {}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
  try
  {
    if(!is_initialized_)
    {
      switch (meas_package.sensor_type_)
      {
        case MeasurementPackage::SensorType::LASER:
          if(use_laser_)
          {
            InitializeFromLidar(&meas_package);
          }
          break;
        case MeasurementPackage::SensorType::RADAR:
          if(use_radar_)
          {
            InitializeFromRadar(&meas_package);
          }
          break;
        default:
          time_us_ = 0;
          std::cerr<<"Undefined Sensor Type...!"<<std::endl;
          return;
      }
      time_us_ = meas_package.timestamp_;
      InitializeNoise();
      is_initialized_ = true;
      return;
    }

    double dt = (meas_package.timestamp_-time_us_)/1000000.0;
    time_us_ = meas_package.timestamp_;
    Prediction(dt);
    switch(meas_package.sensor_type_)
    {
      case MeasurementPackage::SensorType::LASER:
        UpdateLidar(meas_package);
        break;
      case MeasurementPackage::SensorType::RADAR:
        UpdateRadar(meas_package);
        break;
    }
  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }
}

void UKF::Prediction(double delta_t) {
  /**
   * TODO: Complete this function! Estimate the object's location.
   * Modify the state vector, x_. Predict sigma points, the state,
   * and the state covariance matrix.
   */
  try
  {
    VectorXd x_a;
    MatrixXd p_a;
    InitAugmentedStateVect(x_a);
    InitAugmentedCovarMat(p_a);
    MatrixXd L = p_a.llt().matrixL();
    MatrixXd xSig_a;
    Init_X_sig_mat(xSig_a, x_a, L);
    DoPredictSigmaPoints(xSig_a, delta_t);
    PredictMeanAndCovar();
  }
  catch(const std::exception& e)
  {
    std::cerr << e.what() << '\n';
  }
}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
  int n_z = meas_package.raw_measurements_.size();
  int zSig_cols = 2*n_a_ + 1;
  MatrixXd zSig_pred = MatrixXd(n_z, zSig_cols);
  VectorXd l_mean_p = VectorXd(n_z);
  MatrixXd l_cov_p = MatrixXd(n_z,n_z);
  MatrixXd corrMat = MatrixXd(n_x_, n_z);
  zSig_pred.fill(0.0);
  l_mean_p.fill(0.0);
  l_cov_p.fill(0.0);
  corrMat.fill(0.0);
  for(int i = 0; i < zSig_cols; ++i)
  {
    zSig_pred.col(i) = Xsig_pred_.col(i).head(n_z);
    l_mean_p += weights_(i)*zSig_pred.col(i);
  }
  for(int i = 0; i < zSig_cols; ++i)
  {
    VectorXd zDiff = zSig_pred.col(i) - l_mean_p;
    Eigen::Transpose<VectorXd> zDiff_T = zDiff.transpose();
    l_cov_p +=  weights_(i)*zDiff*zDiff_T;
    VectorXd xDiff = Xsig_pred_.col(i) - x_;
    while (xDiff(3) > M_PI) 
      xDiff(3) -= 2. * M_PI;
    while (xDiff(3) < -M_PI) 
      xDiff(3) += 2. * M_PI;
    corrMat += weights_(i)*xDiff*zDiff_T;
  }
  MatrixXd R = MatrixXd(n_z,n_z);
  R.fill(0.0);
  R(0,0) = pow(std_laspx_,2);
  R(1,1) = pow(std_laspy_,2);
  l_cov_p += R;
  MatrixXd kg = corrMat*l_cov_p.inverse();
  VectorXd z_diff = meas_package.raw_measurements_ - l_mean_p;
  x_ += kg*z_diff;
  P_ -= kg*l_cov_p*kg.transpose();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief
   * about the object's position. Modify the state vector, x_, and
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */
  int n_z = meas_package.raw_measurements_.size();
  int zSig_cols = 2*n_a_ + 1;
  MatrixXd zSig_pred = MatrixXd(n_z, zSig_cols);
  VectorXd r_mean_p = VectorXd(n_z);
  MatrixXd r_cov_p = MatrixXd(n_z,n_z);
  MatrixXd corrMat = MatrixXd(n_x_, n_z);
  zSig_pred.fill(0.0);
  r_mean_p.fill(0.0);
  r_cov_p.fill(0.0);
  corrMat.fill(0.0);
  for (size_t i = 0; i < zSig_cols; i++)
  {
    double px = Xsig_pred_(0,i);
    double py = Xsig_pred_(1,i);
    double v = Xsig_pred_(2,i);
    double yaw = Xsig_pred_(3,i);

    zSig_pred(0,i) = sqrt(px*px + py*py);
    zSig_pred(1,i) = atan2(py,px);
    zSig_pred(2,i) = (px*cos(yaw)*v + py*sin(yaw)*v) / zSig_pred(0,i);

    r_mean_p += weights_(i)*zSig_pred.col(i);
  }
  //Predicting the covariance and correlation matrix
  for (size_t i = 0; i < zSig_cols; i++)
  {
    VectorXd zDiff = zSig_pred.col(i) - r_mean_p;
    Eigen::Transpose<VectorXd> zDiff_T = zDiff.transpose();
    while (zDiff(1)> M_PI) zDiff(1)-=2.*M_PI;
    while (zDiff(1)<-M_PI) zDiff(1)+=2.*M_PI;
    r_cov_p += weights_(i)*zDiff*zDiff_T;

    VectorXd xDiff = Xsig_pred_.col(i) - x_;
    while (xDiff(3)> M_PI) xDiff(3)-=2.*M_PI;
    while (xDiff(3)<-M_PI) xDiff(3)+=2.*M_PI;
    corrMat += weights_(i)*xDiff*zDiff_T;
  }
  MatrixXd R = MatrixXd(n_z,n_z);
  R.fill(0.0);
  R(0,0) = pow(std_radr_,2);
  R(1,1) = pow(std_radphi_,2);
  R(2,2) = pow(std_radrd_,2);
  r_cov_p += R;
  MatrixXd kg = corrMat*r_cov_p.inverse();
  VectorXd z_diff = meas_package.raw_measurements_ - r_mean_p;
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  x_ += kg*z_diff;
  P_ -= kg*r_cov_p*kg.transpose();
}

void UKF::InitializeFromRadar(MeasurementPackage* rdr_meas_package){
  double rho = rdr_meas_package->raw_measurements_[0];
  double phi = rdr_meas_package->raw_measurements_[1];
  double rhoDot = rdr_meas_package->raw_measurements_[2];
  double px = rho*cos(phi);
  double py = rho*sin(phi);
  double vx = rhoDot*cos(phi);
  double vy = rhoDot*sin(phi);
  double v = sqrt(pow(vx,2)+pow(vy,2));
  x_<< px, py, v, rho, rhoDot;
  P_<<pow(std_radr_,2),0,0,0,0,
      0,pow(std_radr_,2),0,0,0,
      0, 0, pow(std_radrd_,2), 0, 0,
      0, 0, 0, std_radphi_, 0,
      0, 0, 0, 0, std_radphi_;
}

void UKF::InitializeFromLidar(MeasurementPackage* ldr_meas_package){
  double px = ldr_meas_package->raw_measurements_[0];
  double py = ldr_meas_package->raw_measurements_[1];
  x_<< px, py, 0, 0, 0;
  P_<<pow(std_laspx_,2),0,0,0,0,
      0,pow(std_laspy_,2),0,0,0,
      0, 0, 1, 0, 0,
      0, 0, 0, 1, 0,
      0, 0, 0, 0, 1;
}

void UKF::InitializeNoise(){
  sig_noise_(0,0) = std_a_*std_a_;
  sig_noise_(0,1) = 0;
  sig_noise_(1,0) = 0;
  sig_noise_(1,1) = std_yawdd_*std_yawdd_;
}

void UKF::InitAugmentedStateVect(VectorXd& x_aug){
  x_aug = VectorXd(n_a_);
  x_aug.head(n_x_) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
}

void UKF::InitAugmentedCovarMat(MatrixXd& p_aug){
  p_aug = MatrixXd(n_a_, n_a_);
  p_aug.fill(0.0);
  p_aug.topLeftCorner(n_x_,n_x_) = P_;
  p_aug.bottomRightCorner(2,2) = sig_noise_;
}

void UKF::Init_X_sig_mat(MatrixXd& Xsig_aug, const VectorXd& x_aug, const MatrixXd& sqrt_p_aug){
  Xsig_aug = MatrixXd(n_a_, 2*n_a_ + 1);
  Xsig_aug.col(0) = x_aug;
  for(int i = 0; i < n_a_; ++i)
  {
    Xsig_aug.col(i+1)      = x_aug + sqrt(lambda_+n_a_)*sqrt_p_aug.col(i);
    Xsig_aug.col(i+1+n_a_) = x_aug - sqrt(lambda_+n_a_)*sqrt_p_aug.col(i);
  }
}

void UKF::DoPredictSigmaPoints(const Eigen::MatrixXd& xSigMat, double dt){
  
  int cols = 2*n_a_ + 1;
  double px, py, v, yaw, yawd, nu_a, nu_yawdd;
  for(int i = 0; i < cols; ++i)
  {
    //Getting variable for each calculated sigma point in each column
    px = xSigMat(0, i);
    py = xSigMat(1, i);
    v = xSigMat(2, i);
    yaw = xSigMat(3, i);
    yawd = xSigMat(4, i);
    nu_a = xSigMat(5, i);
    nu_yawdd = xSigMat(6, i);
    //defining and caclulating the predicted sigma points X(k+1|k)
    double px_p, py_p, v_p, yaw_p, yawd_p;
    double sin_yaw = sin(yaw);
    double cos_yaw = cos(yaw);
    if(fabs(yawd)>0.001)
    {
      px_p = px + (v/yawd)*(sin(yaw + yawd*dt) - sin_yaw);
      py_p = py + (v/yawd)*(-cos(yaw + yawd*dt) + cos_yaw);
    }
    else
    {
      px_p = px + v*dt*cos_yaw;
      py_p = py + v*dt*sin_yaw;
    }
    v_p = v;
    yaw_p = yaw + yawd*dt;
    yawd_p = yawd;
    //Adding noise
    px_p += 0.5*dt*dt*cos_yaw*nu_a;
    py_p += 0.5*dt*dt*sin_yaw*nu_a;
    v_p += dt*nu_a;
    yaw_p += 0.5*dt*dt*nu_yawdd;
    yawd_p += dt*nu_yawdd;
    //assign the variables to the matrix column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }
}

void UKF::PredictMeanAndCovar(){
  x_.fill(0.0);
  P_.fill(0.0);
  int colsNum = 2*n_a_ + 1;
  for(int i = 0; i < colsNum; ++i)
  {
    x_ += weights_(i)*Xsig_pred_.col(i);
  }
  for(int i = 0; i < colsNum; ++i)
  {
    VectorXd xDiff = Xsig_pred_.col(i) - x_;
    while(xDiff(3) > M_PI) xDiff(3) -= 2.*M_PI;
    while(xDiff(3) < -M_PI) xDiff(3) += 2.*M_PI;

    P_ += weights_(i)*xDiff*xDiff.transpose();
  }
}