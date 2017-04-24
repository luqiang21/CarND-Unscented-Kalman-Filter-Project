#include "ukf.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 5;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 5;

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
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
 //  time_us_ =
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_aug_;
  NIS_radar_ = 0.;
  NIS_laser_ = 0.;
  weights_ = VectorXd(2*n_aug_+1);                 // weights_ of sigma points

  is_initialized_ = false;
  

  // weights_ for sigma points
  // set vector for weights_
  VectorXd weights_ = VectorXd(2*n_aug_+1);
  double weight_0 = lambda_/(lambda_+n_aug_);
  weights_(0) = weight_0;
  for (int i=1; i<2*n_aug_+1; i++) {
    double weight = 0.5/(n_aug_+lambda_);
    weights_(i) = weight;
  }
  
  
//  x_ << 1, 1, 1, 1, 1;
  
 
  n_z_radar_ = 3;
  
  n_z_lidar_ = 2;
  
  Xsig_pred_ = MatrixXd(n_x_, 2*n_aug_ + 1);          // Predicted sigma points
  Xsig_aug_  = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  
  z_pred_radar_ = VectorXd(n_z_radar_);
  S_radar_ = MatrixXd(n_z_radar_,n_z_radar_);
  Zsig_radar_ = MatrixXd(n_z_radar_, 2 * n_aug_ + 1);
  
  z_pred_lidar_ = VectorXd(n_z_lidar_);
  S_lidar_ = MatrixXd(n_z_lidar_,n_z_lidar_);
  Zsig_lidar_ = MatrixXd(n_z_lidar_, 2 * n_aug_ + 1);

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    // first measurement
    cout << "UKF: 1st" << endl;
    P_ << 1, 0, 0, 0,0,
    0, 1, 0, 0, 0,
    0, 0, 1, 0, 0,
    0, 0, 0, 1, 0,
    0, 0, 0, 0, 1;

    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      /**
       Convert radar from polar to cartesian coordinates and initialize state.
       */
      double px = meas_package.raw_measurements_[0]*cos(meas_package.raw_measurements_[1]);
      double py = meas_package.raw_measurements_[0]*sin(meas_package.raw_measurements_[1]);
      
      x_ << px, py, abs(meas_package.raw_measurements_[2]), 0., 0.;
      
      cout << x_<<endl;
      
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
       Initialize state.
       */
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0., 0., 0.;
      cout << x_<<endl;

    }
    
    if(x_[0] < 0.000001 && x_[1] < 0.000001){
      cout<<"zero input"<<endl;
      x_ << 0.001, 0.001, 0., 0., 0.;
    }
    // done initializing, no need to predict or update
    is_initialized_ = true;
    time_us_ = meas_package.timestamp_;
    
    return;
  }
  
  float delta_t = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
  
  Prediction(delta_t);
  
  if (use_radar_ && meas_package.sensor_type_ == MeasurementPackage::RADAR) {
    UpdateRadar(meas_package);
  }
  else if(use_laser_ && meas_package.sensor_type_ == MeasurementPackage::LASER){
    UpdateLidar(meas_package);

  }else{
    cout<< " sensor unknown or not used" << endl;
  }
  
}
 
  

void UKF::AugmentedSigmaPoints() {
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);
  
  //create augmented mean state
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;
  
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  
  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug_.col(0)  = x_aug;
  for (int i = 0; i< n_aug_; i++)
  {
    Xsig_aug_.col(i+1)       = x_aug + sqrt(lambda_+n_aug_) * L.col(i);
    Xsig_aug_.col(i+1+n_aug_) = x_aug - sqrt(lambda_+n_aug_) * L.col(i);
  }
  
  
  //print result
  std::cout << "Xsig_aug = " << std::endl << Xsig_aug_ << std::endl;
  
  //write result
//  *Xsig_out = Xsig_aug_;
  
}
  
  
void UKF::SigmaPointPrediction(float delta_t) {
//void UKF::SigmaPointPrediction(float delta_t) {
  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  
  for(int i=0; i < (2*n_aug_ + 1); i ++ ){
    VectorXd x_k = Xsig_aug_.col(i);
    
    float px = x_k(0);
    float py = x_k(1);
    float v = x_k(2);
    float psi = x_k(3);
    float psi_dot = x_k(4);
    float noise_a = x_k(5);
    float noise_yawdd = x_k(6);
    float delta_t_2 = delta_t * delta_t;
    
    if(psi_dot < 0.0001){
      px = px + v * cos(psi) * delta_t + 1/2. * delta_t_2 * cos(psi) * noise_a;
      py = py + v * sin(psi) * delta_t + 1/2. * delta_t_2 * sin(psi) * noise_a;
    }else{
      px = px + v / psi_dot * (sin(psi + psi_dot * delta_t) - sin(psi)) + 1/2. * delta_t_2 * cos(psi) * noise_a;
      py = py + v / psi_dot * (-cos(psi + psi_dot * delta_t) + cos(psi)) + 1/2. * delta_t_2 * sin(psi) * noise_a;
      
    }
    
    v = v + delta_t * noise_a;
    psi = psi + psi_dot * delta_t + 1/2. * delta_t_2 * noise_yawdd;
    psi_dot = psi_dot + delta_t * noise_yawdd;
    
    Xsig_pred_.col(i) << px, py, v, psi, psi_dot;
  }
  
  
  //print result
  std::cout << "Xsig_pred = " << std::endl << Xsig_pred_ << std::endl;
  
  //write result
//  *Xsig_out = Xsig_pred_;
  
}
 

void UKF::PredictMeanAndCovariance() {
//  //create vector for predicted state
//  VectorXd x = VectorXd(n_x_);
//  
//  //create covariance matrix for prediction
//  MatrixXd P = MatrixXd(n_x_, n_x_);
  
  x_.fill(0.0);
  P_.fill(0.0);
    //predict state mean
  for(int i=0; i < (2*n_aug_+1); i++){
    VectorXd x_k = Xsig_pred_.col(i);
    x_ += weights_[i] * x_k;
  }
  cout << "I am here"<<Xsig_pred_<<endl;
  
  //predict state covariance matrix
  for(int i=0; i < (2*n_aug_+1); i++){
    VectorXd x_k = Xsig_pred_.col(i);
    x_k -= x_;
    
    while (x_k(3)> M_PI) x_k(3)-=2.*M_PI;
    while (x_k(3)<-M_PI) x_k(3)+=2.*M_PI;
    
    P_ += weights_[i] * x_k * x_k.transpose();
    
  }
  cout << "I am here1"<<endl;

  //print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x_ << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P_ << std::endl;
  
  //write result
//  *x_out = x;
//  *P_out = P;
}
  
  
void UKF::PredictRadarMeasurement() {
//  //mean predicted measurement
//  VectorXd z_pred = VectorXd(n_z_radar_);
//  
//  //measurement covariance matrix S
//  MatrixXd S = MatrixXd(n_z_radar_,n_z_radar_);
  
  //transform sigma points into measurement space
  for(int i=0; i < (2*n_aug_ + 1); i++){
    double px = Xsig_pred_.col(i)(0);
    double py = Xsig_pred_.col(i)(1);
    double v  = Xsig_pred_.col(i)(2);
    double psi = Xsig_pred_.col(i)(3);
    
    double rou = sqrt(px*px + py*py);
    double phi = atan2(py, px);
    double rou_dot;
    
    if(rou < 0.0000001){
      rou_dot = 0.0;
    }else{
      rou_dot = 1/rou * (px*cos(psi)*v + py*sin(psi)*v);
    }
    
    Zsig_radar_.col(i) << rou , phi , rou_dot;
    
  }
  
  //calculate mean predicted measurement
  z_pred_radar_.fill(0.);
  for(int i=0; i < (2*n_aug_ + 1); i++){
    float weight = weights_(i);
    z_pred_radar_ += weight * Zsig_radar_.col(i);
  }
  
  S_radar_.fill(0);
  //calculate measurement covariance matrix S
  for(int i=0; i < (2*n_aug_ + 1); i++){
    float weight = weights_(i);
    VectorXd z_diff = Zsig_radar_.col(i) - z_pred_radar_;
    S_radar_ += weight * z_diff * z_diff.transpose();
  }
  
  MatrixXd R = MatrixXd(n_z_radar_, n_z_radar_);
  R.fill(0);
  R(0,0) = std_radr_ * std_radr_;
  R(1,1) = std_radphi_ * std_radphi_;
  R(2,2) = std_radrd_ * std_radrd_;
  
  S_radar_ += R;
  
    //print result
  std::cout << "z_pred: " << std::endl << z_pred_radar_ << std::endl;
  std::cout << "S: " << std::endl << S_radar_ << std::endl;
  
  //write result
//  *z_out = z_pred;
//  *S_out = S;
}


//void UKF::UpdateRadarState(VectorXd* x_out, MatrixXd* P_out, MeasurementPackage meas_package) {
void UKF::UpdateRadarState (MeasurementPackage meas_package) {
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z_radar_);
  z << meas_package.raw_measurements_;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_radar_);
  
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0; i < 2*n_aug_ + 1; i ++){
    
    VectorXd z_diff = Zsig_radar_.col(i) - z_pred_radar_;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    
    double weight = weights_(i);
    Tc += weight * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S_radar_.inverse();
  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred_radar_;
  
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S_radar_ * K.transpose();
  
  //print result
  std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
  
  //calculate NIS
  NIS_radar_ = z_diff.transpose() * S_radar_.inverse() * z_diff;
  
//  //write result
//  *x_out = x_;
//  *P_out = P_;
}


void UKF::PredictLidarMeasurement() {
  //  //mean predicted measurement
  //  VectorXd z_pred = VectorXd(n_z_radar_);
  //
  //  //measurement covariance matrix S
  //  MatrixXd S = MatrixXd(n_z_radar_,n_z_radar_);
  
  
  //transform sigma points into measurement space
  for(int i=0; i < (2*n_aug_ + 1); i++){
    double px = Xsig_pred_.col(i)(0);
    double py = Xsig_pred_.col(i)(1);
    
    Zsig_lidar_.col(i) << px, py;
    
  }
  
  //calculate mean predicted measurement
  z_pred_lidar_.fill(0.);
  for(int i=0; i < (2*n_aug_ + 1); i++){
    float weight = weights_(i);
    z_pred_lidar_ += weight * Zsig_lidar_.col(i);
  }
  
  S_lidar_.fill(0);
  //calculate measurement covariance matrix S
  for(int i=0; i < (2*n_aug_ + 1); i++){
    float weight = weights_(i);
    VectorXd z_diff = Zsig_lidar_.col(i) - z_pred_lidar_;
    S_lidar_ += weight * z_diff * z_diff.transpose();
  }
  
  MatrixXd R = MatrixXd(n_z_lidar_, n_z_lidar_);
  R.fill(0);
  R(0,0) = std_laspx_ * std_laspx_;
  R(1,1) = std_laspy_ * std_laspy_;
  
  S_lidar_ += R;
  
  //print result
  std::cout << "z_pred: " << std::endl << z_pred_lidar_ << std::endl;
  std::cout << "S: " << std::endl << S_lidar_ << std::endl;
  
  //write result
  //  *z_out = z_pred;
  //  *S_out = S;
}


void UKF::UpdateLidarState (MeasurementPackage meas_package) {
  
  
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z_lidar_);
  z << meas_package.raw_measurements_;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z_lidar_);
  
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0; i < 2*n_aug_ + 1; i ++){
    
    VectorXd z_diff = Zsig_lidar_.col(i) - z_pred_lidar_;
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    
    double weight = weights_(i);
    Tc += weight * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S_lidar_.inverse();
  
  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred_lidar_;
  
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S_lidar_ * K.transpose();
  
  //print result
  std::cout << "Updated state x: " << std::endl << x_ << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P_ << std::endl;
  
  //calculate NIS
  NIS_laser_ = z_diff.transpose() * S_lidar_.inverse() * z_diff;
  
  //  //write result
  //  *x_out = x_;
  //  *P_out = P_;
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  
  AugmentedSigmaPoints();
  SigmaPointPrediction(delta_t);
  PredictMeanAndCovariance();
  
#if UKF_DEBUG
  cout << "prediction completed!" << endl;
#endif
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  
  PredictLidarMeasurement();
  UpdateLidarState(meas_package);

  
  //calculate the NIS
//  NIS_laser_ = z_diff_mean.transpose() * S_inv * z_diff_mean;
#if UKF_DEBUG
  std::cout << "NIS lidar: " << std::endl << NIS_laser_ << std::endl;
#endif
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  
  PredictRadarMeasurement();
  UpdateRadarState( meas_package);
  
  //display the NIS
  
#if UKF_DEBUG
  std::cout << "NIS radar: " << std::endl << NIS_radar_ << std::endl;
#endif
}
