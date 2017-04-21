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
  is_initialized_ = false;
  Xsig_pred = MatrixXd(n_x, 2*n_aug+1);          // Predicted sigma points
//  time_us_ = 
  weights = VectorXd(2*n_aug+1);                 // Weights of sigma points
  n_x_ = 5;
  n_aug_ = 7;
  lambda_ = 3 - n_x;
  NIS_radar_ = 0.;
  NIS_laser_ = 0.;
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
    cout << "EKF: " << endl;
    ekf_.x_ = VectorXd(4);
    ekf_.x_ << 1, 1, 1, 1;
    
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
       Convert radar from polar to cartesian coordinates and initialize state.
       */
      ekf_.x_ << measurement_pack.raw_measurements_[0]*cos(measurement_pack.raw_measurements_[1]), measurement_pack.raw_measurements_[0]*sin(measurement_pack.raw_measurements_[1]),
      measurement_pack.raw_measurements_[2]*cos(measurement_pack.raw_measurements_[1]),
      measurement_pack.raw_measurements_[2]*sin(measurement_pack.raw_measurements_[1]);
      
      cout << ekf_.x_<<endl;
      
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
       Initialize state.
       */
      ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
    
}
  
  
void UKF::AugmentedSigmaPoints(MatrixXd* Xsig_out) {
  
  //set state dimension
  int n_x = 5;
  
  //set augmented dimension
  int n_aug = 7;
  
  //Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a = 0.2;
  
  //Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd = 0.2;
  
  //define spreading parameter
  double lambda = 3 - n_aug;
  
  //set example state
  VectorXd x = VectorXd(n_x);
  x <<   5.7441,
  1.3800,
  2.2049,
  0.5015,
  0.3528;
  
  //create example covariance matrix
  MatrixXd P = MatrixXd(n_x, n_x);
  P <<     0.0043,   -0.0013,    0.0030,   -0.0022,   -0.0020,
  -0.0013,    0.0077,    0.0011,    0.0071,    0.0060,
  0.0030,    0.0011,    0.0054,    0.0007,    0.0008,
  -0.0022,    0.0071,    0.0007,    0.0098,    0.0100,
  -0.0020,    0.0060,    0.0008,    0.0100,    0.0123;
  
  //create augmented mean vector
  VectorXd x_aug = VectorXd(7);
  
  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(7, 7);
  
  //create sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  
  
  //create augmented mean state
  x_aug.head(5) = x;
  x_aug(5) = 0;
  x_aug(6) = 0;
  
  //create augmented covariance matrix
  P_aug.fill(0.0);
  P_aug.topLeftCorner(5,5) = P;
  P_aug(5,5) = std_a*std_a;
  P_aug(6,6) = std_yawdd*std_yawdd;
  
  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();
  
  //create augmented sigma points
  Xsig_aug.col(0)  = x_aug;
  for (int i = 0; i< n_aug; i++)
  {
    Xsig_aug.col(i+1)       = x_aug + sqrt(lambda+n_aug) * L.col(i);
    Xsig_aug.col(i+1+n_aug) = x_aug - sqrt(lambda+n_aug) * L.col(i);
  }
  
  
  //print result
  std::cout << "Xsig_aug = " << std::endl << Xsig_aug << std::endl;
  
  //write result
  *Xsig_out = Xsig_aug;
  
  
}
  
  
void UKF::SigmaPointPrediction(MatrixXd* Xsig_out) {
  
  //set state dimension
  int n_x = 5;
  
  //set augmented dimension
  int n_aug = 7;
  
  //create example sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug, 2 * n_aug + 1);
  Xsig_aug <<
  5.7441,  5.85768,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.63052,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,   5.7441,
  1.38,  1.34566,  1.52806,     1.38,     1.38,     1.38,     1.38,     1.38,   1.41434,  1.23194,     1.38,     1.38,     1.38,     1.38,     1.38,
  2.2049,  2.28414,  2.24557,  2.29582,   2.2049,   2.2049,   2.2049,   2.2049,   2.12566,  2.16423,  2.11398,   2.2049,   2.2049,   2.2049,   2.2049,
  0.5015,  0.44339, 0.631886, 0.516923, 0.595227,   0.5015,   0.5015,   0.5015,   0.55961, 0.371114, 0.486077, 0.407773,   0.5015,   0.5015,   0.5015,
  0.3528, 0.299973, 0.462123, 0.376339,  0.48417, 0.418721,   0.3528,   0.3528,  0.405627, 0.243477, 0.329261,  0.22143, 0.286879,   0.3528,   0.3528,
  0,        0,        0,        0,        0,        0,  0.34641,        0,         0,        0,        0,        0,        0, -0.34641,        0,
  0,        0,        0,        0,        0,        0,        0,  0.34641,         0,        0,        0,        0,        0,        0, -0.34641;
  
  //create matrix with predicted sigma points as columns
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  
  double delta_t = 0.1; //time diff in sec
  
  //predict sigma points
  //avoid division by zero
  //write predicted sigma points into right column
  
  for(int i=0; i < (2*n_aug + 1); i ++ ){
    VectorXd x_k = Xsig_aug.col(i);
    
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
    
    Xsig_pred.col(i) << px, py, v, psi, psi_dot;
  }
  
  
  //print result
  std::cout << "Xsig_pred = " << std::endl << Xsig_pred << std::endl;
  
  //write result
  *Xsig_out = Xsig_pred;
  
}
 

  
void UKF::PredictMeanAndCovariance(VectorXd* x_out, MatrixXd* P_out) {
  
  //set state dimension
  int n_x = 5;
  
  //set augmented dimension
  int n_aug = 7;
  
  //define spreading parameter
  double lambda = 3 - n_aug;
  
  //create example matrix with predicted sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred <<
  5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
  1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
  2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
  0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
  0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
  
  //create vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  
  //create vector for predicted state
  VectorXd x = VectorXd(n_x);
  
  //create covariance matrix for prediction
  MatrixXd P = MatrixXd(n_x, n_x);
  
  //set weights
  double weight;
  for(int i=0; i < (2*n_aug+1); i++){
    if(i == 0){
      weight = lambda / (lambda + n_aug);
    }else{
      weight = 1 / 2. / (lambda + n_aug);
    }
    weights[i] = weight;
  }
  //predict state mean
  for(int i=0; i < (2*n_aug+1); i++){
    VectorXd x_k = Xsig_pred.col(i);
    x += weights[i] * x_k;
  }
  //predict state covariance matrix
  for(int i=0; i < (2*n_aug+1); i++){
    VectorXd x_k = Xsig_pred.col(i);
    x_k -= x;
    
    while (x_k(3)> M_PI) x_k(3)-=2.*M_PI;
    while (x_k(3)<-M_PI) x_k(3)+=2.*M_PI;
    
    P += weights[i] * x_k * x_k.transpose();
    
  }
  
  
  //print result
  std::cout << "Predicted state" << std::endl;
  std::cout << x << std::endl;
  std::cout << "Predicted covariance matrix" << std::endl;
  std::cout << P << std::endl;
  
  //write result
  *x_out = x;
  *P_out = P;
}
  
  
void UKF::PredictRadarMeasurement(VectorXd* z_out, MatrixXd* S_out) {
  
  //set state dimension
  int n_x = 5;
  
  //set augmented dimension
  int n_aug = 7;
  
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  
  //define spreading parameter
  double lambda = 3 - n_aug;
  
  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }
  
  //radar measurement noise standard deviation radius in m
  double std_radr = 0.3;
  
  //radar measurement noise standard deviation angle in rad
  double std_radphi = 0.0175;
  
  //radar measurement noise standard deviation radius change in m/s
  double std_radrd = 0.1;
  
  //create example matrix with predicted sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred <<
  5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
  1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
  2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
  0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
  0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
  
  //create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
  
  //mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  
  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z,n_z);
  
  
  //transform sigma points into measurement space
  for(int i=0; i < (2*n_aug + 1); i++){
    double px = Xsig_pred.col(i)(0);
    double py = Xsig_pred.col(i)(1);
    double v  = Xsig_pred.col(i)(2);
    double psi = Xsig_pred.col(i)(3);
    
    double rou = sqrt(px*px + py*py);
    double phi = atan2(py, px);
    double rou_dot;
    if(rou < 0.0000001){
      rou_dot = 0.0;
    }else{
      rou_dot = 1/rou * (px*cos(psi)*v + py*sin(psi)*v);
    }
    
    Zsig.col(i) << rou , phi , rou_dot;
    
  }
  //calculate mean predicted measurement
  z_pred.fill(0.);
  for(int i=0; i < (2*n_aug + 1); i++){
    float weight = weights(i);
    z_pred += weight * Zsig.col(i);
  }
  
  S.fill(0);
  //calculate measurement covariance matrix S
  for(int i=0; i < (2*n_aug + 1); i++){
    float weight = weights(i);
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weight * z_diff * z_diff.transpose();
  }
  
  MatrixXd R = MatrixXd(n_z, n_z);
  R.fill(0);
  R(0,0) = std_radr * std_radr;
  R(1,1) = std_radphi * std_radphi;
  R(2,2) = std_radrd * std_radrd;
  
  S += R;
  
    //print result
  std::cout << "z_pred: " << std::endl << z_pred << std::endl;
  std::cout << "S: " << std::endl << S << std::endl;
  
  //write result
  *z_out = z_pred;
  *S_out = S;
}


void UKF::UpdateState(VectorXd* x_out, MatrixXd* P_out) {
  
  //set state dimension
  int n_x = 5;
  
  //set augmented dimension
  int n_aug = 7;
  
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  
  //define spreading parameter
  double lambda = 3 - n_aug;
  
  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }
  
  //create example matrix with predicted sigma points in state space
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred <<
  5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
  1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
  2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
  0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
  0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
  
  //create example vector for predicted state mean
  VectorXd x = VectorXd(n_x);
  x <<
  5.93637,
  1.49035,
  2.20528,
  0.536853,
  0.353577;
  
  //create example matrix for predicted state covariance
  MatrixXd P = MatrixXd(n_x,n_x);
  P <<
  0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
  -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
  0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
  -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
  -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;
  
  //create example matrix with sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
  Zsig <<
  6.1190,  6.2334,  6.1531,  6.1283,  6.1143,  6.1190,  6.1221,  6.1190,  6.0079,  6.0883,  6.1125,  6.1248,  6.1190,  6.1188,  6.12057,
  0.24428,  0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
  2.1104,  2.2188,  2.0639,   2.187,  2.0341,  2.1061,  2.1450,  2.1092,  2.0016,   2.129,  2.0346,  2.1651,  2.1145,  2.0786,  2.11295;
  
  //create example vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred <<
  6.12155,
  0.245993,
  2.10313;
  
  //create example matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z,n_z);
  S <<
  0.0946171, -0.000139448,   0.00407016,
  -0.000139448,  0.000617548, -0.000770652,
  0.00407016, -0.000770652,    0.0180917;
  
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z <<
  5.9214,   //rho in m
  0.2187,   //phi in rad
  2.0062;   //rho_dot in m/s
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);
    
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for(int i=0; i < 2*n_aug + 1; i ++){
    
    VectorXd z_diff = Zsig.col(i) - z_pred;
    VectorXd x_diff = Xsig_pred.col(i) - x;
    
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    //angle normalization
    while (x_diff(1)> M_PI) x_diff(1)-=2.*M_PI;
    while (x_diff(1)<-M_PI) x_diff(1)+=2.*M_PI;
    
    double weight = weights(i);
    Tc += weight * x_diff * z_diff.transpose();
  }
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  //update state mean and covariance matrix
  VectorXd z_diff = z - z_pred;
  
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  x = x + K * z_diff;
  P = P - K * S * K.transpose();
  
  //print result
  std::cout << "Updated state x: " << std::endl << x << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P << std::endl;
  
  //write result
  *x_out = x;
  *P_out = P;
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
  
  GenerateSigmaPoints(MatrixXd* Xsig_out);
  AugmentedSigmaPoints(MatrixXd* Xsig_out);
  SigmaPointPrediction(MatrixXd* Xsig_out);
  PredictMeanAndCovariance(VectorXd* x_pred, MatrixXd* P_pred);

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
  
  //set state dimension
  int n_x = 5;
  
  //set augmented dimension
  int n_aug = 7;
  
  //set measurement dimension, radar can measure r, phi, and r_dot
  int n_z = 3;
  
  //define spreading parameter
  double lambda = 3 - n_aug;
  
  //set vector for weights
  VectorXd weights = VectorXd(2*n_aug+1);
  double weight_0 = lambda/(lambda+n_aug);
  weights(0) = weight_0;
  for (int i=1; i<2*n_aug+1; i++) {  //2n+1 weights
    double weight = 0.5/(n_aug+lambda);
    weights(i) = weight;
  }
  
  //create example matrix with predicted sigma points
  MatrixXd Xsig_pred = MatrixXd(n_x, 2 * n_aug + 1);
  Xsig_pred <<
  5.9374,  6.0640,   5.925,  5.9436,  5.9266,  5.9374,  5.9389,  5.9374,  5.8106,  5.9457,  5.9310,  5.9465,  5.9374,  5.9359,  5.93744,
  1.48,  1.4436,   1.660,  1.4934,  1.5036,    1.48,  1.4868,    1.48,  1.5271,  1.3104,  1.4787,  1.4674,    1.48,  1.4851,    1.486,
  2.204,  2.2841,  2.2455,  2.2958,   2.204,   2.204,  2.2395,   2.204,  2.1256,  2.1642,  2.1139,   2.204,   2.204,  2.1702,   2.2049,
  0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337,  0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188,  0.5367, 0.535048,
  0.352, 0.29997, 0.46212, 0.37633,  0.4841, 0.41872,   0.352, 0.38744, 0.40562, 0.24347, 0.32926,  0.2214, 0.28687,   0.352, 0.318159;
  
  //create example vector for predicted state mean
  VectorXd x = VectorXd(n_x);
  x <<
  5.93637,
  1.49035,
  2.20528,
  0.536853,
  0.353577;
  
  //create example matrix for predicted state covariance
  MatrixXd P = MatrixXd(n_x,n_x);
  P <<
  0.0054342,  -0.002405,  0.0034157, -0.0034819, -0.00299378,
  -0.002405,    0.01084,   0.001492,  0.0098018,  0.00791091,
  0.0034157,   0.001492,  0.0058012, 0.00077863, 0.000792973,
  -0.0034819,  0.0098018, 0.00077863,   0.011923,   0.0112491,
  -0.0029937,  0.0079109, 0.00079297,   0.011249,   0.0126972;
  
  //create example matrix with sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug + 1);
  Zsig <<
  6.1190,  6.2334,  6.1531,  6.1283,  6.1143,  6.1190,  6.1221,  6.1190,  6.0079,  6.0883,  6.1125,  6.1248,  6.1190,  6.1188,  6.12057,
  0.24428,  0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,
  2.1104,  2.2188,  2.0639,   2.187,  2.0341,  2.1061,  2.1450,  2.1092,  2.0016,   2.129,  2.0346,  2.1651,  2.1145,  2.0786,  2.11295;
  
  //create example vector for mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred <<
  6.12155,
  0.245993,
  2.10313;
  
  //create example matrix for predicted measurement covariance
  MatrixXd S = MatrixXd(n_z,n_z);
  S <<
  0.0946171, -0.000139448,   0.00407016,
  -0.000139448,  0.000617548, -0.000770652,
  0.00407016, -0.000770652,    0.0180917;
  
  //create example vector for incoming radar measurement
  VectorXd z = VectorXd(n_z);
  z <<
  5.9214,
  0.2187,
  2.0062;
  
  //create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x, n_z);
  
  /*******************************************************************************
   * Student part begin
   ******************************************************************************/
  
  //calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < 2 * n_aug + 1; i++) {  //2n+1 simga points
    
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    //angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
    
    // state difference
    VectorXd x_diff = Xsig_pred.col(i) - x;
    //angle normalization
    while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
    while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;
    
    Tc = Tc + weights(i) * x_diff * z_diff.transpose();
  }
  
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  
  //residual
  VectorXd z_diff = z - z_pred;
  
  //angle normalization
  while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
  while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;
  
  //update state mean and covariance matrix
  x = x + K * z_diff;
  P = P - K*S*K.transpose();
  
  /*******************************************************************************
   * Student part end
   ******************************************************************************/
  
  //print result
  std::cout << "Updated state x: " << std::endl << x << std::endl;
  std::cout << "Updated state covariance P: " << std::endl << P << std::endl;
  
  //write result
  *x_out = x;
  *P_out = P;
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
}
