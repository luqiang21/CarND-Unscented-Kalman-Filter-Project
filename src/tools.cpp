#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  /**
  TODO:
    * Calculate the RMSE here.
  */
  VectorXd RMSE(4);
  RMSE << 0,0,0,0;
  
  // check validity of the inputs;
  if(estimations.size() != ground_truth.size() || estimations.size() < 1){
    std::cout << "invalid ground truth, please check" << std::endl;
    return RMSE;
  }
  
  // accumulating squared residuals
  for(int i=0; i < estimations.size(); ++i){
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    RMSE += residual;
  }
  
  // mean and squared root
  RMSE = RMSE / estimations.size();
  RMSE = RMSE.array().sqrt();
  
  return RMSE;

}
