#include <iostream>
#include "ukf.h"
#include "Eigen/Dense"

using Eigen::MatrixXd;
using Eigen::VectorXd;

/**
 * Initializes Unscented Kalman filter
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  n_x_ = 5;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd(n_x_, n_x_);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 0.2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.2;
  
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
   n_aug_ = 7;

   lambda_ = 3 - n_aug_;


   weights_ = Eigen::VectorXd(2*n_aug_+1);

    weights_(0) = lambda_/(lambda_+n_aug_);
    for (int i=1; i < 2 *n_aug_+1; ++i){
        weights_(i) = 0.5/(n_aug_+lambda_); //0.5379
    }  // 2n+1 weights

    Xsig_pred = MatrixXd(n_x_, 2 * n_aug_ + 1);





}

UKF::~UKF() {}

Eigen::MatrixXd UKF::AugmentedSigmaPoints() {

    // create augmented mean vector
    VectorXd x_aug = VectorXd(UKF::n_aug_);

    // create augmented state covariance
    MatrixXd P_aug = MatrixXd(UKF::n_aug_, UKF::n_aug_);

    // create sigma point matrix
    MatrixXd Xsig_aug = MatrixXd(UKF::n_aug_, 2 * UKF::n_aug_ + 1);


    // create augmented mean state
    x_aug.head(UKF::n_x_ ) = UKF::x_;
    x_aug(UKF::n_x_ ) = 0;
    x_aug(UKF::n_x_ + 1) = 0;

    // create augmented covariance matrix
    P_aug.fill(0.0);
    P_aug.topLeftCorner(UKF::n_x_,UKF::n_x_) = UKF::P_;
    P_aug(UKF::n_x_,UKF::n_x_) = UKF::std_a_*UKF::std_a_;
    P_aug(UKF::n_x_+1,UKF::n_x_+1) = UKF::std_yawdd_*UKF::std_yawdd_;

    // create square root matrix
    MatrixXd L = P_aug.llt().matrixL();

    // create augmented sigma points
    Xsig_aug.col(0)  = x_aug;
    for (int i = 0; i< UKF::n_aug_; ++i) {
        Xsig_aug.col(i+1)             = x_aug + sqrt(UKF::lambda_+UKF::n_aug_) * L.col(i);
        Xsig_aug.col(i+1+UKF::n_aug_) = x_aug - sqrt(UKF::lambda_+UKF::n_aug_) * L.col(i);
    }
    // write result
    return Xsig_aug;
}

void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Make sure you switch between lidar and radar
   * measurements.
   */
   if(!is_initialized_){
       if(meas_package.sensor_type_== MeasurementPackage::RADAR && use_radar_){
           double rho = meas_package.raw_measurements_(0);
           double phi = meas_package.raw_measurements_(1);
           double rho_dot = meas_package.raw_measurements_(2);
           x_.fill(0.0);
           x_(0) = rho * cos(phi);
           x_(1) = rho * sin(phi);
           x_(2) = rho_dot;
           P_.fill(0);
           time_us_ = meas_package.timestamp_;
           is_initialized_ = true;

       }
       if(meas_package.sensor_type_== MeasurementPackage::LASER && use_laser_){
           double X = meas_package.raw_measurements_(0);
           double Y = meas_package.raw_measurements_(1);
           x_.fill(0.0);
           x_(0) = X;
           x_(1) = Y;
           P_.fill(0);
           time_us_ = meas_package.timestamp_;
           is_initialized_ = true;

       }
   }

   if(meas_package.sensor_type_== MeasurementPackage::LASER && use_laser_ && is_initialized_){
       Prediction(meas_package.timestamp_ - time_us_);
       UpdateLidar(meas_package);
       time_us_ = meas_package.timestamp_;
   }
    if(meas_package.sensor_type_== MeasurementPackage::RADAR && use_radar_ && is_initialized_){
        Prediction(meas_package.timestamp_ - time_us_);
        UpdateRadar(meas_package);
        time_us_ = meas_package.timestamp_;
    }
}

void UKF::Prediction(double delta_t) {

    delta_t /= 1000000;
    //std::cout << delta_t << std::endl;

    MatrixXd Xsig_aug = AugmentedSigmaPoints();
    //std::cout << x_ << std::endl;
    //std::cout << weights_.sum() << std::endl;



    for (int i = 0; i< 2*n_aug_+1; ++i) {
        // extract values for better readability
        double p_x = Xsig_aug(0,i);
        double p_y = Xsig_aug(1,i);
        double v = Xsig_aug(2,i);
        double yaw = Xsig_aug(3,i);
        double yawd = Xsig_aug(4,i);
        double nu_a = Xsig_aug(5,i);
        double nu_yawdd = Xsig_aug(6,i);

        // predicted state values
        double px_p, py_p;

        // avoid division by zero
        if (fabs(yawd) > 0.001) {
            px_p = p_x + v/yawd * ( sin (yaw + yawd*delta_t) - sin(yaw));
            py_p = p_y + v/yawd * ( cos(yaw) - cos(yaw+yawd*delta_t) );
        } else {
            px_p = p_x + v*delta_t*cos(yaw);
            py_p = p_y + v*delta_t*sin(yaw);
        }

        double v_p = v;
        double yaw_p = yaw + yawd*delta_t;
        double yawd_p = yawd;

        // add noise
        px_p = px_p + 0.5*nu_a*delta_t*delta_t * cos(yaw);
        py_p = py_p + 0.5*nu_a*delta_t*delta_t * sin(yaw);
        v_p = v_p + nu_a*delta_t;

        yaw_p = yaw_p + 0.5*nu_yawdd*delta_t*delta_t;
        yawd_p = yawd_p + nu_yawdd*delta_t;

        // write predicted sigma point into right column
        Xsig_pred(0,i) = px_p;
        Xsig_pred(1,i) = py_p;
        Xsig_pred(2,i) = v_p;
        Xsig_pred(3,i) = yaw_p;
        Xsig_pred(4,i) = yawd_p;
    }

    x_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
        x_ = x_ + weights_(i) * Xsig_pred.col(i);
    }
    P_.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // iterate over sigma points
        // state difference
        VectorXd x_diff = Xsig_pred.col(i) - x_;
        // angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        P_ = P_ + weights_(i) * x_diff * x_diff.transpose() ;
    }

}

void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use lidar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the lidar NIS, if desired.
   */
    int n_z = 2;

    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);

    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
        // measurement model
        Zsig(0,i) = Xsig_pred(0,i);                   // x
        Zsig(1,i) = Xsig_pred(1,i);                   // y

    }

    // mean predicted measurement
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; ++i) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    // innovation covariance matrix S
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        // angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    // add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<  std_laspx_*std_laspx_, 0,
            0, std_laspy_*std_laspy_;
    S = S + R;

    MatrixXd Tc = MatrixXd(n_x_, n_z);

    // calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        // angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        // state difference
        VectorXd x_diff = Xsig_pred.col(i) - x_;
        // angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    // Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    // residual
    VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // update state mean and covariance matrix

    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();
}

void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
   * TODO: Complete this function! Use radar data to update the belief 
   * about the object's position. Modify the state vector, x_, and 
   * covariance, P_.
   * You can also calculate the radar NIS, if desired.
   */

   int n_z = 3;

    MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

    // mean predicted measurement
    VectorXd z_pred = VectorXd(n_z);

    // measurement covariance matrix S
    MatrixXd S = MatrixXd(n_z,n_z);

    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
        // extract values for better readability
        double p_x = Xsig_pred(0,i);
        double p_y = Xsig_pred(1,i);
        double v  = Xsig_pred(2,i);
        double yaw = Xsig_pred(3,i);

        double v1 = cos(yaw)*v;
        double v2 = sin(yaw)*v;

        // measurement model
        Zsig(0,i) = sqrt(p_x*p_x + p_y*p_y);                       // r
        Zsig(1,i) = atan2(p_y,p_x);                                   // phi
        Zsig(2,i) = (p_x*v1 + p_y*v2) / sqrt(p_x*p_x + p_y*p_y);   // r_dot
    }

    // mean predicted measurement
    z_pred.fill(0.0);
    for (int i=0; i < 2*n_aug_+1; ++i) {
        z_pred = z_pred + weights_(i) * Zsig.col(i);
    }

    // innovation covariance matrix S
    S.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;

        // angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        S = S + weights_(i) * z_diff * z_diff.transpose();
    }

    // add measurement noise covariance matrix
    MatrixXd R = MatrixXd(n_z,n_z);
    R <<  std_radr_*std_radr_, 0, 0,
            0, std_radphi_*std_radphi_, 0,
            0, 0,std_radrd_*std_radrd_;
    S = S + R;

    MatrixXd Tc = MatrixXd(n_x_, n_z);

    // calculate cross correlation matrix
    Tc.fill(0.0);
    for (int i = 0; i < 2 * n_aug_ + 1; ++i) {  // 2n+1 simga points
        // residual
        VectorXd z_diff = Zsig.col(i) - z_pred;
        // angle normalization
        while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
        while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

        // state difference
        VectorXd x_diff = Xsig_pred.col(i) - x_;
        // angle normalization
        while (x_diff(3)> M_PI) x_diff(3)-=2.*M_PI;
        while (x_diff(3)<-M_PI) x_diff(3)+=2.*M_PI;

        Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
    }

    // Kalman gain K;
    MatrixXd K = Tc * S.inverse();

    // residual
    VectorXd z_diff = meas_package.raw_measurements_ - z_pred;

    // angle normalization
    while (z_diff(1)> M_PI) z_diff(1)-=2.*M_PI;
    while (z_diff(1)<-M_PI) z_diff(1)+=2.*M_PI;

    // update state mean and covariance matrix

    x_ = x_ + K * z_diff;
    P_ = P_ - K*S*K.transpose();

}