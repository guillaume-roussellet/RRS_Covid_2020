#include <Rcpp.h>
#include <RcppEigen.h>
#include <iostream>
#include <math.h>

using namespace Rcpp;
using namespace Eigen;

MatrixXd Kronecker_matmat(const MatrixXd& A, const MatrixXd& B){
  
  MatrixXd total = MatrixXd::Zero(A.rows()*B.rows(), A.cols()*B.cols());
  
  for (int i=0; i<A.rows(); i++){
    for (int j=0; j<A.cols(); j++){
      total.block(i*B.rows(),j*B.cols(),B.rows(),B.cols()) = A(i,j) * B;
    }
  }
  
  return total;
}


// Function C to include in the variance-covariance computation
Eigen::MatrixXd C_function(const Eigen::MatrixXd & W,
                           const Eigen::VectorXd & Z,
                           const double & a){
  
  Eigen::VectorXd Ones_vec      = Eigen::VectorXd::Ones(Z.size());
  Eigen::MatrixXd Diag_Z        = Z.asDiagonal();
  
  Eigen::VectorXd To_Be_Diaged  = W.transpose() * Z + Diag_Z * (W * Ones_vec);
  Eigen::MatrixXd Diag_Mat      = To_Be_Diaged.asDiagonal();
  
  Eigen::MatrixXd result        = a * a * (- Diag_Z * W - W.transpose() * Diag_Z + Diag_Mat);
  
  return(result);
  
}



// Function D to include in the vairance-covariance computation
Eigen::MatrixXd D_function(const Eigen::MatrixXd & W_1,
                           const Eigen::MatrixXd & W_2,
                           const Eigen::VectorXd & S, 
                           const Eigen::VectorXd & I, 
                           const double & tau_1, 
                           const double & tau_2, 
                           const Eigen::MatrixXd & Omega){
  
  int n_states = S.size();
  Eigen::MatrixXd Id = Eigen::MatrixXd::Identity(n_states, n_states);
  
  // Eigen::MatrixXd First_term_half = C_function(W, I, tau, 1-tau);
  Eigen::MatrixXd First_term_half1  = C_function(W_1, I, tau_1) + C_function(W_2, I, tau_2);
  Eigen::MatrixXd First_term_half2  = C_function(W_1, S, tau_1) + C_function(W_2, S, tau_2);
  Eigen::MatrixXd First_term        = First_term_half1.array() * First_term_half2.array();
  
  Eigen::MatrixXd Second_term       = First_term_half1.array() * ((Id + Omega) * S * S.transpose() * (Id + Omega).transpose()).array();
  Eigen::MatrixXd Third_term        = First_term_half1.array() * ((Id + Omega) * I * I.transpose() * (Id + Omega).transpose()).array();
 
 
  Eigen::MatrixXd result            = First_term + Second_term + Third_term;
  
  return(result);
  
}




// Main function to compute the variance-covariance matrix of the model. 
Eigen::MatrixXd Condi_var(const Rcpp::List & Parameters,
                          const Eigen::VectorXd & X_past,
                          const Eigen::VectorXd & I_past,
                          const Eigen::MatrixXd & Variance_AR_full,
                          const Eigen::MatrixXd & Theta_TV_var){

  // Retrieving some parameters
  double kappa                    = Parameters("kappa");
  double beta_bar                 = Parameters("beta");
  Eigen::MatrixXd Sigma_mat       = Parameters("Sigma_mat");
  Eigen::MatrixXd Diag_I_past     = I_past.asDiagonal();
  
  
  
  // Taking sizes
  int n_states          = X_past.size();
  int n_W               = 5 * n_states;
  Eigen::MatrixXd Id5   = Eigen::MatrixXd::Identity(5,5);
  Eigen::VectorXd Ones  = Eigen::VectorXd::Ones(n_states);
  
  Eigen::MatrixXd Id              = Eigen::MatrixXd::Identity(n_states, n_states);
  
  Eigen::VectorXd Intercept_Betas = beta_bar * kappa * Ones;
  Eigen::MatrixXd Intercept_diag  = Intercept_Betas.asDiagonal();
  Eigen::MatrixXd AR_Betas        = (1 - kappa) * Id;
  
  double timestep = 1.0/365.0;
  
 
  
  // Computing the first matrix in the conditional variance
  //-------------------------------------------------------
  Eigen::MatrixXd var_premult_TV  = Eigen::MatrixXd::Zero(5, 5);
  var_premult_TV(1,1)             = 1;
  var_premult_TV(2,2)             = 1;
  var_premult_TV(1,2)             = -1;
  var_premult_TV(2,1)             = -1;
  
  Eigen::MatrixXd First_matrix_variance = (Variance_AR_full * Kronecker_matmat(Id5, Diag_I_past)  
    + Kronecker_matmat(var_premult_TV, Theta_TV_var)); 



  // Initializing the Matrix
  Eigen::MatrixXd Final_vcov = Eigen::MatrixXd::Zero(n_W, n_W);

  // Filling up the matrix
  //----------------------
  Eigen::VectorXd sqrt_X = X_past.array().sqrt();
  Eigen::MatrixXd sqrt_X_diag = sqrt_X.asDiagonal();
  Final_vcov.block(4 * n_states, 4 * n_states, n_states, n_states) = 2 * timestep * 
    sqrt_X_diag * Sigma_mat * sqrt_X_diag;
  
  // Last Sum
  //---------
  Final_vcov += First_matrix_variance;

  return(Final_vcov);
  
}




Rcpp::List Build_state(const Rcpp::List & Parameters,
                       const Eigen::MatrixXd & Death_diff_data){


  // Taking the parameters
  double delta = Parameters("delta");
  double gamma = Parameters("gamma");
  double kappa = Parameters("kappa");
  double beta  = Parameters("beta");
  double nu    = ((1 - delta - gamma) * (1 - delta - gamma) * delta + gamma * (1 - delta - gamma))/(1 - delta);

  
  // Taking sizes
  int n_states  = Death_diff_data.cols();
  int n_W       = 5 * n_states;


  // Building objects to retrieve
  Eigen::VectorXd Intercept_full    = Eigen::VectorXd::Zero(n_W);
  Eigen::MatrixXd AR_full           = Eigen::MatrixXd::Zero(n_W, n_W);
  Eigen::MatrixXd Mean_mult_full    = Eigen::MatrixXd::Zero(n_W, n_states);
  Eigen::MatrixXd Variance_AR_full  = Eigen::MatrixXd::Zero(n_W, n_W);

  // Additional matrices and vectors to use
  //---------------------------------------
  Eigen::MatrixXd Id    = Eigen::MatrixXd::Identity(n_states, n_states);
  Eigen::VectorXd Ones  = Eigen::VectorXd::Ones(n_states);

  // Mean model: constructing the matrices without the time-varying elements
  //========================================================================
  // AR matrix
  //----------
  // First row block
  AR_full.block(0, 2 * n_states, n_states, n_states)              = delta * Id;

  // Second row block
  AR_full.block(n_states, n_states, n_states, n_states)           = Id;

  // Third row block
  AR_full.block(2 * n_states, 2 * n_states, n_states, n_states)   = (1 - delta - gamma) * Id;

  // Fourth row block
  AR_full.block(3 * n_states, 2 * n_states, n_states, n_states)   = gamma * Id;
  AR_full.block(3 * n_states, 3 * n_states, n_states, n_states)   = Id;

  // Fifth row block
  AR_full.block(4 * n_states, 4 * n_states, n_states, n_states)   = (1 - kappa) * Id;
  
  // Intercept
  Intercept_full.segment(4 * n_states, n_states)  = beta * kappa * Ones;


  // Term multiplying the nonlinear term in the mean
  //------------------------------------------------
  Mean_mult_full.block(n_states, 0, n_states, n_states)      = - Id;
  Mean_mult_full.block(2 * n_states, 0, n_states, n_states)  = Id;


  // Constructing the terms for the stars
  //-------------------------------------
  Eigen::MatrixXd AR_star         = AR_full.block(n_states, n_states, 3 * n_states, 3 * n_states);
  Eigen::MatrixXd Mean_mult_star  = Mean_mult_full.block(n_states, 0, 3 * n_states, n_states);


  // Variance Model
  //=================
  // First row block
  Variance_AR_full.block(0, 0, n_states, n_states)                        = delta * (1 - delta) * Id;
  Variance_AR_full.block(0, 2 * n_states, n_states, n_states)             = - delta * (1 - delta - gamma) * Id;
  Variance_AR_full.block(0, 3 * n_states, n_states, n_states)             = - delta * gamma * Id;

  // Third row block
  Variance_AR_full.block(2 * n_states, 0, n_states, n_states)             = - delta * (1 - delta - gamma) * Id;
  Variance_AR_full.block(2 * n_states, 2 * n_states, n_states, n_states)  = nu * Id;
  Variance_AR_full.block(2 * n_states, 3 * n_states, n_states, n_states)  = - gamma * (1 - delta - gamma) * Id;

  // Fourth row block
  Variance_AR_full.block(3 * n_states, 0, n_states, n_states)             = - delta * gamma * Id;
  Variance_AR_full.block(3 * n_states, 2 * n_states, n_states, n_states)  = - gamma * (1 - delta - gamma) * Id;
  Variance_AR_full.block(3 * n_states, 3 * n_states, n_states, n_states)  = gamma * (1 - gamma) * Id;

  
  // Send back results
  //==================
  return(Rcpp::List::create(
      Named("Intercept_full")   = Intercept_full,
      Named("AR_full")          = AR_full,
      Named("Variance_AR_full") = Variance_AR_full,
      Named("Mean_mult_full")   = Mean_mult_full,
      Named("AR_star")          = AR_star,
      Named("Mean_mult_star")   = Mean_mult_star));
}


// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::export]]
Rcpp::List EKF_filter_cpp(const Rcpp::List      & Parameters,
                          const Eigen::VectorXd & W0,
                          const Eigen::MatrixXd & P0,
                          const Eigen::MatrixXd & Death_diff_data,
                          const Eigen::MatrixXd & Dummy_theta_F,
                          const Eigen::MatrixXd & Dummy_theta_I,
                          const Eigen::MatrixXd & Dummy_masks,
                          const Eigen::MatrixXd & Weight_commuting,
                          const Eigen::MatrixXd & Weight_travel,
                          const Eigen::ArrayXd  & Pop_vector){
                            

  // Weight matrices are (# [states] * # [states])


  // Retrieve parameters
  Eigen::VectorXd theta_F_high  = Parameters("theta_F_high");
  Eigen::VectorXd theta_F_low   = Parameters("theta_F_low");
  Eigen::VectorXd theta_I_high  = Parameters("theta_I_high");
  Eigen::VectorXd theta_I_low   = Parameters("theta_I_low");
  Eigen::MatrixXd Premult_Beta  = Parameters("Premult_Beta");
  double mask_reduction         = Parameters("Mask_reduc");

  Eigen::MatrixXd Diag_theta_F_high   = theta_F_high.asDiagonal();
  Eigen::MatrixXd Diag_theta_F_low    = theta_F_low.asDiagonal();
  Eigen::MatrixXd Diag_theta_I_high   = theta_I_high.asDiagonal();
  Eigen::MatrixXd Diag_theta_I_low    = theta_I_low.asDiagonal();
  
  // Taking parameter sizes
  int n_states      = Death_diff_data.cols();
  int time          = Death_diff_data.rows();
  int n_W           = 5 * n_states;
  double tau_com    = Parameters("tau_com");
  double tau_trav   = Parameters("tau_trav");

  
  
  // Constructing Theta matrices
  Eigen::MatrixXd Ones_Matrix     = Eigen::MatrixXd::Ones(time, n_states);
  Eigen::VectorXd Ones_vector     = Eigen::VectorXd::Ones(n_states);
  Eigen::MatrixXd Id              = Eigen::MatrixXd::Identity(n_states, n_states);
  Eigen::MatrixXd Id5x5           = Eigen::MatrixXd::Identity(5, 5);
  Eigen::MatrixXd Id5xStates      = Eigen::MatrixXd::Identity(5 * n_states, 5 * n_states);
  Eigen::VectorXd Ones4x1         = Eigen::VectorXd::Ones(4);
  

  
  Eigen::MatrixXd theta_masks     = Ones_Matrix * mask_reduction + Dummy_masks * (1 - mask_reduction);
  Eigen::MatrixXd theta_I_matrix  = Ones_Matrix * Diag_theta_I_low + Dummy_theta_I * (Diag_theta_I_high - Diag_theta_I_low);
  Eigen::MatrixXd theta_F_matrix  = Ones_Matrix * Diag_theta_F_low + Dummy_theta_F * (Diag_theta_F_high - Diag_theta_F_low);
  
  
  // Constructing State-Space
  Rcpp::List State_Space              = Build_state(Parameters, Death_diff_data);

  Eigen::VectorXd Intercept_full      = State_Space("Intercept_full");
  Eigen::MatrixXd AR_full             = State_Space("AR_full");
  Eigen::MatrixXd Mean_mult_full      = State_Space("Mean_mult_full");
  Eigen::MatrixXd Variance_AR_full    = State_Space("Variance_AR_full");
  
  // Measurement Equations stdev
  double Sd_Death                           = Parameters("meas_std_death");
  Eigen::MatrixXd Var_Meas                  = Eigen::MatrixXd::Zero(1 * n_states, 1 * n_states);
  Var_Meas.block(0, 0, n_states, n_states)  = Sd_Death * Sd_Death * Id;
  
  

  // Retrieve the Measurement Equation parameters
  Eigen::MatrixXd Mult_Mat_Meas_Eq                  = Eigen::MatrixXd::Zero(5 * n_states, 1 * n_states);
  Mult_Mat_Meas_Eq.block(0, 0, n_states, n_states)  = Id;


  // Building objects to retrieve
  //-----------------------------
  double pi = 3.141592653589793;

  Eigen::MatrixXd W_t_t           = Eigen::MatrixXd::Zero(n_W, time);
  Eigen::MatrixXd W_t_t_nomodif   = Eigen::MatrixXd::Zero(n_W, time);
  Eigen::MatrixXd W_t_t1          = Eigen::MatrixXd::Zero(n_W, time);
  Eigen::MatrixXd Obs_t_t1        = Eigen::MatrixXd::Zero(1 * n_states, time);
  Eigen::MatrixXd residuals_t_t1  = Eigen::MatrixXd::Zero(1 * n_states, time);

  Rcpp::List P_t_t                = Rcpp::List(time);
  Rcpp::List P_t_t1               = Rcpp::List(time);
  Rcpp::List M_t_t1               = Rcpp::List(time);
  Rcpp::List Total_Gradients      = Rcpp::List(time);

  Eigen::VectorXd loglik          = Eigen::VectorXd::Zero(time);


  // Initialize the Filter
  //----------------------
  Eigen::VectorXd W = W0;
  Eigen::MatrixXd P = P0;

  Eigen::VectorXd W_pred          = Eigen::VectorXd::Zero(n_W);
  Eigen::VectorXd Obs_pred        = Eigen::VectorXd::Zero(1 * n_states);
  Eigen::VectorXd Observables_day = Eigen::VectorXd::Zero(1 * n_states);

  Eigen::MatrixXd P_pred          = Eigen::MatrixXd::Zero(n_W, n_W);
  Eigen::MatrixXd M_pred          = Eigen::MatrixXd::Zero(1 * n_states, 1 * n_states);

  // Jacobian used for the conditional covariance matrix
  Eigen::MatrixXd Jacobian_vcov   = Eigen::MatrixXd::Zero(n_states, n_W);
  Eigen::MatrixXd Total_jacobian  = Eigen::MatrixXd::Zero(n_W, n_W);
  Eigen::MatrixXd Condi_P_pred    = Eigen::MatrixXd::Zero(n_W, n_W);

  
  // Elements that will be updated in the recursions
  //------------------------------------------------
  Eigen::VectorXd theta_F_vector    = theta_F_matrix.row(0);
  Eigen::VectorXd theta_I_vector    = theta_I_matrix.row(0);
  Eigen::VectorXd theta_mask_vector = theta_masks.row(0);
  // COMMUTING MATRICES
  Eigen::MatrixXd Weight_commut_t   = Weight_commuting * Diag_theta_I_high;
  Eigen::MatrixXd Diag_mat_commut   = (Weight_commut_t * Ones_vector).asDiagonal();
  Eigen::MatrixXd Omega_commut      = tau_com * Weight_commut_t.transpose() - tau_com * Diag_mat_commut;
  // TRAVEL MATRICES
  Eigen::MatrixXd Weight_travel_t   = tau_trav * Weight_travel * Diag_theta_I_high * Diag_theta_F_high;
  Eigen::MatrixXd Diag_mat_travel   = (Weight_travel_t * Ones_vector).asDiagonal();
  Eigen::MatrixXd Omega_travel      = Weight_travel_t.transpose() - Diag_mat_travel;
  // DEFINE OMEGA(t)
  Eigen::MatrixXd Omega_t           = Omega_commut + Omega_travel ;


  // Separating the different components
  Eigen::ArrayXd D      = W.segment(0, n_states);
  Eigen::VectorXd S     = W.segment(1 * n_states, n_states);
  Eigen::VectorXd I     = W.segment(2 * n_states, n_states);
  Eigen::ArrayXd R      = W.segment(3 * n_states, n_states);
  Eigen::VectorXd X     = W.segment(4 * n_states, n_states);
  Eigen::ArrayXd Betas  = (Premult_Beta * X).array();

  Eigen::ArrayXd Normalized_Betas = Betas/Pop_vector;

  // Building the nonlinear part of the mean vector
  Eigen::ArrayXd Commut_Travel_S  = ((Id + Omega_t) * S).array();
  Eigen::ArrayXd Commut_Travel_I  = ((Id + Omega_t) * I).array();
  Eigen::ArrayXd Nonlinear_arr    = theta_I_high.array() * Normalized_Betas * Commut_Travel_S * Commut_Travel_I;
  Eigen::VectorXd Nonlinear_vec   = Nonlinear_arr;



  // Building the nonlinear part in the variance-covariance matrix
  Eigen::VectorXd First_vec_Theta_TV      = theta_I_vector.array() * theta_mask_vector.array() * Normalized_Betas;
  Eigen::MatrixXd First_Mat_Theta_TV      = First_vec_Theta_TV * First_vec_Theta_TV.transpose();
  Eigen::MatrixXd Second_Mat_Theta_TV     = D_function(Weight_commut_t, Weight_travel_t, S, I, tau_com, 1, Omega_t);
  Eigen::MatrixXd First_row_Theta_TV      = First_Mat_Theta_TV.array() * Second_Mat_Theta_TV.array();

  Eigen::MatrixXd Diag_mult_Theta_TV_var  = Nonlinear_vec.asDiagonal();
  Eigen::MatrixXd Theta_TV_var            = First_row_Theta_TV + Diag_mult_Theta_TV_var;
  //std::cout << Theta_TV_var << std::endl;




  // LAUNCH THE RECURSIONS
  //======================
  for(int t = 0; t < time; t++){

  
    //======================================================================
    // Elements to be be updated in the recursions
    //------------------------------------------------
    theta_F_vector    = theta_F_matrix.row(t);
    theta_I_vector    = theta_I_matrix.row(t);
    theta_mask_vector = theta_masks.row(t);
    
    // COMMUTING MATRICES
    Weight_commut_t   = Weight_commuting * theta_I_vector.asDiagonal();
    Diag_mat_commut   = (Weight_commut_t * Ones_vector).asDiagonal();
    Omega_commut      = tau_com * Weight_commut_t.transpose() - tau_com * Diag_mat_commut;
    // TRAVEL MATRICES
    Weight_travel_t   = tau_trav * Weight_travel * theta_F_vector.asDiagonal() * theta_I_vector.asDiagonal();
    Diag_mat_travel   = (Weight_travel_t * Ones_vector).asDiagonal();
    Omega_travel      = Weight_travel_t.transpose() - Diag_mat_travel;
    // OMEGA(t)
    Omega_t           = Omega_commut + Omega_travel;
    
    
    // Building the nonlinear part of the mean vector
    //-----------------------------------------------
    Commut_Travel_S                   = ((Id + Omega_t) * S).array();
    Commut_Travel_I                   = ((Id + Omega_t) * I).array();
    Nonlinear_arr                     = theta_I_vector.array() * theta_mask_vector.array() * Normalized_Betas * Commut_Travel_I * Commut_Travel_S;
    Nonlinear_vec                     = Nonlinear_arr;
    
    
    // Building the nonlinear part in the variance-covariance matrix
    //--------------------------------------------------------------
    First_vec_Theta_TV      = theta_I_vector.array() * theta_mask_vector.array() * Normalized_Betas;
    First_Mat_Theta_TV      = First_vec_Theta_TV * First_vec_Theta_TV.transpose();
    Second_Mat_Theta_TV     = D_function(Weight_commut_t, Weight_travel_t, S, I, tau_com, 1, Omega_t);
    First_row_Theta_TV      = First_Mat_Theta_TV.array() * Second_Mat_Theta_TV.array();
    
    Diag_mult_Theta_TV_var  = Nonlinear_vec.asDiagonal();
    Theta_TV_var            = First_row_Theta_TV + Diag_mult_Theta_TV_var;
    
    //===================================================================
    //===================================================================
    // PREDICTION STEP
    //----------------
    // Mean of latent factors
    W_pred = Intercept_full + AR_full * W + Mean_mult_full * Nonlinear_vec;

   

    // Conditional variance matrix
    // Jacobian
    Eigen::VectorXd Jacobian_Block_1  = theta_mask_vector.array() * theta_I_vector.array() * Normalized_Betas * Commut_Travel_I;
    Eigen::VectorXd Jacobian_Block_2  = theta_mask_vector.array() * theta_I_vector.array() * Normalized_Betas * Commut_Travel_S;
    Eigen::VectorXd Jacobian_Block_3  = theta_mask_vector.array() * theta_I_vector.array()/Pop_vector * Commut_Travel_I * Commut_Travel_S;

    Eigen::MatrixXd Diag_Jacobian_Block_1 = Jacobian_Block_1.asDiagonal();

    Jacobian_vcov                                            = 0 * Jacobian_vcov;
    Jacobian_vcov.block(0, n_states, n_states, n_states)     = (Id + Omega_t) * Jacobian_Block_1.asDiagonal();
    Jacobian_vcov.block(0, 2 * n_states, n_states, n_states) = (Id + Omega_t) * Jacobian_Block_2.asDiagonal();
    Jacobian_vcov.block(0, 4 * n_states, n_states, n_states) = Jacobian_Block_3.asDiagonal() * Premult_Beta;

    Total_jacobian = AR_full + Mean_mult_full * Jacobian_vcov;

    // Predict the variance-covariance matrix
    P_pred = Total_jacobian * P * Total_jacobian.transpose() +
      Condi_var(Parameters, X, I, Variance_AR_full, Theta_TV_var);


    // Predicting the observables
    Obs_pred  = Mult_Mat_Meas_Eq.transpose() * W_pred;
    M_pred    = Mult_Mat_Meas_Eq.transpose() * P_pred * Mult_Mat_Meas_Eq + Var_Meas;

    //std::cout << Obs_pred <<std::endl;

    // UPDATING STEP
    //--------------
    Observables_day.segment(0, n_states)  = Death_diff_data.row(t);
    Eigen::VectorXd residual              = Observables_day - Obs_pred ;



    // Determining which components are observed
    Array<bool,Dynamic,1> bool_na(Observables_day.array() == Observables_day.array());

    int n_obs_modif = bool_na.count();

    //std::cout << n_obs_modif <<std::endl;

    // Fill the new components
    Eigen::VectorXd resid_modif             = Eigen::VectorXd::Zero(n_obs_modif);
    Eigen::MatrixXd M_modif                 = Eigen::MatrixXd::Zero(n_obs_modif, n_obs_modif);
    Eigen::MatrixXd Mult_Mat_Meas_Eq_modif  = Eigen::MatrixXd::Zero(n_W, n_obs_modif);


    // Control for unobserved components
    //----------------------------------
    int counter_rows = 0;
    for (int i = 0; i < 1 * n_states; i++){

      if(bool_na(i)==TRUE){

        resid_modif(counter_rows)                 = residual(i);
        Mult_Mat_Meas_Eq_modif.col(counter_rows)  = Mult_Mat_Meas_Eq.col(i);

        int counter_cols = 0 ;
        for (int j = 0; j < 1 * n_states; j++){
          if (bool_na(j) == TRUE){
            M_modif(counter_rows, counter_cols) = M_pred(i,j);
            counter_cols += 1;
          }
        }
        counter_rows += 1;
      }
    }// End of the NA loop




    // Perform the update step and the loglik
    //---------------------------------------
    double det_M_pred      = M_modif.determinant();



    Eigen::MatrixXd M_inv = M_modif.inverse();
    Eigen::MatrixXd Gain  = P_pred * Mult_Mat_Meas_Eq_modif * M_inv;

    Eigen::VectorXd W_up  = W_pred + Gain * resid_modif;
    Eigen::MatrixXd P_up  = P_pred - Gain * Mult_Mat_Meas_Eq_modif.transpose() * P_pred;

    //std::cout << W_up << std::endl;

    double lik_val         = -.5*(n_obs_modif * log(2*pi) + log(det_M_pred) +
                                  resid_modif.transpose() * M_inv * resid_modif);

    W_t_t_nomodif.col(t)  = W_up;



    // Applying the correction on the updated values
    //----------------------------------------------
    // Correcting the evolution
    Eigen::VectorXd I_ToModif = W_up.segment(2 * n_states, n_states);


    // Changing the W_up
    W_up.segment(2 * n_states, n_states)    = I_ToModif;
    for(int i = 0; i < n_W; i++){
      if(W_up(i) < 0){W_up(i) = 0;}
    }



    // Update the values before looping back
    //--------------------------------------
    W = W_up;
    P = P_up;



    // Separating the different components
    D       = W.segment(0, n_states);
    S       = W.segment(1 * n_states, n_states);
    I       = W.segment(2 * n_states, n_states);
    R       = W.segment(3 * n_states, n_states);
    X       = W.segment(4 * n_states, n_states);
    Betas   = Premult_Beta * X;
    
    Normalized_Betas = Betas/Pop_vector;


    // STORING THE VALUES
    //---------------------
    W_t_t.col(t)          = W_up;
    W_t_t1.col(t)         = W_pred;
    Obs_t_t1.col(t)       = Obs_pred;
    residuals_t_t1.col(t) = residual;
    Total_Gradients(t)    = Total_jacobian;

    P_t_t(t)  = P_up;
    P_t_t1(t) = P_pred;
    M_t_t1(t) = M_pred;

    loglik(t) = lik_val;

    


  }
  
  
  Eigen::MatrixXd Beta_updated = Premult_Beta * W_t_t.block(4 * n_states, 0, n_states, time);

  // Sending results back
  //---------------------
  return List::create(
    Named("loglik.vector")  = loglik,
    Named("W.updated")      = W_t_t,
    Named("P.updated")      = P_t_t,
    Named("W.predicted")    = W_t_t1,
    Named("P.predicted")    = P_t_t1,
    Named("Obs.predicted")  = Obs_t_t1,
    Named("M.predicted")    = M_t_t1,
    Named("W.nomodif")      = W_t_t_nomodif,
    Named("W.0")            = W0,
    Named("P.0")            = P0,
    Named("Grads")          = Total_Gradients,
    Named("Mult")           = Mult_Mat_Meas_Eq,
    Named("AR_full")        = AR_full,
    Named("Intercept_full") = Intercept_full,
    Named("residual")       = residuals_t_t1,
    Named("Betas.updated")  = Beta_updated
  );
  

}
















