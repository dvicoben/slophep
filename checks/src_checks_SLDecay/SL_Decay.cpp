#include "SL_Decay.h"
#include "SL_Decay_Config.cpp"
#include <iostream>
#include <stdlib.h>
#include <cmath>
#include <map>
#include <iterator>
#include <string>

//   *************************************   //
//                                           //
//          PUBLIC MEMBER FUNCTION           //
//                                           //
//   *************************************   //


SL_Decay::SL_Decay(int p1, int p2, int p3, int p4, int p5, std::string parametrization)
{
  std::string Decay_chain[5];
  double par[4];
  bool exit_program = true;
  BtoDstMuNu = false;
  BtoDMuNu = false;

  // Check that the initialization paraters correspond to an existing decay.
  for(int i=0; i<npd; i++){
    for(int j=0; j<4; j++){
      par[j] = constr_param_list[i][j];
    } // end j-loop
    if (p1==par[0] && p2==par[1] && p3==par[2] && p4==par[3]){
      exit_program = false;
      break;   
    }
  } // end i-loop
  
  if(exit_program) {
    std::cout << "You inserted a non physical set of parameters in the constructor. Please try again." << std::endl;
    exit (EXIT_FAILURE);
  }

  // Initialize the object by selecting the particles of the decay and their masses. 
  // Set the B meson.
  if (p1==1) {
    B_mass_ = m_B_0;
    Decay_chain[0] = "B_0";
    normailization_index_ = 1;
  }
  else
    if (p1==2) {
      B_mass_ = m_B_pm;
      Decay_chain[0] = "B_pm";
      normailization_index_ = 2;
    }
    else
      if (p1==3) {
	B_mass_ = m_B_s;
	Decay_chain[0] = "B_s";
	normailization_index_ = 3;
      }

  // Set the D* meson.                                   
  if (p2==0){         // For the B -> D lnu decays there are no D* mesons
    D_st_mass_ = 0.;
    Decay_chain[1] = "NO_D_st";
  }
  else
    if (p2==1){
      D_st_mass_ = m_D_st_0;    
      Decay_chain[1] = "D_st_0";
    }
    else
      if (p2==2) {
	D_st_mass_ = m_D_st_pm;
	Decay_chain[1] = "D_st_pm";
      }
      else
	if (p2==3) {
	  D_st_mass_ = m_D_st_s_pm;
	  Decay_chain[1] = "D_st_s_pm";
	}

  // Set the D meson. 
  if (p3==1){
    D_mass_ = m_D_0;
    Decay_chain[2] = "D_0";
  }
  if (p3==2) {
    D_mass_ = m_D_pm;
    Decay_chain[2] = "D_pm";
  }
  if (p3==3) {
    D_mass_ = m_D_s_pm;
    Decay_chain[2] = "D_s_pm";
  }

  // Set the D partner in the D* decay.
  if (p4==0){            // For the B -> D lnu decays there are no D* mesons decaying into D + X
    pi_gamma_mass_ = 0.;
    Decay_chain[3] = "NO_pi_or_gamma";
  }
  if(p4==1){
    pi_gamma_mass_ = m_pi_0;
    Decay_chain[3] = "pi_0";
    is_pi = true;
  }
  if (p4==2){
    pi_gamma_mass_ = m_pi_pm;
    Decay_chain[3] = "pi_pm";
    is_pi = true; 
  }
  if(p4==3){
    pi_gamma_mass_ = m_gamma;
    Decay_chain[3] = "gamma";
    is_pi = false; 
  }

  // Set the lepton.
  if (p5==1){
    lepton_mass_ = m_e_pm;
    Decay_chain[4] = "e_pm";
    lepton_index_ = 1;
  } 
  if (p5==2){
    lepton_mass_ = m_mu_pm;
    Decay_chain[4] = "mu_pm";
    lepton_index_ = 2;
  }
  if (p5==3){
    lepton_mass_ = m_tau_pm;
    Decay_chain[4] = "tau_pm";
    lepton_index_ = 3;
  }


  // Set the boolean flag to distinguish B -> D* lnu from B -> D lnu decays.
  if(p2==0 && p4==0) {BtoDMuNu = true;}
  else {BtoDstMuNu = true;}

 
  // Set the branching ratio of the process D* -> D + something
  if(BtoDstMuNu){
    std::string BR_tag = "BR_" + Decay_chain[1] + "_to_" + Decay_chain[2] + "_" +Decay_chain[3];
    exit_program = true;
    for(int i=0; i<npd-3; i++){
      if(BR_map[i].get_BR_tag() == BR_tag) {
	BR_ = BR_map[i].get_BR();
	exit_program = false;
      }
    }
  }
  else
    if(BtoDMuNu) {exit_program = false;} // In this case nothing has to be checked further.
  
  // Check for correct initialization of the object.
  if(exit_program) {
    std::cout << "No D* BR has been selected. Please try again." << std::endl;
    exit (EXIT_FAILURE);
  }
  
  
  // Initialize Vcb to the most up to date PDG value.
  Vcb_ = Vcb;

  // ********************************** // 
  //       SELECT PARAMETRIZATION       //
  // ********************************** // 
  if(parametrization == "CLN"){
    CLN_ = true;
    BGL_ = false;
    BCL_ = false;
    CLN_MOD_ = false;
    CLN_ORI_ = false;
    CLN_HYB_ = false;
    // Set CLN parameters to default values. 
    if(BtoDstMuNu) {
      rho_squared_ = rho_squared;
      R1_ = R1;
      R2_ = R2;
      h_A1_ = h_A1;
      R0_ = R0;
    }
    else 
      if(BtoDMuNu){
	rho_squared_ = rho_squared_BtoD;
	h_A1_ = h_A1_BtoD;
      }
  }
  else
    if(parametrization == "CLN_ORI"){
      CLN_ = false;
      BGL_ = false;
      BCL_ = false;
      CLN_MOD_ = false;
      CLN_ORI_ = true;
      CLN_HYB_ = false;
      // Set lepton masses to zero and exit if a semitauonic decay is selected.
      if(lepton_index_ == 3) {   // tau case.
	std::cout << "The CLN_ORI parametrization set the lepton's mass to zero, as the original CLN paper does." << std::endl;
	std::cout << "You can not use this parametrization with a semituaonic decay." << std::endl;
	std::cout << "Please, change lepton or use a different version of CLN (CLN, CLN_MOD)." << std::endl;
	exit (EXIT_FAILURE);
      }
      else{ lepton_mass_ = 0.;} // el and mu case.
      // Set CLN parameters to default values.
      if(BtoDstMuNu) { 
	rho_squared_ = rho_squared;
	R1_ = R1;
	R2_ = R2;
	h_A1_ = h_A1;
	R0_ = R0;
      }
      else
	if(BtoDMuNu){
	  rho_squared_ = rho_squared_BtoD;
	  h_A1_ =h_A1_BtoD;
	}
    }
    else
      if(parametrization == "CLN_MOD"){
	CLN_ = false;
	BGL_ = false;
	BCL_ = false;
	CLN_MOD_ = true;
	CLN_ORI_ = false;
	CLN_HYB_ = false;
	// Set CLN parameters to default values.
	if(BtoDstMuNu) {
	  rho_squared_ = rho_squared;
	  R1_ = R1;
	  R2_ = R2;
	  h_A1_ = h_A1;
	  R0_ = R0;
	  // Set the two additional parameters to the standard values of the CLN parametrization.
	  c_z2_ = c_z2;
	  c_z3_ = c_z3;
	}
	else
	  if(BtoDMuNu){
	    rho_squared_ = rho_squared_BtoD;
	    h_A1_ =h_A1_BtoD;
	    // c_z2_ = c_z2;
	    // c_z3_ = c_z3;
	  }
      }
      else
	if(parametrization == "CLN_HYB"){
	  CLN_ = false;
	  BGL_ = false;
	  BCL_ = false;
	  CLN_MOD_ = false;
	  CLN_ORI_ = false;
	  CLN_HYB_ = true;
	  if(BtoDstMuNu) {
	    // Set CLN parameters to default values.
	    R1_ = R1;
	    R2_ = R2;
	    R0_ = R0;
	    // Set the BGL parameters of the f form factors.
	    for(int i=0; i<4; i++){ a_f_BGL_[i] = a_f_BGL[i]; }
	    chi_T_1p_ = chi_T_1p;
	  }
	  else
	    if(BtoDMuNu){
	      std::cout << "The parametrization CLN_HYB is meant to be uesed only for the B -> D* l nu decays." << std::endl;
	      std::cout << "Please, select another parametrization for your B -> D l nu decay." << std::endl;
	      exit (EXIT_FAILURE);
	    }
	}
	else
	  if(parametrization == "BGL"){
	    CLN_ = false;
	    BGL_ = true;
	    BCL_ = false;
	    CLN_MOD_ = false;
	    CLN_ORI_ = false;
	    CLN_HYB_ = false;
	    // Set BGL parameters to default values.
	    if(BtoDstMuNu) {
	      for(int i=0; i<4; i++) {
		a_f_BGL_[i] = a_f_BGL[i];
		a_F1_BGL_[i] = a_F1_BGL[i];
		a_g_BGL_[i] = a_g_BGL[i];
		a_F2_BGL_[i] = a_F2_BGL[i];
	      }
	      // Set a_F1[0] to the correct value.
	      a_F1_BGL_[0] = a_F1_0_BGL();
	      
	      chi_T_1p_ = chi_T_1p;
	      chi_T_1m_ = chi_T_1m;
	      chi_L_1p_ = chi_L_1p;
	    }
	    else
	      if(BtoDMuNu){
		for(int i=0; i<4; i++) {
		  a_BGL_[i] = a_BGL[i];
		  b_BGL_[i] = b_BGL[i];
		}
		chi_L_0_ = chi_L_0;
		chi_Tt_0_ = chi_Tt_0;
	      }
	  }
	  else
	    if(parametrization == "BCL"){
	      CLN_ = false;
	      BGL_ = false;
	      BCL_ = true;
	      CLN_MOD_ = false;
	      CLN_ORI_ = false;
	      CLN_HYB_ = false;
	      // Set BCL parameters to default values.
	      if(BtoDstMuNu) {
	      for(int i=0; i<4; i++) {
		a_f_BCL_[i] = a_f_BCL[i];
		a_F1_BCL_[i] = a_F1_BCL[i];
		a_g_BCL_[i] = a_g_BCL[i];
		a_F2_BCL_[i] = a_F2_BCL[i];
	      }
	      
	      N_f_ = N_f;
	      N_F1_ = N_F1;
	      N_g_ = N_g;
	      N_F2_ = N_F2;
	      }
	      else
		if(BtoDMuNu){
		  for(int i=0; i<4; i++) {
                    a_BCL_[i] = a_BCL[i];
                    b_BCL_[i] = b_BCL[i];
                  }
		}
	    }
	    else{
	      std::cout << "You inserted a wrong type of parametrization." << std::endl;
	      std::cout << "Please type CLN, CLN_MOD, CLN_ORI, CLN_HYB, BGL or BCL." << std::endl;
	      exit (EXIT_FAILURE);
	    }
  
  
} // End of constructor.


// Second constructor. USE ONLY FOR B -> D* lnu DECAYS.
SL_Decay::SL_Decay(double m1, double m2, double m3, double m4, double m5, double BR, std::string pi_or_gamma, std::string parametrization){
  
  B_mass_ = m1;
  D_st_mass_ = m2;
  D_mass_ = m3;
  pi_gamma_mass_ = m4;
  lepton_mass_ = m5;
  BR_ = BR;
  if(pi_or_gamma == "pi") is_pi = true;
  else 
    if(pi_or_gamma == "gamma") is_pi = false;
    else 
      {   
	std::cout << "You inserted a wrong tag. Please type pi or gamma." << std::endl;
	exit (EXIT_FAILURE);
      } 
     
  // Initialize Vcb to the most up to date PDG value.
  Vcb_ = Vcb;
  
  // ********************************** //
  //       SELECT PARAMETRIZATION       //
  // ********************************** //
  if (parametrization == "CLN"){
    CLN_ = true;
    BGL_ = false;
    BCL_ = false;
    CLN_MOD_ = false;
    CLN_ORI_ = false;
    CLN_HYB_ = false;
    // Set CLN parameters to default values.
    if(BtoDstMuNu) {
      rho_squared_ = rho_squared;
      h_A1_ = h_A1;
      R1_ = R1;
      R2_ = R2;
      R0_ = R0;
    }
    else 
      if(BtoDMuNu){
	rho_squared_ = rho_squared_BtoD;
	h_A1_ = h_A1_BtoD;
      }
  }
  else
    if(parametrization == "CLN_MOD"){
      CLN_ = false;
      BGL_ = false;
      BCL_ = false;
      CLN_MOD_ = true;
      CLN_ORI_ = false;
      CLN_HYB_ = false;
      // Set CLN parameters to default values.
      if(BtoDstMuNu) {
	rho_squared_ = rho_squared;
	R1_ = R1;
	R2_ = R2;
	h_A1_ = h_A1;
	R0_ = R0;
	// Set the two additional parameters to the standard values of the CLN parametrization.
	c_z2_ = c_z2;
	c_z3_ = c_z3;
      }
      else
	if(BtoDMuNu){
	  rho_squared_ = rho_squared_BtoD;
	  h_A1_ =h_A1_BtoD;
	  // c_z2_ = c_z2;                                                                                                                                                                                
	  // c_z3_ = c_z3;
	}
    }
    else
      if(parametrization == "CLN_ORI"){
	CLN_ = false;
	BGL_ = false;
	BCL_ = false;
	CLN_MOD_ = false;
	CLN_ORI_ = true;
	CLN_HYB_ = false;
	// Set lepton masses to zero and exit if a semitauonic decay is selected.
	if(lepton_index_ == 3) {   // tau case.
	  std::cout << "The CLN_ORI parametrization set the lepton's mass to zero, as the original CLN paper does." << std::endl;
	  std::cout << "You can not use this parametrization with a semituaonic decay." << std::endl;
	  std::cout << "Please, change lepton or use a different version of CLN (CLN, CLN_MOD)." << std::endl;
	  std::exit (EXIT_FAILURE);
	}
	else{ lepton_mass_ = 0.;} // el and mu case.
	// Set CLN parameters to default values.
	if(BtoDstMuNu) {
	  rho_squared_ = rho_squared;
	  R1_ = R1;
	  R2_ = R2;
	  h_A1_ = h_A1;
	  R0_ = R0;
	}
	else
	  if(BtoDMuNu) {
	    rho_squared_ = rho_squared_BtoD;
	    h_A1_ =h_A1_BtoD;
	  }
      }
      else
	if(parametrization == "CLN_HYB"){
          CLN_ = false;
          BGL_ = false;
          BCL_ = false;
          CLN_MOD_ = false;
          CLN_ORI_ = false;
          CLN_HYB_ = true;
          // Set CLN parameters to default values.
	  if(BtoDstMuNu) {
	    R1_ = R1;
	    R2_ = R2;
	    R0_ = R0;
	    // Set the BGL parameters of the f form factors.
	    for(int i=0; i<4; i++) { a_f_BGL_[i] = a_f_BGL[i]; }
	  }
	  else
	    if(BtoDMuNu) {
	      std::cout << "The parametrization CLN_HYB is meant to be uesed only for the B -> D* l nu decays." << std::endl;
	      std::cout << "Please, select another parametrization for your B -> D l nu decay." << std::endl;
              exit (EXIT_FAILURE);
	    }
	}
	else
	  if(parametrization == "BGL"){
	    CLN_ = false;
	    BGL_ = true;
	    BCL_ = false;
	    CLN_MOD_ = false;
	    CLN_ORI_ = false;
	    CLN_HYB_ = false;
	    // Set BGL parameters to default values.
	    if(BtoDstMuNu) {
	    for(int i=0; i<4; i++) {
	      a_f_BGL_[i] = a_f_BGL[i];
	      a_F1_BGL_[i] = a_F1_BGL[i];
	      a_g_BGL_[i] = a_g_BGL[i];
	      a_F2_BGL_[i] = a_F2_BGL[i];
	    }
	    // Set a_F1_[0] and to the correct value.
	    a_F1_BGL_[0] = a_F1_0_BGL();
	    
	    chi_T_1p_ = chi_T_1p;
	    chi_T_1m_ = chi_T_1m;
	    chi_L_1p_ = chi_L_1p;
	    }
	    else
	      if(BtoDMuNu) {
		for(int i=0; i<4; i++) {
                  a_BGL_[i] = a_BGL[i];
                  b_BGL_[i] = b_BGL[i];
		}
                chi_L_0_ = chi_L_0;
                chi_Tt_0_ = chi_Tt_0;
	      }
	  }
	  else
	    if(parametrization == "BCL"){
	      CLN_ = false;
	      BGL_ = false;
	      BCL_ = true;
	      CLN_MOD_ = false;
	      CLN_ORI_ = false;
	      CLN_HYB_ = false;
	      // Set BCL parameters to default values.
	      if(BtoDstMuNu) {
	      for(int i=0; i<4; i++) {
		a_f_BCL_[i] = a_f_BCL[i];
		a_F1_BCL_[i] = a_F1_BCL[i];
		a_g_BCL_[i] = a_g_BCL[i];
		a_F2_BCL_[i] = a_F2_BCL[i];
	      }
	      N_f_ = N_f;
	      N_F1_ =N_F1;
	      N_g_ = N_g;
	      N_F2_ =N_F2;
	      }
	      else
		if(BtoDMuNu){
		  for(int i=0; i<4; i++) {
		    a_BCL_[i] = a_BCL[i];
		    b_BCL_[i] = b_BCL[i];
		  }
		} 
	    }
	    else{
	      std::cout << "You inserted a wrong type of parametrization." << std::endl;
	      std::cout << "Please type CLN, CLN_MOD, CLN_ORI, CLN_HYB, BGL or BCL." << std::endl;
	      exit (EXIT_FAILURE);
	    }
  
  
} // End of constructor.

SL_Decay::~SL_Decay() {};  // Destructor.

  
//  ******************************  //
//     DIFFERENTIAL DECAY WIDTH     //
//  ******************************  //
// Gamma function immplementation. Valid only for B -> D* lnu decays.
double SL_Decay::Gamma_Diff_4d(double Q2, double cos_theta_l, double cos_theta_V, double k) {
  
  // If this function is called for a B -> D lnu decay print an error message and exit.
  if(BtoDMuNu){
    std::cout << "The 4D differential decay width is available ONLY for B -> D* lnu decays." << std::endl;
    std::cout << "For B -> D lnu decays, please, use the Gamma_Diff_1d function." << std::endl;
    exit (EXIT_FAILURE);
  }
  
  Q2_ = Q2;
  cos_theta_l_ = cos_theta_l;
  cos_theta_V_ = cos_theta_V;
  k_ = k;

  //Define middle-level variables
  double func_Diff_Gamma_4d;
  double constant, A, B, C, D, E, F, G;  // Gamma_func = constant*(A+B+C+D+E+F+G)

  if(is_pi){   // Gamma func. for pions.
  
    constant = Nf()*D_st_P_module()*pow( (1 - ( pow(lepton_mass_, 2)/Q2_ ) ), 2); 
    A = I_1s()*( 1 - pow(cos_theta_V_, 2) ) + I_1c()*pow(cos_theta_V_, 2);
    B = ( I_2s()*(1 - pow(cos_theta_V_, 2)) + I_2c()*pow(cos_theta_V_, 2) )*( 2*pow(cos_theta_l_, 2) -1 );
    C = I_3()*( 1 - pow(cos_theta_V_, 2) )*( 1 - pow(cos_theta_l_, 2) )*cos(2*k_);
    D = 4*I_4()*sqrt( 1 - pow(cos_theta_V_, 2) )*cos_theta_V_*sqrt( 1 - pow(cos_theta_l_, 2) )*cos_theta_l_*cos(k_);
    E = 2*I_5()*sqrt( 1 - pow(cos_theta_V_, 2) )*cos_theta_V_*sqrt( 1 - pow(cos_theta_l_, 2) )*cos(k_);
    F = ( I_6s()*( 1 - pow(cos_theta_V_, 2)  ) + I_6c()*pow(cos_theta_V_ ,2) )*cos_theta_l_;
    G = 2*I_7()*sqrt( 1 - pow(cos_theta_V_, 2) )*cos_theta_V_*sqrt( 1 - pow(cos_theta_l_, 2) )*sin(k_);

  }  // End of if scope.

  else {   // Gamma func. for photons.

    constant = Nf()*D_st_P_module()*pow( (1 - ( pow(lepton_mass_, 2)/Q2_ ) ), 2);
    A = I_1s()*( 1 - pow(cos_theta_V_, 2) ) + 2*I_1c()*( 1 + pow(cos_theta_V_, 2) );
    B = ( I_2s()*( 1 - pow(cos_theta_V_, 2) ) + 2*I_2c()*( 1 + pow(cos_theta_V_, 2) )  )*( 2*pow(cos_theta_l_, 2) - 1 );
    C = I_3()*( 1 - pow(cos_theta_V_, 2) )*( 1 - pow(cos_theta_l_, 2) )*cos(2*k_);
    D = 4*I_4()*sqrt( 1 - pow(cos_theta_V_, 2) )*cos_theta_V_*sqrt( 1 - pow(cos_theta_l_, 2) )*cos_theta_l_*cos(k_);
    E = 2*I_5()*sqrt( 1 - pow(cos_theta_V_, 2) )*cos_theta_V_*sqrt( 1 - pow(cos_theta_l_, 2) )*cos(k_);
    F = (  I_6s()*( 1 - pow(cos_theta_V_, 2) ) + 2*I_6c()*( 1 + pow(cos_theta_V_, 2) )  )*cos_theta_l_;
    G = 2*I_7()*sqrt( 1 - pow(cos_theta_V_, 2) )*cos_theta_V_*sqrt( 1 - pow(cos_theta_l_, 2) )*sin(k_);

  }  // End of else scope.
    
  func_Diff_Gamma_4d = constant*( A + B + C + D + E + F + G);
  return func_Diff_Gamma_4d;
  
}  // End of Gamma_Diff_4d

double SL_Decay::Gamma_Diff_1d(double Q2) {

  Q2_ = Q2;
  double func_Gamma_Diff_1d = 0.;
  double constant = 0.;
  if(BtoDstMuNu){    // B -> D* lnu decays
    if (is_pi) {
      constant = Nf()*D_st_P_module()*pow( (1 - ( pow(lepton_mass_, 2)/Q2_ ) ), 2);
      func_Gamma_Diff_1d = (8./9.)*M_PI*( 6*I_1s() + 3*I_1c() -2*I_2s() - I_2c() );
    }
    else {
      constant = Nf()*D_st_P_module()*pow( (1 - ( pow(lepton_mass_, 2)/Q2_ ) ), 2);
      func_Gamma_Diff_1d = (16./9.)*M_PI*( 3*I_1s() + 12*I_1c() - I_2s() -4*I_2c() );
    }
  }
  else
    if(BtoDMuNu){   // B -> D lnu decays
      constant = pow(eta_EW ,2)*pow(G_Fermi ,2)*pow(Vcb ,2)*B_mass_*pow(lambda_func_prime(), 1/2)/(192*pow(M_PI ,3));
      func_Gamma_Diff_1d = pow( 1 - (pow(lepton_mass_, 2)/Q2_), 2) *(c_l_p()*pow(f_p() ,2) + c_l_0()*pow(f_0() ,2));
    }
    
  func_Gamma_Diff_1d *= constant;
  return func_Gamma_Diff_1d;
}


void SL_Decay::Set_Vcb(double x) { Vcb_ = x; };


// -------------------------------- //
// CLN and CLN_ORI setter functions //
// -------------------------------- //
void SL_Decay::Set_CLN(double rho_squared, double R1, double R2, double h_A1){
  if(CLN_ || CLN_ORI_){
  rho_squared_ = rho_squared;
  R1_ = R1;
  R2_ = R2;
  h_A1_ = h_A1;
  // R0 remain initializated to the default value set by the constructor.
  }
  else
    if(CLN_MOD_){
      std::cout << "You can not use the CLN setter method for an SL_Decay object with the CLN_MOD parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(CLN_HYB_){
	std::cout << "You can not use the CLN setter method for an SL_Decay object with the CLN_HYB parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
	exit (EXIT_FAILURE);
      }
      else
	if(BGL_){
	  std::cout << "You can not use the CLN setter method for an SL_Decay object with the BGL parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
	  exit (EXIT_FAILURE);
	}
	else
	  if(BCL_){
	    std::cout << "You can not use the CLN setter method for an SL_Decay object with the BCL parameterization." << std::endl;
	    std::cout << "Please try again." << std::endl;
	    exit (EXIT_FAILURE);
	  }
}

void SL_Decay::Set_h_A1(double h_A1) { 
  if(CLN_ || CLN_MOD_ || CLN_ORI_){ h_A1_ = h_A1; }
  else
    if(CLN_HYB_){
      std::cout << "You can not use a CLN setter method for an SL_Decay object with the CLN_HYB parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(BGL_){
	std::cout << "You can not use a CLN setter function for an SL_Decay object with the BGL parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
	exit (EXIT_FAILURE);
      }
      else
	if(BCL_){
	  std::cout << "You can not use a CLN setter function for an SL_Decay object with the BCL parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
	  exit (EXIT_FAILURE);
	}
}

void SL_Decay::Set_R1(double R1) { 
  if(CLN_ || CLN_MOD_ || CLN_ORI_ || CLN_HYB_) { R1_ = R1; }
  else
    if(BGL_){
      std::cout << "You can not use a CLN setter function for an SL_Decay object with the BGL parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(BCL_){
	std::cout << "You can not use a CLN setter function for an SL_Decay object with the BCL parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
        exit (EXIT_FAILURE);
      }
}

void SL_Decay::Set_R2(double R2) { 
  if(CLN_ || CLN_MOD_ || CLN_ORI_ || CLN_HYB_) { R2_ = R2; }
  else
    if(BGL_){
      std::cout << "You can not use a CLN setter function for an SL_Decay object with the BGL parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(BCL_){
	std::cout << "You can not use a CLN setter function for an SL_Decay object with the BCL parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
        exit (EXIT_FAILURE);
      }
}

void SL_Decay::Set_R0(double R0) { 
  if(CLN_ || CLN_MOD_ || CLN_ORI_ || CLN_HYB_) { R0_ = R0; }
  else
    if(BGL_){
      std::cout << "You can not use a CLN setter function for an SL_Decay object with the BGL parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(BCL_){
	std::cout << "You can not use a CLN setter function for an SL_Decay object with the BCL parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
        exit (EXIT_FAILURE);
      }
}

void SL_Decay::Set_rho_squared(double rho_squared) { 
  if(CLN_ || CLN_MOD_ || CLN_ORI_) { rho_squared_ = rho_squared; } 
  else
    if(CLN_HYB_){
      std::cout << "You can not use a CLN setter method for an SL_Decay object with the CLN_HYB parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(BGL_){
	std::cout << "You can not use a CLN setter function for an SL_Decay object with the BGL parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
	exit (EXIT_FAILURE);
      }
      else
	if(BCL_){
	  std::cout << "You can not use a CLN setter function for an SL_Decay object with the BCL parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
	  exit (EXIT_FAILURE);
	}
}

// --------------------------------- //
// CLN_MOD setter and getter methods //
// --------------------------------- //
void SL_Decay::Set_CLN_MOD(double rho_squared, double R1, double R2, double h_A1, double coeff_z2, double coeff_z3) {
  if(CLN_MOD_){
    rho_squared_ = rho_squared;
    R1_ = R1;
    R2_ = R2;
    h_A1_ = h_A1;
    // R0 remain initializated to the default value set by the constructor.
    c_z2_ = coeff_z2;
    c_z3_ = coeff_z3;
  }
  else
    if(CLN_){
      std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the CLN parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(CLN_ORI_){
	std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the CLN_ORI parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
	exit (EXIT_FAILURE);
      }
      else
	if(CLN_HYB_){
	  std::cout << "You can not use the CLN_MOD setter method for an SL_Decay object with the CLN_HYB parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
	  exit (EXIT_FAILURE);
	}
	else
	  if(BGL_){
	    std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the BGL parameterization." << std::endl;
	    std::cout << "Please try again." << std::endl;
	    exit (EXIT_FAILURE);
	  } 
	  else
	    if(BCL_){
	      std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the BCL parameterization." << std::endl;
	      std::cout << "Please try again." << std::endl;
	      exit (EXIT_FAILURE);
	    }
}

void SL_Decay::Set_cz2(double coeff_z2){
  if(CLN_MOD_){ c_z2_ = coeff_z2; }
  else
    if(CLN_){
      std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the CLN parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(CLN_ORI_){
	std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the CLN_ORI parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
        exit (EXIT_FAILURE);
      }
      else
	if(CLN_HYB_){
	  std::cout << "You can not use a CLN_MOD setter method for an SL_Decay object with the CLN_HYB parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
	  exit (EXIT_FAILURE);
	}
	else
	  if(BGL_){
	    std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the BGL parameterization." << std::endl;
	    std::cout << "Please try again." << std::endl;
	    exit (EXIT_FAILURE);
	  }
	  else
	    if(BCL_){
	      std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the BCL parameterization." << std::endl;
	      std::cout << "Please try again." << std::endl;
	      exit (EXIT_FAILURE);
	    }
}

void SL_Decay::Set_cz3(double coeff_z3){
  if(CLN_MOD_){ c_z3_ = coeff_z3; }
  else
    if(CLN_){
      std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the CLN parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(CLN_ORI_){
	std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the CLN_ORI parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
        exit (EXIT_FAILURE);
      }
      else
	if(CLN_HYB_){
	  std::cout << "You can not use a CLN_MOD setter method for an SL_Decay object with the CLN_HYB parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
	  exit (EXIT_FAILURE);
	}
	else
	  if(BGL_){
	    std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the BGL parameterization." << std::endl;
	    std::cout << "Please try again." << std::endl;
	    exit (EXIT_FAILURE);
	  }
	  else
	    if(BCL_){
	      std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the BCL parameterization." << std::endl;
	      std::cout << "Please try again." << std::endl;
	      exit (EXIT_FAILURE);
	    }
}

double SL_Decay::Get_cz2() {
  double c_z2 = 0.; 
  if(CLN_MOD_) { c_z2 = c_z2_; }
  else
    if(CLN_){
      std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the CLN parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      c_z2 = 1.;
      exit (EXIT_FAILURE);
    }
    else
      if(CLN_ORI_){
	std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the CLN_ORI parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
	c_z2 = 1.;
	exit (EXIT_FAILURE);
      }
      else
        if(CLN_HYB_){
	  std::cout << "You can not use a CLN_MOD setter method for an SL_Decay object with the CLN_HYB parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
	  c_z2 = 1.;
	  exit (EXIT_FAILURE);
        }
        else
          if(BGL_){
	    std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the BGL parameterization." << std::endl;
	    std::cout << "Please try again." << std::endl;
	    c_z2 = 1.;
	    exit (EXIT_FAILURE);
          }
          else
            if(BCL_){
	      std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the BCL parameterization." << std::endl;
	      std::cout << "Please try again." << std::endl;
	      c_z2 = 1.;
	      exit (EXIT_FAILURE);
            }

  return c_z2; 
}

double SL_Decay::Get_cz3() {
  double c_z3 = 0.; 
  if(CLN_MOD_) { c_z3 = c_z3_; } 
  else
    if(CLN_){
      std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the CLN parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      c_z3 = 1.;
      exit (EXIT_FAILURE);
    }
    else
      if(CLN_ORI_){
	std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the CLN_ORI parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
	c_z3 = 1.;
        exit (EXIT_FAILURE);
      }
      else
        if(CLN_HYB_){
	  std::cout << "You can not use a CLN_MOD setter method for an SL_Decay object with the CLN_HYB parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
	  c_z3 = 1.;
          exit (EXIT_FAILURE);
        }
        else
          if(BGL_){
	    std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the BGL parameterization." << std::endl;
	    std::cout << "Please try again." << std::endl;
	    c_z3 = 1.;
            exit (EXIT_FAILURE);
          }
          else
            if(BCL_){
	      std::cout << "You can not use a CLN_MOD setter function for an SL_Decay object with the BCL parameterization." << std::endl;
	      std::cout << "Please try again." << std::endl;
	      c_z3 = 1.;
              exit (EXIT_FAILURE);
            }

  return c_z3;
}

// ------------------------------- //
// BGL parameters setter functions //
// ------------------------------- //
void SL_Decay::Set_BGL_f(double a_f_0, double a_f_1, double a_f_2, double a_f_3){
  if(BGL_ || CLN_HYB_){
    a_f_BGL_[0] = a_f_0;
    a_f_BGL_[1] = a_f_1; 
    a_f_BGL_[2] = a_f_2;
    a_f_BGL_[3] = a_f_3;
    a_F1_BGL_[0] = a_F1_0_BGL(); // a_F1_[0] is a function of a_f_[0]. 
  }
  else
    if(CLN_){
      std::cout << "You can not use a BGL setter function for an SL_Decay object with the CLN parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(CLN_MOD_){
	std::cout << "You can not use a BGL setter method for an SL_Decay object with the CLN_MOD parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
	exit (EXIT_FAILURE);
      }
      else
	if(CLN_ORI_){
	  std::cout << "You can not use a BGL setter function for an SL_Decay object with the CLN_ORI parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
	  exit (EXIT_FAILURE);
	}
  
	else
	  if(BCL_){
	    std::cout << "You can not use a BGL setter function for an SL_Decay object with the BCL parameterization." << std::endl;
	    std::cout << "Please try again." << std::endl;
	    exit (EXIT_FAILURE);
	  }
  
}

void SL_Decay::Set_BGL_g(double a_g_0, double a_g_1, double a_g_2, double a_g_3){
  if(BGL_){
  a_g_BGL_[0] = a_g_0;
  a_g_BGL_[1] = a_g_1;
  a_g_BGL_[2] = a_g_2;
  a_g_BGL_[3] = a_g_3;
  }
  else
    if(CLN_){
      std::cout << "You can not use a BGL setter function for an SL_Decay object with the CLN parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(CLN_MOD_){
	std::cout << "You can not use a BGL setter method for an SL_Decay object with the CLN_MOD parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
        exit (EXIT_FAILURE);
      }
      else
        if(CLN_ORI_){
	  std::cout << "You can not use a BGL setter function for an SL_Decay object with the CLN_ORI parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
          exit (EXIT_FAILURE);
        }
	else
          if(CLN_HYB_){
	    std::cout << "You can not use a BGL setter function for an SL_Decay object with the CLN_HYB parameterization." << std::endl;
	    std::cout << "Please try again." << std::endl;
            exit (EXIT_FAILURE);
          }
	  else
	    if(BCL_){
	      std::cout << "You can not use a BGL setter function for an SL_Decay object with the BCL parameterization." << std::endl;
	      std::cout << "Please try again." << std::endl;
	      exit (EXIT_FAILURE);
	    } 
  
} 

void SL_Decay::Set_BGL_F1(double a_F1_0, double a_F1_1, double a_F1_2, double a_F1_3){
  if(BGL_){
    a_F1_BGL_[0] = a_F1_0;
    a_F1_BGL_[1] = a_F1_1;
    a_F1_BGL_[2] = a_F1_2;
    a_F1_BGL_[3] = a_F1_3;
    a_F1_BGL_[0] = a_F1_0_BGL(); // a_F1_[0] is a function of a_f_[0].
  }
  else
    if(CLN_){
      std::cout << "You can not use a BGL setter function for an SL_Decay object with the CLN parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(CLN_MOD_){
	std::cout << "You can not use a BGL setter method for an SL_Decay object with the CLN_MOD parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
        exit (EXIT_FAILURE);
      }
      else
        if(CLN_ORI_){
	  std::cout << "You can not use a BGL setter function for an SL_Decay object with the CLN_ORI parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
          exit (EXIT_FAILURE);
        }
	else
          if(CLN_HYB_){
	    std::cout << "You can not use a BGL setter function for an SL_Decay object with the CLN_HYB parameterization." << std::endl;
	    std::cout << "Please try again." << std::endl;
            exit (EXIT_FAILURE);
          }
	  else
	    if(BCL_){
	      std::cout << "You can not use a BGL setter function for an SL_Decay object with the BCL parameterization." << std::endl;
	      std::cout << "Please try again." << std::endl;
	      exit (EXIT_FAILURE);
	    }
  
}

void SL_Decay::Set_BGL_F2(double a_F2_0, double a_F2_1, double a_F2_2, double a_F2_3){
  if(BGL_){
    a_F2_BGL_[0] = a_F2_0;
    a_F2_BGL_[1] = a_F2_1;
    a_F2_BGL_[2] = a_F2_2;
    a_F2_BGL_[3] = a_F2_3;
  }
  else
    if(CLN_){
      std::cout << "You can not use a BGL setter functions for an SL_Decay object with the CLN parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(CLN_MOD_){
	std::cout << "You can not use a BGL setter method for an SL_Decay object with the CLN_MOD parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
        exit (EXIT_FAILURE);
      }
      else
        if(CLN_ORI_){
	  std::cout << "You can not use a BGL setter function for an SL_Decay object with the CLN_ORI parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
          exit (EXIT_FAILURE);
        }
	else
          if(CLN_HYB_){
	    std::cout << "You can not use a BGL setter function for an SL_Decay object with the CLN_HYB parameterization." << std::endl;
	    std::cout << "Please try again." << std::endl;
            exit (EXIT_FAILURE);
          }
	  else
	    if(BCL_){
	      std::cout << "You can not use a BGL setter functions for an SL_Decay object with the BCL parameterization." << std::endl;
	      std::cout << "Please try again." << std::endl;
	      exit (EXIT_FAILURE);
	    }
  
}

// void SL_Decay::Set_BGL_chi_T_1p(double x) { chi_T_1p_ = x; }
// void SL_Decay::Set_BGL_chi_T_1m(double x) { chi_T_1m_ = x; }
// void SL_Decay::Set_BGL_chi_t_0m(double x) { chi_t_0m_ = x;  }


// -------------------- //
// BCL Setter functions //
// -------------------- //
void SL_Decay::Set_BCL_f(double a_f_0, double a_f_1, double a_f_2, double a_f_3){
  if(BCL_){
    a_f_BCL_[0] = a_f_0;
    a_f_BCL_[1] = a_f_1;
    a_f_BCL_[2] = a_f_2;
    a_f_BCL_[3] = a_f_3;
}
  else
    if(CLN_){
      std::cout << "You can not use a BCL setter function for an SL_Decay object with the CLN parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(CLN_MOD_){
	std::cout << "You can not use a BCL stter method for an SL_Decay object with the CLN_MOD parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
        exit (EXIT_FAILURE);
      }
      else
        if(CLN_ORI_){
	  std::cout << "You can not use a BCL setter function for an SL_Decay object with the CLN_ORI parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
          exit (EXIT_FAILURE);
        }
	else
          if(CLN_HYB_){
	    std::cout << "You can not use a BCL setter function for an SL_Decay object with the CLN_HYB parameterization." << std::endl;
	    std::cout << "Please try again." << std::endl;
            exit (EXIT_FAILURE);
          }
	  else
	    if(BGL_){
	      std::cout << "You can not use a BCL setter function for an SL_Decay object with the BGL parameterization." << std::endl;
	      std::cout << "Please try again." << std::endl;
	      exit (EXIT_FAILURE);
	    }
  
}

void SL_Decay::Set_BCL_g(double a_g_0, double a_g_1, double a_g_2, double a_g_3){
  if(BCL_){
    a_g_BCL_[0] = a_g_0;
    a_g_BCL_[1] = a_g_1;
    a_g_BCL_[2] = a_g_2;
    a_g_BCL_[3] = a_g_3;
  }
  else
    if(CLN_){
      std::cout << "You can not use a BCL setter function for an SL_Decay object with the CLN parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(CLN_MOD_){
	std::cout << "You can not use a BCL setter method for an SL_Decay object with the CLN_MOD parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
        exit (EXIT_FAILURE);
      }
      else
        if(CLN_ORI_){
	  std::cout << "You can not use a BCL setter function for an SL_Decay object with the CLN_ORI parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
          exit (EXIT_FAILURE);
        }
	else
	  if(CLN_HYB_){
	    std::cout << "You can not use a BCL setter function for an SL_Decay object with the CLN_HYB parameterization." << std::endl;
	    std::cout << "Please try again." << std::endl;
	    exit (EXIT_FAILURE);
	  }
	  else
	    if(BGL_){
	      std::cout << "You can not use a BCL setter function for an SL_Decay object with the BGL parameterization." << std::endl;
	      std::cout << "Please try again." << std::endl;
	      exit (EXIT_FAILURE);
	    }
  
}

void SL_Decay::Set_BCL_F1(double a_F1_0, double a_F1_1, double a_F1_2, double a_F1_3){
  if(BCL_){
    a_F1_BCL_[0] = a_F1_0;
    a_F1_BCL_[1] = a_F1_1;
    a_F1_BCL_[2] = a_F1_2;
    a_F1_BCL_[3] = a_F1_3;
  }
  else
    if(CLN_){
      std::cout << "You can not use a BCL setter function for an SL_Decay object with the CLN parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(CLN_MOD_){
	std::cout << "You can not use a BCL setter method for an SL_Decay object with the CLN_MOD parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
        exit (EXIT_FAILURE);
      }
      else
        if(CLN_ORI_){
	  std::cout << "You can not use a BCL setter function for an SL_Decay object with the CLN_ORI parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
          exit (EXIT_FAILURE);
        }
	else
          if(CLN_HYB_){
	    std::cout << "You can not use a BCL setter function for an SL_Decay object with the CLN_HYB parameterization." << std::endl;
	    std::cout << "Please try again." << std::endl;
            exit (EXIT_FAILURE);
          }
	  else
	    if(BGL_){
	      std::cout << "You can not use a BCL setter function for an SL_Decay object with the BGL parameterization." << std::endl;
	      std::cout << "Please try again." << std::endl;
	      exit (EXIT_FAILURE);
	    }
  
}

void SL_Decay::Set_BCL_F2(double a_F2_0, double a_F2_1, double a_F2_2, double a_F2_3){
  if(BCL_){
    a_F2_BCL_[0] = a_F2_0;
    a_F2_BCL_[1] = a_F2_1;
    a_F2_BCL_[2] = a_F2_2;
    a_F2_BCL_[3] = a_F2_3;
  }
  else
    if(CLN_){
      std::cout << "You can not use a BCL setter function for an SL_Decay object with the CLN parameterization." << std::endl;
      std::cout << "Please try again." << std::endl;
      exit (EXIT_FAILURE);
    }
    else
      if(CLN_MOD_){
	std::cout << "You can not use a BCL setter method for an SL_Decay object with the CLN_MOD parameterization." << std::endl;
	std::cout << "Please try again." << std::endl;
        exit (EXIT_FAILURE);
      }
      else
        if(CLN_ORI_){
	  std::cout << "You can not use a BCL setter function for an SL_Decay object with the CLN_ORI parameterization." << std::endl;
	  std::cout << "Please try again." << std::endl;
          exit (EXIT_FAILURE);
        }
	else
          if(CLN_HYB_){
	    std::cout << "You can not use a BCL setter function for an SL_Decay object with the CLN_HYB parameterization." << std::endl;
	    std::cout << "Please try again." << std::endl;
            exit (EXIT_FAILURE);
          }
	  else
	    if(BGL_){
	      std::cout << "You can not use a BCL setter function for an SL_Decay object with the BGL parameterization." << std::endl;
	      std::cout << "Please try again." << std::endl;
	      exit (EXIT_FAILURE);
	    }
  
}

// Setter and Getter functions for the class objects.
void SL_Decay::Set_Q2(double Q2) { Q2_ = Q2; }
void SL_Decay::Set_cos_theta_l(double cos_theta_l) { cos_theta_l_ = cos_theta_l; }
void SL_Decay::Set_cos_theta_V(double cos_theta_V) { cos_theta_V_ = cos_theta_V; }
void SL_Decay::Set_k(double k) { k_ = k; } 
double SL_Decay::Get_Q2() { return Q2_;  }
double SL_Decay::Get_cos_theta_l() { return cos_theta_l_; }
double SL_Decay::Get_cos_theta_V() { return cos_theta_V_; }
double SL_Decay::Get_k() { return k_; }

double SL_Decay::Get_Q2_max() {
  double Q2_max = pow( (B_mass_ - D_st_mass_), 2); 
  return Q2_max;
}

double SL_Decay::Get_Q2_min() {
  double Q2_min = 0.;
  if(CLN_ORI_){                             // CLN_ORI paramtrization 
    if(lepton_index_ == 1) { Q2_min = pow( (0.510999 * pow(10, -3)), 2); } // electron
    else
      if(lepton_index_ == 2) { Q2_min = pow( 0.10566, 2); }                // muon
  } 
  else { Q2_min = pow (lepton_mass_, 2); }  // All other parametrizations
  
  return Q2_min;
}

std::string SL_Decay::Get_parameterization(){
  std::string param;
  if(CLN_) { param = "CLN"; }
  else
    if(BGL_) { param = "BGL"; }
    else
      if(BCL_) {param = "BCL"; }
      else
	if(CLN_MOD_) {param = "CLN_MOD";}
	else
	  if(CLN_ORI_) {param = "CLN_ORI";}
	  else
	    if(CLN_HYB_) {param = "CLN_HYB";}
  
  return param;
}

double SL_Decay::Get_BR(){ return BR_; }

//   ********************************************   //
//                                                  //
//            PRIVATE MEMBER FUNCTIONS              //
//                                                  // 
//   ********************************************   //



// I1s, I1c, I2s, I2c, I3, I4, I5, I6s, I6c, I7 function implementation.
double SL_Decay::I_1s() {
  double func_I_1s;
  if(is_pi) {      // I1s for pions.
    func_I_1s = ( pow(f_H_p(), 2) + pow(f_H_m(), 2) )*( pow(lepton_mass_, 2) + 3*Q2_ )/2;
  }
  else{            // I1s for photons.
    func_I_1s = 2*pow(lepton_mass_, 2)*pow(f_H_t(), 2) + pow(f_H_0(), 2)*( pow(lepton_mass_, 2) + Q2_ ); 
  }

  return func_I_1s;
} 

double SL_Decay::I_1c() {
  double func_I_1c;
  if(is_pi) {      // Ics for pions.
    func_I_1c = 2*( 2*pow(lepton_mass_, 2)*pow(f_H_t(), 2) + pow(f_H_0(), 2)*( pow(lepton_mass_, 2) + Q2_ ) );
  }
  else{            // Ics for photons.
    func_I_1c = ( pow(f_H_p(), 2) + pow(f_H_m(), 2) )*( pow(lepton_mass_, 2) + 3*Q2_ )/8;
  }

  return func_I_1c;
} 

double SL_Decay::I_2s() {
  double func_I_2s;
  if(is_pi) {      // I2s for pions.
    func_I_2s = ( pow(f_H_p(), 2) + pow(f_H_m(), 2)  )*( Q2_ - pow(lepton_mass_, 2) )/2;
  }
  else{            // I2s for photons.
    func_I_2s = pow(f_H_0(), 2)*( pow(lepton_mass_, 2) - Q2_ );
  }

  return func_I_2s;
} 

double SL_Decay::I_2c() {
  double func_I_2c;
  if(is_pi) {      // I2c for pions.
    func_I_2c = 2*pow(f_H_0() ,2)*( pow(lepton_mass_, 2) - Q2_ );
  }
  else{            // I2c for photons.
    func_I_2c = ( pow(f_H_p(), 2) + pow(f_H_m(), 2) )*( Q2_ - pow(lepton_mass_, 2) )/8;
  }

  return func_I_2c;
} 

double SL_Decay::I_3() {
  double func_I_3;
  if(is_pi) {      // I3 for pions.
    func_I_3 = 2*f_H_p()*f_H_m()*( pow(lepton_mass_, 2) - Q2_ );
  }
  else{            // I3 for photons.
    func_I_3 = -f_H_p()*f_H_m()*( pow(lepton_mass_, 2) - Q2_ );
  }

  return func_I_3;
} 

double SL_Decay::I_4() {
  double func_I_4;
  if(is_pi) {      // I4 for pions.
    func_I_4 = f_H_0()*( f_H_p() + f_H_m() )*( pow(lepton_mass_, 2) - Q2_ );
  }
  else{            // I4 for photons.
    func_I_4 = -f_H_0()*( f_H_p() + f_H_m() )*( pow(lepton_mass_, 2) - Q2_ )/2;
  }

  return func_I_4;
} 

double SL_Decay::I_5() {
  double func_I_5;
  double add1, add2;
  if(is_pi) {      // I5 for pions.
    add1 = -2*( f_H_p() + f_H_m() )*f_H_t()*pow(lepton_mass_, 2);
    add2 = -2*f_H_0()*( f_H_p() - f_H_m() )*Q2_;
  }
  else{            // I5 for photons.
    add1 = ( f_H_p() + f_H_m() )*f_H_t()*pow(lepton_mass_, 2);
    add2 = f_H_0()*( f_H_p() - f_H_m() )*Q2_;
  }

  func_I_5 = add1 + add2;
  return func_I_5;
} 

double SL_Decay::I_6s() {
  double func_I_6s;
  if(is_pi) {      // I6s for pions.
    func_I_6s = 2*( pow(f_H_p(), 2) - pow(f_H_m(), 2) )*Q2_;
  }
  else{            // I6s for photons.
    func_I_6s = -4*f_H_0()*f_H_t()*pow(lepton_mass_, 2);
  }

  return func_I_6s;
} 

double SL_Decay::I_6c() {
  double func_I_6c;
  if(is_pi) {      // I6c for pions.
    func_I_6c = -8*f_H_0()*f_H_t()*pow(lepton_mass_, 2);
  }
  else{            // I6c for photons.
    func_I_6c = ( pow(f_H_p(), 2) - pow(f_H_m(), 2) )*Q2_/2;
  }

  return func_I_6c;
} 

double SL_Decay::I_7() {
  double func_I_7;
  if(is_pi) {      // I7 for pions.
    func_I_7 = 0.;
  }
  else{            // I7 for photons.
    func_I_7 = 0.;
  }

  return func_I_7;
} 



// H0, H+, H- and Ht functions implementation.
double SL_Decay::f_H_0() {

  double func_H0 = 0;
  // H_0 definition for CLN parameterization.
  if(CLN_ || CLN_MOD_ || CLN_ORI_ || CLN_HYB_){
    double num = pow( (B_mass_ + D_st_mass_), 2)*( pow(B_mass_, 2) - pow(D_st_mass_, 2) - Q2_  )*f_A1() - lambda_func( pow(B_mass_, 2), pow(D_st_mass_, 2), Q2_ )*f_A2();
    double den = 2*D_st_mass_*(B_mass_ + D_st_mass_)*sqrt(Q2_);
    func_H0 = num/den;
  }
  else
    if(BGL_){  // H_0 definition for BGL parameterization. 
      func_H0 = F1()/sqrt(Q2_);
    }
    else{      // H_0 definition for BCL parameterization.
      func_H0 = F1()/sqrt(Q2_);
    }

  return func_H0;
}

double SL_Decay::f_H_p() {

  double func_Hp = 0;
  // H+ definition for CLN parameterization.
  if(CLN_ || CLN_MOD_ || CLN_ORI_ || CLN_HYB_){
    double num = pow( (B_mass_ + D_st_mass_), 2)*f_A1() - sqrt( lambda_func( pow(B_mass_, 2), pow(D_st_mass_, 2), Q2_ ) )*f_V();
    double den = B_mass_ + D_st_mass_;
    func_Hp = num/den;
  }
  else
    if(BGL_){  // H+ definition for BGL parameterization.
      func_Hp = f() - B_mass_*D_st_mass_*sqrt( pow(w(), 2) -1 )*g();
    }
    else{      // H+ definition for BCL parameterization.
      func_Hp = f() - B_mass_*D_st_mass_*sqrt( pow(w(), 2) -1 )*g();
    }

  return func_Hp;
}

double SL_Decay::f_H_m() {

  double func_Hm = 0;
  if(CLN_ || CLN_MOD_ || CLN_ORI_ || CLN_HYB_){    // H- definition for CLN parametrization.
    double num = pow( (B_mass_ + D_st_mass_), 2)*f_A1() + sqrt( lambda_func( pow(B_mass_, 2), pow(D_st_mass_, 2), Q2_ ) )*f_V();
    double den = B_mass_ + D_st_mass_;
    func_Hm = num/den;
  }
  else
    if(BGL_){  // H- definition for BGL parametrization. 
      func_Hm = f() + B_mass_*D_st_mass_*sqrt( pow(w(), 2) -1 )*g();
    }
    else{      // H- definition for BCL parametrization.
      func_Hm = f() + B_mass_*D_st_mass_*sqrt( pow(w(), 2) -1 )*g();
    }

  return func_Hm;
}

double SL_Decay::f_H_t() {

  double func_Ht = 0;
  if(CLN_ || CLN_MOD_ || CLN_ORI_ || CLN_HYB_){   // H_t definition for CLN parametrization.
    func_Ht = -sqrt( lambda_func( pow(B_mass_, 2), pow(D_st_mass_, 2), Q2_ ) )*f_A0() / sqrt(Q2_);
  }
  else
    if(BGL_){  // H_t definition for BGL parametrization.
      func_Ht = -F2()*B_mass_*( sqrt(r())*(1+r())*sqrt(pow(w(), 2) - 1) )/ sqrt(1 + pow(r(), 2) - 2*w()*r());
    }
    else{   // H_t definition for BCL parametrization.
      func_Ht = F2();
    }

  return func_Ht;
}



//  ****************************************************  //
//                                                        //                                                                                                                                              
//      Form factors used only for B -> D lnu decays      //
//                                                        //                                                                                                                                              
//  ****************************************************  //

double SL_Decay::f_p(){
  double f_p = 0.;
  
  if(CLN_ || CLN_MOD_ || CLN_ORI_ || CLN_HYB_){
    f_p = h_A1_ * ( 1 - 8*rho_squared_*z() + (51*rho_squared_ -10)*pow(z(), 2) - (252*rho_squared_ -84)*pow(z(), 3) ); 
  }
  else
    if(CLN_MOD_){
      f_p = h_A1_ * ( 1 - 8*rho_squared_*z() + c_z2_*pow(z(), 2) + c_z3_*pow(z(), 3)  );
    }
    else
      if(BGL_){
	f_p = 1/(P_p()*phi_p())*( a_BGL_[0] + a_BGL_[1]*z() + a_BGL_[2]*pow(z(), 2) - a_BGL_[3]*pow(z(), 3) );
      }
      else
	if(BCL_){
	  int N = 4; // Truncation index of the series.
	  for(int k=0; k<N; k++) {
	    f_p += ( 1/(1 - Q2_/pow(m_B_st_c[1][0], 2)) )*( a_BCL_[k]*( pow(z(), k) - pow(-1, k-N-1)*(k/(N+1))*pow(z(), N+1) )  );
	  }
	}
  
  return f_p;
}


double SL_Decay::f_0(){
  double f_0 = 0.;

  if(CLN_ || CLN_MOD_ || CLN_ORI_ || CLN_HYB_){
    f_0 = f_p() * pow( 2*sqrt(r())/(1+r()) ,2) * (1 + w())/2 * (1 - 0.0068*(1-w()) + 0.0017*pow(1-w(), 2) - 0.0013*pow(1- w(), 3) );
  }
  else
    if(BGL_){
      f_0 = 1/(P_0()*phi_0())*( b_BGL_[0] + b_BGL_[1]*z() + b_BGL_[2]*pow(z(), 2) + b_BGL_[3]*pow(z(), 3) );
    }
    else
      if(BCL_){
	int N =4; // Truncation index of the series.
	for(int k=0; k<N; k++) {
	  f_0 += ( 1/(1 - Q2_/pow(m_B_st_c_0p[0], 2)) ) * b_BCL_[k]*pow(z(), k);
	}

      }

  return f_0;
}

// A1, A2, A0 and V function implementation.
double SL_Decay::R_star() {
  double func_R_star = 2*sqrt( B_mass_*D_st_mass_ ) / (B_mass_ + D_st_mass_);
  return func_R_star;
}


double SL_Decay::f_A1() {
  double func_A1 = 0.;
  if(CLN_ || CLN_ORI_ || CLN_HYB_) { func_A1 = (w()+1)*R_star()*f_h_A1()/2; }
  else
    if(CLN_MOD_) {func_A1 = (w()+1)*R_star()*f_h_A1_mod()/2; }
  return func_A1;
}

double SL_Decay::f_A2() {
  double func_A2 = f_R2()*f_h_A1()/R_star();
  return func_A2;
}

double SL_Decay::f_A0() {
  double func_A0 = f_R0()*f_h_A1()/R_star();
  return func_A0;
}

double SL_Decay::f_V() {
  double func_V = f_R1()*f_h_A1()/R_star();
  return func_V;
}



// CLN Paramatrizazion implemantation
double SL_Decay::w(){
  double w = 0.; 
  if(BtoDstMuNu) { w = ( pow(B_mass_, 2) + pow(D_st_mass_, 2) - Q2_ ) / ( 2*B_mass_*D_st_mass_ ); }
  else
    if(BtoDMuNu) { w = ( pow(B_mass_, 2) + pow(D_mass_, 2) - Q2_ ) / ( 2*B_mass_*D_mass_ ); }
  return w;
}

double SL_Decay::z(){
  double z = ( sqrt(w()+1) - sqrt(2) ) / ( sqrt(w()+1) + sqrt(2)  );
  return z;
}

double SL_Decay::f_h_A1(){
  double func_h_A1 = 0.;
  if(CLN_ || CLN_ORI_ || CLN_MOD_) {
    func_h_A1 = h_A1_*( 1 - 8*rho_squared_*z() + (53*rho_squared_ - 15)*pow(z(), 2) - (231*rho_squared_ - 91)*pow(z(), 3)  );
  }
  else 
    if(CLN_HYB_){ func_h_A1 = 1/(sqrt(B_mass_*D_st_mass_)*(1+w())) * f(); }
  return func_h_A1;
}

double SL_Decay::f_R1() {
  double func_R1 = R1_ - 0.12*(w() - 1) + 0.05*pow( (w()-1), 2);
  return func_R1;
}

double SL_Decay::f_R2() {
  double func_R2 = R2_ + 0.11*(w() - 1) - 0.06*pow( (w()-1), 2);
  return func_R2;
}

double SL_Decay::f_R0() {
  double func_R0 = R0_ - 0.11*(w() - 1) + 0.01*pow( (w()-1), 2);
  return func_R0;
}

// CLN_MOD Parametrization
double SL_Decay::f_h_A1_mod(){
  double func_h_A1_mod = h_A1_*( 1 - 8*rho_squared_*z() + c_z2_*pow(z(), 2) - c_z3_*pow(z(), 3)  );
  return func_h_A1_mod;
}

// Function Miscellaneous.
double SL_Decay::Nf()
{
  double Nf = 3*pow(G_Fermi, 2)*pow(Vcb_, 2)*BR_ / ( 128*pow(2*M_PI, 4)*pow(B_mass_, 2) ); // M_PI def. in cmath.
  return Nf;
}

double SL_Decay::lambda_func(double x, double y, double z)
{
  double lambda;
  lambda = x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z;
  return lambda;
}

double SL_Decay::D_st_P_module() {

  double D_st_mom_module = D_st_mass_ * sqrt( w()*w() -1 );
  if ( Q2_ > pow(lepton_mass_, 2) ) { return D_st_mom_module; }
  else { return 0; };

}


//  **********************************  //
//                                      //
//     BGL and BCL Member Function      //
//                                      //
//  **********************************  //


double SL_Decay::f() {
  double func_f = 0.;
  if(BGL_ || CLN_HYB_){  // BGL parameterization. 
    func_f = ( 1/(P1_p()*phi_f()) )*( a_f_BGL_[0] + a_f_BGL_[1]*z() + a_f_BGL_[2]*pow(z(), 2) + a_f_BGL_[3]*pow(z(), 3) );
  }
  else
    if(BCL_){  // BCL parameterization.
      for(int k=0; k<4; k++){  // in the BCL formula K = 4. https://arxiv.org/abs/0807.2722v3
	func_f += a_f_BCL_[k] * pow(z(), k);
	// func_f += a_f_BCL_[k]*( pow(z(), k) - pow(-1, k-4)*(k/4)*pow(z(), 4) );
      }
      func_f *= Q_f();
    }

  return func_f;
}

double SL_Decay::F1() {
  double func_F1 = 0.;
  if(BGL_){  // BGL parameterization.
    func_F1 = ( 1/(P1_p()*phi_F1()) )*( a_F1_BGL_[0] + a_F1_BGL_[1]*z() + a_F1_BGL_[2]*pow(z(), 2) + a_F1_BGL_[3]*pow(z(), 3) );
  }
  else
    if(BCL_){  // BCL parameterization.
      for(int k=0; k<4; k++){  // in the BCL formula K = 4. https://arxiv.org/abs/0807.2722v3
	func_F1 += a_F1_BCL_[k] * pow(z(), k);
	// func_F1 += a_F1_BCL_[k]*( pow(z(), k) - pow(-1, k-4)*(k/4)*pow(z(), 4) );
      }
      func_F1 *= Q_F1();
    }

  return func_F1;
}

double SL_Decay::g() {
  double func_g = 0.;
  if(BGL_){  // BGL parameterization.
    func_g = ( 1/(P1_m()*phi_g()) )*( a_g_BGL_[0] + a_g_BGL_[1]*z() + a_g_BGL_[2]*pow(z(), 2) + a_g_BGL_[3]*pow(z(), 3) );
  }
  else
    if(BCL_){  // BCL parameterization.
      for(int k=0; k<4; k++){  // in the BCL formula K = 4. https://arxiv.org/abs/0807.2722v3 
	func_g += a_g_BCL_[k] * pow(z(), k);
	// func_g += a_g_BCL_[k]*( pow(z(), k) - pow(-1, k-4)*(k/4)*pow(z(), 4) );
      }
      func_g *= Q_g();
    }
  
  return func_g;
}

double SL_Decay::F2() {
  double func_F2 = 0.;
  if(BGL_){  // BGL parameterization.
    func_F2 = ( sqrt(r())/((1+r())*P2()*phi_F2()) )*( a_F2_BGL_[0] + a_F2_BGL_[1]*z() + a_F2_BGL_[2]*pow(z(), 2) + a_F2_BGL_[3]*pow(z(), 3) );
  } 
  else
    if(BCL_){  // BCL parameterization.
      for(int k=0; k<4; k++){  // in the BCL formula K = 4. https://arxiv.org/abs/0807.2722v3
	func_F2 += a_F2_BCL_[k] * pow(z(), k);
	// func_F2 += a_F2_BCL_[k]*( pow(z(), k) - pow(-1, k-4)*(k/4)*pow(z(), 4) ); 
      }
      func_F2 *= Q_F2();
    }

  return func_F2;
}

// BCL Member function
double SL_Decay::Q_f() {
  double func_Qf = N_f_ / ( 1 - Q2_/pow(m_B_st_c[0][0], 2) );
  return func_Qf;
}

double SL_Decay::Q_F1() {
  double func_QF1 = N_F1_ / ( 1 - Q2_/pow(m_B_st_c[0][0], 2) );
  return func_QF1;
}

double SL_Decay::Q_g() {
  double func_Qg = N_g_ / ( 1 - Q2_/pow(m_B_st_c[1][0], 2) );
  return func_Qg;
}

double SL_Decay::Q_F2() {
  double func_QF2 = N_F2_ / ( 1 - Q2_/pow(m_B_st_c_0m[0], 2) );
  return func_QF2;
}
// End of BCL Member function implemetation


//  *******************************************  //
//               Outer Functions                 //
//  *******************************************  //

double SL_Decay::phi_g() {
  double func_phi_g;
  double constant = sqrt( n_I/( 3*M_PI*chi_T_1m_ ) );
  double num = 16*pow(r(), 2)*pow(1+z(), 2)*pow(1-z(), -1./2.);
  double den = pow( ( (1+r())*(1-z()) + 2*sqrt(r())*(1+z()) ), 4);
  func_phi_g = constant*num/den;
  
  return func_phi_g;
}

double SL_Decay::phi_f() {
  double func_phi_f;
  double constant = ( 4*r()/pow(B_mass_, 2) )*sqrt( n_I/( 3*M_PI*chi_T_1p_ ) );
  double num = (1+z())*pow(1-z(), 3./2.);
  double den = pow( ( (1+r())*(1-z()) + 2*sqrt(r())*(1+z()) ), 4);
  func_phi_f = constant*num/den;

  return func_phi_f;
}

double SL_Decay::phi_F1() {
  double func_phi_F1;
  double constant = ( 4*r()/pow(B_mass_, 3) )*sqrt( n_I/( 6*M_PI*chi_T_1p_ ) );
  double num = (1+z())*pow(1-z(), 5./2.);
  double den = pow( ( (1+r())*(1-z()) + 2*sqrt(r())*(1+z()) ), 5);
  func_phi_F1 = constant*num/den;

  return func_phi_F1;
}

double SL_Decay::phi_F2() {
  double func_phi_F2;
  double constant = sqrt( n_I/(M_PI*chi_L_1p_) );
  double num = 8*sqrt(2)*pow(r(), 2)*pow(1+z(), 2);
  double den = sqrt(1-z())*pow( ( (1+r())*(1-z()) + 2*sqrt(r())*(1+z()) ),  4);
  func_phi_F2 = constant*num/den;

  return func_phi_F2;
}

double SL_Decay::phi_p(){
  double phi_p_func;
  double k_p = (8*pow(r(), 2)/B_mass_)*sqrt(8*n_I/(3*M_PI*chi_Tt_0_));
  double num = (1 + pow(z(), 2))*sqrt(1-z());
  double den = pow( ( (1+r())*(1-z()) + 2*sqrt(r())*(1+z()) ), 5);
  phi_p_func = k_p*num/den;
  
  return phi_p_func;
}

double SL_Decay::phi_0(){
  double phi_0_func;
  double k_0 = r()*(1-pow(r(), 2))*sqrt(8*n_I/(M_PI*chi_L_0_));
  double num = (1 - pow(z(), 2))*sqrt(1-z());
  double den = pow( ( (1+r())*(1-z()) + 2*sqrt(r())*(1+z()) ), 4);
  phi_0_func = k_0*num/den;

  return phi_0_func;
}

//  *******************************************  //
//              Blaschke Factors                 //
//  *******************************************  //

double SL_Decay::P1_p() {

  int i = 0;  // index to select mass states with +1 quantum number.
  double func_P1_p = 1.;
  double norm;
  for(int j=0; j<n_reson_1; j++) {
    func_P1_p *= ( z() - z_P(i,j) )/( 1 - z()*z_P(i,j) );
  }
  if(normailization_index_ == 1 || normailization_index_ == 2 ) { norm = 1.; } // No changes of the Blaschke factor for B0 and B+-.
  else{ norm = P_1p_coeff; } // Scale Bs Blascke factor. See documentation ... Value set in the SL_Decay_Config.cpp file.
  func_P1_p *= norm;  // func_P1_p *= norm;

  return func_P1_p;
}

double SL_Decay::P1_m() {

  int i = 1;  // index to select mass states with -1 quantum number.
  double func_P1_m = 1.;
  double norm;
  int n_resonances;
  if(normailization_index_ == 1 || normailization_index_ == 2 ) { n_resonances = n_reson_1 -1; } // The heaviest resonance is eliminated for the B0 and B+- case because is too close to threshold.
  else{ n_resonances = n_reson_1; } // The Bs has a higher threshold and therefore also the heaviest resonance in included.

  for(int j=0; j<n_resonances; j++) { 
    func_P1_m *= ( z() - z_P(i,j) )/( 1 - z()*z_P(i,j) );
  }
  if(normailization_index_ == 1 || normailization_index_ == 2 ) { norm = 1.; } // No changes of the Blaschke factor for B0 and B+-.
  else{ norm = P_1m_coeff; } // Scale Bs Blascke factor. See documentation ... Value set in the SL_Decay_Config.cpp file.
  func_P1_m *= norm;  // func_P1_m *= norm;

  return func_P1_m;
}

double SL_Decay::P2() {
  double func_P2 = 1.;
  double norm;
  for(int i=0; i<n_reson_0m; i++) {  // All resonance included. Maybe last one should be eliminated.
    func_P2 *= ( z() - z_P2(i) )/( 1 - z()*z_P2(i) );
  }
  if(normailization_index_ == 1 || normailization_index_ == 2 ) { norm = 1.; } // No changes of the Blaschke factor for B0 and B+-.
  else{ norm = P_2_coeff; } // Scale Bs Blascke factor. See documentation ... Value set in the SL_Decay_Config.cpp file.
  func_P2 *= norm;  // func_P2 *= norm;

  return func_P2;
}


double SL_Decay::P_p(){
  double func_Pp = 1.;
  for(int i=0; i<n_reson_1 -1; i++) {  // Last resonance in the Config file is too heavy. 
    func_Pp *= ( z() - z_Pp(i) )/( 1 - z()*z_Pp(i) );
  }
  return func_Pp;
}

double SL_Decay::P_0(){
  double func_P0 = 1.;
  for(int i=0; i<n_reson_0p; i++) {  // Last resonance in the Config file is too heavy.
    func_P0 *= ( z() - z_P0(i) )/( 1 - z()*z_P0(i) );
  }
  return func_P0;
}


double SL_Decay::z_P(int i, int j){
  double func_z_P = ( sqrt( t_p() - pow(m_B_st_c[i][j], 2) ) - sqrt( t_p() - t_m() ) ) / ( sqrt( t_p() - pow(m_B_st_c[i][j], 2) ) + sqrt( t_p() - t_m() ) );
  return func_z_P;
}

double SL_Decay::z_P2(int i) {
  double func_z_P2 = ( sqrt( t_p() - pow(m_B_st_c_0m[i], 2) ) - sqrt( t_p() - t_m() ) ) / ( sqrt( t_p() - pow(m_B_st_c_0m[i], 2) ) + sqrt( t_p() - t_m() ) );
  return func_z_P2;
}

double SL_Decay::z_Pp(int i){
  double z_Pp_func = ( sqrt( t_p() - pow(m_B_st_c[1][i], 2) ) - sqrt( t_p() - t_m() ) ) / ( sqrt( t_p() - pow(m_B_st_c[1][i], 2) ) + sqrt( t_p() - t_m() ) );  // Only 1- resonances.
  return z_Pp_func;
}

double SL_Decay::z_P0(int i){
  double z_P0_func = ( sqrt( t_p() - pow(m_B_st_c_0p[i], 2) ) - sqrt( t_p() - t_m() ) ) / ( sqrt( t_p() - pow(m_B_st_c_0p[i], 2) ) + sqrt( t_p() - t_m() ) );
  return z_P0_func;
}


double SL_Decay::r() { 
  double r =0.;
  if(BtoDstMuNu) { r = D_st_mass_/B_mass_; }
  else
    if(BtoDMuNu) { r = D_mass_/B_mass_; }
  return r;
}

double SL_Decay::a_F1_0_BGL() {
  double x = a_f_BGL_[0]*(B_mass_ - D_st_mass_)/( sqrt(2)*B_mass_*( 1 + r() + 2*sqrt(r()) ) );
  return x;
}

double SL_Decay::t_p() { 
  double t_p = 0.;
  if(BtoDstMuNu) { t_p = pow(B_mass_ + D_st_mass_, 2);}
  else
    if(BtoDMuNu) { t_p = pow(B_mass_ + D_mass_, 2);}
  return t_p;
}

double SL_Decay::t_m() {
  double t_m = 0.;
  if(BtoDstMuNu) { t_m = pow(B_mass_ - D_st_mass_, 2);}
  else
    if(BtoDMuNu) { t_m = pow(B_mass_ - D_mass_, 2);}
  return t_m;
}


// functions useful only for B -> D lnu decays.
double SL_Decay::lambda_func_prime() { 
  double lambda = pow( (Q2_ - pow(B_mass_, 2) - pow(D_mass_, 2)), 2) -4*pow(B_mass_, 2)*pow(D_mass_, 2);
  return lambda;
} 

double SL_Decay::c_l_p() { 
  double clp = lambda_func_prime()/pow(B_mass_, 4)*(1 + pow(lepton_mass_, 2)/2*Q2_);
  return clp;
}


double SL_Decay::c_l_0() {
  double cl0 = pow( 1 - pow(r(), 2), 2)*3*pow(lepton_mass_, 2)/(2*Q2_);
  return cl0;
}

