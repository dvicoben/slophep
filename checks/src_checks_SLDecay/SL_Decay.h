#ifndef SL_Decay_h
#define SL_Decay_h

#include<vector>
#include<string>

class SL_Decay
{
  //   ********************************************   //
  //                    DECAY SCHEME                  //
  //                                                  //
  //                 M --> A +  B                     //
  //                       |    |                     // 
  //                       |    -> B1 + B2            //
  //                       |                          // 
  //                       ->  A1 + A2                //
  //                                                  //
  //   ********************************************   //

public:
// private : 
  
  // B is always a W decaying leptonically. B1 is the charged lepton and B2 the associated neutrino.
  double B_mass_;         // B meson.
  double D_st_mass_;      // D* meson.
  double D_mass_;         // D meson.
  double pi_gamma_mass_;  // pi or gamma.
  double lepton_mass_;    // charged lepton.
  double BR_;
  int normailization_index_; // Index introduced to scale the Bs Blaschke factors.
  int lepton_index_;         // 1=el, 2=mu, 3=tau. Used in CLN_ORG parametrization.
  bool is_pi;                // boolean variable to distinguish cases: D*->D+pi, D*->D+gamma
  bool CLN_;                 // CLN parametrization with explict lepton mass dependence (usable also with semitauonic decays).
  bool BGL_;                 // BGL parametrization with explict lepton mass dependence (usable also with semitauonic decays).
  bool BCL_;                 // BCL parametrization. Not fully debugged.
  bool CLN_MOD_;             // CLN parametrization with explict lepton mass dependence + two more parameters in hA1(usable also with semitauonic decays).
  bool CLN_ORI_;             // CLN parametrization accordig to the original paper. Lepton mass set to zero everywhere except for the phase space.
  bool CLN_HYB_;             // CLN parametrization with explict lepton mass dependence with hA1 form factor expressed in function of the corresponding BGL form factor (HYBRID).
  bool BtoDstMuNu;           // Variable to distinguish B -> D* l nu from B -> D l nu decays.
  bool BtoDMuNu;             // Variable to distinguish B -> D* l nu from B -> D l nu decays.
 
  // Differential amplitude parameters.
  double Q2_;
  double cos_theta_l_;
  double cos_theta_V_;
  double k_;

  // Different parametrizations use different value of Vcb. 
  double Vcb_;

  // CLN Parametrization parameters.
  double rho_squared_;  
  double R1_;
  double R2_;
  double h_A1_;
  double R0_;

  // CLN_MOD Parametrization parameters.
  // All of the CLN parameters plus two.
  double c_z2_;
  double c_z3_;

  // CLN Parameters for the B -> D lnu   

  // BGL Parameterization parameters.
  double a_f_BGL_[4];
  double a_F1_BGL_[4];
  double a_g_BGL_[4];
  double a_F2_BGL_[4];
  double chi_T_1p_;
  double chi_T_1m_;
  double chi_L_1p_;
  // BGL Parameters for the B -> D lnu. 
  double a_BGL_[4];
  double b_BGL_[4];
  double chi_L_0_;
  double chi_Tt_0_;


  // BCL Parameterization parameters.
  double a_f_BCL_[4];
  double a_F1_BCL_[4];
  double a_g_BCL_[4];
  double a_F2_BCL_[4];
  double N_f_;
  double N_F1_;
  double N_g_;
  double N_F2_;
  // BCL Parameters for the B -> D lnu.
  double a_BCL_[4];
  double b_BCL_[4];

// public:

  // Member function
  double Nf();                                         // Compute the constant factor Nf.
  double lambda_func(double x, double y, double z);    // Triangular function.
  double D_st_P_module();
  double r();
  // Used only for B -> D lnu decays.
  double lambda_func_prime();
  double c_l_0();
  double c_l_p();
  

  //  *****************************  //
  //            I functions          //
  //  *****************************  //
  // I1s, I1c, I2s, I2c, I3, I4, I5, I6s, I6c, I7 functions declarations.
  double I_1s();
  double I_1c();
  double I_2s();
  double I_2c();
  double I_3();
  double I_4();
  double I_5();
  double I_6s();
  double I_6c();
  double I_7();

  //  *****************************  //
  //       Helicity Amplitudes       //
  //  *****************************  //
  // H0, H+, H- and H_t functions.
  double f_H_0();
  double f_H_p();
  double f_H_m();
  double f_H_t();

  //  ****************************************************  //
  //      Form factors used only for B -> D lnu decays      //
  //  ****************************************************  //
  double f_p();
  double f_0();

  // A1, A2, A0 and V functions.
  double R_star();                                                                                                                                     
  double f_A1();
  double f_A2();
  double f_A0();
  double f_V();

  //  *****************************  //
  //       CLN member functions      //
  //  *****************************  //
  double f_h_A1();
  double f_R1();
  double f_R2();
  double f_R0();

  // CLN_MOD member function
  double f_h_A1_mod();
  


  //  *****************************  //
  //       BGL member functions      //
  //  *****************************  //
  double a_F1_0_BGL();

  double t_p();
  double t_m();
  double z_P(int i, int j);    // Blaschke factor.
  double z_P2(int i);          // Blaschke factor.
  
  double P1_p();
  double P1_m();
  double P2();

  double phi_g();              // Outer function.
  double phi_f();              // Outer function.
  double phi_F1();             // Outer function.
  double phi_F2();             // Outer function.

  double f();
  double F1();
  double g();
  double F2();

  // BGL functions ONLY for B -> D lnu
  double P_p();                // Blaschke factor.
  double P_0();                // Blaschke factor.
  double phi_p();              // Outer function.
  double phi_0();              // Outer function.
  double z_Pp(int i);
  double z_P0(int i);


  //  *****************************  //
  //       BCL member functions      //
  //  *****************************  //
  double Q_f();
  double Q_F1();
  double Q_g();
  double Q_F2();

  // BCL functions ONLY for B -> D lnu 



  // public:
    
  // Basic functions
  SL_Decay(int p1, int p2, int p3, int p4, int p5, std::string parametrization);   // Standard constructor.
  SL_Decay(double m1, double m2, double m3, double m4, double m5, double BR, std::string pi_or_gamma, std::string parametrization);  // Second constructor.
  ~SL_Decay();

  // Gamma function definition.
  double Gamma_Diff_4d(double Q2, double cos_theta_l, double cos_theta_V, double k);
  // Gamma function integrated in cos_theta_l, cos_theta_V and k.
  double Gamma_Diff_1d(double Q2);
  
  // Setter and Getter functions for the class object.
  void Set_Q2 (double Q2);
  void Set_cos_theta_l (double cos_theta_l);
  void Set_cos_theta_V (double cos_theta_V);
  void Set_k (double k);
  double Get_Q2();
  double Get_cos_theta_l();
  double Get_cos_theta_V();
  double Get_k();
  double Get_Q2_max();
  double Get_Q2_min();
  std::string Get_parameterization(); 
  double Get_BR();

  // Vcb is a parameter common to all the parametrizations.
  void Set_Vcb(double x);
  // Useful Q^2 reparametrizations.
  double w();  // function of Q^2.
  double z();  // function of w -> z(w(Q^2)).

  // CLN Setter Functions. 
  void Set_CLN(double rho_squared, double R1, double R2, double h_A1);
  void Set_h_A1 (double h_A1);
  void Set_R1 (double R1);
  void Set_R2 (double R2);
  void Set_R0 (double R0);
  void Set_rho_squared (double rho_squared);

  // CLN_MOD setter and getter methods.
  // All CLN setter methods except for Set_CLN() plus three more.
  void Set_CLN_MOD(double rho_squared, double R1, double R2, double h_A1, double coeff_z2, double coeff_z3);
  void Set_cz2(double coeff_z2);
  void Set_cz3(double coeff_z3);
  double Get_cz2();
  double Get_cz3();

  //                         BGL Setter Functions
  //    *************************************************************   //
  //                        !!!! WARNING !!!!                           //
  //  BGL parameters initialization is critical due to the dependency   //
  //  of a_F1[0] from a_f[0] and vice versa. In case of multiple        //
  //  initializations, keep in mind this dependency.  //
  //  More information about that can be found on Sec. 6 of the         //
  //  Internal Note at: https://cds.cern.ch/record/2313977              //
  //                        !!!! WARNING !!!!                           //
  //    *************************************************************   //
  void Set_BGL_f(double a_f_0, double a_f_1, double a_f_2, double a_f_3);
  void Set_BGL_g(double a_g_0, double a_g_1, double a_g_2, double a_g_3);
  void Set_BGL_F1(double a_F1_0, double a_F1_1, double a_F1_2, double a_F1_3);
  void Set_BGL_F2(double a_F2_0, double a_F2_1, double a_F2_2, double a_F2_3);
  // void Set_BGL_chi_T_1p(double x);
  // void Set_BGL_chi_T_1m(double x);
  // void Set_BGL_chi_t_0m(double x);

  // BCL Setter Function
  void Set_BCL_f(double a_f_0, double a_f_1, double a_f_2, double a_f_3);
  void Set_BCL_g(double a_g_0, double a_g_1, double a_g_2, double a_g_3);
  void Set_BCL_F1(double a_F1_0, double a_F1_1, double a_F1_2, double a_F1_3);
  void Set_BCL_F2(double a_F2_0, double a_F2_1, double a_F2_2, double a_F2_3);

};


#endif // ifndef SL_Decay_h
