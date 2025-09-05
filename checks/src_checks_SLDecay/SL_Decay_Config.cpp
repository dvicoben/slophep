#include "SL_Decay_Config.h"
#include <iostream>
#include <cmath>
#include <map>
#include <algorithm>

BR_mapping::BR_mapping(std::string s, double val) {
    dacay_tag = s;
    BR_value = val;
}
BR_mapping::~BR_mapping() {};
std::string BR_mapping::get_BR_tag(void) const { return dacay_tag; }
double BR_mapping::get_BR(void) const { return BR_value; }


// Particles Masses [GeV]
// B Mesons
const double m_B_0 = 5.27963;
const double m_B_pm = 5.27932;
const double m_B_s = 5.36689;
// D* Mesons
const double m_D_st_0 = 2.00685;
const double m_D_st_pm = 2.01026;
const double m_D_st_s_pm = 2.1121;
// D Mesons
const double m_D_0 = 1.86483;
const double m_D_pm = 1.86959;
const double m_D_s_pm = 1.96828;
//Possible D meson partners in the D* decay
const double m_pi_0 = 0.13497;
const double m_pi_pm = 0.13957;
const double m_gamma = 0.;
// Leptons
const double m_e_pm = 0.510999 * pow(10, -3);
const double m_mu_pm = 0.10566;
const double m_tau_pm = 1.77686;


// List of all physical constructor's parameters combination.
// The 5th paramaters, regarding the lepton, has been neglected becuase all their combination are ok.
const int npd = 10;   // Number of possible decay.
const int constr_param_list[npd][4] = {
    {1, 2, 1, 2},
    {1, 2, 2, 1},
    {1, 2, 2, 3},
    {2, 1, 1, 1},
    {2, 1, 1, 3},
    {3, 3, 3, 1},
    {3, 3, 3, 3},
    {1, 0, 2, 0}, // parameters relative to the B->Dlnu decays
    {2, 0, 1, 0}, // parameters relative to the B->Dlnu decays
    {3, 0, 3, 0}  // parameters relative to the B->Dlnu decays
};

/*
// Branching Ratios
std::map<std::string, double> BR_map = {
    {"BR_D_st_0_to_D_0_pi_0", 0.},
    {"BR_D_st_0_to_D_0_gamma", 0.},
    {"BR_D_st_pm_to_D_0_pi_pm", 0.},
    {"BR_D_st_pm_to_D_pm_pi_0", 0.},
    {"BR_D_st_pm_to_D_pm_gamma", 0.},
    {"BR_D_st_s_pm_to_D_s_pm_pi_0", 0.},
    {"BR_D_st_s_pm_to_D_s_pm_gamma", 0.}
};
*/

// Branching Ratios
const BR_mapping BR_map[npd-3] = {   // The BR value of D* mesons is only needed for B -> D*lnu decays. 
    BR_mapping("BR_D_st_0_to_D_0_pi_0", 0.647),
    BR_mapping("BR_D_st_0_to_D_0_gamma", 0.353),
    BR_mapping("BR_D_st_pm_to_D_0_pi_pm", 0.677),
    BR_mapping("BR_D_st_pm_to_D_pm_pi_0", 0.307),
    BR_mapping("BR_D_st_pm_to_D_pm_gamma", 0.016),
    BR_mapping("BR_D_st_s_pm_to_D_s_pm_pi_0", 0.058),
    BR_mapping("BR_D_st_s_pm_to_D_s_pm_gamma", 0.935),
};

 
// Miscellaneous
const double G_Fermi = 1.166378 * pow(10, -5);  // [(GeV)^-2]
const double Vcb = 41.1 * pow(10, -3);          // PDG value (2018).
// const double Vcb = 37.4 * pow(10, -3);       // Theoretical paper value.  

//  ************************************  //
//                                        //
//            CLN PARAMETERS              //
//                                        //
//  ***********************************   //

// Theoric Paper with a set of 'default' CLN parameters: http://arxiv-export-lb.library.cornell.edu/abs/1801.10468
// The following are the CLN parameters used by LHCb
const double rho_squared = 1.205;
const double R1 = 1.404;
const double R2 = 0.854;
const double h_A1 = 0.921;
const double R0 = 1.;

// Parameter only for the B -> D lnu decays
const double h_A1_BtoD = 1.;
const double rho_squared_BtoD = 1.;


//  ************************************  //
//                                        //
//           CLN_MOD PARAMETERS           //
//                                        //
//  ***********************************   //

// Coeffcients of z^2 and z^3 in the h_A1(w) form factor of the CLN parametrization.
const double c_z2 = 53*rho_squared - 15;  // z^2 original coeff. Set to the default CLN value. Modifialble through setter method.
const double c_z3 = 231*rho_squared - 91; // z^3 original coeff. Set to the default CLN value. Modifialble through setter method.


//  ************************************  //
//                                        //
//     BGL Parameters and Constants       //
//                                        //
//  ************************************  //

const double chi_T_1p = 3.894 * pow(10, -4);   // [(GeV)^-2]
const double chi_T_1m = 5.131 * pow(10, -4);   // [(GeV)^-2]
const double chi_L_1p = 1.9421 * pow(10, -2);  // Linked to chi_0-. Pure number.
const double n_I = 2.6; 

// B_st_c resonances masses with quantum number +-1 
const int n_reson_1 = 4;
const double m_B_st_c[2][n_reson_1] = {  // [+1, -1][excited states] 
  { 6.739, 6.750, 7.145, 7.150 },    // +1 excited states [GeV] 
  { 6.329, 6.920, 7.020, 7.280 }     // -1 excited states [GeV]
};

// B_st_c resonances masses with quantum number 0-
const int n_reson_0m = 3;
const double m_B_st_c_0m[n_reson_0m] = { 6.275, 6.842, 7.250 };  // [GeV]
const double dec_const = 0.427;  // [GeV]

// Series parameters
const double a_f_BGL[4] = { 0.01224, -0.052, 1.0, 0. };
const double a_F1_BGL[4] = { 0., -0.0070, 0.089, 0. };   // a_F1[0] set to 0 for simplicity reasons.
                                                         // The correct value is set in the constructor using the member function a_0_F1.
const double a_g_BGL[4] = { 0.0289, 0.08, -1.0, 0. };
const double a_F2_BGL[4] = { 0.0595, -0.318, 0., 0. };     // a_F2[2] = 0.

// Blascke Factors normalization coefficients for Bs SL decays.
const double P_1p_coeff = 2.02159;  // Obtained including the last 1^- resonance and calculating the avarage of 10 scaling coefficients along the Q2 range.
const double P_1m_coeff = 2.52733;
const double P_2_coeff = 1.73835;
// const double P_1p_coeff = 2.02159; // Obtained calculating the avarage of 10 scaling coefficients along the Q2 range. 
// const double P_1m_coeff = 1.54435;
// const double P_2_coeff = 1.73835;

//  ***************************************************  //
//     BGL parameters ONLY for the B -> D lnu decays     //
//  ***************************************************  //

// Series parameters
const double a_BGL[4] = {0.01565, -0.0353, -0.043, 0.194};
const double b_BGL[4] = {0.07932, -0.214, 0.17, -0.958};

const double chi_L_0 = 6.204 * pow(10, -3); 
const double chi_Tt_0 = 5.131 * pow(10, -4);

// B_st_c resonances masses with quantum number 0+
const int n_reson_0p = 2;
const double m_B_st_c_0p[n_reson_0p] = { 6.716, 7.121};  // [GeV]


//  ************************************  //
//                                        //
//     BCL Parameters and Constants       //
//                                        //
//  ************************************  //

const double N_f = 300.;
const double N_F1 = 7000.;
const double N_g = 5.;
const double N_F2 = 10.; // Unknown parameter. Set randomly.


// See paper https://arxiv.org/abs/1711.11013, table XII. 
const double a_f_BCL[4] = { 0.01502, -0.311, -0.02, -0.007 };
const double a_F1_BCL[4] = { 0.002946, -0.0152, 0.21, 0.03 };
const double a_g_BCL[4] = { 0.109, -0.29, -0.02, -0.00009 };
const double a_F2_BCL[4] = { 0.01678, -0.1246, 0.17, 0.01 }; // Unknown parameters. Set randomly, but with the same magnitude order of the others. 

//  ***************************************************  //
//     BCL parameters ONLY for the B -> D lnu decays     //
//  ***************************************************  //

// Series parameters
const double a_BCL[4] = {0.421, -0.469, -0.178, -0.825};
const double b_BCL[4] = {0.0, 0.0, 0.0, 0.0};                // Set to zero for now. Need to find a paper from where I can take these parameters.


//  **************************************  //
//                                          //
//     Parameters for B -> D lnu decays     //
//                                          //
//  **************************************  //

const double eta_EW = 1.0066;          // electroweak correction factor. Upadate if a more accurate estimation exists.

