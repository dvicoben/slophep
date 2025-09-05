//    **********************************************************    //
//                                                                  //
//                    CONSTRUCTOR PARAMETERS LIST                   //
//                                                                  //
//                                --                                //
//              B0   -----------> 1 |                               // 
//              B+-  -----------> 2 | ---> param 1                  //
//              Bs   -----------> 3 |                               //
//                                --                                //
//                                                                  //
//                                --                                //
//              B -> D lnu  ----> 0 |                               //
//              D*0  -----------> 1 |                               //
//              D*+- -----------> 2 | ---> param 2                  //
//              D*s  -----------> 3 |                               //
//                                --                                //
//                                                                  //
//                                --                                //
//              D0   -----------> 1 |                               //
//              D+-  -----------> 2 | ---> param 3                  //
//              Ds   -----------> 3 |                               //
//                                --                                //
//                                                                  //
//                                --                                //
//              B -> D lnu  ----> 0 |                               //
//              pi0  -----------> 1 |                               //
//              pi+- -----------> 2 | ---> param 4                  //
//              gamma ----------> 3 |                               //
//                                --                                // 
//                                                                  //
//                                --                                //
//              e+-  -----------> 1 |                               //
//              mu+- -----------> 2 | ---> param 5                  //
//              tau+- ----------> 3 |                               //
//                                --                                //
//                                                                  //
//                           -------                                //
//                       "CLN"      |                               //
//                       "CLN_ORI"  |                               //
//                       "CLN_MOD"  |                               //
//                       "CLN_HYB"  |                               //
//                       "BGL"      | ---> param 6                  //
//                       "BCL"      |                               //
//                           -------                                //
//                                                                  //
//    **********************************************************    //


#include "SL_Decay.h"
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

int main() {
  int b{1};
  int dst{2};
  int d0{1};
  int pi{2};
  int lep{2};
  std::string ffpar{"BGL"};
  SL_Decay sl{b, dst, d0, pi, lep, ffpar};

  std::cout << "BGL coefficients" << '\n';
  for (unsigned int i=0; i<4; ++i)
  {
    std::cout << "a_f["  << i << "]= " << sl.a_f_BGL_[i] << '\n';
    std::cout << "a_g["  << i << "]= " << sl.a_g_BGL_[i] << '\n';
    std::cout << "a_F1[" << i << "]= " << sl.a_F1_BGL_[i] << '\n';
    std::cout << "a_F2[" << i << "]= " << sl.a_F2_BGL_[i] << '\n';
  }

  std::cout << "BGL params" << '\n';
  std::cout << "chiT1m " << sl.chi_T_1m_ << '\n';
  std::cout << "chiT1p " << sl.chi_T_1p_ << '\n';
  std::cout << "chiL1p " << sl.chi_L_1p_ << '\n';
  std::cout << "r " << sl.r() << '\n';
  
  double q2min{sl.Get_Q2_min()};
  double q2max{sl.Get_Q2_max()};
  double q2delta{q2max - q2min};
  unsigned int nq2{100};
  double q2step{q2delta/static_cast<double>(nq2)};
  
  std::cout << "q2 z g f F1 F2 Bsck1p Bsck1m, Bsck0 Phig Phif PhiF1 PhiF2 \n";
  for (unsigned int i{0}; i <= nq2; ++i) 
  {
    double iq2{q2min + q2step*i};
    sl.Set_Q2(iq2);
    std::cout << iq2     << " "
              << sl.z() << " "
              << sl.g()  << " "
              << sl.f()  << " "
              << sl.F1() << " "
              << sl.F2() << " "
              << sl.P1_p() << " "
              << sl.P1_m() << " "
              << sl.P2()   << " "
              << sl.phi_g()   << " "
              << sl.phi_f()   << " "
              << sl.phi_F1()  << " "
              << sl.phi_F2()  << " "
              <<'\n' ;
  }

  
  // SL_Decay x1(1, 2, 1, 2, 2, "CLN"); 
  // SL_Decay x2(1, 2, 2, 1, 2, "CLN");
  // SL_Decay x3(1, 2, 2, 3, 2, "CLN");
  // SL_Decay x4(2, 1, 1, 1, 2, "BGL"); 
  // SL_Decay x5(2, 1, 1, 3, 2, "BGL");
  // SL_Decay x6(3, 3, 3, 1, 2, "CLN");
  // SL_Decay x7(3, 3, 3, 3, 2, "BGL");
  // SL_Decay x8(1, 2, 2, 1, 2, "CLN_MOD");
 
  return 0;
}
