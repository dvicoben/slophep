#include <iostream>
#include <vector>
#include <string>
#include "Hammer/FormFactors/BCL/FFBtoPiBCL.hh"
#include "Hammer/SettingsHandler.hh"
#include "Hammer/Hammer.hh"
#include "Hammer/Math/Units.hh"


const std::vector<double> BToPi_MASSES = {5.27966, 0.13957039000000002};

class BToPiFFBCLWrapper: public Hammer::FFBtoPiBCL {
public:
  BToPiFFBCLWrapper(Hammer::SettingsHandler& sh)
    : Hammer::FFBtoPiBCL()
  {
    setSettingsHandler(sh);
    initSettings();
    defineSettings();
  }
  void performEvalAtQ2(double qsq) {
    evalAtPSPoint({qsq}, BToPi_MASSES);
  }
};



int main() {
  Hammer::Hammer ham;
  ham.setUnits("GeV");
  

  Hammer::SettingsHandler sh;
  BToPiFFBCLWrapper FF(sh);
  double q2min{0.012};
  double q2max{pow(BToPi_MASSES[0]-BToPi_MASSES[1], 2)};
  size_t N{100};
  double q2step = (q2max-q2min)/N;
  
  std::cout << "q2 F0 Fp\n";
  for (size_t i{0}; i < N; ++i) {
    double iq2 = q2min + q2step*i;
    FF.performEvalAtQ2(iq2);
    auto tensor = FF.getTensor();
    std::cout << iq2 << " "
              << tensor.element({1}).real() << " "
              << tensor.element({2}).real() << "\n";
  }
}