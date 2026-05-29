#include <iostream>
#include <vector>
#include <string>
#include "Hammer/FormFactors/BCL/FFBstoKBCL.hh"
#include "Hammer/SettingsHandler.hh"
#include "Hammer/Hammer.hh"
#include "Hammer/Math/Units.hh"


const std::vector<double> BsToK_MASSES = {5.36692, 0.493677};

class BsToKFFBCLWrapper: public Hammer::FFBstoKBCL {
public:
  BsToKFFBCLWrapper(Hammer::SettingsHandler& sh)
    : Hammer::FFBstoKBCL()
  {
    setSettingsHandler(sh);
    initSettings();
    defineSettings();
  }
  void performEvalAtQ2(double qsq) {
    evalAtPSPoint({qsq}, BsToK_MASSES);
  }
};



int main() {
  Hammer::Hammer ham;
  ham.setUnits("GeV");
  

  Hammer::SettingsHandler sh;
  BsToKFFBCLWrapper FF(sh);
  double q2min{0.012};
  double q2max{pow(BsToK_MASSES[0]-BsToK_MASSES[1], 2)};
  size_t N{100};
  double q2step = (q2max-q2min)/N;
  double etaVcb = 1.0066*41.5e-3;
  
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