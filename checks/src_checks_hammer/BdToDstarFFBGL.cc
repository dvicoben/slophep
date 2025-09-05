#include <iostream>
#include <vector>
#include "Hammer/FormFactors/FFBtoDstarBGL.hh"
#include "Hammer/SettingsHandler.hh"
#include "Hammer/Hammer.hh"
#include "Hammer/Math/Units.hh"
#include <string>

// A wrapper class to help initialise and give access to protected member functions
class FFWrapper: public Hammer::FFBtoDstarBGL {

public:
  FFWrapper(Hammer::SettingsHandler& sh) 
    : Hammer::FFBtoDstarBGL()
  {  
    setSettingsHandler(sh);
    initSettings();
    defineSettings();
  }

  void performEvalAtQ2(double qsq) {
    std::vector<double> masses = {5.27963, 2.0126};
    evalAtPSPoint({qsq}, masses);
  }
};

int main() {
  Hammer::Hammer ham;
  ham.setUnits("GeV");

  Hammer::SettingsHandler sh;
  FFWrapper FF(sh);

  double q2min{0.012};
  double q2max{10.68};
  size_t N{100};
  double q2step = (q2max-q2min)/N;
  double etaVcb = 1.0066*41.5e-3;

  // std::cout << "Fps" << 
  for (size_t i{0}; i < N; ++i) {
    double iq2 = q2min + q2step*i;
    FF.performEvalAtQ2(iq2);
    auto tensor = FF.getTensor();
    std::cout << iq2 << " "
              << tensor.element({0}).real()*etaVcb << " "
              << tensor.element({1}).real()*etaVcb << " "
              << tensor.element({2}).real()*etaVcb << " "
              << tensor.element({3}).real()*etaVcb << " "
              << tensor.element({4}).real()*etaVcb << "\n";
  }

  return 0;
}