#include <iostream>
#include <vector>
#include "Hammer/FormFactors/FFBtoDstarBGL.hh"
#include "Hammer/SettingsHandler.hh"
#include "Hammer/Hammer.hh"
#include "Hammer/Math/Units.hh"
#include <string>

#include "HammerFFWrappers.h"

int main() {
  Hammer::Hammer ham;
  // std::vector<std::string> include = {"BD*MuNu", "D*DPi"};
  // ham.includeDecay(include);
  ham.setUnits("GeV");

  Hammer::SettingsHandler sh;
  BdToDstFFWrapper FF(sh);

  double q2min{0.012};
  double q2max{10.68};
  size_t N{100};
  double q2step = (q2max-q2min)/N;
  double etaVcb = 1.0066*41.5e-3;

  std::cout << "q2 Ff Fg Fm Fp\n";
  for (size_t i{0}; i < N; ++i) {
    double iq2 = q2min + q2step*i;
    FF.performEvalAtQ2(iq2);
    auto tensor = FF.getTensor();
    std::cout << iq2 << " "
              << tensor.element({1}).real() << " "
              << tensor.element({2}).real() << " "
              << tensor.element({3}).real() << " "
              << tensor.element({4}).real() << "\n";
  }

  return 0;
}