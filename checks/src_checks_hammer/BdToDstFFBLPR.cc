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
  BdToDstFFBLPRWrapper FF(sh);

  double q2min{0.012};
  double q2max{pow(BdToDst_MASSES[0]-BdToDst_MASSES[1], 2)};
  size_t N{100};
  double q2step = (q2max-q2min)/N;

  // // Fs
  // result.element({0}) = -((Hps * Mc) / sqMbMc);
  // // Ff
  // result.element({1}) = (Ha1 * ((Mb + Mc)*(Mb + Mc) - Sqq)) / (2. * sqMbMc);
  // // Fg
  // result.element({2}) = Hv / (2. * sqMbMc);
  // // Fm
  // result.element({3}) = -(Ha2 * sqrt(Mc / Mb3)) / 2. + Ha3 / (2. * sqMbMc);
  // // Fp
  // result.element({4}) = -(Ha2 * sqrt(Mc / Mb3)) / 2. - Ha3 / (2. * sqMbMc);
  // // Fzt
  // result.element({5}) = (Ht3 * Mc) / (2. * pow(Mb * Mc, 1.5));
  // // Fmt
  // result.element({6}) = (Ht1 * (Mb - Mc)) / (2. * sqMbMc) - (Ht2 * (Mb + Mc)) / (2. * sqMbMc);
  // // Fpt
  // result.element({7}) = (Ht2 * (Mb - Mc)) / (2. * sqMbMc) - (Ht1 * (Mb + Mc)) / (2. * sqMbMc);


  std::cout << "q2 Fs Ff Fg Fm Fp Fzt Fmt Fpt\n";
  for (size_t i{0}; i < N; ++i) {
    double iq2 = q2min + q2step*i;
    FF.performEvalAtQ2(iq2);
    auto tensor = FF.getTensor();
    std::cout << iq2 << " "
              << tensor.element({0}).real() << " "
              << tensor.element({1}).real() << " "
              << tensor.element({2}).real() << " "
              << tensor.element({3}).real() << " "
              << tensor.element({4}).real() << " "
              << tensor.element({5}).real() << " "
              << tensor.element({6}).real() << " "
              << tensor.element({7}).real() << "\n";
  }

  return 0;
}