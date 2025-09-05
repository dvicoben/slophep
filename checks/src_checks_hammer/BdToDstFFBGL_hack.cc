#include <iostream>
#include <vector>
#include "Hammer/FormFactors/FFBtoDstarBGL.hh"
#include "Hammer/SettingsHandler.hh"
#include "Hammer/Hammer.hh"
#include "Hammer/Math/Units.hh"
#include <string>

#include "HammerFFWrappers.h"
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
    std::vector<double> masses = {5.27966, 2.01026};
    evalAtPSPoint({qsq}, masses);
  }

  void evalAtPSPoint(const std::vector<double>& point, const std::vector<double>& masses) {
    auto& result = getTensor();
    result.clearData();

    double Mb = 0.;
    double Mc = 0.;
    double unitres = 1.;
    if(masses.size() >= 2) {
        Mb = masses[0];
        Mc = masses[1];
        unitres = _units;
    }
    else {
        Mb = this->masses()[0];
        Mc = this->masses()[1];
    }
    const double Sqq = point[0];
    const double Mb2 = Mb*Mb;
    const double Mb3 = Mb2*Mb;
    const double Mc2 = Mc*Mc;
    const double rC = Mc/Mb;
    const double rC2 = rC*rC;
    const double sqrC = sqrt(rC);
    const double sqrt2 = sqrt(2.);
    double pi = M_PI;

    double w = (Mb2 + Mc2 - Sqq) / (2. * Mb * Mc);
    //safety measure if w==1.0
    //if(isZero(w - 1.0)) w += 1e-6;
    const double w2 = w*w;

    const size_t nmax = 4;
    const double z = (sqrt(w+1) - sqrt2)/(sqrt(w+1) + sqrt2);
    std::vector<double> zpow{1.};
    for (size_t n = 1; n < nmax; ++n){
        zpow.push_back(zpow[n-1]*z);
    }

    const std::vector<double>& ag=(*getSetting<std::vector<double>>("avec"));
    const std::vector<double>& af=(*getSetting<std::vector<double>>("bvec"));
    const std::vector<double>& aF1=(*getSetting<std::vector<double>>("cvec"));
    const std::vector<double>& aP1=(*getSetting<std::vector<double>>("dvec"));

    const double nc=2.6;
    const double etaEWVcb = 1.0066*(*getSetting<double>("Vcb"));
    const double chim = (*getSetting<double>("Chim"))/(unitres*unitres);
    const double chip = (*getSetting<double>("Chip"))/(unitres*unitres);
    const double chimL=(*getSetting<double>("ChimL"));

    const double tp=(Mb+Mc)*(Mb+Mc)/Mb2;
    const double tm=(Mb-Mc)*(Mb-Mc)/Mb2;
    const double sqtptm = sqrt(tp - tm);

    std::vector<double>& BcStatesf = (*getSetting<std::vector<double>>("BcStatesf"));
    std::vector<double>& BcStatesg = (*getSetting<std::vector<double>>("BcStatesg"));
    std::vector<double>& BcStatesP1 = (*getSetting<std::vector<double>>("BcStatesP1"));

    double Pf = 1.;
    for(size_t n = 0; n< BcStatesf.size(); ++n){
        double sqtpBc = sqrt(tp-pow(BcStatesf[n]*unitres/Mb,2));
        Pf *= ((z-((sqtpBc-sqtptm)/(sqtpBc+sqtptm)))/(1.-z*((sqtpBc-sqtptm)/(sqtpBc+sqtptm))));
    }

    double PF1 = Pf;

    double Pg = 1.;
    for(size_t n = 0; n< BcStatesg.size(); ++n){
        double sqtpBc = sqrt(tp-pow(BcStatesg[n]*unitres/Mb,2));
        Pg *= ((z-((sqtpBc-sqtptm)/(sqtpBc+sqtptm)))/(1.-z*((sqtpBc-sqtptm)/(sqtpBc+sqtptm))));
    }

    double PP1 = 1.;
    for(size_t n = 0; n< BcStatesP1.size(); ++n){
        double sqtpBc = sqrt(tp-pow(BcStatesP1[n]*unitres/Mb,2));
        PP1 *= ((z-((sqtpBc-sqtptm)/(sqtpBc+sqtptm)))/(1.-z*((sqtpBc-sqtptm)/(sqtpBc+sqtptm))));
    }

    const double phig=sqrt(256.*nc/(3*pi*chip))*((rC2*pow(1+z,2)*pow(1-z,-0.5))/pow(((1+rC)*(1-z)+2*sqrC*(1+z)),4));
    const double phif=(1./Mb2)*sqrt(16.*nc/(3*pi*chim))*((rC*(1+z)*pow(1-z,1.5))/pow(((1+rC)*(1-z)+2*sqrC*(1+z)),4));
    const double phiF1=(1./Mb3)*sqrt(8.*nc/(3*pi*chim))*((rC*(1+z)*pow(1-z,2.5))/pow(((1+rC)*(1-z)+2*sqrC*(1+z)),5));

    const double phif_0 = 4.*rC*sqrt(nc/chim)/(Mb2*sqrt(3*pi)*pow(1+2*sqrC+rC,4));
    const double phiF1_0 = 2.*sqrt(2/(3*pi))*rC*sqrt(nc/chim)/(Mb3*pow(1+2*sqrC+rC,5));

    const double phiP1=sqrt(nc/(pi*chimL))*((8.*sqrt2*rC2*pow(1+z,2)*pow(1-z,-0.5)))/pow(((1+rC)*(1-z)+2*sqrC*(1+z)),4);

    double g=0;
    for(size_t n = 0; n< ag.size(); ++n) g += ag[n] * zpow[n];
    g /= (Pg*phig);

    double f=0;
    for(size_t n = 0; n< af.size(); ++n) f += af[n] * zpow[n];
    f /= (Pf*phif);

    double F1=af[0]*(Mb-Mc)*phiF1_0/phif_0;
    for(size_t n = 0; n< aF1.size(); ++n) F1 += aF1[n] * zpow[n+1];
    F1 /= (PF1*phiF1);

    //From 1707.09509 (sqrC/(1+rC) maps from F2)
    double P1=0;
    for(size_t n = 0; n< aP1.size(); ++n) P1 += aP1[n]*zpow[n];
    // P1 *= sqrC/((1+rC)*PP1*phiP1);
    P1 *= 1/(PP1*phiP1);

    //Mapping to amplitude FF basis
    const double Fpf = (rC-w)/(2.*rC*Mb2*(w2 - 1));
    const double FpF1 = 1./(2.*rC*Mb3*(w2 - 1));
    const double Fmf = (rC+w)/(2.*rC*Mb2*(w2 - 1));
    const double FmF1 = 1./(2.*rC*Mb3*(w2 - 1))*(rC2-1)/(1 + rC2 - 2*rC*w);
    const double FmP1 = sqrC*(rC+1)/(Mb*(1 + rC2 - 2*rC*w));

    // set elements, mapping to amplitude FF basis
    // Fs (dim 0)
    // result.element({0}) = 0;
    // Ff (dim 1)
    result.element({1}) = f;
    // Fg (dim -1) Definition in BGL (hep-ph/9705252) differs from Manohar, Wise or 1610.02045 by factor of 2
    result.element({2}) = g;
    // Fm (dim -1)
    result.element({3}) = F1;
    // Fp (dim -1)
    result.element({4}) = P1;
    // Fzt (dim -2)
    // result.element({5}) = 0;
    // Fmt (dim 0)
    // result.element({6}) = 0;
    // Fpt (dim 0)
    // result.element({7}) = 0;

    std::cout << Sqq << " "
              << z << " "
              << f/etaEWVcb << " "
              << g/etaEWVcb << " "
              << F1/etaEWVcb << " "
              << P1/etaEWVcb << " "
              << Pf << " "
              << Pg << " "
              << PP1 << " "
              << phif << " "
              << phig << " "
              << phiF1 << " "
              << phiP1 << " "
              << chim << " "
              << chip << " "
              << chimL << " "
              << rC << "\n";

    result*=(1./etaEWVcb);
  }
};


int main() {
  Hammer::Hammer ham;
  // std::vector<std::string> include = {"BD*MuNu", "D*DPi"};
  // ham.includeDecay(include);
  ham.setUnits("GeV");

  Hammer::SettingsHandler sh;
  FFWrapper FF(sh);

  double q2min{0.012};
  double q2max{10.68};
  size_t N{100};
  double q2step = (q2max-q2min)/N;
  double etaVcb = 1.0066*41.5e-3;

  std::cout << "q2 z f g F1 F2 Blf Blg BlP1 phif phig phiF1 phiP1\n";
  for (size_t i{0}; i < N; ++i) {
    double iq2 = q2min + q2step*i;
    FF.performEvalAtQ2(iq2);
    // auto tensor = FF.getTensor();
    // std::cout << iq2 << " "
    //           << tensor.element({1}).real() << " "
    //           << tensor.element({2}).real() << " "
    //           << tensor.element({3}).real() << " "
    //           << tensor.element({4}).real() << "\n";
  }

  return 0;
}