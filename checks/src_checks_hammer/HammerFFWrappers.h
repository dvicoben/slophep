#ifndef CHECKS_HAMMER_FFWRAPPERS_H
#define CHECKS_HAMMER_FFWRAPPERS_H

// A series of wrapper classes to help initialise and give access to protected member functions

#include "Hammer/FormFactors/FFBtoDstarBGL.hh"
#include "Hammer/SettingsHandler.hh"
#include <vector>

const std::vector<double> BdToDst_MASSES = {5.27966, 2.01026};

class BdToDstFFWrapper: public Hammer::FFBtoDstarBGL {
public:
  BdToDstFFWrapper(Hammer::SettingsHandler& sh) 
    : Hammer::FFBtoDstarBGL()
  {  
    setSettingsHandler(sh);
    initSettings();
    defineSettings();
  }

  void performEvalAtQ2(double qsq) {
    evalAtPSPoint({qsq}, BdToDst_MASSES);
  }
};



#endif