#ifndef CHECKS_HAMMER_FFWRAPPERS_H
#define CHECKS_HAMMER_FFWRAPPERS_H

// A series of wrapper classes to help initialise and give access to protected member functions

#include "Hammer/FormFactors/FFBtoDstarBGL.hh"
#include "Hammer/FormFactors/FFBtoDstarBLPR.hh"
#include "Hammer/SettingsHandler.hh"
#include <vector>

const std::vector<double> BdToDst_MASSES = {5.27966, 2.01026};

class BdToDstFFBGLWrapper: public Hammer::FFBtoDstarBGL {
public:
  BdToDstFFBGLWrapper(Hammer::SettingsHandler& sh) 
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


class BdToDstFFBLPRWrapper: public Hammer::FFBtoDstarBLPR {
public:
  BdToDstFFBLPRWrapper(Hammer::SettingsHandler& sh) 
    : Hammer::FFBtoDstarBLPR()
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