#ifndef CHECKS_HAMMER_FFWRAPPERS_H
#define CHECKS_HAMMER_FFWRAPPERS_H

// A series of wrapper classes to help initialise and give access to protected member functions

#include "Hammer/FormFactors/FFBtoDstarBGL.hh"
#include "Hammer/FormFactors/FFBtoDstarBLPR.hh"
#include "Hammer/FormFactors/FFBtoDBLPR.hh"
#include "Hammer/SettingsHandler.hh"
#include <vector>

const std::vector<double> BdToDst_MASSES = {5.27966, 2.01026};
const std::vector<double> BdToD_MASSES = {5.27966, 1.86966};
const std::vector<double> BuToD_MASSES = {5.27934, 1.86484};

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


class BdToDFFBLPRWrapper: public Hammer::FFBtoDBLPR {
public:
  BdToDFFBLPRWrapper(Hammer::SettingsHandler& sh) 
    : Hammer::FFBtoDBLPR()
  {  
    setSettingsHandler(sh);
    initSettings();
    defineSettings();
  }

  void performEvalAtQ2(double qsq) {
    evalAtPSPoint({qsq}, BdToD_MASSES);
  }
};


class BuToDFFBLPRWrapper: public Hammer::FFBtoDBLPR {
public:
  BuToDFFBLPRWrapper(Hammer::SettingsHandler& sh) 
    : Hammer::FFBtoDBLPR()
  {  
    setSettingsHandler(sh);
    initSettings();
    defineSettings();
  }

  void performEvalAtQ2(double qsq) {
    evalAtPSPoint({qsq}, BuToD_MASSES);
  }
};


#endif