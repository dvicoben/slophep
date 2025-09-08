#ifndef CHECKS_HAMMER_FFWRAPPERS_H
#define CHECKS_HAMMER_FFWRAPPERS_H

// A series of wrapper classes to help initialise and give access to protected member functions

#include "Hammer/FormFactors/FFBtoDstarBGL.hh"
#include "Hammer/FormFactors/FFBtoDstarCLN.hh"
#include "Hammer/FormFactors/FFBtoDstarBLPR.hh"
#include "Hammer/FormFactors/FFBtoDBGL.hh"
#include "Hammer/FormFactors/FFBtoDCLN.hh"
#include "Hammer/FormFactors/FFBtoDBLPR.hh"
#include "Hammer/SettingsHandler.hh"
#include <vector>

const std::vector<double> BdToDst_MASSES = {5.27966, 2.01026};
const std::vector<double> BdToD_MASSES = {5.27966, 1.86966};
const std::vector<double> BuToD_MASSES = {5.27934, 1.86484};

// BdToDst BGL
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
// BdToDst CLN
class BdToDstFFCLNWrapper: public Hammer::FFBtoDstarCLN {
public:
  BdToDstFFCLNWrapper(Hammer::SettingsHandler& sh) 
    : Hammer::FFBtoDstarCLN()
  {  
    setSettingsHandler(sh);
    initSettings();
    defineSettings();
  }

  void performEvalAtQ2(double qsq) {
    evalAtPSPoint({qsq}, BdToDst_MASSES);
  }
};
// BdToDst BLPR
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


// BdToD BGL
class BdToDFFBGLWrapper: public Hammer::FFBtoDBGL {
public:
  BdToDFFBGLWrapper(Hammer::SettingsHandler& sh) 
    : Hammer::FFBtoDBGL()
  {  
    setSettingsHandler(sh);
    initSettings();
    defineSettings();
  }

  void performEvalAtQ2(double qsq) {
    evalAtPSPoint({qsq}, BdToD_MASSES);
  }
};
// BdToD CLN
class BdToDFFCLNWrapper: public Hammer::FFBtoDCLN {
public:
  BdToDFFCLNWrapper(Hammer::SettingsHandler& sh) 
    : Hammer::FFBtoDCLN()
  {  
    setSettingsHandler(sh);
    initSettings();
    defineSettings();
  }

  void performEvalAtQ2(double qsq) {
    evalAtPSPoint({qsq}, BdToD_MASSES);
  }
};
// BdToD BLPR
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


// BuToD BGL
class BuToDFFBGLWrapper: public Hammer::FFBtoDBGL {
public:
  BuToDFFBGLWrapper(Hammer::SettingsHandler& sh) 
    : Hammer::FFBtoDBGL()
  {  
    setSettingsHandler(sh);
    initSettings();
    defineSettings();
  }

  void performEvalAtQ2(double qsq) {
    evalAtPSPoint({qsq}, BuToD_MASSES);
  }
};
// BuToD CLN
class BuToDFFCLNWrapper: public Hammer::FFBtoDCLN {
public:
  BuToDFFCLNWrapper(Hammer::SettingsHandler& sh) 
    : Hammer::FFBtoDCLN()
  {  
    setSettingsHandler(sh);
    initSettings();
    defineSettings();
  }

  void performEvalAtQ2(double qsq) {
    evalAtPSPoint({qsq}, BuToD_MASSES);
  }
};
// BuToD BLPR
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