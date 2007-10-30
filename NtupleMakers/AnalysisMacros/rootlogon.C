
//------------------------------------------------------------------------------
//  rootlogon.C: a sample ROOT logon macro allowing use of ROOT script
//               compiler in CDF RunII environment. The name of this macro file
//               is defined by the .rootrc file
//
//  USESHLIBS variable has to be set to build MyAna libraries locally:
//
//  setenv USESHLIBS 1
//
//  July 7 2003 A.Taffard
//------------------------------------------------------------------------------

{

#include <iomanip.h>
#include <time.h>
#include <TObjectTable.h>

  const char* exec_name = gApplication->Argv(0);
  // the line below tells ROOT script compiler
  // where to look for the include files

  // load in ROOT physics vectors and event
  // generator libraries

  gSystem->Load("$ROOTSYS/lib/libPhysics.so");
  
  gSystem->Load("$ROOTSYS/lib/libEG.so");
 
  gSystem->Load("$ROOTSYS/lib/libMatrix.so");
  gSystem->Load("$ROOTSYS/lib/libHist.so");
  gSystem->Load("$ROOTSYS/lib/libMathCore.so");
  
  gSystem->Load("$ROOTSYS/lib/libCore.so");

// load Gremlins ana Utils
//  gSystem->CompileMacro("~yagil/ana/anaUtils.cc","k");

  // This line reports the process ID which simplifies
  // debugging
//  gInterpreter->ProcessLine(".! ps | grep root");

  //ADD HERE any of your own libs to load automatically when starting root.


//  gInterpreter->ProcessLine(".x /Users/yagil/ana/settings.C");
//gInterpreter->ProcessLine(".x /Users/yagil/ana/RootStyle.C");
  gROOT->ProcessLine(".L ./RootStyle.cc");
    gROOT->ProcessLine("Style = RootStyle()");

}

