#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "NtupleMakers/ElectronAnalyzer/interface/NewElectrons.h"
#include "NtupleMakers/ElectronAnalyzer/interface/NewElectronsB.h"
#include "NtupleMakers/ElectronAnalyzer/interface/Conversion.h"
#include "NtupleMakers/ElectronAnalyzer/interface/SeedEfficiency.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(NewElectrons);  // for efficiency
DEFINE_ANOTHER_FWK_MODULE(NewElectronsB); // for fake rate
DEFINE_ANOTHER_FWK_MODULE(Conversion);    // to study conversion
DEFINE_ANOTHER_FWK_MODULE(SeedEfficiency);
