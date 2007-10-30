#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "NtupleMakers/ElectronAnalyzer/interface/NewElectrons.h"
#include "NtupleMakers/ElectronAnalyzer/interface/NewElectronsB.h"
#include "NtupleMakers/ElectronAnalyzer/interface/SeedEfficiency.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(NewElectrons);
DEFINE_ANOTHER_FWK_MODULE(NewElectronsB);
DEFINE_ANOTHER_FWK_MODULE(SeedEfficiency);
