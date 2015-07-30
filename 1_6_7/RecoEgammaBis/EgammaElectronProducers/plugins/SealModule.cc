#include "PluginManager/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/eventsetupdata_registration_macro.h"

//#include "GlobalCtfElectronProducer.h"
#include "GlobalGsfElectronProducer.h"
#include "SubSeedsCollectionProducer.h"


DEFINE_SEAL_MODULE();
/*DEFINE_ANOTHER_FWK_MODULE(GlobalCtfElectronProducer);*/
DEFINE_ANOTHER_FWK_MODULE(GlobalGsfElectronProducer);
DEFINE_ANOTHER_FWK_MODULE(SubSeedsCollectionProducer);
