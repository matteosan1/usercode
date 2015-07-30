#include "FWCore/PluginManager/interface/ModuleDef.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "HLTrigger/Egamma/interface/HLTEgammaEtFilter.h"
#include "HLTrigger/Egamma/interface/HLTEgammaDoubleEtFilter.h"
#include "HLTrigger/Egamma/interface/HLTEgammaEcalIsolFilter.h"
#include "HLTrigger/Egamma/interface/HLTEgammaHcalIsolFilter.h"
#include "HLTrigger/Egamma/interface/HLTEgammaHcalDBCFilter.h"
#include "HLTrigger/Egamma/interface/HLTEgammaHOEFilter.h"
#include "HLTrigger/Egamma/interface/HLTPhotonTrackIsolFilter.h"
#include "HLTrigger/Egamma/interface/HLTElectronPixelMatchFilter.h"
#include "HLTrigger/Egamma/interface/HLTPMMassFilter.h"
#include "HLTrigger/Egamma/interface/HLTPMDocaFilter.h"


#include "HLTrigger/Egamma/interface/HLTEgammaL1MatchFilterRegional.h"
#include "HLTrigger/Egamma/interface/HLTElectronEoverpFilterRegional.h"
#include "HLTrigger/Egamma/interface/HLTElectronTrackIsolFilterRegional.h"

#include "HLTrigger/Egamma/interface/HLTEgammaDoubleEtPhiFilter.h"
#include "HLTrigger/Egamma/interface/HLTElectronOneOEMinusOneOPFilterRegional.h"
#include "HLTrigger/Egamma/interface/HLTElectronRegionalSeedFilter.h"

DEFINE_SEAL_MODULE();
DEFINE_ANOTHER_FWK_MODULE(HLTEgammaEtFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTEgammaDoubleEtFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTEgammaEcalIsolFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTEgammaHcalIsolFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTEgammaHcalDBCFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTEgammaHOEFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTPhotonTrackIsolFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTElectronPixelMatchFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTPMMassFilter);
DEFINE_ANOTHER_FWK_MODULE(HLTPMDocaFilter);


DEFINE_ANOTHER_FWK_MODULE(HLTEgammaL1MatchFilterRegional);
DEFINE_ANOTHER_FWK_MODULE(HLTElectronEoverpFilterRegional);
DEFINE_ANOTHER_FWK_MODULE(HLTElectronTrackIsolFilterRegional);

DEFINE_ANOTHER_FWK_MODULE(HLTEgammaDoubleEtPhiFilter);

DEFINE_ANOTHER_FWK_MODULE(HLTElectronOneOEMinusOneOPFilterRegional);

DEFINE_ANOTHER_FWK_MODULE(HLTElectronRegionalSeedFilter);


