
//#include <FWCore/Framework/interface/Handle.h>
#include <FWCore/Framework/interface/Event.h>
#include <FWCore/Framework/interface/EventSetup.h>
#include <FWCore/Framework/interface/ESHandle.h>

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoEgamma/EgammaHLTProducers/interface/EcalRecHitsMerger.h"

#include "FWCore/Utilities/interface/Exception.h"


using namespace edm;
using namespace std;


EcalRecHitsMerger::EcalRecHitsMerger(const edm::ParameterSet& pset) {

 debug_ = pset.getUntrackedParameter<bool>("debug");

 EgammaSourceEB_ = pset.getUntrackedParameter<edm::InputTag>("EgammaSource_EB");
 MuonsSourceEB_  = pset.getUntrackedParameter<edm::InputTag>("MuonsSource_EB");
 TausSourceEB_  = pset.getUntrackedParameter<edm::InputTag>("TausSource_EB");
 JetsSourceEB_   = pset.getUntrackedParameter<edm::InputTag>("JetsSource_EB");
 RestSourceEB_   = pset.getUntrackedParameter<edm::InputTag>("RestSource_EB");

 EgammaSourceEE_ = pset.getUntrackedParameter<edm::InputTag>("EgammaSource_EE");
 MuonsSourceEE_  = pset.getUntrackedParameter<edm::InputTag>("MuonsSource_EE");
 TausSourceEE_  = pset.getUntrackedParameter<edm::InputTag>("TausSource_EE");
 JetsSourceEE_   = pset.getUntrackedParameter<edm::InputTag>("JetsSource_EE");
 RestSourceEE_   = pset.getUntrackedParameter<edm::InputTag>("RestSource_EE");

 OutputLabelEB_ = pset.getUntrackedParameter<std::string>("OutputLabel_EB");
 OutputLabelEE_ = pset.getUntrackedParameter<std::string>("OutputLabel_EE");

 produces<EcalRecHitCollection>(OutputLabelEB_);
 produces<EcalRecHitCollection>(OutputLabelEE_);

}



EcalRecHitsMerger::~EcalRecHitsMerger() {
}


void EcalRecHitsMerger::beginJob(const edm::EventSetup& c){
}

void EcalRecHitsMerger::endJob(){
}

void EcalRecHitsMerger::produce(edm::Event & e, const edm::EventSetup& iSetup){

 if (debug_) cout << " EcalRecHitMerger : Run " << e.id().run() << " Event " << e.id().event() << endl;

 std::vector< edm::Handle<EcalRecHitCollection> > EcalRecHits_done;
 e.getManyByType(EcalRecHits_done);

 std::auto_ptr<EcalRecHitCollection> EBMergedRecHits(new EcalRecHitCollection);
 std::auto_ptr<EcalRecHitCollection> EEMergedRecHits(new EcalRecHitCollection);

 unsigned int nColl = EcalRecHits_done.size();

 int nEB = 0;
 int nEE = 0;


 for (unsigned int i=0; i < nColl; i++) {

   std::string instance = EcalRecHits_done[i].provenance()->productInstanceName();
   std::string module_label = EcalRecHits_done[i].provenance()->moduleLabel();

   if ( module_label != "ecalRegionalEgammaRecHitTmp" && 
	module_label != "ecalRegionalMuonsRecHitTmp" &&
	module_label != "ecalRegionalJetsRecHitTmp" &&
 	module_label != "ecalRegionalTausRecHitTmp" &&
 	module_label != "ecalRegionalRestRecHitTmp" ) continue;

   if (instance == "EcalRecHitsEB")  {
	nEB += EcalRecHits_done[i] -> size();
   }
   else if (instance == "EcalRecHitsEE") {
	nEE += EcalRecHits_done[i] -> size();
   }

 }

 EBMergedRecHits -> reserve(nEB);
 EEMergedRecHits -> reserve(nEE);
 if (debug_) cout << " Number of EB Rechits to merge  = " << nEB << endl;
 if (debug_) cout << " Number of EE Rechits to merge  = " << nEE << endl;

 for (unsigned int i=0; i < nColl; i++) {
   std::string instance = EcalRecHits_done[i].provenance()->productInstanceName(); 

   std::string module_label = EcalRecHits_done[i].provenance()->moduleLabel();

   if ( module_label != "ecalRegionalEgammaRecHitTmp" &&
        module_label != "ecalRegionalMuonsRecHitTmp" &&
        module_label != "ecalRegionalJetsRecHitTmp" &&
        module_label != "ecalRegionalTausRecHitTmp" &&
        module_label != "ecalRegionalRestRecHitTmp" ) continue;

   if (instance == "EcalRecHitsEB") {
	for (EcalRecHitCollection::const_iterator it=EcalRecHits_done[i]->begin(); it !=EcalRecHits_done[i]->end(); it++) {
		EBMergedRecHits -> push_back(*it);
  	}
   }
   else if (instance == "EcalRecHitsEE") {
	for (EcalRecHitCollection::const_iterator it=EcalRecHits_done[i]->begin(); it !=EcalRecHits_done[i]->end(); it++) {
		EEMergedRecHits -> push_back(*it);
	}
   }

 }


 // cout << " avant le put " << endl;
 e.put(EBMergedRecHits,OutputLabelEB_);
 e.put(EEMergedRecHits,OutputLabelEE_);
 // cout << " apres le put " << endl;

}

