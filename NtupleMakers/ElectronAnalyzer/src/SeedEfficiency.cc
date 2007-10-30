#include "NtupleMakers/ElectronAnalyzer/interface/SeedEfficiency.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
#include "DataFormats/TrajectorySeed/interface/BasicTrajectorySeed.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"

#include "SimTracker/TrackerHitAssociation/interface/TrackerHitAssociator.h"

#include "DataFormats/Common/interface/OwnVector.h"
//#include "DataFormats/Common/interface/RangeMap.h"

#include "TFile.h"
#include "TTree.h"

using namespace std;
using namespace edm;
//using namespace reco;

typedef edm::OwnVector<TrackingRecHit> recHitContainer;
typedef recHitContainer::const_iterator const_iterator;
typedef std::pair<const_iterator,const_iterator> range;
   
SeedEfficiency::SeedEfficiency(const ParameterSet& pset) {
  fileName = pset.getParameter<std::string>("RootFileName");
  trajectorySeedName =  pset.getParameter<std::string>("TrajectorySeedName");
  simHitName   =  pset.getParameter<std::string>("SimHitName");
  for(int i=0; i<3; ++i)
    nAss[i] = 0;

  nSeed = 0;
}

SeedEfficiency::~SeedEfficiency() {}

void SeedEfficiency::beginJob(const EventSetup& eventSetup) {

  file = new TFile(fileName.c_str(), "recreate");
  tree = new TTree("event","Event data");
}

void SeedEfficiency::endJob() {

  cout << "FINAL REPORT: " << endl;
  cout << "Total number of seeds: " << nSeed << endl;
  cout << "Seed with 0 ass hits : " << nAss[0] << endl;
  cout << "Seed with 1 ass hits : " << nAss[1] << endl;
  cout << "Seed with 2 ass hits : " << nAss[2] << endl;
  
  file->Write();
  file->Close();
}

void SeedEfficiency::analyze(const Event & event, const EventSetup& eventSetup) {
  
  cout << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;
  
  Handle<TrajectorySeedCollection> htjs;
  event.getByLabel(trajectorySeedName, htjs);
  const TrajectorySeedCollection* tjSeed = htjs.product();
  
  TrackerHitAssociator *associator = new TrackerHitAssociator(event);

  nSeed += tjSeed->size();
  int el;

  for (TrajectorySeedCollection::const_iterator itSeed = tjSeed->begin(); itSeed != tjSeed->end(); ++itSeed) { 
    
    //cout << itSeed->nHits() << endl;
    
    vector<PSimHit> assSimHits;
    range rhRange = itSeed->recHits();   
    el = 0;
    for(recHitContainer::const_iterator itRH = rhRange.first; itRH != rhRange.second; ++itRH) {
      //cout << (*itRH).localPosition() << endl;
      assSimHits = associator->associateHit(*itRH);
      for(unsigned int i=0; i<assSimHits.size(); ++i) {
        if (abs(assSimHits[i].particleType()) == 11) {
          cout << assSimHits[i].particleType() << endl;
          el++;
        }
      }
    }  
    
    nAss[el]++;

  }
}

