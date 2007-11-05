// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "RecoEgammaBis/EgammaElectronProducers/plugins/GlobalGsfElectronProducer.h"
#include "RecoEgammaBis/EgammaElectronAlgos/interface/ElectronAlgoB.h"
#include "DataFormats/EgammaReco/interface/ElectronPixelSeed.h"
#include "DataFormats/EgammaReco/interface/SeedSuperClusterAssociation.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"

#include <iostream>

using namespace reco;
 
GlobalGsfElectronProducer::GlobalGsfElectronProducer(const edm::ParameterSet& iConfig) : conf_(iConfig)
{
  //register your products
  produces<PixelMatchGsfElectronCollection>();

  //create algo
  algo_ = new
    ElectronAlgoB(iConfig.getParameter<double>("maxEOverPBarrel"),
		  iConfig.getParameter<double>("maxEOverPEndcaps"),
                           iConfig.getParameter<double>("minEOverPBarrel"),
			   iConfig.getParameter<double>("minEOverPEndcaps"),
			   iConfig.getParameter<double>("hOverEConeSize"),
			   iConfig.getParameter<double>("maxHOverE"),
			   iConfig.getParameter<double>("maxDeltaEta"),
			   iConfig.getParameter<double>("maxDeltaPhi"),
			   iConfig.getParameter<double>("ptCut"),
			   iConfig.getParameter<bool>("highPtPreselection"),
			   iConfig.getParameter<double>("highPtMin"));
}


GlobalGsfElectronProducer::~GlobalGsfElectronProducer()
{
  delete algo_;
}

void GlobalGsfElectronProducer::beginJob(edm::EventSetup const&iSetup) 
{     
  algo_->setupES(iSetup,conf_);  
}

// ------------ method called to produce the data  ------------
void GlobalGsfElectronProducer::produce(edm::Event& e, const edm::EventSetup& iSetup) 
{

  // Create the output collections   
  std::auto_ptr<PixelMatchGsfElectronCollection> pOutEle(new PixelMatchGsfElectronCollection);
  
  // invoke algorithm
    algo_->run(e,*pOutEle);

  // put result into the Event
    e.put(pOutEle);
  
}


