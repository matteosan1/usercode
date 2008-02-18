// -*- C++ -*-
//
// Package:    EgammaHLTProducers
// Class:      EgammaHLTPixelMatchElectronProducers
// 
/**\class EgammaHLTPixelMatchElectronProducers RecoEgamma/ElectronProducers/src/EgammaHLTPixelMatchElectronProducers.cc

 Description: EDProducer of HLT Electron objects

*/
//
// Original Author: Monica Vazquez Acosta (CERN)
// $Id: EgammaHLTPixelMatchElectronProducers.cc,v 1.2 2007/10/16 09:13:48 ghezzi Exp $
//
//

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTPixelMatchElectronProducers.h"
#include "RecoEgamma/EgammaHLTAlgos/interface/EgammaHLTPixelMatchElectronAlgo.h"
#include "DataFormats/EgammaReco/interface/ElectronPixelSeedFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronPixelSeed.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"

// OUR ALGO
#include "RecoEgamma/EgammaHLTAlgos/interface/EgammaHLTGlobalElectronAlgo.h"

#include <iostream>

using namespace reco;
 
EgammaHLTPixelMatchElectronProducers::EgammaHLTPixelMatchElectronProducers(const edm::ParameterSet& iConfig) : conf_(iConfig) {
  //register your products
  produces<ElectronCollection>();

  std::string algoType_ = conf_.getParameter<std::string>("AlgoType");

  //create algo
  
  if (algoType_ == "std") {
    algoStd_ = new EgammaHLTPixelMatchElectronAlgo();
  } else {
    algoGge_ = new EgammaHLTGlobalElectronAlgo();
  }
}


EgammaHLTPixelMatchElectronProducers::~EgammaHLTPixelMatchElectronProducers() {
  delete algoStd_;
  delete algoGge_;
}

void EgammaHLTPixelMatchElectronProducers::beginJob(edm::EventSetup const&iSetup) {
  if (algoType_ == "std") {
    algoStd_->setupES(iSetup,conf_);  
  } else {
    algoGge_->setupES(iSetup,conf_);  
  } 
}

// ------------ method called to produce the data  ------------
void EgammaHLTPixelMatchElectronProducers::produce(edm::Event& e, const edm::EventSetup& iSetup) {

  // Create the output collections   
  std::auto_ptr<ElectronCollection> pOutEle(new ElectronCollection);
  
  // invoke algorithm
  if (algoType_ == "std") {
    algoStd_->run(e,*pOutEle);
  } else {
    algoGge_->run(e,*pOutEle);
  } 
    

  // put result into the Event
    e.put(pOutEle);
  
}


