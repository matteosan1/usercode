/** \class EgammaHLTElectronTrackIsolationProducers
 *
 *  \author Monica Vazquez Acosta (CERN)
 * 
 * $Id: EgammaHLTElectronTrackIsolationProducers.cc,v 1.3 2007/03/07 09:22:03 monicava Exp $
 *
 */

#include "RecoEgamma/EgammaHLTProducers/interface/EgammaHLTElectronTrackIsolationProducers.h"

// Framework
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Utilities/interface/Exception.h"

#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronIsolationAssociation.h"

#include "DataFormats/Common/interface/RefToBase.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/RefProd.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"




EgammaHLTElectronTrackIsolationProducers::EgammaHLTElectronTrackIsolationProducers(const edm::ParameterSet& config) : conf_(config)
{

  electronProducer_         = conf_.getParameter<edm::InputTag>("electronProducer");
  trackProducer_                = conf_.getParameter<edm::InputTag>("trackProducer");

  egTrkIsoPtMin_                = conf_.getParameter<double>("egTrkIsoPtMin");
  egTrkIsoConeSize_             = conf_.getParameter<double>("egTrkIsoConeSize");
  egTrkIsoZSpan_                = conf_.getParameter<double>("egTrkIsoZSpan");
  egTrkIsoRSpan_                = conf_.getParameter<double>("egTrkIsoRSpan");
  egTrkIsoVetoConeSize_         = conf_.getParameter<double>("egTrkIsoVetoConeSize");


  test_ = new EgammaHLTTrackIsolation(egTrkIsoPtMin_,egTrkIsoConeSize_,
				      egTrkIsoZSpan_,egTrkIsoRSpan_,egTrkIsoVetoConeSize_);


  //register your products
  produces < reco::ElectronIsolationMap >();

}


EgammaHLTElectronTrackIsolationProducers::~EgammaHLTElectronTrackIsolationProducers(){delete test_;}


//
// member functions
//

// ------------ method called to produce the data  ------------
void
EgammaHLTElectronTrackIsolationProducers::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  // Get the HLT filtered objects
  edm::Handle<reco::ElectronCollection> electronHandle;
  iEvent.getByLabel(electronProducer_,electronHandle);

 // Get the general tracks
  edm::Handle<reco::TrackCollection> trackHandle;
  iEvent.getByLabel(trackProducer_, trackHandle);
  const reco::TrackCollection* trackCollection = trackHandle.product();

  reco::ElectronIsolationMap isoMap;

  for(reco::ElectronCollection::const_iterator iElectron = electronHandle->begin(); iElectron != electronHandle->end(); iElectron++){


    reco::ElectronRef electronref(reco::ElectronRef(electronHandle,iElectron - electronHandle->begin()));
    reco::TrackRef electrontrackref = iElectron->track();

    float isol = test_->electronPtSum(&(*electrontrackref),trackCollection);
    if(electrontrackref->pt() != 0. ) isol = isol/electrontrackref->pt();

    isoMap.insert(electronref, isol);


  }

  std::auto_ptr<reco::ElectronIsolationMap> isolMap(new reco::ElectronIsolationMap(isoMap));
  iEvent.put(isolMap);

}

//define this as a plug-in
//DEFINE_FWK_MODULE(EgammaHLTTrackIsolationProducers);
