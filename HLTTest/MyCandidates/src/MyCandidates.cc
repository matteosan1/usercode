// -*- C++ -*-
//
// Package:    MyCandidates
// Class:      MyCandidates
//
// Original Author:  Matteo Sani,40 3-A02,+41227671577,
//         Created:  Thu Feb 14 14:06:52 CET 2013
// $Id$
//
//


#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/ElectronIsolationAssociation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaHcalIsolation.h"
#include "RecoEgamma/EgammaHLTAlgos/interface/EgammaHLTHcalIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/ElectronTkIsolation.h"
#include "RecoEgamma/EgammaIsolationAlgos/interface/EgammaRecHitIsolation.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputer.h"
#include "RecoLocalCalo/HcalRecAlgos/interface/HcalSeverityLevelComputerRcd.h"
#include "CondFormats/HcalObjects/interface/HcalChannelQuality.h"
#include "CondFormats/DataRecord/interface/HcalChannelQualityRcd.h"

#include "Geometry/Records/interface/CaloGeometryRecord.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"

#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"

#include <string>
#include <iostream>

class MyCandidates : public edm::EDAnalyzer {
public:
  explicit MyCandidates(const edm::ParameterSet&);
  ~MyCandidates();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  int matchElectronToSC(trigger::TriggerObject electronRef, edm::Handle<reco::RecoEcalCandidateCollection> candsH);
  int findEleRef(reco::RecoEcalCandidateRef ref, edm::Handle<reco::ElectronCollection> hltEleH); 
  int findEleRef(reco::RecoEcalCandidateRef ref, edm::Handle<reco::GsfElectronCollection> eleH);

  Int_t hlt_n0;
  Float_t hlt_et0[8];
  Float_t hlt_eta0[8];
  Float_t hlt_phi0[8];
  Float_t hlt_ecaliso0[8];
  Float_t hlt_sieie0[8];
  Float_t hlt_tkiso0[8];
  Float_t hlt_hcaliso0[8];
  Float_t hlt_hoe0[8];
  Float_t hlt_deta0[8];
  Float_t hlt_dphi0[8];
  Int_t hlt_mishits0[8];

  Int_t hlt_n1;
  Float_t hlt_et1[8];
  Float_t hlt_eta1[8];
  Float_t hlt_phi1[8];
  Float_t hlt_ecaliso1[8];
  Float_t hlt_sieie1[8];
  Float_t hlt_tkiso1[8];
  Float_t hlt_hcaliso1[8];
  Float_t hlt_hoe1[8];
  Float_t hlt_deta1[8];
  Float_t hlt_dphi1[8];
  
  Int_t el_n;
  Float_t el_eta[10];
  Float_t el_et[10];
  Float_t el_phi[10];
  Float_t el_iso[10];
  Float_t el_isostd[10];
  Float_t el_sieie[10];
  Float_t el_tkiso[10];
  Float_t el_hcaliso[10];
  Float_t el_hoe[10];
  Float_t el_hcalisostd[10];
  Float_t el_hoestd[10];
  Float_t el_deta[10];
  Float_t el_dphi[10];
  Int_t el_mishits[10];
  
  Float_t rhoReco1, rhoReco2, rhoHlt;
  Int_t nVtxReco, nVtxHlt;

  std::string fileName;
  TFile* file;
  TTree* tree;

  edm::InputTag tagCandCollection, isoMapTag, hltEleCandidates;
  edm::InputTag probeCandCollection, isoMapProbe;
  edm::InputTag sieieMapProbe, hoeMapProbe, hcalMapProbe;
  edm::InputTag detaMapProbe, dphiMapProbe, tkisoMapProbe;
  
  EgammaHLTHcalIsolation* hcalIsolationAlgoHE;
  EgammaHLTHcalIsolation* hcalIsolationAlgoHcal;
};

MyCandidates::MyCandidates(const edm::ParameterSet& iConfig) {
  
  fileName = iConfig.getParameter<std::string>("RootFileName");
  probeCandCollection = iConfig.getParameter<edm::InputTag>("ProbeCandCollection");
  isoMapProbe = iConfig.getParameter<edm::InputTag>("IsoMapProbeCollection");
  sieieMapProbe = iConfig.getParameter<edm::InputTag>("SieieMapProbeCollection");
  hcalMapProbe = iConfig.getParameter<edm::InputTag>("HcalMapProbeCollection");
  hoeMapProbe = iConfig.getParameter<edm::InputTag>("HoeMapProbeCollection");
  detaMapProbe = iConfig.getParameter<edm::InputTag>("DetaMapProbeCollection");
  dphiMapProbe = iConfig.getParameter<edm::InputTag>("DphiMapProbeCollection");
  tkisoMapProbe = iConfig.getParameter<edm::InputTag>("TkisoMapProbeCollection");
  tagCandCollection = iConfig.getParameter<edm::InputTag>("TagCandCollection");
  isoMapTag = iConfig.getParameter<edm::InputTag>("IsoMapTagCollection");
  hltEleCandidates = iConfig.getParameter<edm::InputTag>("HLTEleCollection");
}

MyCandidates::~MyCandidates() 
{}

int MyCandidates::findEleRef(reco::RecoEcalCandidateRef ref, edm::Handle<reco::ElectronCollection> hltEleH) {

  int index = -1;
  for (unsigned int i=0; i<hltEleH->size(); i++) {
    reco::ElectronRef cand(hltEleH, i);
    if (cand->superCluster() == ref->superCluster()) 
      return i;
  }

  return index;
}


int MyCandidates::findEleRef(reco::RecoEcalCandidateRef ref, edm::Handle<reco::GsfElectronCollection> eleH) {

  int index = -1;
  float minDR = 0.2;
  for (unsigned int i=0; i<eleH->size(); i++) {
    reco::GsfElectronRef cand(eleH, i);
    float dR = deltaR(ref->p4(), cand->p4());
    if (dR < minDR) {
      index = i;
      minDR = dR;
    }
  }

  return index;
}

int MyCandidates::matchElectronToSC(trigger::TriggerObject electronRef, edm::Handle<reco::RecoEcalCandidateCollection> candsH) {
  
  int index = -1;
  float minDR = 0.2;
  for (unsigned int i=0; i<candsH->size(); i++) {
    reco::RecoEcalCandidateRef cand(candsH, i);
    float dR = deltaR(cand->p4(), electronRef.particle().p4());
    if (dR < minDR) {
      index = i;
      minDR = dR;
    }
  }
  
  return index;
}

void MyCandidates::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<EcalRecHitCollection> ecalBarrelRecHitHandle; //EcalRecHitCollection is a typedef to
  iEvent.getByLabel("ecalRecHit", "EcalRecHitsEB", ecalBarrelRecHitHandle);

  edm::Handle<EcalRecHitCollection> ecalEndcapRecHitHandle;
  iEvent.getByLabel("ecalRecHit", "EcalRecHitsEE", ecalEndcapRecHitHandle);

  EcalRecHitMetaCollection ecalBarrelHits(*ecalBarrelRecHitHandle);
  EcalRecHitMetaCollection ecalEndcapHits(*ecalEndcapRecHitHandle);
  
  edm::ESHandle<CaloGeometry> caloGeomHandle;
  iSetup.get<CaloGeometryRecord>().get(caloGeomHandle);
  const CaloGeometry* caloGeom = caloGeomHandle.product();
  //if (caloGeomCacheId_!=iSetup.get<CaloGeometryRecord>().cacheIdentifier()) {
  //  iSetup.get<CaloGeometryRecord>().get(caloGeom_);
  //  caloGeomCacheId_=iSetup.get<CaloGeometryRecord>().cacheIdentifier();
  //}
  
  //edm::ESHandle<EcalSeverityLevelAlgo> sevlv;
  //iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);
  const EcalSeverityLevelAlgo* sevLevel = 0;//sevlv.product();

  EgammaRecHitIsolation ecalBarrelIsol(0.3, 3.0, 1.5, 0., 0.095, caloGeomHandle, &ecalBarrelHits, sevLevel, DetId::Ecal);
  ecalBarrelIsol.setUseNumCrystals(true);
  EgammaRecHitIsolation ecalEndcapIsol(0.3, 3.0, 1.5, 0.110, 0., caloGeomHandle, &ecalEndcapHits, sevLevel, DetId::Ecal);
  ecalEndcapIsol.setUseNumCrystals(true);

  //edm::Handle<reco::VertexCollection> vertexH;
  //iEvent.getByLabel("offlinePrimaryVerticesWithBS", vertexH);
  //nVtxReco = vertexH->size();
  edm::Handle<reco::VertexCollection> vertexH;
  iEvent.getByLabel("pixelVertices", vertexH);
  nVtxReco = vertexH->size();
  
  edm::Handle<reco::VertexCollection> hltVertexH;
  iEvent.getByLabel("hltPixelVertices", hltVertexH);
  nVtxHlt = hltVertexH->size();
  
  edm::Handle<double> rho1H;
  iEvent.getByLabel("kt6PFJets", "rho", rho1H);
  rhoReco1 = *(rho1H.product());
  
  edm::Handle<double> rho2H;
  iEvent.getByLabel("kt6CaloJets", "rho", rho2H);
  rhoReco2 = *(rho2H.product());
  
  edm::Handle<double> hltRhoH;
  iEvent.getByLabel("hltKT6CaloJets", "rho", hltRhoH);
  rhoHlt = *(hltRhoH.product());
 
  // Get the barrel hcal hits
  //edm::Handle<HBHERecHitCollection>  hbheRecHitHandle;
  //iEvent.getByLabel("hbhereco", hbheRecHitHandle);
  //const HBHERecHitCollection* hbheRecHitCollection = hbheRecHitHandle.product();
  ////HBHERecHitMetaCollection* hbhe= new HBHERecHitMetaCollection(*hbheRecHitHandle);
  //
  //edm::Handle<HBHERecHitCollection>  hbheRecHitHandle2;
  //iEvent.getByLabel("hltHbhereco", hbheRecHitHandle2);
  //const HBHERecHitCollection* HLThbheRecHitCollection = hbheRecHitHandle2.product();
  //
  //edm::ESHandle<HcalChannelQuality> hcalChStatus;
  ////iSetup.get<HcalChannelQualityRcd>().get(hcalChStatus);
  //edm::ESHandle<HcalSeverityLevelComputer> hcalSevLvlComp;
  ////iSetup.get<HcalSeverityLevelComputerRcd>().get(hcalSevLvlComp);
  //
  //hcalIsolationAlgo = new EgammaHcalIsolation(0.3, 0.15, 0.7, 0.8, -1., -1., caloGeomHandle, hbhe);
  hcalIsolationAlgoHE = new EgammaHLTHcalIsolation(0.7, 0.8, -1., -1., 0.0, 0.15, -1);
  hcalIsolationAlgoHcal = new EgammaHLTHcalIsolation(0.7, 0.8, -1., -1., 0.16, 0.29, -1);

  //edm::Handle<reco::TrackCollection> trackHandle;
  //iEvent.getByLabel("hltEcalActivityEgammaRegionalCTFFinalFitWithMaterial", trackHandle);
  //const reco::TrackCollection* trackCollection = trackHandle.product();
  //
  //edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  //iEvent.getByLabel("hltOnlineBeamSpot", recoBeamSpotHandle);
  //const reco::BeamSpot::Point& beamSpotPosition = recoBeamSpotHandle->position(); 
  //ElectronTkIsolation tkIsoAlgo(0.3, 0.015, 0.015, 0.015, 0.015, 1.0, 0.2, 999999., trackCollection, beamSpotPosition);

  edm::Handle<reco::GsfElectronCollection> eleH;
  iEvent.getByLabel("gsfElectrons", eleH);

  edm::Handle<edm::ValueMap<double> > isoH;   
  iEvent.getByLabel("eleIsoFromDepsEcalFromHitsByCrystalFull03", isoH);

  edm::Handle<edm::ValueMap<double> > tkIsoH;   
  iEvent.getByLabel("eleIsoFromDepsTk03", tkIsoH);

  edm::Handle<reco::RecoEcalCandidateCollection> tagCandsH;
  iEvent.getByLabel(tagCandCollection, tagCandsH);

  edm::Handle<reco::RecoEcalCandidateCollection> probeCandsH;
  iEvent.getByLabel(probeCandCollection, probeCandsH);

  edm::Handle<reco::ElectronCollection> hltEleH;
  iEvent.getByLabel(hltEleCandidates, hltEleH);

  /*
  // GSF
  edm::Handle<reco::RecoEcalCandidateIsolationMap> vmProbe4H;
  const reco::RecoEcalCandidateIsolationMap* mapDeta = 0;
  iEvent.getByLabel(detaMapProbe, vmProbe4H);
  if (!vmProbe4H.failedToGet())
    mapDeta = vmProbe4H.product();

  edm::Handle<reco::RecoEcalCandidateIsolationMap> vmProbe5H;
  const reco::RecoEcalCandidateIsolationMap* mapDphi = 0;
  iEvent.getByLabel(dphiMapProbe, vmProbe5H);
  if (!vmProbe5H.failedToGet())
    mapDphi = vmProbe5H.product();
  */

  edm::Handle<reco::ElectronIsolationMap> vmProbe4H;
  const reco::ElectronIsolationMap* mapDeta = 0;
  iEvent.getByLabel(detaMapProbe, vmProbe4H);
  if (!vmProbe4H.failedToGet())
    mapDeta = vmProbe4H.product();

  edm::Handle<reco::ElectronIsolationMap> vmProbe5H;
  const reco::ElectronIsolationMap* mapDphi = 0;
  iEvent.getByLabel(dphiMapProbe, vmProbe5H);
  if (!vmProbe5H.failedToGet())
    mapDphi = vmProbe5H.product();

  edm::Handle<reco::ElectronIsolationMap> vmProbe6H;
  const reco::ElectronIsolationMap* mapTkiso = 0;
  iEvent.getByLabel(tkisoMapProbe, vmProbe6H);
  if (!vmProbe6H.failedToGet())
    mapTkiso = vmProbe6H.product();

  edm::Handle<reco::RecoEcalCandidateIsolationMap> vmTagH;
  const reco::RecoEcalCandidateIsolationMap* mapTag = 0;
  iEvent.getByLabel(isoMapTag, vmTagH);
  if (!vmTagH.failedToGet())
    mapTag = vmTagH.product();
 
  edm::Handle<reco::RecoEcalCandidateIsolationMap> vmProbe1H;
  const reco::RecoEcalCandidateIsolationMap* mapSieie = 0;
  iEvent.getByLabel(sieieMapProbe, vmProbe1H);
  if (!vmProbe1H.failedToGet())
    mapSieie = vmProbe1H.product();

  edm::Handle<reco::RecoEcalCandidateIsolationMap> vmProbe2H;
  const reco::RecoEcalCandidateIsolationMap* mapHoe = 0;
  iEvent.getByLabel(hoeMapProbe, vmProbe2H);
  if (!vmProbe2H.failedToGet())
    mapHoe = vmProbe2H.product();

  edm::Handle<reco::RecoEcalCandidateIsolationMap> vmProbe3H;
  const reco::RecoEcalCandidateIsolationMap* mapHcal = 0;
  iEvent.getByLabel(hcalMapProbe, vmProbe3H);
  if (!vmProbe3H.failedToGet())
    mapHcal = vmProbe3H.product();

  edm::Handle<reco::RecoEcalCandidateIsolationMap> vmProbeH;
  const reco::RecoEcalCandidateIsolationMap* mapProbe = 0;
  iEvent.getByLabel(isoMapProbe, vmProbeH);
  if (!vmProbeH.failedToGet())
    mapProbe = vmProbeH.product();

  //added Trigger Objects
  edm::InputTag trigEventTag("hltTriggerSummaryAOD","","TEST"); 
  edm::Handle<trigger::TriggerEvent> trigEvent; 
  iEvent.getByLabel(trigEventTag,trigEvent);

  std::vector<std::string> temp_names;
  temp_names.clear();
  
  temp_names.push_back("hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter");//hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17PMMassFilter");
  temp_names.push_back("hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter");//hltEle32CaloIdTCaloIsoTTrkIdTTrkIsoTSC17TrackIsoFilter");

  std::vector<std::string>::iterator filter_it;

  std::vector<trigger::TriggerObject> ElectronRefs0;
  std::vector<trigger::TriggerObject> ElectronRefs1;

  if (trigEvent.isValid()) {
    for (filter_it = temp_names.begin(); filter_it != temp_names.end(); ++filter_it){
      const trigger::TriggerObjectCollection & triggerObjects = trigEvent->getObjects();
      trigger::size_type filter1_idx = trigEvent -> filterIndex (edm::InputTag(*filter_it,"","TEST") ) ;   
      trigger::size_type n_filters    = trigEvent -> sizeFilters();
      if ( filter1_idx < n_filters ) {
	const trigger::Keys & triggerKeys ( trigEvent -> filterKeys ( filter1_idx ) );
	const int nkeys = triggerKeys.size();
	for (int ikey = 0; ikey < nkeys; ++ikey ) {
	  if (*filter_it == "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter")
	    ElectronRefs0.push_back(triggerObjects[ triggerKeys [ikey] ]);
	  else if (*filter_it == "hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter")
	    ElectronRefs1.push_back(triggerObjects[ triggerKeys [ikey]]);	  
	}
      }
    }
  }

  hlt_n0 = 0;
  for (unsigned int i=0; i<ElectronRefs0.size(); i++) {
    if (hlt_n0 >= 8) 
      break;
    trigger::TriggerObject pho = ElectronRefs0[i];
    hlt_eta0[i]       = pho.eta();
    hlt_phi0[i]       = pho.phi();
    hlt_et0[i]        = pho.et();
    
    int index = matchElectronToSC(pho, probeCandsH); 
    hlt_ecaliso0[i] = 99999.;
    hlt_hcaliso0[i] = 99999.;
    hlt_hoe0[i]     = 99999.;
    hlt_sieie0[i]   = 99999.;
    hlt_deta0[i]   = 99999.;
    hlt_dphi0[i]   = 99999.;
    hlt_tkiso0[i]   = 99999.;

    //std::cout << index << std::endl;    
    if (index != -1) {
      reco::RecoEcalCandidateRef cand(probeCandsH, index);
      //int myEle = findEleRef(cand, eleH);
      //if (myEle != -1) {
      //	Float_t temp1 = hcalIsolationAlgo->getESum(cand->superCluster()->eta(), cand->superCluster()->phi(), HLThbheRecHitCollection, caloGeom, hcalSevLvlComp.product(), hcalChStatus.product()); 
      //	reco::GsfElectronRef ele(eleH, myEle);
      //	Float_t temp2 = hcalIsolationAlgo->getESum(ele->superCluster()->eta(), ele->superCluster()->phi(), hbheRecHitCollection, caloGeom, hcalSevLvlComp.product(), hcalChStatus.product()); 
      //	if (temp1 != temp2) {
      //	  std::cout << iEvent.id().run() << " " << iEvent.id().event() << " - " << deltaR(cand->p4(), ele->p4()) << " - ";
      //	  std::cout << "ONLINE: " << temp1 << " - ";
      //	  std::cout << "OFFLINE: " << temp2 << std::endl;
      //	}
      //}

      if (mapProbe)
	hlt_ecaliso0[i] = (*mapProbe)[cand];

      if (mapHcal)
	hlt_hcaliso0[i] = (*mapHcal)[cand];

      if (mapHoe) 
	hlt_hoe0[i] = (*mapHoe)[cand];

      if (mapSieie)
	hlt_sieie0[i] = (*mapSieie)[cand];
      
      int eleInd = findEleRef(cand, hltEleH);
      //std::cout << eleInd << std::endl;
      if (eleInd != -1) {
	reco::ElectronRef e(hltEleH, eleInd);
	hlt_mishits0[i] = e->track()->trackerExpectedHitsInner().numberOfLostHits();

	if (mapDeta)
	  hlt_deta0[i] = (*mapDeta)[e];
	if (mapDphi)
	  hlt_dphi0[i] = (*mapDphi)[e];
	if (mapTkiso)
	  hlt_tkiso0[i] = (*mapTkiso)[e];
      }
    }

    hlt_n0++;
  }
  
  hlt_n1 = 0;
  for (unsigned int i=0; i<ElectronRefs1.size(); i++) {
    if (hlt_n1 >= 8) 
      break;
    trigger::TriggerObject pho=ElectronRefs1[i];
    hlt_eta1[i] = pho.eta();
    hlt_phi1[i] = pho.phi();
    hlt_et1[i]  = pho.et();

    hlt_ecaliso1[i]   = 99999.;
    if (mapTag) {
      int index = matchElectronToSC(pho, tagCandsH); 
      if (index != -1) {
	reco::RecoEcalCandidateRef cand(tagCandsH, index);
	hlt_ecaliso1[i]   = (*mapTag)[cand];
      }
    }    
    
    hlt_n1++;
  }

  el_n = 0;
  for (unsigned int i=0; i<eleH->size(); i++) {
    reco::GsfElectronRef ele(eleH, i);
    el_et[el_n]  = ele->superCluster()->energy()*sin((2*atan(exp(-ele->superCluster()->eta()))));
    el_eta[el_n] = ele->superCluster()->eta();
    el_phi[el_n] = ele->superCluster()->phi();
    el_iso[el_n] = ecalBarrelIsol.getEtSum(&(*ele)) + ecalEndcapIsol.getEtSum(&(*ele));
    //el_iso[el_n] = (*isoH)[ele];
    el_isostd[el_n] = ele->dr03EcalRecHitSumEt();
    el_sieie[el_n] = ele->sigmaIetaIeta();
    //const reco::Track* eleTrk = &*ele->gsfTrack();
    //el_tkiso[el_n] = tkIsoAlgo.getIso(eleTrk).second;//(*tkIsoH)[ele];//ele->dr03TkSumPt();

    el_hoestd[el_n] = ele->hcalOverEcal();
    el_hcalisostd[el_n] = ele->dr03HcalTowerSumEt();  
    //el_hoe[el_n] = hcalIsolationAlgoHE->getESum(ele->superCluster()->eta(), ele->superCluster()->phi(), hbheRecHitCollection,caloGeom, hcalSevLvlComp.product(), hcalChStatus.product()); 
    //el_hcaliso[el_n] = hcalIsolationAlgoHcal->getEtSum(ele->superCluster()->eta(), ele->superCluster()->phi(),hbheRecHitCollection, caloGeom, hcalSevLvlComp.product(), hcalChStatus.product()); 
    el_deta[el_n] = ele->deltaEtaSuperClusterTrackAtVtx();
    el_dphi[el_n] = ele->deltaPhiSuperClusterTrackAtVtx();
    el_mishits[el_n] = ele->gsfTrack()->trackerExpectedHitsInner().numberOfLostHits();
    
    el_n++;
  }

  tree->Fill();
}

void  MyCandidates::beginJob() {

  file = new TFile(fileName.c_str(), "recreate");
  tree = new TTree("event", "Event data");
  
  tree->Branch("rhoReco1", &rhoReco1, "rhoReco1/F");
  tree->Branch("rhoReco2", &rhoReco2, "rhoReco2/F");
  tree->Branch("rhoHlt",   &rhoHlt,   "rhoHlt/F");
  tree->Branch("nvtxHlt",  &nVtxHlt,  "nvtxHlt/I");
  tree->Branch("nvtxReco", &nVtxReco, "nvtxReco/I");

  tree->Branch("n0",       &hlt_n0,       "n0/I");
  tree->Branch("eta0",     &hlt_eta0,     "eta0[n0]/F");
  tree->Branch("pt0",      &hlt_et0,      "pt0[n0]/F");
  tree->Branch("phi0",     &hlt_phi0,     "phi0[n0]/F");
  tree->Branch("ecaliso0", &hlt_ecaliso0, "ecaliso0[n0]/F");
  tree->Branch("tkiso0",   &hlt_tkiso0,   "tkiso0[n0]/F");
  tree->Branch("hcaliso0", &hlt_hcaliso0, "hcaliso0[n0]/F");
  tree->Branch("hoe0",     &hlt_hoe0,     "hoe0[n0]/F");
  tree->Branch("sieie0",   &hlt_sieie0,   "sieie0[n0]/F");
  tree->Branch("deta0",    &hlt_deta0,    "deta0[n0]/F");
  tree->Branch("dphi0",    &hlt_dphi0,    "dphi0[n0]/F");
  tree->Branch("mishits0", &hlt_mishits0, "mishits0[n0]/I");

  tree->Branch("n1"   ,    &hlt_n1,       "n1/I");
  tree->Branch("eta1",     &hlt_eta1,     "eta1[n1]/F");
  tree->Branch("pt1",      &hlt_et1,      "pt1[n1]/F");
  tree->Branch("phi1",     &hlt_phi1,     "phi1[n1]/F");
  tree->Branch("ecaliso1", &hlt_ecaliso1, "ecaliso1[n1]/F");
  tree->Branch("tkiso1",   &hlt_tkiso1,   "tkiso1[n1]/F");
  tree->Branch("hcaliso1", &hlt_hcaliso1, "hcaliso1[n1]/F");
  tree->Branch("hoe1",     &hlt_hoe1,     "hoe1[n1]/F");
  tree->Branch("sieie1",   &hlt_sieie1,   "sieie1[n1]/F");
  tree->Branch("deta1",    &hlt_deta1,    "deta1[n1]/F");
  tree->Branch("dphi1",    &hlt_dphi1,    "dphi1[n1]/F");

  tree->Branch("el_n",          &el_n,          "el_n/I");
  tree->Branch("el_eta",        &el_eta,        "el_eta[el_n]/F");
  tree->Branch("el_et",         &el_et,         "el_et [el_n]/F");
  tree->Branch("el_phi",        &el_phi,        "el_phi[el_n]/F");
  tree->Branch("el_iso",        &el_iso,        "el_iso[el_n]/F");
  tree->Branch("el_isostd",     &el_isostd,     "el_isostd[el_n]/F");
  tree->Branch("el_sieie",      &el_sieie,      "el_sieie[el_n]/F");
  tree->Branch("el_tkiso",      &el_tkiso,      "el_tkiso[el_n]/F");
  tree->Branch("el_hoe",        &el_hoe,        "el_hoe[el_n]/F");
  tree->Branch("el_hcaliso",    &el_hcaliso,    "el_hcaliso[el_n]/F");
  tree->Branch("el_hoestd",     &el_hoestd,     "el_hoestd[el_n]/F");
  tree->Branch("el_hcalisostd", &el_hcalisostd, "el_hcalisostd[el_n]/F");
  tree->Branch("el_deta",       &el_deta,       "el_deta[el_n]/F");
  tree->Branch("el_dphi",       &el_dphi,       "el_dphi[el_n]/F");
  tree->Branch("el_mishits",    &el_mishits,    "el_mishits[n0]/I");
}

void MyCandidates::endJob() {
  file->cd();
  tree->Write(0, TObject::kWriteDelete);
  file->Close();
}

void MyCandidates::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(MyCandidates);
