// -*- C++ -*-
//
// Package:    MyCandidates
// Class:      TurnOn
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

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "DataFormats/HLTReco/interface/TriggerTypeDefs.h"
#include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateIsolation.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "CommonTools/Utils/interface/StringCutObjectSelector.h"

#include "EGamma/EGammaAnalysisTools/interface/EGammaCutBasedEleId.h"

#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Common/interface/ValueMap.h"

#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"

#include <string>
#include <iostream>

class TurnOn : public edm::EDAnalyzer {
public:
  typedef std::vector< edm::Handle< edm::ValueMap<double> > > IsoDepositVals;
  
  explicit TurnOn(const edm::ParameterSet&);
  ~TurnOn();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  template <class T1, class T2>  
  bool matchByDeltaR(const T1 &, 
		     const std::vector<T2> &,
		     const double);

  trigger::TriggerObjectCollection selectTriggerObjects(const trigger::TriggerObjectCollection &,
							const trigger::TriggerEvent &, const edm::InputTag&,
							const std::string);

  std::vector<bool> selectElectrons(edm::Handle<reco::GsfElectronCollection>, 
				    const StringCutObjectSelector<reco::GsfElectron> &,
				    edm::Handle<reco::VertexCollection>,
				    edm::Handle<reco::ConversionCollection>,
				    double, reco::BeamSpot, IsoDepositVals, EgammaCutBasedEleId::WorkingPoint wp);

  std::vector<std::pair<trigger::TriggerObject, trigger::TriggerObject>> removeTags(const trigger::TriggerObjectCollection & all,
										    const trigger::TriggerObjectCollection & tags);

  std::vector<reco::GsfElectron> selectOfflineProbe(std::vector<std::pair<trigger::TriggerObject, trigger::TriggerObject>> hltPairs,
						    std::vector<std::pair<reco::GsfElectron, reco::GsfElectron>> offlinePairs);
  
  trigger::TriggerObjectCollection purgeTrailingLeg(const trigger::TriggerObjectCollection & all,
						    const trigger::TriggerObjectCollection & leading);
private:
  virtual void beginJob() ;
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob() ;
  
  std::string fileName;
  TFile* file;
  TTree* tree;

  edm::InputTag inputCollection;
  edm::InputTag triggerResultsLabel;
  edm::InputTag triggerSummaryLabel;
  //  std::vector<std::string> modules;
  //std::string recoCuts;
  std::string hltCuts;
  bool doMatching;
  float dR;

  Float_t pt1, pt2, nvtx1, nvtx2, eta1, eta2;
  Float_t ptden, nvtxden, etaden;

  //std::vector<edm::InputTag> moduleLabels;
  HLTConfigProvider hltConfig;
  std::vector<std::string> realHltPaths, triggerPaths;
  std::vector<edm::InputTag> inputTagIsoValElectronsPFId_;
};

TurnOn::TurnOn(const edm::ParameterSet& iConfig) {
  
  fileName = iConfig.getParameter<std::string>("RootFileName");
  inputTagIsoValElectronsPFId_   = iConfig.getParameter< std::vector<edm::InputTag> >("IsoValElectronPF");
  
  dR = 0.2;
  hltCuts = "";
}

TurnOn::~TurnOn() 
{}

std::vector<reco::GsfElectron> TurnOn::selectOfflineProbe(std::vector<std::pair<trigger::TriggerObject, trigger::TriggerObject>> hltPairs,
							  std::vector<std::pair<reco::GsfElectron, reco::GsfElectron>> offlinePairs) {
  
  float dM = 50;
  std::vector<reco::GsfElectron> offlineProbes;
  //reco::GsfElectron offlineTag;
  
  for (unsigned int i=0; i<hltPairs.size(); i++) {
    for (unsigned int j=0; j<offlinePairs.size(); j++) {
      float dr11 = deltaR(offlinePairs[j].first, hltPairs[i].first);
      float dr22 = deltaR(offlinePairs[j].second, hltPairs[i].second);

      float mass = (offlinePairs[j].first.p4() + offlinePairs[j].second.p4()).M();      
      if (fabs(mass - 91.186) < dM) { 
	if (dr11 < 0.05 && dr22<0.05) {
	  dM = fabs(mass - 91.186);
	  offlineProbes.push_back(offlinePairs[j].second);
	  //offlineTag = offlinePairs[j].first;
	} 
      }
    }
  }
  
  return offlineProbes;
}

trigger::TriggerObjectCollection TurnOn::purgeTrailingLeg(const trigger::TriggerObjectCollection & all,
							  const trigger::TriggerObjectCollection & leading) {
  
  trigger::TriggerObjectCollection result;
  for (unsigned int i=0; i<all.size(); i+=2) {
    //bool isTag = false;
    //std::cout << "all: " << all[i].pt() << " " << all[i+1].pt() << std::endl;
    for (unsigned int j=0; j<leading.size(); j++) {
      float dR1 = deltaR(all[i], leading[j]);
      float dR2 = deltaR(all[i+1], leading[j]);

      if (dR1<0.05)
	result.push_back(all[i+1]);

      if (dR2<0.05)
	result.push_back(all[i]);
    }
  }

  return result;
}

std::vector<std::pair<trigger::TriggerObject, trigger::TriggerObject>> TurnOn::removeTags(const trigger::TriggerObjectCollection & all,
											  const trigger::TriggerObjectCollection & tags) {
  
  std::vector<std::pair<trigger::TriggerObject, trigger::TriggerObject>> pairs;

  for (unsigned int j=0; j<all.size(); j+=2) {  
    //std::cout << "Probes: " << all[j].pt() << " " << all[j].eta() << std::endl;
    //std::cout << "Probes: " << all[j+1].pt() << " " << all[j+1].eta() << std::endl;
    for (unsigned int i=0; i<tags.size(); i++) {
      //std::cout << "TAG: " << tags[i].pt() << std::endl;
      
      float dR = deltaR(tags[i], all[j]);
      float dE = fabs(tags[i].energy() - all[j].energy());
      
      if (dR < 0.05 && dE < 1.0) {
	pairs.push_back(std::make_pair(all[j], all[j+1]));
	//std::cout << all[j].pt() << " " << all[j+1].pt() << std::endl;
	break;
      }
    }
  }
    
  return pairs;     
}

trigger::TriggerObjectCollection TurnOn::selectTriggerObjects(const trigger::TriggerObjectCollection & triggerObjects,
							      const trigger::TriggerEvent & triggerSummary,
							      const edm::InputTag& moduleLabel, const std::string hltcuts) {
  
  trigger::TriggerObjectCollection selectedObjects;
  StringCutObjectSelector<trigger::TriggerObject> selector(hltcuts);

  size_t filterIndex = triggerSummary.filterIndex(moduleLabel);
  
  if (filterIndex < triggerSummary.sizeFilters()) {
    const trigger::Keys &keys = triggerSummary.filterKeys(filterIndex);
    
    for (size_t j = 0; j < keys.size(); j++) {
      trigger::TriggerObject foundObject = triggerObjects[keys[j]];
      if (selector(foundObject)) {
	selectedObjects.push_back(foundObject);
      } 
    }
  }

  return selectedObjects;
}

std::vector<bool> TurnOn::selectElectrons(edm::Handle<reco::GsfElectronCollection> allElectrons, 
					  const StringCutObjectSelector<reco::GsfElectron> &selector,
					  edm::Handle<reco::VertexCollection> vtxs,
					  edm::Handle<reco::ConversionCollection> conversions,
					  double rho,
					  reco::BeamSpot beamspot, IsoDepositVals electronIsoVals, 
					  EgammaCutBasedEleId::WorkingPoint wp) {

  std::vector<bool> result(allElectrons->size(), false);
  
  for (unsigned int i=0; i<allElectrons->size(); i++) {
    reco::GsfElectronRef ref(allElectrons, i);
    float el_pfiso_charged = (*(electronIsoVals[0].product()))[ref];
    float el_pfiso_photon  = (*(electronIsoVals[1].product()))[ref];
    float el_pfiso_neutral = (*(electronIsoVals[2].product()))[ref];
    
    if (selector(*ref) && PassWP(wp, *ref,
				 conversions, beamspot, vtxs,
				 el_pfiso_charged, el_pfiso_photon, el_pfiso_neutral, rho))
      result[i] = true;
  }

  return result;
}

template <class T1, class T2> 
bool TurnOn::matchByDeltaR(const T1 & element1, 
			   const std::vector<T2> & collection2,
			   const double maxDeltaR) {
  
  const size_t n = collection2.size();

  float dr = 0.1;
  for (size_t j = 0; j < n; j++) {
    //std::cout << collection2[j].pt() << " " << collection2[j].eta() << std::endl;
    float r = deltaR(element1, collection2[j]);
    if (r < dr)
      return true;
  }

  return false;
}
      
void TurnOn::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  IsoDepositVals electronIsoVals(3);
  for (size_t j = 0; j<inputTagIsoValElectronsPFId_.size(); ++j) {
    iEvent.getByLabel(inputTagIsoValElectronsPFId_[j], electronIsoVals[j]);
  }

  edm::Handle<reco::VertexCollection> vertexH;
  iEvent.getByLabel("offlinePrimaryVerticesWithBS", vertexH);
  int nVertices = vertexH->size();

  edm::Handle<reco::ConversionCollection> conversionsH;
  iEvent.getByLabel("allConversions", conversionsH);
  
  edm::Handle<double> rhoH;
  iEvent.getByLabel("kt6PFJets", "rho", rhoH);
  double rho = *(rhoH.product());

  edm::Handle<reco::GsfElectronCollection> eleH;
  iEvent.getByLabel("gsfElectrons", eleH);
  reco::GsfElectronCollection electrons = *(eleH.product());

  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  iEvent.getByLabel("offlineBeamSpot", recoBeamSpotHandle);
  reco::BeamSpot bs = *(recoBeamSpotHandle.product()); 

  std::auto_ptr<reco::GsfElectronCollection> offlineProbeElectrons(new reco::GsfElectronCollection());

  //added Trigger Objects
  edm::InputTag trigEventTag("hltTriggerSummaryAOD", "", "HLT"); 
  edm::Handle<trigger::TriggerEvent> trigEvent; 
  iEvent.getByLabel(trigEventTag, trigEvent);

  if(!trigEvent.isValid()) {
    edm::LogError("MyCandidates/TurnOn") << "Missing triggerSummary with label " << triggerSummaryLabel <<std::endl;
    return;
  }

  edm::Handle<edm::TriggerResults> triggerResults;
  edm::InputTag tResults("TriggerResults", "", "HLT");
  iEvent.getByLabel(tResults, triggerResults);
  if(!triggerResults.isValid()) {
    edm::LogError("MyCandidates/TurnOn") << "Missing triggerResults with label " << triggerResultsLabel <<std::endl;
    return;
  }

  // Update when necessary the trigger table
  bool changedConfig = false;
  //if (!hltConfig.init(iEvent.getRun(), iSetup, triggerResultsLabel.process(), changedConfig)) {
  if (!hltConfig.init(iEvent.getRun(), iSetup, "HLT", changedConfig)) {
    edm::LogError("HLTMatchingFilter") << "Initialization of HLTConfigProvider failed!!"; 
    return;
  }

  if (changedConfig) {
    triggerPaths.clear();
    realHltPaths.clear();
    triggerPaths.push_back("HLT_Ele17_CaloIdVT_CaloIsoVT_TrkIdT_TrkIsoVT_Ele8_Mass50_v*");
    
    
    for (size_t i = 0; i < triggerPaths.size(); i++) {
      realHltPaths = hltConfig.matched(hltConfig.triggerNames(), triggerPaths[i]);
      //TRegexp pattern(triggerPaths[i]);
      //for (size_t j = 0; j < hltConfig.triggerNames().size(); j++) {
      //	if (TString(hltConfig.triggerNames()[j]).Contains(pattern)) {
      //	  realHltPaths.push_back(hltConfig.triggerNames()[j]);
      //	}
      // }
    }
  }

  // Check if the event passes the specified trigger
  bool passTrigger = false;
  for (size_t j = 0; j < hltConfig.size(); j++) {
    for (size_t i = 0; i < realHltPaths.size(); i++) {
      if (hltConfig.triggerName(j) == realHltPaths[i]) {
	if (triggerResults->accept(j)) { 
	  passTrigger = true;
	  break;
	}
      }
    }
  }

  if (!passTrigger)
    return;

  // Define Offline T&P pairs
  std::vector<bool> probeOfflineSelection = selectElectrons(eleH, std::string(""), vertexH, conversionsH, rho, bs, electronIsoVals, EgammaCutBasedEleId::WorkingPoint::MEDIUM);  
  std::vector<bool> tagOfflineSelection   = selectElectrons(eleH, std::string(""), vertexH, conversionsH, rho, bs, electronIsoVals, EgammaCutBasedEleId::WorkingPoint::TIGHT);

  std::vector<std::pair<reco::GsfElectron, reco::GsfElectron>> offlinePairs;
  for (unsigned int i=0; i<probeOfflineSelection.size(); i++) {
    if (!probeOfflineSelection[i])
      continue;
    for (unsigned int j=0; j<tagOfflineSelection.size(); j++) {
      if (i==j)
	continue;

      if (!tagOfflineSelection[j])
	continue;

      offlinePairs.push_back(std::make_pair(electrons[j], electrons[i]));
    }
  }

  //std::cout << "OFFLINE" << std::endl;
  //for(unsigned int i=0; i<offlinePairs.size(); i++) 
  //  std::cout << offlinePairs[i].first.pt() << " " << offlinePairs[i].second.pt() << std::endl;

  // Select HLT probes
  trigger::TriggerObjectCollection allTriggerObjects = trigEvent->getObjects();

  edm::InputTag moduleLabel("hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8PMMassFilter", "", "HLT");
  trigger::TriggerObjectCollection hltTagAndProbes = selectTriggerObjects(allTriggerObjects, *trigEvent, moduleLabel, hltCuts);

  moduleLabel = edm::InputTag("hltEle17CaloIdVTCaloIsoVTTrkIdTTrkIsoVTEle8TrackIsoFilter", "", "HLT");
  trigger::TriggerObjectCollection hltTags = trigEvent->getObjects();

  auto hltPairs = removeTags(hltTagAndProbes, hltTags);
  
  //std::cout << "HLT" << std::endl;
  //for(unsigned int i=0; i<hltPairs.size(); i++) 
  //  std::cout << hltPairs[i].first.pt() << " " << hltPairs[i].second.pt() << std::endl;
  
  // Denominator
  auto offlineProbes = selectOfflineProbe(hltPairs, offlinePairs);
  
  //std::cout << "Denominator" << std::endl;
  //for (unsigned int i=0; i<offlineProbes.size(); i++) 
    //std::cout << offlineProbes[i].pt() << std::endl;

  // ProtoNumerator
  //std::cout << "_________" << std::endl;
  moduleLabel = edm::InputTag("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoDZ", "", "HLT");
  trigger::TriggerObjectCollection hltAllPairs = selectTriggerObjects(allTriggerObjects, *trigEvent, moduleLabel, hltCuts);
  //for (unsigned int i=0; i<hltAllPairs.size(); i+=2) 
  //  std::cout << hltAllPairs[i].pt() << " " << hltAllPairs[i+1].pt() << std::endl;

  //std::cout << "_________" << std::endl;
  moduleLabel = edm::InputTag("hltEle17TightIdLooseIsoEle8TightIdLooseIsoTrackIsoFilter", "", "HLT");
  trigger::TriggerObjectCollection hltLeadingLeg  = selectTriggerObjects(allTriggerObjects, *trigEvent, moduleLabel, hltCuts);
  //for (unsigned int i=0; i<hltLeadingLeg.size(); i++) 
  //  std::cout << hltLeadingLeg[i].pt() << std::endl;

  trigger::TriggerObjectCollection hltTrailingLeg = purgeTrailingLeg(hltAllPairs, hltLeadingLeg);
  //std::cout << "_________" << std::endl;
  //for (unsigned int i=0; i<hltTrailingLeg.size(); i++) 
  //  std::cout << hltTrailingLeg[i].pt() << std::endl;

  for (unsigned int i=0; i<offlineProbes.size(); i++) {
    bool leadingMatches  = matchByDeltaR(offlineProbes[i], hltLeadingLeg, 0.05);
    bool trailingMatches = matchByDeltaR(offlineProbes[i], hltTrailingLeg, 0.05);
    
    float e     = offlineProbes[i].superCluster()->energy();
    float eta   = offlineProbes[i].superCluster()->position().eta();
    float theta = offlineProbes[i].superCluster()->position().Theta();
    
    //std::cout << "OFFLINE PROBE " << e*sin(theta) << " " << eta << std::endl;
    
    etaden  = eta;
    ptden   = e*sin(theta);
    nvtxden = nVertices;
    
    eta1  = 999.;
    pt1   = 999.;
    nvtx1 = 999.;
    eta2  = 999.;
    pt2   = 999.;
    nvtx2 = 999.;
    
    if (leadingMatches) {
      eta1  = eta;
      pt1   = e*sin(theta);
      nvtx1 = nVertices;
    } 
    
    if (trailingMatches) {
      eta2  = eta;
      pt2   = e*sin(theta);
      nvtx2 = nVertices;
    } 
    
    tree->Fill();
  }
}

void  TurnOn::beginJob() {
  file = new TFile(fileName.c_str(), "recreate");
  tree = new TTree("tree", "tree");

  tree->Branch("pt1",    &pt1,    "pt1/F");
  tree->Branch("pt2",    &pt2,    "pt2/F");
  tree->Branch("ptden",  &ptden,  "ptden/F");
  tree->Branch("eta1",   &eta1,   "eta1/F");
  tree->Branch("eta2",   &eta2,   "eta2/F");
  tree->Branch("etaden", &etaden, "etaden/F");
  tree->Branch("nvtx1",   &nvtx1,   "nvtx1/F");
  tree->Branch("nvtx2",   &nvtx2,   "nvtx2/F");
  tree->Branch("nvtxden", &nvtxden, "nvtxden/F");

}

void TurnOn::endJob() {
  file->cd();
  tree->Write();
  file->Close();
}

void TurnOn::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(TurnOn);
