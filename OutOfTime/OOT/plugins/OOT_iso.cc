// -*- C++ -*-
//
// Package:    Test/OOT
// Class:      OOT
// 
/**\class OOT OOT.cc Test/OOT/plugins/OOT.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matteo Sani
//         Created:  Tue, 08 Apr 2014 19:49:02 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "TFile.h"
#include "TTree.h"

class OOT_iso : public edm::EDAnalyzer {
public:
  explicit OOT_iso(const edm::ParameterSet&);
  ~OOT_iso();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  float pu_n_true[16];

  TTree* tree;
  TFile* file;
};


OOT_iso::OOT_iso(const edm::ParameterSet& iConfig) {

}


OOT_iso::~OOT_iso() {
 
}

void OOT_iso::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<std::vector< PileupSummaryInfo> > PupInfo;
  iEvent.getByLabel("addPileupInfo", PupInfo);
  
  std::vector<PileupSummaryInfo>::const_iterator PVI;
  int n = 0;
  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
    pu_n_true[n] = PVI->getTrueNumInteractions();
    n++;
  }

  edm::Handle<double> rhoH;
  iEvent.getByLabel("fixedGridRhoFastjetAllCalo", rhoH);
  if (rhoH.failedToGet())
    return;
  const double* rhoV = rhoH.product();
  
  edm::Handle<std::vector<reco::Vertex>> vtxH;
  iEvent.getByLabel("hltPixelVertices", vtxH);
  if (vtxH.failedToGet())
    return;
  const std::vector<reco::Vertex>* vtxC = vtxH.product();

  edm::Handle<LumiScalersCollection> lumiH;;
  iEvent.getByLabel("hltScalersRawToDigi", lumiH);
  const LumiScalersCollection* lumiC = lumiH.product();

  rho = *rhoV;
  lumi = lumiC->front().instantLumi();
  vtx = vtxC->size();
  pu = lumiC->front().pileup();

  edm::Handle<reco::GsfElectronCollection> elH;
  iEvent.getByLabel("gsfElectron", elH);

  for (unsigned int i=0; i<elH->size(); i++) {
    reco::GsfElectronRef ref(elH, i);
    iso[n] = ref->dr03EcalRecHitSumEt();
    et[n] = ref->et();
    n++;
  }
  
  tree->Fill();
}

void  OOT_iso::beginJob() {

  file = new TFile("pu_test.root", "recreate");
  tree = new TTree("tree", "tree");
  tree->Branch("pu_n", &pu_n, "pu_n[16]/I");
  tree->Branch("pu_n_true", &pu_n_true, "pu_n_true[16]/F");
  tree->Branch("pu_bx", &pu_bunchcrossing, "pu_bx[16]/I");
}

void OOT_iso::endJob() {
  tree->Write();
  file->Close();
}

void OOT_iso::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(OOT_iso);
