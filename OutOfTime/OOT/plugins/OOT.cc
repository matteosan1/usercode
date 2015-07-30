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

#include "TFile.h"
#include "TTree.h"

class OOT : public edm::EDAnalyzer {
public:
  explicit OOT(const edm::ParameterSet&);
  ~OOT();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  
  int pu_n[16];
  float pu_n_true[16];
  int pu_bunchcrossing[16];

  TTree* tree;
  TFile* file;
};


OOT::OOT(const edm::ParameterSet& iConfig) {

}


OOT::~OOT() {
 
}

void OOT::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  edm::Handle<std::vector< PileupSummaryInfo> > PupInfo;
  iEvent.getByLabel("addPileupInfo", PupInfo);
  
  std::vector<PileupSummaryInfo>::const_iterator PVI;

  //const int getPU_NumInteractions() const { return num_PU_vertices_; }
  //const std::vector<float>& getPU_instLumi() const { return instLumi_; }
  //const std::vector<edm::EventID>& getPU_EventID() const { return eventInfo_; }
  //const int getBunchCrossing() const { return bunchCrossing_;}
  //const float getTrueNumInteractions() const { return TrueNumInteractions_;}

  int n = 0;
  for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
    pu_bunchcrossing[n] = PVI->getBunchCrossing();
    pu_n[n] = PVI->getPU_NumInteractions();
    pu_n_true[n] = PVI->getTrueNumInteractions();
    //lumi = PVI->getPU_instLumi();
    n++;
  }
  
  tree->Fill();
}

void  OOT::beginJob() {

  file = new TFile("pu_test.root", "recreate");
  tree = new TTree("tree", "tree");
  tree->Branch("pu_n", &pu_n, "pu_n[16]/I");
  tree->Branch("pu_n_true", &pu_n_true, "pu_n_true[16]/F");
  tree->Branch("pu_bx", &pu_bunchcrossing, "pu_bx[16]/I");
}

void OOT::endJob() {
  tree->Write();
  file->Close();
}

void OOT::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
DEFINE_FWK_MODULE(OOT);
