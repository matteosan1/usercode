// -*- C++ -*-
//
// Package:    MiniAOD/Validation
// Class:      Validation
// 
/**\class Validation Validation.cc MiniAOD/Validation/plugins/Validation.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Matteo Sani
//         Created:  Thu, 30 Oct 2014 12:46:57 GMT
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

#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

#include "TFile.h"
#include "TTree.h"

#include <string>
#include <map>

class Validation : public edm::EDAnalyzer {
public:
  explicit Validation(const edm::ParameterSet&);
  ~Validation();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;


  edm::EDGetTokenT<edm::View<reco::GsfElectron> > candSrcToken_;
  TFile* file;
  TTree* tree;
  std::string outputFileName;
  Long_t run, event;
  Int_t n;
  Float_t pt[10], et[10], eta[10], sieie[10], phi[10];
  Float_t fbrem[10], deta[10], dphi[10], eop[10];
};

Validation::Validation(const edm::ParameterSet& iConfig) : 
  candSrcToken_(consumes<edm::View<reco::GsfElectron> >(iConfig.getParameter<edm::InputTag>("electrons"))) {

  outputFileName = iConfig.getParameter<std::string>("outputFileName");
}

Validation::~Validation() {
  
  file->cd();
  tree->Write();
  file->Close();
}

void Validation::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  run = iEvent.id().run();
  event = iEvent.id().event();

  edm::Handle<edm::View<reco::GsfElectron> > elH;
  iEvent.getByToken(candSrcToken_, elH);

  n = 0;
  for(size_t i=0; i<elH->size(); ++i) {
    edm::Ptr<reco::GsfElectron> e = elH->ptrAt(i);
    
    if (i==10)
      continue;
    
    pt[i] = e->pt();
    et[i] = e->superCluster()->energy()/cosh(e->superCluster()->eta());
    eta[i] = e->superCluster()->eta();
    phi[i] = e->superCluster()->phi();
    sieie[i] = e->full5x5_sigmaIetaIeta();
    fbrem[i] = e->fbrem();
    deta[i] = e->deltaEtaSuperClusterTrackAtVtx();
    dphi[i] = e->deltaPhiSuperClusterTrackAtVtx();
    eop[i] = e->eSuperClusterOverP();

    n++;
  }
  
  tree->Fill();
}


void Validation::beginJob() {
  file = new TFile(outputFileName.c_str(), "recreate");
  tree = new TTree("validation", "");
  tree->Branch("run",   &run,   "run/L");
  tree->Branch("event", &event, "event/L");
  tree->Branch("n",     &n,     "n/I");
  tree->Branch("pt",    &pt,    "pt[n]/F");
  tree->Branch("et",    &et,    "et[n]/F");
  tree->Branch("eta",   &eta,   "eta[n]/F");
  tree->Branch("phi",   &phi,   "phi[n]/F");
  tree->Branch("sieie", &sieie, "sieie[n]/F");
  tree->Branch("fbrem", &fbrem, "fbrem[n]/F");
  tree->Branch("deta",  &deta,  "deta[n]/F");
  tree->Branch("dphi",  &dphi,  "dphi[n]/F");
  tree->Branch("eop",   &eop,   "eop[n]/F");
}

void Validation::endJob() 
{}

void Validation::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(Validation);
