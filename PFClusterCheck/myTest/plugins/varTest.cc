// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "RecoEcal/EgammaCoreTools/interface/EcalClusterTools.h"

#include "DataFormats/ParticleFlowReco/interface/PFCluster.h"
#include "DataFormats/ParticleFlowReco/interface/PFClusterFwd.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHitFwd.h"
#include "DataFormats/ParticleFlowReco/interface/PFRecHit.h"


#include "Geometry/CaloGeometry/interface/CaloGeometry.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"
#include "Geometry/EcalAlgo/interface/EcalBarrelGeometry.h"
#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"

#include "TFile.h"
#include "TTree.h"

#include "DataFormats/Math/interface/deltaR.h"

class varTest : public edm::EDAnalyzer {
public:
  explicit varTest(const edm::ParameterSet&);
  //explicit varTest(const edm::ParameterSet&,  edm::ConsumesCollector &&cc);
  ~varTest();
  
  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  int findPho(reco::PhotonRef ref, edm::Handle<reco::GenParticleCollection> gpH);
  
private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  TFile* file;
  TTree* tree;
  Float_t ele;
  Float_t eleta;
  Float_t elphi;
  Float_t elecorr;
  Float_t e1x3;
  Float_t e3x3;
  Float_t e2x5max;
  Float_t e2x2;
  Float_t e5x5;
  Int_t iseb, size;

  Int_t truePU;

  // edm::EDGetTokenT<EcalRecHitCollection> _recHitsEB;
  //edm::EDGetTokenT<EcalRecHitCollection> _recHitsEE;  
  edm::InputTag _recHitsEB;
  edm::InputTag _recHitsEE;  
  std::string outfilename;
  edm::InputTag eleLabel;
};

varTest::varTest(const edm::ParameterSet& iConfig) {

  outfilename = iConfig.getParameter<std::string>("OutputFileName");
  eleLabel = iConfig.getParameter<edm::InputTag>("label");

  //_recHitsEB = cc.consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("recHitsEBLabel"));
  //_recHitsEE = cc.consumes<EcalRecHitCollection>(iConfig.getParameter<edm::InputTag>("recHitsEELabel"));
  _recHitsEB = iConfig.getParameter<edm::InputTag>("recHitsEBLabel");
  _recHitsEE = iConfig.getParameter<edm::InputTag>("recHitsEELabel");

  file = new TFile(outfilename.c_str(), "recreate");
  tree = new TTree("tree", "tree");
  tree->Branch("ele", &ele, "ele/F");
  tree->Branch("eleta", &eleta, "eleta/F");
  tree->Branch("elphi", &elphi, "elephi/F");
  tree->Branch("size", &size, "size/I");
  tree->Branch("iseb", &iseb, "iseb/I");
  tree->Branch("e1x3", &e1x3, "e1x3/F");
  tree->Branch("e2x2", &e2x2, "e2x2/F");
  tree->Branch("e2x5max", &e2x5max, "e2x5max/F");
  tree->Branch("e3x3", &e3x3, "e3x3/F");
  tree->Branch("e5x5", &e5x5, "e5x5/F");
  tree->Branch("elecorr", &elecorr, "elecorr/F");
}


varTest::~varTest() {
  file->cd();
  tree->Write();
  file->Close();
}
 
int varTest::findPho(reco::PhotonRef ref, edm::Handle<reco::GenParticleCollection> gpH) {
  
  int index = -1;
  float dRMin = 0.3;
  for (unsigned int i=0; i<gpH->size(); i++) {
    reco::GenParticleRef gp(gpH, i);
    if (gp->pt() > 4.) {
      //if (abs(gp->pdgId()) == 11) {
      if (abs(gp->pdgId()) == 22) {
	//std::cout << gp->et() << " " << gp->status() << std::endl;
	//if (gp->status() == 1 or gp->status() == 23) 
	{
	  float dr = deltaR(*ref, *gp);
	  if (dr < dRMin) {
	    dr = dRMin;
	    index = i;
	  }
	}
      }
    }
  }

  return index;
}

void varTest::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {

  const CaloTopology *topology;
  edm::ESHandle<CaloTopology> pTopology;
  iSetup.get<CaloTopologyRecord>().get(pTopology);
  topology = pTopology.product();

  //edm::Handle<reco::GenParticleCollection> gpH;
  //iEvent.getByLabel("genParticles", gpH);
  
  edm::Handle<std::vector<PileupSummaryInfo>> puH;
  iEvent.getByLabel("addPileupInfo", puH);
  truePU = (*puH)[0].getTrueNumInteractions();
  
  edm::Handle<reco::PhotonCollection> pEleH;
  iEvent.getByLabel(eleLabel, pEleH);
    
  for (size_t i=0; i<pEleH->size(); i++) {
      
    reco::PhotonRef el(pEleH, i);

    //elrawe = el->superCluster()->rawEnergy();
    for (reco::CaloCluster_iterator bit = el->superCluster()->clustersBegin(); bit!=el->superCluster()->clustersEnd(); ++bit) {
      const reco::CaloClusterPtr bc = *bit;
      
      DetId id = bc->hitsAndFractions()[0].first;
      bool isBarrel=(id.subdetId() == EcalBarrel);
      
      // Rech-Hits related
      edm::Handle<EcalRecHitCollection> prechits;
      iEvent.getByLabel((isBarrel ? _recHitsEB : _recHitsEE) ,prechits);
    
      ele = bc->energy();
      elecorr = bc->correctedEnergy();
      size = 0;//lazyTool.n5x5(*bc);
      eleta = bc->eta();
      elphi = bc->phi();
      e1x3    = EcalClusterTools::e1x3(*bc, &(*prechits), topology);
      e2x2    = EcalClusterTools::e2x2(*bc, &(*prechits), topology);
      e2x5max = EcalClusterTools::e2x5Max(*bc, &(*prechits), topology);
      e3x3    = EcalClusterTools::e3x3(*bc, &(*prechits), topology);
      e5x5    = EcalClusterTools::e5x5(*bc, &(*prechits), topology);

      tree->Fill();
    }
  }
}

void varTest::beginJob() {}

void varTest::endJob() {}

void varTest::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(varTest);
