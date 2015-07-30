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

namespace {
  typedef reco::PFCluster::EEtoPSAssociation::value_type EEPSPair;
  bool sortByKey(const EEPSPair& a, const EEPSPair& b) {
    return a.first < b.first;
  } 
}

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
  Float_t nozsele;
  Float_t eleta;
  Float_t elphi;
  Float_t elecorr;
  Float_t e1x3;
  Float_t e3x3;
  Float_t e2x5max;
  Float_t e2x2;
  Float_t e5x5;
  Float_t e2nd;
  Float_t emax;
  Float_t pre1;
  Float_t pre2;
  Float_t nozse1x3;
  Float_t nozse3x3;
  Float_t nozse2x5max;
  Float_t nozse2x2;
  Float_t nozse5x5;
  Float_t nozse2nd;
  Float_t nozsemax;
  Int_t iseb, size;

  Int_t truePU;
  edm::InputTag eetopsSrc_;
  edm::InputTag _recHitsEB;
  edm::InputTag _recHitsEE;  
  std::string outfilename;
  edm::InputTag eleLabel;
};

varTest::varTest(const edm::ParameterSet& iConfig) {

  outfilename = iConfig.getParameter<std::string>("OutputFileName");
  eleLabel = iConfig.getParameter<edm::InputTag>("label");
  eetopsSrc_ = iConfig.getParameter<edm::InputTag>("eetops");

  _recHitsEB = iConfig.getParameter<edm::InputTag>("recHitsEBLabel");
  _recHitsEE = iConfig.getParameter<edm::InputTag>("recHitsEELabel");

  file = new TFile(outfilename.c_str(), "recreate");
  tree = new TTree("tree", "tree");
  tree->Branch("ele", &ele, "ele/F");
  tree->Branch("nozsele", &nozsele, "nozsele/F");
  tree->Branch("eleta", &eleta, "eleta/F");
  tree->Branch("elphi", &elphi, "elephi/F");
  tree->Branch("size", &size, "size/I");
  tree->Branch("iseb", &iseb, "iseb/I");
  tree->Branch("e1x3", &e1x3, "e1x3/F");
  tree->Branch("e2x2", &e2x2, "e2x2/F");
  tree->Branch("e2x5max", &e2x5max, "e2x5max/F");
  tree->Branch("e3x3", &e3x3, "e3x3/F");
  tree->Branch("e5x5", &e5x5, "e5x5/F");
  tree->Branch("emax", &emax, "emax/F");
  tree->Branch("e2nd", &e2nd, "e2nd/F");
  tree->Branch("pre1", &pre1, "pre1/F");
  tree->Branch("pre2", &pre2, "pre2/F");
  tree->Branch("nozse1x3", &nozse1x3, "nozse1x3/F");
  tree->Branch("nozse2x2", &nozse2x2, "nozse2x2/F");
  tree->Branch("nozse2x5max", &nozse2x5max, "nozse2x5max/F");
  tree->Branch("nozse3x3", &nozse3x3, "nozse3x3/F");
  tree->Branch("nozse5x5", &nozse5x5, "nozse5x5/F");
  tree->Branch("nozsemax", &nozsemax, "nozsemax/F");
  tree->Branch("nozse2nd", &nozse2nd, "nozse2nd/F");
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


  edm::Handle<reco::PFCluster::EEtoPSAssociation> eetops;
  iEvent.getByLabel(eetopsSrc_,eetops);
  const reco::PFCluster::EEtoPSAssociation assoc = *(eetops.product());

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
      
      nozsele = 0;
      // Rech-Hits related
      edm::Handle<EcalRecHitCollection> prechits;
      iEvent.getByLabel((isBarrel ? _recHitsEB : _recHitsEE) ,prechits);
      const EcalRecHitCollection* cH = prechits.product();

      const std::vector< std::pair<DetId, float>> v_id = bc->hitsAndFractions();
      std::vector< std::pair<DetId, float>>::const_iterator it;
      //for (it = v_id.begin(); it != v_id.end(); ++it) {
      for (unsigned int i=0; i<v_id.size(); i++) {
	//energy2 += (EcalClusterTools::recHitEnergy(v_id[i].first, cH))*v_id[i].second;
	nozsele += noZS::EcalClusterTools::recHitEnergy(v_id[i].first, cH);
	//std::cout << "ID: " << cH->find(v_id[i].first)->energy() << " " << v_id[i].second << std::endl;
        //energy += cH->find(v_id[i].first)->energy();
      }

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
      e2nd    = EcalClusterTools::e2nd(*bc, &(*prechits));//, topology);
      emax    = EcalClusterTools::eMax(*bc, &(*prechits));//, topology);

      pre1 = 0;
      pre2 = 0;
      if(!isBarrel) {
        auto ee_key_val = std::make_pair(bc.key(), edm::Ptr<reco::PFCluster>());
        const auto clustops = std::equal_range(assoc.begin(),
					       assoc.end(),
					       ee_key_val,
					       sortByKey);
        for( auto i_ps = clustops.first; i_ps != clustops.second; ++i_ps) {
	  edm::Ptr<reco::PFCluster> psclus(i_ps->second);
          switch( psclus->layer() ) {
          case PFLayer::PS1:
            pre1 += psclus->energy();
            break;
          case PFLayer::PS2:
            pre2 += psclus->energy();
            break;
          default:
            break;
          }
        }
      }


      nozse1x3    = noZS::EcalClusterTools::e1x3(*bc, &(*prechits), topology);
      nozse2x2    = noZS::EcalClusterTools::e2x2(*bc, &(*prechits), topology);
      nozse2x5max = noZS::EcalClusterTools::e2x5Max(*bc, &(*prechits), topology);
      nozse3x3    = noZS::EcalClusterTools::e3x3(*bc, cH, topology);
      nozse5x5    = noZS::EcalClusterTools::e5x5(*bc, cH, topology);
      nozse2nd    = noZS::EcalClusterTools::e2nd(*bc, cH);//, topology);
      nozsemax    = noZS::EcalClusterTools::eMax(*bc, cH);//, topology);
      //std::cout << "--------" << std::endl;
      //std::cout << e5x5 << " " << e5x5/ele << " " << ele << " " << nozse5x5 << " " << nozse5x5/nozsele << " " << nozsele << std::endl;

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
