#include "NtupleMakers/ElectronAnalyzer/interface/NewElectronsC.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "CLHEP/HepMC/GenParticle.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#include "TFile.h"
#include "TTree.h"

using namespace std;
using namespace edm;
using namespace reco;

/// Constructor
NewElectronsC::NewElectronsC(const ParameterSet& pset) {
  fileName = pset.getParameter<std::string>("RootFileName");
  baselineEleCollName =  pset.getParameter<std::string>("BaselineEleCollName");
  customEleCollName   =  pset.getParameter<std::string>("CustomEleCollName");
}

NewElectronsC::~NewElectronsC() {}

void NewElectronsC::beginJob(const EventSetup& eventSetup) {

  file = new TFile(fileName.c_str(), "recreate");
  tree = new TTree("event","Event data");

  //tree->Branch("run", &run, "run/I");
  //tree->Branch("id", &id, "id/I");

  tree->Branch("mc_pt", &mc_pt, "mc_pt/F");
  tree->Branch("mc_eta", &mc_eta, "mc_eta/F");
  tree->Branch("mc_phi", &mc_phi, "mc_phi/F");
  tree->Branch("mc_id", &mc_id, "mc_id/I");
  tree->Branch("mc_e", &mc_e, "mc_e/F");
  tree->Branch("mc_mother", &mc_mother, "mc_mother/I");
  tree->Branch("mc_crack", &mc_crack, "mc_crack/I"); // 0 no crack, 1 crack
  tree->Branch("mc_dr", &mc_dr, "mc_dr/F");

  tree->Branch("mc1_pt", &mc1_pt, "mc1_pt/F");
  tree->Branch("mc1_eta", &mc1_eta, "mc1_eta/F");
  tree->Branch("mc1_phi", &mc1_phi, "mc1_phi/F");
  tree->Branch("mc1_id", &mc1_id, "mc1_id/I");
  tree->Branch("mc1_e", &mc1_e, "mc1_e/F");
  tree->Branch("mc1_mother", &mc1_mother, "mc1_mother/I");
  tree->Branch("mc1_crack", &mc1_crack, "mc1_crack/I"); // 0 no crack, 1 crack
  tree->Branch("mc1_dr", &mc1_dr, "mc1_dr/F");
  /*
  tree->Branch("sc_e", &sc_e, "sc_e/F");
  tree->Branch("sc_rawe", &sc_rawe, "sc_rawe/F");
  tree->Branch("sc_et", &sc_et, "sc_et/F");
  tree->Branch("sc_eta", &sc_eta, "sc_eta/F");
  tree->Branch("sc_phi", &sc_phi, "sc_phi/F");
  tree->Branch("sc_dr", &sc_dr, "sc_dr/F");
  tree->Branch("sc_type", &sc_type, "sc_type/I"); // 0 barrel, 1 endcap

  tree->Branch("tk_pt", &tk_pt, "tk_e/F");
  tree->Branch("tk_eta", &tk_eta, "tk_eta/F");
  tree->Branch("tk_phi", &tk_phi, "tk_phi/F");
  tree->Branch("tk_dr", &tk_dr, "tk_dr/F");
  tree->Branch("tk_nhit", &tk_nhit, "tk_nhit/I");
  //tree->Branch("tk_layer", &tk_layer, "tk_layer/I");
  //tree->Branch("tk_subdet", &tk_subdet, "tk_subdet/I");
  */
  /*
  tree->Branch("el_pt", &el_pt, "el_pt/F");
  tree->Branch("el_e", &el_e, "el_e/F");
  tree->Branch("el_eta", &el_eta, "el_eta/F");
  tree->Branch("el_phi", &el_phi, "el_phi/F");
  tree->Branch("el_eopin", &el_eopin, "el_eopin/F");
  tree->Branch("el_eopout", &el_eopout, "el_eopout/F");
  tree->Branch("el_pout", &el_pout, "el_pout/F");
  tree->Branch("el_fbrem", &el_fbrem, "el_fbrem/F");
  tree->Branch("el_hoe", &el_hoe, "el_hoe/F");
  tree->Branch("el_detain", &el_detain, "el_detain/F");
  tree->Branch("el_dphiin", &el_dphiin, "el_dphiin/F");
  tree->Branch("el_detaout", &el_detaout, "el_detaout/F");
  tree->Branch("el_dphiout", &el_dphiout, "el_dphiout/F");
  tree->Branch("el_e3x3", &el_e3x3, "el_e3x3/F");
  tree->Branch("el_e5x5", &el_e5x5, "el_e5x5/F");
  tree->Branch("el_eseed", &el_eseed, "el_eseed/F");
  tree->Branch("el_spp", &el_spp, "el_spp/F");
  tree->Branch("el_see", &el_see, "el_see/F");
  tree->Branch("el_class", &el_class, "el_class/I");
  tree->Branch("el_nsihit", &el_nsihits, "el_nsihit/I");
  tree->Branch("el_npxhit", &el_npxhits, "el_npxhit/I"); 
  tree->Branch("el_detinnerhit", &el_detinnerhit, "el_detinnerhit/I");
  tree->Branch("el_rinnerhit", &el_rinnerhit, "el_rinnerhit/F");
  tree->Branch("el_z0", &el_z0, "el_z0/F");
  tree->Branch("el_tkiso", &el_tkiso, "el_tkiso/F");
  */
  tree->Branch("el1_pt", &el1_pt, "el1_pt/F");
  tree->Branch("el1_e", &el1_e, "el1_e/F");
  tree->Branch("el1_eta", &el1_eta, "el1_eta/F");
  tree->Branch("el1_phi", &el1_phi, "el1_phi/F");
  tree->Branch("el1_eopin", &el1_eopin, "el1_eopin/F");
  tree->Branch("el1_eopout", &el1_eopout, "el1_eopout/F");
  tree->Branch("el1_pout", &el1_pout, "el1_pout/F");
  tree->Branch("el1_fbrem", &el1_fbrem, "el1_fbrem/F");
  tree->Branch("el1_hoe", &el1_hoe, "el1_hoe/F");
  tree->Branch("el1_detain", &el1_detain, "el1_detain/F");
  tree->Branch("el1_dphiin", &el1_dphiin, "el1_dphiin/F");
  tree->Branch("el1_detaout", &el1_detaout, "el1_detaout/F");
  tree->Branch("el1_dphiout", &el1_dphiout, "el1_dphiout/F");
  tree->Branch("el1_e3x3", &el1_e3x3, "el1_e3x3/F");
  tree->Branch("el1_e5x5", &el1_e5x5, "el1_e5x5/F");
  tree->Branch("el1_eseed", &el1_eseed, "el1_eseed/F");
  tree->Branch("el1_spp", &el1_spp, "el1_spp/F");
  tree->Branch("el1_see", &el1_see, "el1_see/F");
  tree->Branch("el1_class", &el1_class, "el1_class/I");
  tree->Branch("el1_nsihit", &el1_nsihits, "el1_nsihit/I");
  tree->Branch("el1_npxhit", &el1_npxhits, "el1_npxhit/I"); 
  tree->Branch("el1_detinnerhit", &el1_detinnerhit, "el1_detinnerhit/I");
  tree->Branch("el1_rinnerhit", &el1_rinnerhit, "el1_rinnerhit/F");
  tree->Branch("el1_z0", &el1_z0, "el1_z0/F");
  tree->Branch("el1_tkiso", &el1_tkiso, "el1_tkiso/F");
}

void NewElectronsC::endJob() {

  file->Write();
  file->Close();
}

void NewElectronsC::analyze(const Event & event, const EventSetup& eventSetup) {

  cout << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;
 
  // access the tracker
  edm::ESHandle<TrackerGeometry> theTrackerGeometry;
  eventSetup.get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);
  const TrackerGeometry& theTracker(*theTrackerGeometry);

  Handle<HepMCProduct> evt;
  event.getByLabel("source", evt);
  myGenEvent = new HepMC::GenEvent(*(evt->GetEvent()));

  Handle<TrackCollection> tkh;
  event.getByLabel("ctfWithMaterialTracks", tkh);
  const TrackCollection* tracks = tkh.product();
  TrackCollection::const_iterator itt;

  Handle<PixelMatchGsfElectronCollection> elh;
  event.getByLabel(baselineEleCollName, elh);
  const PixelMatchGsfElectronCollection* copy = elh.product();
  PixelMatchGsfElectronCollection::const_iterator ite; 

  //Handle<SuperClusterCollection> sch1;
  //event.getByLabel("hybridSuperClusters", sch1);
  //const SuperClusterCollection* scb = sch1.product();

  //Handle<SuperClusterCollection> sch2;
  //event.getByLabel("islandSuperClusters", "islandEndcapSuperClusters", sch2);
  //const SuperClusterCollection* sce = sch2.product();
  //SuperClusterCollection::const_iterator itscb, itsce;

  Handle<PixelMatchGsfElectronCollection> elh1;
  event.getByLabel(customEleCollName, elh1);
  const PixelMatchGsfElectronCollection*  electrons1 = elh1.product();
  PixelMatchGsfElectronCollection::const_iterator ite1;

  for(ite1 = electrons1->begin(); ite1 != electrons1->end(); ++ite1) {
    el1_pt = ite1->pt(); 
    el1_eta = ite1->eta(); 
    el1_e = ite1->energy();
    el1_phi = ite1->phi(); 
    el1_eopin = ite1->eSuperClusterOverP();
    el1_eopout = ite1->eSeedClusterOverPout();
    el1_hoe = ite1->hadronicOverEm();
    el1_dphiin = ite1->deltaPhiSuperClusterTrackAtVtx();
    el1_detain = ite1->deltaEtaSuperClusterTrackAtVtx();
    el1_dphiout = ite1->deltaPhiSeedClusterTrackAtCalo();
    el1_detaout = ite1->deltaEtaSeedClusterTrackAtCalo();
    float pin  = ite1->trackMomentumAtVtx().R();
    float pout = ite1->trackMomentumOut().R();
    el1_pout = pout;
    el1_fbrem = (pin-pout)/pin;
    el1_class = ite1->classification();
    R9_25_gsf(event, &(*ite1), el1_eseed, el1_e3x3, el1_e5x5, el1_spp, el1_see);
    int a, b;
    nHits(ite1->gsfTrack(), a, b);
    el1_npxhits = a;
    el1_nsihits = b;

    int index = 1;
    while(1) {
      TrackingRecHitRef hit = ite1->gsfTrack()->recHit(ite1->gsfTrack()->recHitsSize()-index);
      
      if (hit->isValid()) {
        GlobalPoint hitPosition = theTracker.idToDet(hit->geographicalId())->surface().toGlobal(hit->localPosition());
        GlobalPoint pos(hitPosition.x()-ite1->gsfTrack()->vx(), hitPosition.y()-ite1->gsfTrack()->vy(), hitPosition.z()-ite1->gsfTrack()->vz());
        el1_rinnerhit = sqrt(pow(pos.perp(),2) + pow(pos.z(),2));
        subDetector(hit, a, b);
        el1_detinnerhit = a;
        break;
      }
      index++;
    }
    
    el1_z0 = ite1->gsfTrack()->vz();
    el1_tkiso = trackIsolation(ite1->trackMomentumAtVtx(), ite1->vertex(), tracks);
    
    int type = 0;
    float theta = 1.;
    double dR, dRmin = 0.05;
    HepMC::GenEvent::particle_const_iterator nearMC;
    for (HepMC::GenEvent::particle_const_iterator itmc1 = myGenEvent->particles_begin(); itmc1 != myGenEvent->particles_end(); ++itmc1) { 
      
      if ((abs((*itmc1)->pdg_id()) == 11) && ((*itmc1)->status() != 3)) {      
        if (((*itmc1)->momentum().perp() > 5.) && (fabs((*itmc1)->momentum().eta()) < 2.5)) {
          
          math::XYZVector mcv((*itmc1)->momentum().px(), (*itmc1)->momentum().py(), (*itmc1)->momentum().pz());          
          
          dR = ROOT::Math::VectorUtil::DeltaR(ite1->p4(), mcv);
          if (dR < dRmin) {
            dRmin = dR;
            nearMC = itmc1;
          }
          
        }
      }
    }
        
    if (dRmin < 0.05) {
      mc1_dr = dRmin;
      mc1_mother = mother(*nearMC);
      mc1_pt = (*nearMC)->momentum().perp();
      mc1_eta = (*nearMC)->momentum().eta();
      mc1_phi = (*nearMC)->momentum().phi();
      mc1_e = (*nearMC)->momentum().e();
      mc1_id = (*nearMC)->pdg_id();
      
      // check if it is in a crack
      if (inCrack(fabs((*nearMC)->momentum().eta())))
        mc1_crack = 1;
      else
        mc1_crack = 0;
      
    } else {
      mc1_dr = 0.05;
      mc1_mother = 0;
      mc1_pt = 0;
      mc1_eta = 0;
      mc1_phi = 0;
      mc1_e = 0;
      mc1_id = 0;
      mc1_crack = -1;

    }
    tree->Fill();
  }
  /*
  // remove duplicate electrons
  PixelMatchGsfElectronCollection electrons;
  PixelMatchGsfElectronCollection::const_iterator it1, it2;
  
  for(it1=copy->begin(); it1!=copy->end(); ++it1) {
    
    bool isRemoved = false;
    for(it2=copy->begin(); it2!=copy->end(); ++it2) {
      if (it1 == it2)
        continue;
      if (((*it1).superCluster().id() == (*it2).superCluster().id()) &&
          ((*it1).superCluster().index() == (*it2).superCluster().index())) {
        
        float deltaEp1 = fabs((*it1).eSuperClusterOverP() - 1.);
        float deltaEp2 = fabs((*it2).eSuperClusterOverP() - 1.);
        if (deltaEp1 > deltaEp2) {
          isRemoved = true;
          break;
        }
      }
    }
    
    if (!isRemoved)
      electrons.push_back(*it1);
  }
  
  // new electrons collection
  for(ite = electrons.begin(); ite != electrons.end(); ++ite) {
    el_pt = ite->pt(); 
    el_eta = ite->eta(); 
    el_e = ite->energy();
    el_phi = ite->phi(); 
    el_eopin = ite->eSuperClusterOverP();
    el_eopout = ite->eSeedClusterOverPout();
    el_hoe = ite->hadronicOverEm();
    el_dphiin = ite->deltaPhiSuperClusterTrackAtVtx();
    el_detain = ite->deltaEtaSuperClusterTrackAtVtx();
    el_dphiout = ite->deltaPhiSeedClusterTrackAtCalo();
    el_detaout = ite->deltaEtaSeedClusterTrackAtCalo();
    float pin  = ite->trackMomentumAtVtx().R();
    float pout = ite->trackMomentumOut().R();
    el_pout = pout;
    el_fbrem = (pin-pout)/pin;
    el_class = ite->classification();
    R9_25_gsf(event, &(*ite), el_eseed, el_e3x3, el_e5x5, el_spp, el_see);
    int a, b;
    nHits(ite->gsfTrack(), a, b);
    el_npxhits = a;
    el_nsihits = b;

    int index = 1;
    while(1) {
      TrackingRecHitRef hit = ite->gsfTrack()->recHit(ite->gsfTrack()->recHitsSize()-index);
      
      if (hit->isValid()) {
        GlobalPoint hitPosition = theTracker.idToDet(hit->geographicalId())->surface().toGlobal(hit->localPosition());
        GlobalPoint pos(hitPosition.x()-ite->gsfTrack()->vx(), hitPosition.y()-ite->gsfTrack()->vy(), hitPosition.z()-ite->gsfTrack()->vz());
        el_rinnerhit = sqrt(pow(pos.perp(),2) + pow(pos.z(),2));
        subDetector(hit, a, b);
        el_detinnerhit = a;
        break;
      }
      index++;
    }
    
    el_z0 = ite->gsfTrack()->vz();
    el_tkiso = trackIsolation(ite->trackMomentumAtVtx(), ite->vertex(), tracks);
    int type = 0;
    float theta = 1.;
    double dR, dRmin = 0.05;
    HepMC::GenEvent::particle_const_iterator nearMC1;
    for (HepMC::GenEvent::particle_const_iterator itmc = myGenEvent->particles_begin(); itmc != myGenEvent->particles_end(); ++itmc) { 
      
      if ((abs((*itmc)->pdg_id()) == 11) && ((*itmc)->status() != 3)) {      
        if (((*itmc)->momentum().perp() > 5.) && (fabs((*itmc)->momentum().eta()) < 2.5)) {
          
          math::XYZVector mcv((*itmc)->momentum().px(), (*itmc)->momentum().py(), (*itmc)->momentum().pz());          
          
          dR = ROOT::Math::VectorUtil::DeltaR(ite->p4(), mcv);
          if (dR < dRmin) {
            dRmin = dR;
            nearMC1 = itmc;
          }
          
        }
      }
    }
        
    if (dRmin < 0.05) {
      mc_dr = dRmin;
      mc_mother = mother(*nearMC1);
      mc_pt = (*nearMC1)->momentum().perp();
      mc_eta = (*nearMC1)->momentum().eta();
      mc_phi = (*nearMC1)->momentum().phi();
      mc_e = (*nearMC1)->momentum().e();
      mc_id = (*nearMC1)->pdg_id();
      
      // check if it is in a crack
      if (inCrack(fabs((*nearMC1)->momentum().eta())))
        mc_crack = 1;
      else
        mc_crack = 0;
      
    } else {
      mc_dr = 0.05;
      mc_mother = 0;
      mc_pt = 0;
      mc_eta = 0;
      mc_phi = 0;
      mc_e = 0;
      mc_id = 0;
      mc_crack = -1;

    }
  }
  */
}

bool NewElectronsC::inCrack(float eta) {

  return (eta < 0.018 ||
          (eta>0.423 && eta<0.461) ||
          (eta>0.770 && eta<0.806) ||
          (eta>1.127 && eta<1.163) ||
          (eta>1.460 && eta<1.558));
}

int NewElectronsC::mother(HepMC::GenParticle *p) {
  
  while (p->production_vertex()) {
    HepMC::GenVertex* inVertex = p->production_vertex();
    for(std::set<HepMC::GenParticle*>::const_iterator iter = inVertex->particles_in_const_begin();
        iter != inVertex->particles_in_const_end();iter++) {
      if ((*iter)->pdg_id() != p->pdg_id()) {
        return (*iter)->pdg_id();
      } else {
        p = *iter;
        break;
      }
    }
  }
  
  return -1;
}

void NewElectronsC::R9_25_gsf(const Event & event, const reco::PixelMatchGsfElectron* e,
                            float& eseed, float& e3x3, float& e5x5, float& spp, float& see) {
  
  reco::SuperClusterRef sclRef=e->superCluster();

  edm::Handle<reco::BasicClusterShapeAssociationCollection> bH, eH;
  event.getByLabel("hybridSuperClusters", "hybridShapeAssoc", bH);
  const reco::BasicClusterShapeAssociationCollection* barrelClShp = bH.product();
  event.getByLabel("islandBasicClusters", "islandEndcapShapeAssoc", eH);
  const reco::BasicClusterShapeAssociationCollection* endcapClShp = eH.product();

  reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr;
  DetId id = sclRef->seed()->getHitsByDetId()[0];
  if (id.subdetId() == EcalBarrel) {
    seedShpItr = barrelClShp->find(sclRef->seed());
  } else {
    seedShpItr = endcapClShp->find(sclRef->seed());
  }

  // Get the ClusterShapeRef corresponding to the BasicCluster
  const reco::ClusterShapeRef& seedShapeRef = seedShpItr->val;

  eseed = sclRef->seed()->energy();
  e3x3 = seedShapeRef->e3x3();
  e5x5 = seedShapeRef->e5x5();
  spp = sqrt(seedShapeRef->covPhiPhi());
  see = sqrt(seedShapeRef->covEtaEta());
}

void NewElectronsC::nHits(const reco::GsfTrackRef t, int& nPixelHits, int& nSiTkHits) {

  // loop sugli hits e conta il risultato facile no ?
  nPixelHits = 0; 
  nSiTkHits = 0;

  for(size_t i = 0; i < t->recHitsSize(); ++i) {

    TrackingRecHitRef hit = t->recHit(i);
    
    if (hit->isValid()) {
      DetId detid(hit->geographicalId());       
      unsigned int subdetId = static_cast<unsigned int>(detid.subdetId());
      if ((subdetId > 2) && (subdetId < 7))
        nSiTkHits++;
      if ((subdetId == 2) || (subdetId == 1))
        nPixelHits++;

    }
  }
}

//trackRelIsolation(el->trackMomentumAtVtx(), el->vertex(), tracks, 0.3, 0.01, 0.1, 999.9, 0.5, 1.5, 7);
double NewElectronsC::trackIsolation(const math::XYZVector momentum, 
                                    const math::XYZPoint vertex,
                                    const TrackCollection* tracks) {

  double dRConeMax = 0.3;
  double dRConeMin = 0.01;
  double tkVtxDMax = 0.1;
  double vtxDiffDMax = 999.9;
  double vtxDiffZMax = 0.5;
  double ptMin = 1.5;
  unsigned int nHits = 7;
  double isoResult = -10.;

  if ( tracks == 0 ) {
    return isoResult;
  }
  
  double sumPt = 0;
  
  std::vector<Track>::const_iterator iTk;
  for (iTk = tracks->begin(); iTk != tracks->end(); ++iTk){
    double dR = ROOT::Math::VectorUtil::DeltaR(momentum, iTk->momentum());
    //exclude tks in veto cone (set it to small number to 
    //exclude this track
    double dZ = fabs(vertex.z() - iTk->vz());
    double d0 = sqrt(iTk->vertex().perp2());
    double dD0 = sqrt((iTk->vertex() - vertex).perp2());
    
    if (dR < dRConeMin) 
      continue;
    
    if ( dR < dRConeMax 
         && dZ < vtxDiffZMax
         && d0 < tkVtxDMax
         && dD0 < vtxDiffDMax 
         && iTk->pt() >= ptMin
         && iTk->found() > nHits){
      sumPt += iTk->pt();
    }
  }
  
  isoResult = sumPt;

  return isoResult;
}

void NewElectronsC::subDetector(TrackingRecHitRef hit, int& subdet, int& layer) {

  DetId detid(hit->geographicalId());       
  unsigned int subdetId = static_cast<unsigned int>(detid.subdetId());
  switch (subdetId) {
  case 1:
    {
      PXBDetId thePXBDetId(detid.rawId());
      layer = thePXBDetId.layer();
      subdet = 1;
      break;
    }
  case 2:
    {
      PXFDetId thePXFDetId(detid.rawId());
      layer = thePXFDetId.disk();
      subdet = 2;
      break;
    }
  case StripSubdetector::TIB:
    {
      TIBDetId theTIBDetId(detid.rawId());
      layer = theTIBDetId.layer();
      subdet = 3;
      break;
    }
  case StripSubdetector::TID:
    {
      TIDDetId theTIDDetId(detid.rawId());
      layer = theTIDDetId.wheel();
      subdet = 4;
      break;
    }
  case StripSubdetector::TOB:
    {
      TOBDetId theTOBDetId(detid.rawId());
      layer = theTOBDetId.layer();
      subdet = 5;
      break;
    }
  case StripSubdetector::TEC:
    {
      TECDetId theTECDetId(detid.rawId());             
      layer = theTECDetId.wheel();
      subdet = 6;
      break;
    }
  }

}
