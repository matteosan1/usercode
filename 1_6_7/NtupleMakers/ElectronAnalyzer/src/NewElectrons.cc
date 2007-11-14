#include "NtupleMakers/ElectronAnalyzer/interface/NewElectrons.h"

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
NewElectrons::NewElectrons(const ParameterSet& pset) {
  fileName = pset.getParameter<std::string>("RootFileName");
  baselineEleCollName =  pset.getParameter<std::string>("BaselineEleCollName");
  customEleCollName   =  pset.getParameter<std::string>("CustomEleCollName");
}

NewElectrons::~NewElectrons() {}

void NewElectrons::beginJob(const EventSetup& eventSetup) {

  file = new TFile(fileName.c_str(), "recreate");
  tree = new TTree("event","Event data");

  tree->Branch("run", &run, "run/I");
  tree->Branch("id", &id, "id/I");

  //tree->Branch("mc_n", &nMC, "mc_n/I");
  tree->Branch("mc_pt", &mc_pt, "mc_pt/F");
  tree->Branch("mc_eta", &mc_eta, "mc_eta/F");
  tree->Branch("mc_phi", &mc_phi, "mc_phi/F");
  tree->Branch("mc_id", &mc_id, "mc_id/I");
  tree->Branch("mc_e", &mc_e, "mc_e/F");
  tree->Branch("mc_mother", &mc_mother, "mc_mother/I");
  tree->Branch("mc_crack", &mc_crack, "mc_crack/I"); // 0 no crack, 1 crack

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

  tree->Branch("el_pt", &el_pt, "el_pt/F");
  tree->Branch("el_e", &el_e, "el_e/F");
  tree->Branch("el_eta", &el_eta, "el_eta/F");
  tree->Branch("el_phi", &el_phi, "el_phi/F");
  tree->Branch("el_dr", &el_dr, "el_dr/F");
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

  tree->Branch("el1_pt", &el1_pt, "el1_pt/F");
  tree->Branch("el1_e", &el1_e, "el1_e/F");
  tree->Branch("el1_eta", &el1_eta, "el1_eta/F");
  tree->Branch("el1_phi", &el1_phi, "el1_phi/F");
  tree->Branch("el1_dr", &el1_dr, "el1_dr/F");
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

void NewElectrons::endJob() {

  file->Write();
  file->Close();
}

void NewElectrons::analyze(const Event & event, const EventSetup& eventSetup) {

  cout << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;
 
  // access the tracker
  edm::ESHandle<TrackerGeometry> theTrackerGeometry;
  eventSetup.get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);
  const TrackerGeometry& theTracker(*theTrackerGeometry);

  Handle<HepMCProduct> evt;
  event.getByLabel("source", evt);
  myGenEvent = new HepMC::GenEvent(*(evt->GetEvent()));

  Handle<PixelMatchGsfElectronCollection> elh;
  event.getByLabel(baselineEleCollName, elh);
  const PixelMatchGsfElectronCollection* copy = elh.product();
  PixelMatchGsfElectronCollection::const_iterator ite; 

  Handle<SuperClusterCollection> sch1;
  event.getByLabel("correctedHybridSuperClusters", sch1);
  const SuperClusterCollection* scb = sch1.product();

  Handle<SuperClusterCollection> sch2;
  event.getByLabel("correctedEndcapSuperClustersWithPreshower", sch2);
  const SuperClusterCollection* sce = sch2.product();
  SuperClusterCollection::const_iterator itscb, itsce;

  Handle<TrackCollection> tkh;
  event.getByLabel("ctfWithMaterialTracks", tkh);
  const TrackCollection* tracks = tkh.product();
  TrackCollection::const_iterator itt;

  //Handle<GlobalCtfElectronCollection> elh1;
  Handle<PixelMatchGsfElectronCollection> elh1;
  event.getByLabel(customEleCollName, elh1);
  //const GlobalCtfElectronCollection*  electrons1 = elh1.product();
  //GlobalCtfElectronCollection::const_iterator ite1;
  const PixelMatchGsfElectronCollection*  electrons1 = elh1.product();
  PixelMatchGsfElectronCollection::const_iterator ite1;

  for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it) { 

    if ((abs((*it)->pdg_id()) == 11) && ((*it)->status() != 3)) {      
      if (((*it)->momentum().perp() > 5.) && (fabs((*it)->momentum().eta()) < 2.5)) {

        mc_mother = mother(*it);
        mc_pt = (*it)->momentum().perp();
        mc_eta = (*it)->momentum().eta();
        mc_phi = (*it)->momentum().phi();
        mc_e = (*it)->momentum().e();
        mc_id = (*it)->pdg_id();

        // check if it is in a crack
        if (inCrack(fabs((*it)->momentum().eta())))
          mc_crack = 1;
        else
          mc_crack = 0;

        math::XYZVector mcv((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz());
        int type = 0;
        float theta = 1.;
        double dR, dRmin = 0.2;
        SuperClusterCollection::const_iterator nearSC;
        
        // loop over barrel superclusters
        for(itscb = scb->begin(); itscb != scb->end(); ++itscb) {
          math::XYZVector scv(itscb->x(), itscb->y(), itscb->z());
          
          dR = ROOT::Math::VectorUtil::DeltaR(scv, mcv);
          
          if (dR < dRmin) {
            dRmin = dR;
            nearSC = itscb;
            type = 0;
            theta = scv.Theta();
          }
        }
       
        for(itsce = sce->begin(); itsce != sce->end(); ++itsce) {
          math::XYZVector scv(itsce->x(), itsce->y(), itsce->z());
          dR = ROOT::Math::VectorUtil::DeltaR(scv, mcv);
          if (dR < dRmin) {
            dRmin = dR;
            nearSC = itsce;
            type = 1;
            theta = scv.Theta();
          }
        }
        // store info about SC	  
        if (dRmin < 0.2) {
          sc_e = nearSC->energy(); 
          sc_rawe = nearSC->rawEnergy();
          sc_et = sin(theta)*nearSC->energy();
          sc_eta = nearSC->eta();
          sc_phi = nearSC->phi(); 
          sc_dr = dRmin;
          sc_type = type;
        } else {
          sc_e = 0.;
          sc_rawe = 0;
          sc_et = 0.;
          sc_eta = 0.;
          sc_phi = 0.;
          sc_type = -1;
          sc_dr = 0.2; 
        }
  
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

        dRmin = 0.1;
        PixelMatchGsfElectronCollection::const_iterator nearElectron;	
        for(ite = electrons.begin(); ite != electrons.end(); ++ite) {
          dR = ROOT::Math::VectorUtil::DeltaR(ite->p4(), mcv);
          if (dR < dRmin) {
            dRmin = dR;
            nearElectron = ite;
          }
        }
        // strore info about Ele
        if (dRmin < 0.1) {
          el_pt = nearElectron->pt(); 
          el_eta = nearElectron->eta(); 
          el_e = nearElectron->energy();
          el_phi = nearElectron->phi(); 
          el_dr = dRmin; 
          el_eopin = nearElectron->eSuperClusterOverP();
          el_eopout = nearElectron->eSeedClusterOverPout();
          el_hoe = nearElectron->hadronicOverEm();
          el_dphiin = nearElectron->deltaPhiSuperClusterTrackAtVtx();
          el_detain = nearElectron->deltaEtaSuperClusterTrackAtVtx();
          el_dphiout = nearElectron->deltaPhiSeedClusterTrackAtCalo();
          el_detaout = nearElectron->deltaEtaSeedClusterTrackAtCalo();
          float pin  = nearElectron->trackMomentumAtVtx().R();
          float pout = nearElectron->trackMomentumOut().R();
          el_pout = pout;
          el_fbrem = (pin-pout)/pin;
          el_class = nearElectron->classification();
          R9_25_gsf(event, &(*nearElectron), el_eseed, el_e3x3, el_e5x5, el_spp, el_see);
          int a, b;
          nHits(nearElectron->gsfTrack(), a, b);
          el_npxhits = a;
          el_nsihits = b;

          int index = 1;
          while(1) {
            TrackingRecHitRef hit = nearElectron->gsfTrack()->recHit(nearElectron->gsfTrack()->recHitsSize()-index);
            
            if (hit->isValid()) {    

              GlobalPoint hitPosition = theTracker.idToDet(hit->geographicalId())->surface().toGlobal(hit->localPosition());
              GlobalPoint pos(hitPosition.x()-nearElectron->gsfTrack()->vx(), hitPosition.y()-nearElectron->gsfTrack()->vy(),
                              hitPosition.z()-nearElectron->gsfTrack()->vz());
              //std::cout << "Inner: " <<  HitPosition.perp() << "  " << HitPosition.z() << std::endl;
              el_rinnerhit = sqrt(pow(pos.perp(),2) + pow(pos.z(),2));
              //std::cout << "Inner: " << el_rinnerhit << std::endl; 
              subDetector(hit, a, b);
              el_detinnerhit = a;
              break;
            }
            index++;
          }

          el_z0 = nearElectron->gsfTrack()->vz();
          el_tkiso = trackIsolation(nearElectron->trackMomentumAtVtx(), nearElectron->vertex(), tracks);
        } else {
          el_pt = 0.;
          el_e = 0.;
          el_eta = 0.;
          el_phi = 0.;
          el_dr = 0.1;
          el_eopin = 0.;
          el_eopout = 0.;
          el_pout = 0;
          el_hoe = 0.;
          el_dphiin = 0.;
          el_detain = 0.;
          el_dphiout = 0.;
          el_detaout = 0.;
          el_fbrem = 0.;
          el_eseed = 0.;
          el_e3x3 = 0.;
          el_e5x5 = 0.;
          el_spp = 0.;
          el_see = 0.;
          el_class = -1;
          el_npxhits = -1;
          el_nsihits = -1;
          el_rinnerhit = 0.;
          el_detinnerhit = -1;
          el_z0 = -1;
          el_tkiso = -1;
        }

        // new electrons collection
        //GlobalCtfElectronCollection::const_iterator nearElectron1;
        PixelMatchGsfElectronCollection::const_iterator nearElectron1;
        dRmin = 0.1;
        for(ite1 = electrons1->begin(); ite1 != electrons1->end(); ++ite1) {
          dR = ROOT::Math::VectorUtil::DeltaR(ite1->p4(), mcv);
          if (dR < dRmin) {
            dRmin = dR;
            nearElectron1 = ite1;
          }
        }
 
        // strore info about Ele
        if (dRmin < 0.1) {
          el1_pt = nearElectron1->pt(); 
          el1_eta = nearElectron1->eta(); 
          el1_e = nearElectron1->energy();
          el1_phi = nearElectron1->phi(); 
          el1_dr = dRmin; 
          el1_eopin = nearElectron1->eSuperClusterOverP();
          el1_eopout = nearElectron1->eSeedClusterOverPout();
          el1_hoe = nearElectron1->hadronicOverEm();
          el1_dphiin = nearElectron1->deltaPhiSuperClusterTrackAtVtx();
          el1_detain = nearElectron1->deltaEtaSuperClusterTrackAtVtx();
          el1_dphiout = nearElectron1->deltaPhiSeedClusterTrackAtCalo();
          el1_detaout = nearElectron1->deltaEtaSeedClusterTrackAtCalo();
          float pin  = nearElectron1->trackMomentumAtVtx().R();
          float pout = nearElectron1->trackMomentumOut().R();
          el1_pout = pout;
          el1_fbrem = (pin-pout)/pin;
          el1_class = nearElectron1->classification();
          R9_25_gsf(event, &(*nearElectron1), el1_eseed, el1_e3x3, el1_e5x5, el1_spp, el1_see);
          int a, b;
          nHits(nearElectron1->gsfTrack(), a, b);
          el1_npxhits = a;
          el1_nsihits = b;
          int index = 1;
          while(1) {
            TrackingRecHitRef hit = nearElectron1->gsfTrack()->recHit(nearElectron1->gsfTrack()->recHitsSize()-index);
            
            if (hit->isValid()) {
              GlobalPoint hitPosition = theTracker.idToDet(hit->geographicalId())->surface().toGlobal(hit->localPosition());
              GlobalPoint pos(hitPosition.x()-nearElectron1->gsfTrack()->vx(), hitPosition.y()-nearElectron1->gsfTrack()->vy(),
                              hitPosition.z()-nearElectron1->gsfTrack()->vz());
              //std::cout << "Inner: " <<  HitPosition.perp() << "  " << HitPosition.z() << std::endl;
              el1_rinnerhit = sqrt(pow(pos.perp(),2) + pow(pos.z(),2));
              //math::XYZVector clusterGlobalPos(sc->x() - track->vx(), sc->y() - track->vy(), sc->z() - track->vz())
              //el1_rinnerhit = HitPosition.perp();
              subDetector(hit, a, b);
              el1_detinnerhit = a;
              break;
            }
            index++;
          }

          el1_z0 = nearElectron1->gsfTrack()->vz();
          el1_tkiso = trackIsolation(nearElectron1->trackMomentumAtVtx(), nearElectron1->vertex(), tracks);
        } else {
          el1_pt = 0.;
          el1_e = 0.;
          el1_eta = 0.;
          el1_phi = 0.;
          el1_dr = 0.1; 
          el1_eopin = 0.;
          el1_eopout = 0.;
          el1_pout = 0;
          el1_hoe = 0.;
          el1_dphiin = 0.;
          el1_detain = 0.;
          el1_dphiout = 0.;
          el1_detaout = 0.;
          el1_fbrem = 0.;
          el1_eseed = 0.;
          el1_e3x3 = 0.;
          el1_e5x5 = 0.;
          el1_spp = 0.;
          el1_see = 0.;
          el1_class = -1;
          el1_npxhits = -1;
          el1_nsihits = -1;
          el1_rinnerhit = 0.;
          el1_detinnerhit = -1;
          el1_z0 = -1;
          el1_tkiso = -1;
        }
        
        // loop over combinatorial track finder tracks
        dRmin = 0.1;
        TrackCollection::const_iterator nearTk;
        for(itt = tracks->begin(); itt != tracks->end(); ++itt) {
          dR = ROOT::Math::VectorUtil::DeltaR(itt->momentum(), mcv);
          if (dR < dRmin) {
            dRmin = dR;
            nearTk = itt;
          }
        }
         // strore info about Tk
        if (dRmin < 0.1) {
          tk_pt = nearTk->pt();
          tk_nhit = nearTk->found();
          
          tk_eta = nearTk->eta(); 
          tk_phi = nearTk->phi(); 
          tk_dr = dRmin;
        } else {
          tk_pt = 0.;
          tk_nhit = 0;
          //tk_subdet = -1;
          //tk_layer = -1;
          tk_eta = 0.;
          tk_phi = 0.;
          tk_dr = 0.2; 
        }
        run = event.id().run();
        id = event.id().event();
        tree->Fill();
      }
    }
  }
}

bool NewElectrons::inCrack(float eta) {

  return (eta < 0.018 ||
          (eta>0.423 && eta<0.461) ||
          (eta>0.770 && eta<0.806) ||
          (eta>1.127 && eta<1.163) ||
          (eta>1.460 && eta<1.558));
}

int NewElectrons::mother(HepMC::GenParticle *p) {
  
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

void NewElectrons::R9_25_gsf(const Event & event, const reco::PixelMatchGsfElectron* e,
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

/*
void NewElectrons::R9_25_ctf(const Event & event, const reco::GlobalCtfElectron* e,
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
*/
void NewElectrons::nHits(const reco::GsfTrackRef t, int& nPixelHits, int& nSiTkHits) {

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
double NewElectrons::trackIsolation(const math::XYZVector momentum, 
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

void NewElectrons::subDetector(TrackingRecHitRef hit, int& subdet, int& layer) {

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
