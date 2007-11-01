/*
  logical flow of the selection:
  1) save all SC with Et > 5. GeV
  2) for each saved SC search for matched following object (within dR = 0.1):
     a) Standard electron
     b) Our electron
     c) MC particle
     d) CTF track
     e) nearest MC particle to stored CTF track
*/

#include "NtupleMakers/ElectronAnalyzer/interface/Conversion.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"

#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackReco/interface/Track.h"

#include "TFile.h"
#include "TTree.h"

using namespace std;
using namespace edm;
using namespace reco;

/// Constructor
Conversion::Conversion(const ParameterSet& pset) {
  fileName = pset.getParameter<std::string>("RootFileName");
  baselineEleCollName =  pset.getParameter<std::string>("BaselineEleCollName");
  customEleCollName   =  pset.getParameter<std::string>("CustomEleCollName");
}

Conversion::~Conversion() {}

void Conversion::beginJob(const EventSetup& eventSetup) {

  file = new TFile(fileName.c_str(), "recreate");
  tree = new TTree("event","Event data");

  tree->Branch("run", &run, "run/I");
  tree->Branch("id", &id, "id/I");

  // tree->Branch("sc_n", &nSC, "sc_n/I");
  tree->Branch("sc_e", &sc_e, "sc_e/F");
  tree->Branch("sc_rawe", &sc_rawe, "sc_rawe/F");
  tree->Branch("sc_et", &sc_et, "sc_et/F");
  tree->Branch("sc_eta", &sc_eta, "sc_eta/F");
  tree->Branch("sc_phi", &sc_phi, "sc_phi/F");
  tree->Branch("sc_dr", &sc_dr, "sc_dr/F");
  tree->Branch("sc_type", &sc_type, "sc_type/I"); // 0 barrel, 1 endcap

  tree->Branch("mcsc_pt", &mcsc_pt, "mcsc_pt/F");
  tree->Branch("mcsc_dr", &mcsc_dr, "mcsc_dr/F");
  tree->Branch("mcsc_eta", &mcsc_eta, "mcsc_eta/F");
  tree->Branch("mcsc_phi", &mcsc_phi, "mcsc_phi/F");
  tree->Branch("mcsc_id", &mcsc_id, "mcsc_id/I");
  tree->Branch("mcsc_e", &mcsc_e, "mcsc_e/F");
  tree->Branch("mcsc_mother", &mcsc_mother, "mcsc_mother/I");
  tree->Branch("mcsc_crack", &mcsc_crack, "mcsc_crack/I"); // 0 no crack, 1 crack

  tree->Branch("mctk_pt", &mctk_pt, "mctk_pt/F");
  tree->Branch("mctk_dr", &mctk_dr, "mctk_dr/F");
  tree->Branch("mctk_eta", &mctk_eta, "mctk_eta/F");
  tree->Branch("mctk_phi", &mctk_phi, "mctk_phi/F");
  tree->Branch("mctk_id", &mctk_id, "mctk_id/I");
  tree->Branch("mctk_e", &mctk_e, "mctk_e/F");
  tree->Branch("mctk_mother", &mctk_mother, "mctk_mother/I");
  tree->Branch("mctk_crack", &mctk_crack, "mctk_crack/I"); // 0 no crack, 1 crack

  tree->Branch("mctk1_pt", &mctk1_pt, "mctk1_pt/F");
  tree->Branch("mctk1_dr", &mctk1_dr, "mctk1_dr/F");
  tree->Branch("mctk1_eta", &mctk1_eta, "mctk1_eta/F");
  tree->Branch("mctk1_phi", &mctk1_phi, "mctk1_phi/F");
  tree->Branch("mctk1_id", &mctk1_id, "mctk1_id/I");
  tree->Branch("mctk1_e", &mctk1_e, "mctk1_e/F");
  tree->Branch("mctk1_mother", &mctk1_mother, "mctk1_mother/I");
  tree->Branch("mctk1_crack", &mctk1_crack, "mctk1_crack/I"); // 0 no crack, 1 crack

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
  tree->Branch("el_tkpt", &el_tkpt, "el_tkpt/F");
  tree->Branch("el_tketa", &el_tketa, "el_tketa/F");
  tree->Branch("el_tkphi", &el_tkphi, "el_tkphi/F");

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
  tree->Branch("el1_tkpt", &el1_tkpt, "el1_tkpt/F");
  tree->Branch("el1_tketa", &el1_tketa, "el1_tketa/F");
  tree->Branch("el1_tkphi", &el1_tkphi, "el1_tkphi/F");
}

void Conversion::endJob() {

  file->Write();
  file->Close();
}

void Conversion::analyze(const Event & event, const EventSetup& eventSetup) {

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
  event.getByLabel("hybridSuperClusters", sch1);
  const SuperClusterCollection* scb = sch1.product();

  Handle<SuperClusterCollection> sch2;
  event.getByLabel("islandSuperClusters", "islandEndcapSuperClusters", sch2);
  const SuperClusterCollection* sce = sch2.product();
  SuperClusterCollection::const_iterator itsc, itscb, itsce;
  SuperClusterCollection sc;
  sc.insert(sc.end(), scb->begin(), scb->end());
  sc.insert(sc.end(), sce->begin(), sce->end());

  Handle<PixelMatchGsfElectronCollection> elh1;
  event.getByLabel(customEleCollName, elh1);
  const PixelMatchGsfElectronCollection*  electrons1 = elh1.product();
  PixelMatchGsfElectronCollection::const_iterator ite1;

  Handle<TrackCollection> tkh;
  event.getByLabel("ctfWithMaterialTracks", tkh);
  const TrackCollection* tracks = tkh.product();
  TrackCollection::const_iterator itt;

  for (itsc = sc.begin(); itsc != sc.end(); ++itsc) {
    
    math::XYZVector scv(itsc->x(), itsc->y(), itsc->z());
    
    if (sin(scv.Theta())*itsc->energy() > 5.) {
      sc_e = itsc->energy(); 
      sc_rawe = itsc->rawEnergy();
      sc_et = sin(scv.Theta())*itsc->energy();
      sc_eta = itsc->eta();
      sc_phi = itsc->phi(); 
      
      if (fabs(itsc->eta()) < 1.47)
        sc_type = 0;
      else
        sc_type = 1;
      
      double dR, dRmin = 0.1;
      HepMC::GenEvent::particle_const_iterator nearMC;
      for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it) { 
        
        if ((*it)->status() == 1) {      
          math::XYZVector mcv((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz());
          dR = ROOT::Math::VectorUtil::DeltaR(scv, mcv);
          
          if (dR < dRmin) {
            dRmin = dR;
            nearMC = it;
          }
        }
      }
      
      if (dRmin < 0.1) {
        mcsc_dr = dRmin;
        mcsc_mother = mother(*nearMC);
        mcsc_pt = (*nearMC)->momentum().perp();
        mcsc_eta = (*nearMC)->momentum().eta();
        mcsc_phi = (*nearMC)->momentum().phi();
        mcsc_e = (*nearMC)->momentum().e();
        mcsc_id = (*nearMC)->pdg_id();
        
        // check if it is in a crack
        if (inCrack(fabs((*nearMC)->momentum().eta())))
          mcsc_crack = 1;
        else
          mcsc_crack = 0;  
      } else {
        mcsc_dr = 0.1;
        mcsc_mother = -1;
        mcsc_pt = -1;
        mcsc_eta = -1;
        mcsc_phi = -1;
        mcsc_e = -1;
        mcsc_id = 0;
        mcsc_crack = -1;
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
        dR = ROOT::Math::VectorUtil::DeltaR(ite->p4(), scv);
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
        el_eseed = nearElectron->superCluster()->seed()->energy();
        el_e3x3 = nearElectron->seedClusterShape()->e3x3();
        el_e5x5 = nearElectron->seedClusterShape()->e5x5();
        el_spp = sqrt(nearElectron->seedClusterShape()->covPhiPhi());
        el_see = sqrt(nearElectron->seedClusterShape()->covEtaEta());
        
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
            subDetector(hit, a, b);
            el_detinnerhit = a;
            break;
          }
          index++;
        }

        el_z0 = nearElectron->gsfTrack()->vz();
        el_tkiso = trackIsolation(nearElectron->trackMomentumAtVtx(), nearElectron->vertex(), tracks);
        el_tkpt = nearElectron->gsfTrack()->pt();
        el_tketa = nearElectron->gsfTrack()->eta(); 
        el_tkphi = nearElectron->gsfTrack()->phi();  

        double dR, dRmin = 0.1;
        HepMC::GenEvent::particle_const_iterator nearMC;
        for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it) { 
          
          if ((*it)->status() == 1) {      
            math::XYZVector mcv((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz());
            dR = ROOT::Math::VectorUtil::DeltaR(nearElectron->gsfTrack()->innerMomentum(), mcv);
            
            if (dR < dRmin) {
              dRmin = dR;
              nearMC = it;
            }
          }
        }
        
        if (dRmin < 0.1) {
          mctk_dr = dRmin;
          mctk_mother = mother(*nearMC);
          mctk_pt = (*nearMC)->momentum().perp();
          mctk_eta = (*nearMC)->momentum().eta();
          mctk_phi = (*nearMC)->momentum().phi();
          mctk_e = (*nearMC)->momentum().e();
          mctk_id = (*nearMC)->pdg_id();
          
          // check if it is in a crack
          if (inCrack(fabs((*nearMC)->momentum().eta())))
            mctk_crack = 1;
          else
            mctk_crack = 0;  
        } else {
          mctk_dr = 0.1;
          mctk_mother = -1;
          mctk_pt = -1;
          mctk_eta = -1;
          mctk_phi = -1;
          mctk_e = -1;
          mctk_id = 0;
          mctk_crack = -1;
        }
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
        el_tkpt = 0.;
        el_tketa = 0.; 
        el_tkphi = 0.; 
        mctk_dr = 0.1;
        mctk_mother = -1;
        mctk_pt = -1;
        mctk_eta = -1;
        mctk_phi = -1;
        mctk_e = -1;
        mctk_id = -1;
        mctk_crack = -1;
      }
      
      // new electrons collection
      //GlobalCtfElectronCollection::const_iterator nearElectron1;
      PixelMatchGsfElectronCollection::const_iterator nearElectron1;
      dRmin = 0.1;
      for(ite1 = electrons1->begin(); ite1 != electrons1->end(); ++ite1) {
        dR = ROOT::Math::VectorUtil::DeltaR(ite1->p4(), scv);
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
        el1_eseed = nearElectron1->superCluster()->seed()->energy();
        el1_e3x3 = nearElectron1->seedClusterShape()->e3x3();
        el1_e5x5 = nearElectron1->seedClusterShape()->e5x5();
        el1_spp = sqrt(nearElectron1->seedClusterShape()->covPhiPhi());
        el1_see = sqrt(nearElectron1->seedClusterShape()->covEtaEta());

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
            subDetector(hit, a, b);
            el1_detinnerhit = a;
            break;
          }
          index++;
        }
        
        el1_z0 = nearElectron1->gsfTrack()->vz();
        el1_tkiso = trackIsolation(nearElectron1->trackMomentumAtVtx(), nearElectron1->vertex(), tracks);
        el1_tkpt = nearElectron1->gsfTrack()->pt();
        el1_tketa = nearElectron1->gsfTrack()->eta(); 
        el1_tkphi = nearElectron1->gsfTrack()->phi(); 
        
        double dR, dRmin = 0.1;
        HepMC::GenEvent::particle_const_iterator nearMC;
        for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it) { 
          
          if ((*it)->status() == 1) {      
            math::XYZVector mcv((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz());
            dR = ROOT::Math::VectorUtil::DeltaR(nearElectron1->gsfTrack()->innerMomentum(), mcv);
            
            if (dR < dRmin) {
              dRmin = dR;
              nearMC = it;
            }
          }
        }
        
        if (dRmin < 0.1) {
          mctk1_dr = dRmin;
          mctk1_mother = mother(*nearMC);
          mctk1_pt = (*nearMC)->momentum().perp();
          mctk1_eta = (*nearMC)->momentum().eta();
          mctk1_phi = (*nearMC)->momentum().phi();
          mctk1_e = (*nearMC)->momentum().e();
          mctk1_id = (*nearMC)->pdg_id();
          
          // check if it is in a crack
          if (inCrack(fabs((*nearMC)->momentum().eta())))
            mctk1_crack = 1;
          else
            mctk1_crack = 0;  
        } else {
          mctk1_dr = 0.1;
          mctk1_mother = -1;
          mctk1_pt = -1;
          mctk1_eta = -1;
          mctk1_phi = -1;
          mctk1_e = -1;
          mctk1_id = 0;
          mctk1_crack = -1;
        }
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
        el1_rinnerhit = 0.;
        el1_detinnerhit = -1;
        el1_npxhits = -1;
        el1_nsihits = -1;
        el1_z0 = -1;
        el1_tkiso = -1;
        el1_tkpt = 0.;
        el1_tketa = 0.; 
        el1_tkphi = 0.; 
        mctk1_dr = 0.1;
        mctk1_mother = -1;
        mctk1_pt = -1;
        mctk1_eta = -1;
        mctk1_phi = -1;
        mctk1_e = -1;
        mctk1_id = -1;
        mctk1_crack = -1;
      }
      
      run = event.id().run();
      id = event.id().event();
      tree->Fill();
    }
  }
}

bool Conversion::inCrack(float eta) {

  return (eta < 0.018 ||
          (eta>0.423 && eta<0.461) ||
          (eta>0.770 && eta<0.806) ||
          (eta>1.127 && eta<1.163) ||
          (eta>1.460 && eta<1.558));
}

int Conversion::mother(HepMC::GenParticle *p) {
  
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

void Conversion::nHits(const reco::GsfTrackRef t, int& nPixelHits, int& nSiTkHits) {

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
double Conversion::trackIsolation(const math::XYZVector momentum, 
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

void Conversion::subDetector(TrackingRecHitRef hit, int& subdet, int& layer) {

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
