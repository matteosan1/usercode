#include "NtupleMakers/HiggsAnalyzer/interface/hzz.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
//#include "HepMC/GenEvent.h"
//#include "HepMC/GenVertex.h"
//#include "HepMC/GenParticle.h"

//#include "CLHEP/HepMC/GenParticle.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EgammaReco/interface/ClusterShape.h"

#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"

//#include "RecoEgammaBis/EgammaElectronAlgos/interface/ElectronAlgoA.h"

#include "TFile.h"
#include "TTree.h"

using namespace std;
using namespace edm;
using namespace reco;

/// Constructor
hzz::hzz(const ParameterSet& pset) {
  fileName = pset.getParameter<std::string>("RootFileName");
  baselineEleCollName =  pset.getParameter<std::string>("BaselineEleCollName");
  customEleCollName   =  pset.getParameter<std::string>("CustomEleCollName");
}

hzz::~hzz() {}

void hzz::beginJob(const EventSetup& eventSetup) {

  file = new TFile(fileName.c_str(), "recreate");
  tree = new TTree("event","Event data");

  tree->Branch("run", &run, "run/I");
  tree->Branch("id", &id, "id/I");

  tree->Branch("mc_n", &nMC, "mc_n/I");
  tree->Branch("mc_pt", mc_pt, "mc_pt[mc_n]/F");
  tree->Branch("mc_eta", mc_eta, "mc_eta[mc_n]/F");
  tree->Branch("mc_phi", mc_phi, "mc_phi[mc_n]/F");
  tree->Branch("mc_id", mc_id, "mc_id[mc_n]/I");
  tree->Branch("mc_e", mc_e, "mc_e[mc_n]/F");
  tree->Branch("mc_mother", mc_mother, "mc_mother[mc_n]/I");
  tree->Branch("mc_crack", mc_crack, "mc_crack[mc_n]/I"); // 0 no crack, 1 crack

  tree->Branch("sc_e", sc_e, "sc_e[mc_n]/F");
  tree->Branch("sc_rawe", sc_rawe, "sc_rawe[mc_n]/F");
  tree->Branch("sc_et", sc_et, "sc_et[mc_n]/F");
  tree->Branch("sc_eta", sc_eta, "sc_eta[mc_n]/F");
  tree->Branch("sc_phi", sc_phi, "sc_phi[mc_n]/F");
  tree->Branch("sc_dr", sc_dr, "sc_dr[mc_n]/F");
  tree->Branch("sc_type", sc_type, "sc_type[mc_n]/I"); // 0 barrel, 1 endcap

  /*
  tree->Branch("tk_pt", tk_pt, "tk_e[mc_n]/F");
  tree->Branch("tk_eta", tk_eta, "tk_eta[mc_n]/F");
  tree->Branch("tk_phi", tk_phi, "tk_phi[mc_n]/F");
  tree->Branch("tk_dr", tk_dr, "tk_dr[mc_n]/F");
  tree->Branch("tk_nhit", tk_nhit, "tk_nhit[mc_n]/I");
  tree->Branch("tk_layer", tk_layer, "tk_layer[mc_n]/I");
  tree->Branch("tk_subdet", tk_subdet, "tk_subdet[mc_n]/I");
  */

  tree->Branch("el_pt", el_pt, "el_pt[mc_n]/F");
  tree->Branch("el_e", el_e, "el_e[mc_n]/F");
  tree->Branch("el_eta", el_eta, "el_eta[mc_n]/F");
  tree->Branch("el_phi", el_phi, "el_phi[mc_n]/F");
  tree->Branch("el_dr", el_dr, "el_dr[mc_n]/F");
  tree->Branch("el_eopin", el_eopin, "el_eopin[mc_n]/F");
  tree->Branch("el_eopout", el_eopout, "el_eopout[mc_n]/F");
  tree->Branch("el_pout", el_pout, "el_pout[mc_n]/F");
  tree->Branch("el_fbrem", el_fbrem, "el_fbrem[mc_n]/F");
  tree->Branch("el_hoe", el_hoe, "el_hoe[mc_n]/F");
  tree->Branch("el_detain", el_detain, "el_detain[mc_n]/F");
  tree->Branch("el_dphiin", el_dphiin, "el_dphiin[mc_n]/F");
  tree->Branch("el_detaout", el_detaout, "el_detaout[mc_n]/F");
  tree->Branch("el_dphiout", el_dphiout, "el_dphiout[mc_n]/F");
  tree->Branch("el_e3x3", el_e3x3, "el_e3x3[mc_n]/F");
  tree->Branch("el_e5x5", el_e5x5, "el_e5x5[mc_n]/F");
  tree->Branch("el_eseed", el_eseed, "el_eseed[mc_n]/F");
  tree->Branch("el_spp", el_spp, "el_spp[mc_n]/F");
  tree->Branch("el_see", el_see, "el_see[mc_n]/F");
  tree->Branch("el_class", el_class, "el_class[mc_n]/I");
  tree->Branch("el_tkhits", el_tkhits, "el_tkhits[mc_n]/I");
  tree->Branch("el_pixelhits", el_pixelhits, "el_pixelhits[mc_n]/I");
  tree->Branch("el_gsfpt", el_gsfpt, "el_gsfpt[mc_n]/F");
  tree->Branch("el_gsfeta", el_gsfeta, "el_gsfeta[mc_n]/F");
  tree->Branch("el_scpt", el_scpt, "el_scpt[mc_n]/F");
  tree->Branch("el_sceta", el_sceta, "el_sceta[mc_n]/F");
  tree->Branch("el_scphi", el_scphi, "el_scphi[mc_n]/F");



  tree->Branch("el1_pt", el1_pt, "el1_pt[mc_n]/F");
  tree->Branch("el1_e", el1_e, "el1_e[mc_n]/F");
  tree->Branch("el1_eta", el1_eta, "el1_eta[mc_n]/F");
  tree->Branch("el1_phi", el1_phi, "el1_phi[mc_n]/F");
  tree->Branch("el1_dr", el1_dr, "el1_dr[mc_n]/F");
  tree->Branch("el1_eopin", el1_eopin, "el1_eopin[mc_n]/F");
  tree->Branch("el1_eopout", el1_eopout, "el1_eopout[mc_n]/F");
  tree->Branch("el1_pout", el1_pout, "el1_pout[mc_n]/F");
  tree->Branch("el1_fbrem", el1_fbrem, "el1_fbrem[mc_n]/F");
  tree->Branch("el1_hoe", el1_hoe, "el1_hoe[mc_n]/F");
  tree->Branch("el1_detain", el1_detain, "el1_detain[mc_n]/F");
  tree->Branch("el1_dphiin", el1_dphiin, "el1_dphiin[mc_n]/F");
  tree->Branch("el1_detaout", el1_detaout, "el1_detaout[mc_n]/F");
  tree->Branch("el1_dphiout", el1_dphiout, "el1_dphiout[mc_n]/F");
  tree->Branch("el1_e3x3", el1_e3x3, "el1_e3x3[mc_n]/F");
  tree->Branch("el1_e5x5", el1_e5x5, "el1_e5x5[mc_n]/F");
  tree->Branch("el1_eseed", el1_eseed, "el1_eseed[mc_n]/F");
  tree->Branch("el1_spp", el1_spp, "el1_spp[mc_n]/F");
  tree->Branch("el1_see", el1_see, "el1_see[mc_n]/F");
  tree->Branch("el1_class", el1_class, "el1_class[mc_n]/I");
  tree->Branch("el1_tkhits", el1_tkhits, "el1_tkhits[mc_n]/I");
  tree->Branch("el1_pixelhits", el1_pixelhits, "el1_pixelhits[mc_n]/I");
  tree->Branch("el1_scpt", el1_scpt, "el1_scpt[mc_n]/F");
  tree->Branch("el1_sceta", el1_sceta, "el1_sceta[mc_n]/F");
  tree->Branch("el1_scphi", el1_scphi, "el1_scphi[mc_n]/F");

  tree->Branch("nz", &nz, "nz/I");
  tree->Branch("mz1", &mz1, "mz1/F");
  tree->Branch("mz2", &mz2, "mz2/F");
  tree->Branch("mzh", &mzh, "mzh/F");
  tree->Branch("mzl", &mzl, "mzl/F");
  tree->Branch("mh", &mh, "mh/F");

  tree->Branch("rmz1", &rmz1, "rmz1/F");
  tree->Branch("rmz2", &rmz2, "rmz2/F");
  tree->Branch("rmh", &rmh, "rmh/F");
   
  tree->Branch("qmz1", &qmz1, "qmz1/F");
  tree->Branch("qmz2", &qmz2, "qmz2/F");
  tree->Branch("qmh", &qmh, "qmh/F");
                                  
}

void hzz::endJob() {

  file->Write();
  file->Close();
}

void hzz::analyze(const Event & event, const EventSetup& eventSetup) {


  Handle<HepMCProduct> evt;
  event.getByLabel("source", evt);
  myGenEvent = new HepMC::GenEvent(*(evt->GetEvent()));
  HepMC::GenParticle* mce2=0; HepMC::GenParticle* mcp2=0; HepMC::GenParticle* mce1=0; HepMC::GenParticle* mcp1=0; HepMC::GenParticle* it=0;

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
  SuperClusterCollection::const_iterator itscb, itsce;

  // Handle<TrackCollection> tkh;
  //event.getByLabel("ctfWithMaterialTracks", tkh);
  //const TrackCollection* tracks = tkh.product();
  //TrackCollection::const_iterator itt;

  //Handle<GlobalCtfElectronCollection> elh1;
  Handle<PixelMatchGsfElectronCollection> elh1;
  event.getByLabel(customEleCollName, elh1);
  //const GlobalCtfElectronCollection*  electrons1 = elh1.product();
  //GlobalCtfElectronCollection::const_iterator ite1;
  const PixelMatchGsfElectronCollection*  electrons1 = elh1.product();
  PixelMatchGsfElectronCollection::const_iterator ite1;


  math::XYZTLorentzVector rce1(0,0,0,0), rcp1(0,0,0,0), rce2(0,0,0,0), rcp2(0,0,0,0);
  math::XYZTLorentzVector qce1(0,0,0,0), qcp1(0,0,0,0), qce2(0,0,0,0), qcp2(0,0,0,0);
  bool ele1=false, pos1=false, ele2=false, pos2=false;
  bool z1=false, z2=false;
  nz=0; 
  for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it) { 
    
    if ( (abs((*it)->pdg_id()) == 23) && (*it)->end_vertex() ) 
      {
	//cout << "debug 1" << endl;
	//cout << (*it)->pdg_id() << " status: " << (*it)->status() << endl;
	//cout << "ndau: " << (*it)->end_vertex()->particles_out_size() << endl ;
	
	int ndau = (*it)->end_vertex()->particles_out_size() ;
	// First Z:
	if ( ndau == 3 && z1 == false ) 
	  {
	    z1 = true;
	    for(HepMC::GenVertex::particles_out_const_iterator junk = (*it)->end_vertex()->particles_out_const_begin();
		junk != (*it)->end_vertex()->particles_out_const_end(); ++junk)
	      {
		if ( (*junk)->pdg_id()  ==  11 ) { 
		  pos1 = true;
		  mcp1 = *junk; }
		if ( (*junk)->pdg_id()  == -11 ) { 
		  ele1 = true;
		  mce1 = *junk; }
		// cout << "duaghters of z1 " <<  (*junk)->pdg_id() << " " << (*junk)->status() << endl;
		if (ele1 && pos1 && nz==0) 
		  { 
		    nz++; 
		}
	      }
	    //cout << "nz: " << nz << endl;
	    continue;
	  }
	// Second Z:
	if ( ndau == 3 && z1 == true && z2 == false ) 
	  {
	    z2 = true;
	    for(HepMC::GenVertex::particles_out_const_iterator junk = (*it)->end_vertex()->particles_out_const_begin();
		junk != (*it)->end_vertex()->particles_out_const_end(); ++junk)
	      {
		if ( (*junk)->pdg_id()  ==  11 ) { 
		  pos2 = true;
		  mcp2 = (*junk); }
		if ( (*junk)->pdg_id()  == -11 ) { 
		  ele2 = true;
		  mce2 = *junk; }
		// cout << "duaghters of z2 " << (*junk)->pdg_id() << " " << (*junk)->status() << endl;
		if (ele2 && pos2 && nz==1) 
		  { 
		    nz++; 
		  }
	      }
	    //	    cout << "nz: " << nz << endl;
	    continue;
	  }
      }
  }
  if ( ele1 && pos1 ) {
    // cout << " trying to compute mz1 "<< endl;
    math::XYZVector ve1((mce1)->momentum().px(), (mce1)->momentum().py(), (mce1)->momentum().pz());
    math::XYZVector vp1((mcp1)->momentum().px(), (mcp1)->momentum().py(), (mcp1)->momentum().pz());
    mz1 = sqrt ( 2. * (mce1)->momentum().e() * (mcp1)->momentum().e() * ( 1 - ROOT::Math::VectorUtil::CosTheta(ve1, vp1) ) );
    //    cout << " mz1 " << mz1 << " " << endl;
  }
  if ( ele2 && pos2 ){
    //cout << " trying to compute mz2 "<< endl;
    math::XYZVector ve2((mce2)->momentum().px(), (mce2)->momentum().py(), (mce2)->momentum().pz());
    math::XYZVector vp2((mcp2)->momentum().px(), (mcp2)->momentum().py(), (mcp2)->momentum().pz());
    mz2 = sqrt ( 2. * (mce2)->momentum().e() * (mcp2)->momentum().e() * ( 1 - ROOT::Math::VectorUtil::CosTheta(ve2, vp2) ) );
    //    cout << " mz2: " << mz2 << endl;
  }
  
   math::XYZTLorentzVector vz1( mce1->momentum().px()+mcp1->momentum().px(),
				mce1->momentum().py()+mcp1->momentum().py(), 
				mce1->momentum().pz()+mcp1->momentum().pz(),
				mce1->momentum().e() +mcp1->momentum().e() );
   math::XYZTLorentzVector vz2( mce2->momentum().px()+mcp2->momentum().px(),
				mce2->momentum().py()+mcp2->momentum().py(), 
				mce2->momentum().pz()+mcp2->momentum().pz(),
				mce2->momentum().e() +mcp2->momentum().e() );
   /*
  math::XYZTLorentzVector vz1( mce1->p4()+mcp1->p4() );
  math::XYZTLorentzVector vz2( mce2->p4()+mcp2->p4() );
   */
  math::XYZTLorentzVector vh = vz1 + vz2;
   mh = vh.mass();
   //cout << " mH " << mh << endl;
   //cout << " mz1, mz2: " << vz1.mass() << " " << vz2.mass() << endl;

   mzh = mz1;
   mzl = mz2;
   if (mzh < mzl) 
     {
       mzh = mz2;
       mzl = mz1;
     }
   
  
   
   for (nMC=0; nMC<4; ++nMC) { 
     if (nMC == 0) { it = mce1; }
     if (nMC == 1) { it = mcp1; }
     if (nMC == 2) { it = mce2; }
     if (nMC == 3) { it = mcp2; }
     
     //    if ((abs((*it)->pdg_id()) == 11) && ((*it)->status() != 3)) {      
     //      if ( ((it) == mce1) ||  ((it) == mcp1) ||  ((it) == mce2) || ((it) == mcp2) )  {
     
     mc_mother[nMC] = mother(it);
     mc_pt[nMC] = (it)->momentum().perp();
     mc_eta[nMC] = (it)->momentum().eta();
     mc_phi[nMC] = (it)->momentum().phi();
     mc_e[nMC] = (it)->momentum().e();
     mc_id[nMC] = (it)->pdg_id();
    
    // check if it is in a crack
    if (inCrack(fabs((it)->momentum().eta())))
      mc_crack[nMC] = 1;
    else
      mc_crack[nMC] = 0;
    
    math::XYZVector mcv((it)->momentum().px(), (it)->momentum().py(), (it)->momentum().pz());
    int type = 0;
    float theta = 1.;
    double dR, dRmin = 0.2;
    SuperClusterCollection::const_iterator nearSC;
    
    // loop over barrel superclusters
    for(itscb = scb->begin(); itscb != scb->end(); ++itscb) {
      math::XYZVector scv(itscb->x(), itscb->y(), itscb->z());
      /*
      std::cout << "DEBUG_CLUSTER: sc eta,phi, energy: " 
		<< scv.Eta() << " , "  
		<< scv.Phi() << " , " 
		<< itscb->energy() << std::endl;
      */
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
      /*
      std::cout << "DEBUG_CLUSTER: sc eta,phi, energy: " 
		<< scv.Eta() << " , "  
		<< scv.Phi() << " , " 
		<< itsce->energy() << std::endl;
      */
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
      sc_e[nMC] = nearSC->energy(); 
      sc_rawe[nMC] = nearSC->rawEnergy();
      sc_et[nMC] = sin(theta)*nearSC->energy();
      sc_eta[nMC] = nearSC->eta();
      sc_phi[nMC] = nearSC->phi(); 
      sc_dr[nMC] = dRmin;
      sc_type[nMC] = type;
    } else {
      sc_e[nMC] = 0.;
      sc_rawe[nMC] = 0;
      sc_et[nMC] = 0.;
      sc_eta[nMC] = 0.;
      sc_phi[nMC] = 0.;
      sc_type[nMC] = -1;
      sc_dr[nMC] = 0.2; 
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
      if (nMC == 0) { qce1 = nearElectron->p4(); }
      if (nMC == 1) { qcp1 = nearElectron->p4(); }
      if (nMC == 2) { qce2 = nearElectron->p4(); }
      if (nMC == 3) { qcp2 = nearElectron->p4(); }
      /*
      std::cout << "DEBUG, qelectron phi,eta, pt, pz, p, e, mass: " 
		<< nearElectron->p4().Phi() << " , "
		<< nearElectron->p4().Eta() << " , "
		<< nearElectron->p4().Pt() << " , "
		<< nearElectron->p4().Pz() << " , "
		<< nearElectron->p4().P() << " , "
		<< nearElectron->p4().E() << " , "
		<< nearElectron->p4().M()  << std::endl;
      */
      el_pt[nMC] = sqrt(nearElectron->trackMomentumAtVtx().perp2()); 
      el_eta[nMC] = nearElectron->eta(); 
      el_e[nMC] = nearElectron->energy();
      el_phi[nMC] = nearElectron->phi(); 
      el_dr[nMC] = dRmin; 
      el_eopin[nMC] = nearElectron->eSuperClusterOverP();
      el_eopout[nMC] = nearElectron->eSeedClusterOverPout();
      el_hoe[nMC] = nearElectron->hadronicOverEm();
      el_dphiin[nMC] = nearElectron->deltaPhiSuperClusterTrackAtVtx();
      el_detain[nMC] = nearElectron->deltaEtaSuperClusterTrackAtVtx();
      el_dphiout[nMC] = nearElectron->deltaPhiSeedClusterTrackAtCalo();
      el_detaout[nMC] = nearElectron->deltaEtaSeedClusterTrackAtCalo();
      float pin  = nearElectron->trackMomentumAtVtx().R();
      float pout = nearElectron->trackMomentumOut().R();
      el_pout[nMC] = pout;
      el_fbrem[nMC] = (pin-pout)/pin;
      el_class[nMC] = nearElectron->classification();
      el_pixelhits[nMC] = nearElectron->gsfTrack()->hitPattern().numberOfValidPixelHits();
      el_tkhits[nMC] = nearElectron->gsfTrack()->hitPattern().numberOfValidTrackerHits();
      el_gsfpt[nMC] = nearElectron->gsfTrack()->pt();
      el_gsfeta[nMC] = nearElectron->gsfTrack()->eta();
      math::XYZVector trackGlobalDir(nearElectron->gsfTrack()->momentum());
      reco::GsfTrackRef track = nearElectron->gsfTrack();
      reco::SuperClusterRef sc = nearElectron->superCluster(); 
      math::XYZVector clusterDir(sc->x()-track->vx(), sc->y()-track->vy(), sc->z()-track->vz());
      el_scpt[nMC] = sc->energy()*sin(clusterDir.theta());
      el_sceta[nMC] = clusterDir.eta(); 
      el_scphi[nMC] = clusterDir.phi();

      R9_25_gsf(event, &(*nearElectron), el_eseed[nMC], el_e3x3[nMC], el_e5x5[nMC], el_spp[nMC], el_see[nMC]);
    } else {
      if (nMC == 0) { qce1 = math::XYZTLorentzVector(0,0,0,0); }
      if (nMC == 1) { qcp1 = math::XYZTLorentzVector(0,0,0,0); }
      if (nMC == 2) { qce2 = math::XYZTLorentzVector(0,0,0,0); }
      if (nMC == 3) { qcp2 = math::XYZTLorentzVector(0,0,0,0); }
      el_pt[nMC] = 0.;
      el_e[nMC] = 0.;
      el_eta[nMC] = 0.;
      el_phi[nMC] = 0.;
      el_dr[nMC] = 0.1;
      el_eopin[nMC] = 0.;
      el_eopout[nMC] = 0.;
      el_pout[nMC] = 0;
      el_hoe[nMC] = 0.;
      el_dphiin[nMC] = 0.;
      el_detain[nMC] = 0.;
      el_dphiout[nMC] = 0.;
      el_detaout[nMC] = 0.;
      el_fbrem[nMC] = 0.;
      el_eseed[nMC] = 0.;
      el_e3x3[nMC] = 0.;
      el_e5x5[nMC] = 0.;
      el_spp[nMC] = 0.;
      el_see[nMC] = 0.;
      el_class[nMC] = -1;
      el_pixelhits[nMC] = -1;
      el_tkhits[nMC] = -1;
      el_gsfpt[nMC] = 0.;
      el_gsfeta[nMC] = 0.;
      el_scpt[nMC] = 0.;
      el_sceta[nMC] = 0.; 
      el_scphi[nMC] = 0.; 

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
    //cout << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;
    // strore info about Ele
    if (dRmin < 0.1) {
      if (nMC == 0) { rce1 = nearElectron1->p4(); }
      if (nMC == 1) { rcp1 = nearElectron1->p4(); }
      if (nMC == 2) { rce2 = nearElectron1->p4(); }
      if (nMC == 3) { rcp2 = nearElectron1->p4(); }
      /*
      std::cout << "DEBUG, relectron phi,eta, pt,pz,p, e,  mass: " 
		<< nearElectron1->p4().Phi() << " , "
		<< nearElectron1->p4().Eta() << " , "
		<< nearElectron1->p4().Pt() << " , "
		<< nearElectron1->p4().Pz() << " , "
		<< nearElectron1->p4().P() << " , "
		<< nearElectron1->p4().E() << " , "
		<< nearElectron1->p4().M()  << std::endl;
      */

      el1_pt[nMC] = sqrt(nearElectron1->trackMomentumAtVtx().perp2()); 
      el1_eta[nMC] = nearElectron1->eta(); 
      el1_e[nMC] = nearElectron1->energy();
      el1_phi[nMC] = nearElectron1->phi(); 
      el1_dr[nMC] = dRmin; 
      el1_eopin[nMC] = nearElectron1->eSuperClusterOverP();
      el1_eopout[nMC] = nearElectron1->eSeedClusterOverPout();
      el1_hoe[nMC] = nearElectron1->hadronicOverEm();
      el1_dphiin[nMC] = nearElectron1->deltaPhiSuperClusterTrackAtVtx();
      el1_detain[nMC] = nearElectron1->deltaEtaSuperClusterTrackAtVtx();
      el1_dphiout[nMC] = nearElectron1->deltaPhiSeedClusterTrackAtCalo();
      el1_detaout[nMC] = nearElectron1->deltaEtaSeedClusterTrackAtCalo();
      float pin  = nearElectron1->trackMomentumAtVtx().R();
      float pout = nearElectron1->trackMomentumOut().R();
      el1_pout[nMC] = pout;
      el1_fbrem[nMC] = (pin-pout)/pin;
      el1_class[nMC] = nearElectron1->classification();
      el1_pixelhits[nMC] = nearElectron1->gsfTrack()->hitPattern().numberOfValidPixelHits();
      el1_tkhits[nMC] = nearElectron1->gsfTrack()->hitPattern().numberOfValidTrackerHits();
      el1_gsfpt[nMC] = nearElectron1->gsfTrack()->pt();
      el1_gsfeta[nMC] = nearElectron1->gsfTrack()->eta();
      math::XYZVector trackGlobalDir(nearElectron1->gsfTrack()->momentum());
      reco::GsfTrackRef track = nearElectron1->gsfTrack();
      reco::SuperClusterRef sc = nearElectron1->superCluster(); 
      math::XYZVector clusterDir(sc->x()-track->vx(), sc->y()-track->vy(), sc->z()-track->vz());
      el1_scpt[nMC] = sc->energy()*sin(clusterDir.theta());
      el1_sceta[nMC] = clusterDir.eta(); 
      el1_scphi[nMC] = clusterDir.phi();

      R9_25_gsf(event, &(*nearElectron1), el1_eseed[nMC], el1_e3x3[nMC], el1_e5x5[nMC], el1_spp[nMC], el1_see[nMC]);
    } else {
      if (nMC == 0) { rce1 = math::XYZTLorentzVector(0,0,0,0); }
      if (nMC == 1) { rcp1 = math::XYZTLorentzVector(0,0,0,0); }
      if (nMC == 2) { rce2 = math::XYZTLorentzVector(0,0,0,0); }
      if (nMC == 3) { rcp2 = math::XYZTLorentzVector(0,0,0,0); }
      el1_pt[nMC] = 0.;
      el1_e[nMC] = 0.;
      el1_eta[nMC] = 0.;
      el1_phi[nMC] = 0.;
      el1_dr[nMC] = 0.1; 
      el1_eopin[nMC] = 0.;
      el1_eopout[nMC] = 0.;
      el1_pout[nMC] = 0;
      el1_hoe[nMC] = 0.;
      el1_dphiin[nMC] = 0.;
      el1_detain[nMC] = 0.;
      el1_dphiout[nMC] = 0.;
      el1_detaout[nMC] = 0.;
      el1_fbrem[nMC] = 0.;
      el1_eseed[nMC] = 0.;
      el1_e3x3[nMC] = 0.;
      el1_e5x5[nMC] = 0.;
      el1_spp[nMC] = 0.;
      el1_see[nMC] = 0.;
      el1_class[nMC] = -1;
      el1_pixelhits[nMC] = -1;
      el1_tkhits[nMC] = -1;
      el1_gsfpt[nMC] = 0.;
      el1_gsfeta[nMC] = 0.;
      el1_scpt[nMC] = 0.;
      el1_sceta[nMC] = 0.;
      el1_scphi[nMC] = 0.; 

    }
   }
  

   //calc MZ1, MZ@ & MH for CMS baseline:
   math::XYZTLorentzVector vqz1( qce1 + qcp1 );
   math::XYZTLorentzVector vqz2( qce2 + qcp2 );

   qmz1 = 0.0; qmz2 = 0.0;
   qmz1 = vqz1.mass();
   qmz2 = vqz2.mass();

   math::XYZTLorentzVector vqh = vqz1 + vqz2;

   qmh = 0.0;
   if (qmz1 >0. && qmz2>0.) {qmh = vqh.mass(); }
   

   //calc MZ1, MZ@ & MH for S&T:
   math::XYZTLorentzVector vrz1( rce1 + rcp1 );
   math::XYZTLorentzVector vrz2( rce2 + rcp2 );

  
   rmz1 = 0.0; rmz2 = 0.0;
   rmz1 = vrz1.mass();
   rmz2 = vrz2.mass();

   math::XYZTLorentzVector vrh = vrz1 + vrz2;
   rmh = 0.0;
   if (rmz1 >0. && rmz2>0.) {rmh = vrh.mass(); }


   //Fill tree:
   run = event.id().run();
   id = event.id().event();
   tree->Fill();
   
}








bool hzz::inCrack(float eta) {

  return (eta < 0.018 ||
          (eta>0.423 && eta<0.461) ||
          (eta>0.770 && eta<0.806) ||
          (eta>1.127 && eta<1.163) ||
          (eta>1.460 && eta<1.558));
}

int hzz::mother(HepMC::GenParticle *p) {
  
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


void hzz::R9_25_gsf(const Event & event, const reco::PixelMatchGsfElectron* e,
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
  spp = seedShapeRef->covPhiPhi();
  see = seedShapeRef->covEtaEta();
}
