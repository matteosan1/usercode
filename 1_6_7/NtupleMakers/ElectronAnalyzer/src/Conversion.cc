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

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"


//__________________________added by PDK_______________________
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
//_____________________________________________________________


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
  logFileName         =  pset.getParameter<std::string>("LogFileName");
}

Conversion::~Conversion() {}

void Conversion::beginJob(const EventSetup& eventSetup) {

  logfile.open(logFileName.c_str());
  
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
  tree->Branch("el_q", &el_q, "el_q/I");
  tree->Branch("el_eta", &el_eta, "el_eta/F");
  tree->Branch("el_phi", &el_phi, "el_phi/F");
  tree->Branch("el_dr", &el_dr, "el_dr/F");
  tree->Branch("el_sceta", &el_sceta, "el_sceta/F");
  tree->Branch("el_scphi", &el_scphi, "el_scphi/F");
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
  tree->Branch("el_xhit", &el_xhit, "el_xhit/F");
  tree->Branch("el_yhit", &el_yhit, "el_yhit/F");
  tree->Branch("el_zhit", &el_zhit, "el_zhit/F");
  tree->Branch("el_detinnerhit", &el_detinnerhit, "el_detinnerhit/I");
  tree->Branch("el_layerinnerhit", &el_detinnerhit, "el_layerinnerhit/I");
  tree->Branch("el_rinnerhit", &el_rinnerhit, "el_rinnerhit/F");
  tree->Branch("el_z0", &el_z0, "el_z0/F");
  tree->Branch("el_d0", &el_d0, "el_d0/F");
  tree->Branch("el_d0err", &el_d0err, "el_d0err/F");
  tree->Branch("el_tkiso", &el_tkiso, "el_tkiso/F");
  tree->Branch("el_tkpt", &el_tkpt, "el_tkpt/F");
  tree->Branch("el_tketa", &el_tketa, "el_tketa/F");
  tree->Branch("el_tkphi", &el_tkphi, "el_tkphi/F");
  
  tree->Branch("el1_pt", &el1_pt, "el1_pt/F");
  tree->Branch("el1_e", &el1_e, "el1_e/F");
  tree->Branch("el1_q", &el1_q, "el1_q/I");
  tree->Branch("el1_eta", &el1_eta, "el1_eta/F");
  tree->Branch("el1_phi", &el1_phi, "el1_phi/F");
  tree->Branch("el1_dr", &el1_dr, "el1_dr/F");
  tree->Branch("el1_sceta", &el1_sceta, "el1_sceta/F");
  tree->Branch("el1_scphi", &el1_scphi, "el1_scphi/F");
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
  tree->Branch("el1_xhit", &el1_xhit, "el1_xhit/F");
  tree->Branch("el1_yhit", &el1_yhit, "el1_yhit/F");
  tree->Branch("el1_zhit", &el1_zhit, "el1_zhit/F");
  tree->Branch("el1_detinnerhit", &el1_detinnerhit, "el1_detinnerhit/I");
  tree->Branch("el1_layerinnerhit", &el1_detinnerhit, "el1_layerinnerhit/I");
  tree->Branch("el1_rinnerhit", &el1_rinnerhit, "el1_rinnerhit/F");
  tree->Branch("el1_z0", &el1_z0, "el1_z0/F");
  tree->Branch("el1_d0", &el1_d0, "el1_d0/F");
  tree->Branch("el1_d0err", &el1_d0err, "el1_d0err/F");
  tree->Branch("el1_tkiso", &el1_tkiso, "el1_tkiso/F");
  tree->Branch("el1_tkpt", &el1_tkpt, "el1_tkpt/F");
  tree->Branch("el1_tketa", &el1_tketa, "el1_tketa/F");
  tree->Branch("el1_tkphi", &el1_tkphi, "el1_tkphi/F");
    
  //Sim branches
  //px,py,pz and e of the photon closest to the standard electron's track
  tree->Branch("el_sim_gpx", &el_sim_gpx, "el_sim_gpx/F");
  tree->Branch("el_sim_gpy", &el_sim_gpy, "el_sim_gpy/F");
  tree->Branch("el_sim_gpz", &el_sim_gpz, "el_sim_gpz/F");
  tree->Branch("el_sim_gpt", &el_sim_gpt, "el_sim_gpt/F");
  tree->Branch("el_sim_ge", &el_sim_ge, "el_sim_ge/F");
  //position of the vertex where the conversion happened
  tree->Branch("el_sim_gvx", &el_sim_gvx, "el_sim_gvx/F");
  tree->Branch("el_sim_gvy", &el_sim_gvy, "el_sim_gvy/F");
  tree->Branch("el_sim_gvz", &el_sim_gvz, "el_sim_gvz/F");
  tree->Branch("el_sim_gvr", &el_sim_gvr, "el_sim_gvr/F");
  tree->Branch("el_sim_gveta", &el_sim_gveta, "el_sim_gveta/F");
  tree->Branch("el_sim_gvphi", &el_sim_gvphi, "el_sim_gvphi/F");
  //px,py,pz, e of the 2 decay electrons
  tree->Branch("el_sim_decay1px", &el_sim_decay1px, "el_sim_decay1px/F");
  tree->Branch("el_sim_decay1py", &el_sim_decay1py, "el_sim_decay1py/F");
  tree->Branch("el_sim_decay1pz", &el_sim_decay1pz, "el_sim_decay1pz/F");
  tree->Branch("el_sim_decay1pt", &el_sim_decay1pt, "el_sim_decay1pt/F");
  tree->Branch("el_sim_decay1e",  &el_sim_decay1e,  "el_sim_decay1e/F" );
  tree->Branch("el_sim_decay1pid",  &el_sim_decay1pid,  "el_sim_decay1pid/I" );
  tree->Branch("el_sim_decay2px", &el_sim_decay2px, "el_sim_decay2px/F");
  tree->Branch("el_sim_decay2py", &el_sim_decay2py, "el_sim_decay2py/F");
  tree->Branch("el_sim_decay2pz", &el_sim_decay2pz, "el_sim_decay2pz/F");
  tree->Branch("el_sim_decay2pt", &el_sim_decay2pt, "el_sim_decay2pt/F");
  tree->Branch("el_sim_decay2e",  &el_sim_decay2e,  "el_sim_decay2e/F" );
  tree->Branch("el_sim_decay2pid",  &el_sim_decay2pid,  "el_sim_decay2pid/I" );

  //UCSD electron
  tree->Branch("el1_sim_gpx", &el1_sim_gpx, "el1_sim_gpx/F");
  tree->Branch("el1_sim_gpy", &el1_sim_gpy, "el1_sim_gpy/F");
  tree->Branch("el1_sim_gpz", &el1_sim_gpz, "el1_sim_gpz/F");
  tree->Branch("el1_sim_gpt", &el1_sim_gpt, "el1_sim_gpt/F");
  tree->Branch("el1_sim_ge", &el1_sim_ge, "el1_sim_ge/F");
  //position of the vertex where the conversion happened
  tree->Branch("el1_sim_gvx", &el1_sim_gvx, "el1_sim_gvx/F");
  tree->Branch("el1_sim_gvy", &el1_sim_gvy, "el1_sim_gvy/F");
  tree->Branch("el1_sim_gvz", &el1_sim_gvz, "el1_sim_gvz/F");
  tree->Branch("el1_sim_gvr", &el1_sim_gvr, "el1_sim_gvr/F");
  tree->Branch("el1_sim_gveta", &el1_sim_gveta, "el1_sim_gveta/F");
  tree->Branch("el1_sim_gvphi", &el1_sim_gvphi, "el1_sim_gvphi/F");
  //px,py,pz, e of the 2 decay electrons
  tree->Branch("el1_sim_decay1px", &el1_sim_decay1px, "el1_sim_decay1px/F");
  tree->Branch("el1_sim_decay1py", &el1_sim_decay1py, "el1_sim_decay1py/F");
  tree->Branch("el1_sim_decay1pz", &el1_sim_decay1pz, "el1_sim_decay1pz/F");
  tree->Branch("el1_sim_decay1pt", &el1_sim_decay1pt, "el1_sim_decay1pt/F");
  tree->Branch("el1_sim_decay1e",  &el1_sim_decay1e,  "el1_sim_decay1e/F" );
  tree->Branch("el1_sim_decay1pid",  &el1_sim_decay1pid,  "el1_sim_decay1pid/I" );
  tree->Branch("el1_sim_decay2px", &el1_sim_decay2px, "el1_sim_decay2px/F");
  tree->Branch("el1_sim_decay2py", &el1_sim_decay2py, "el1_sim_decay2py/F");
  tree->Branch("el1_sim_decay2pz", &el1_sim_decay2pz, "el1_sim_decay2pz/F");
  tree->Branch("el1_sim_decay2pt", &el1_sim_decay2pt, "el1_sim_decay2pt/F");
  tree->Branch("el1_sim_decay2e",  &el1_sim_decay2e,  "el1_sim_decay2e/F" );
  tree->Branch("el1_sim_decay2pid",  &el1_sim_decay2pid,  "el1_sim_decay2pid/I" );

  //px,py,pz and e of the photon closest to the sc
  tree->Branch("sc_sim_gpx", &sc_sim_gpx, "sc_sim_gpx/F");
  tree->Branch("sc_sim_gpy", &sc_sim_gpy, "sc_sim_gpy/F");
  tree->Branch("sc_sim_gpz", &sc_sim_gpz, "sc_sim_gpz/F");
  tree->Branch("sc_sim_gpt", &sc_sim_gpt, "sc_sim_gpt/F");
  tree->Branch("sc_sim_ge", &sc_sim_ge, "sc_sim_ge/F");
  //position of the vertex where the conversion happened
  tree->Branch("sc_sim_gvx", &sc_sim_gvx, "sc_sim_gvx/F");
  tree->Branch("sc_sim_gvy", &sc_sim_gvy, "sc_sim_gvy/F");
  tree->Branch("sc_sim_gvz", &sc_sim_gvz, "sc_sim_gvz/F");
  tree->Branch("sc_sim_gvr", &sc_sim_gvr, "sc_sim_gvr/F");
  tree->Branch("sc_sim_gveta", &sc_sim_gveta, "sc_sim_gveta/F");
  tree->Branch("sc_sim_gvphi", &sc_sim_gvphi, "sc_sim_gvphi/F");
  //px,py,pz, e of the 2 decay electrons
  tree->Branch("sc_sim_decay1px", &sc_sim_decay1px, "sc_sim_decay1px/F");
  tree->Branch("sc_sim_decay1py", &sc_sim_decay1py, "sc_sim_decay1py/F");
  tree->Branch("sc_sim_decay1pz", &sc_sim_decay1pz, "sc_sim_decay1pz/F");
  tree->Branch("sc_sim_decay1pt", &sc_sim_decay1pt, "sc_sim_decay1pt/F");
  tree->Branch("sc_sim_decay1e",  &sc_sim_decay1e,  "sc_sim_decay1e/F" );
  tree->Branch("sc_sim_decay1pid",  &sc_sim_decay1pid,  "sc_sim_decay1pid/I" );
  tree->Branch("sc_sim_decay2px", &sc_sim_decay2px, "sc_sim_decay2px/F");
  tree->Branch("sc_sim_decay2py", &sc_sim_decay2py, "sc_sim_decay2py/F");
  tree->Branch("sc_sim_decay2pz", &sc_sim_decay2pz, "sc_sim_decay2pz/F");
  tree->Branch("sc_sim_decay2pt", &sc_sim_decay2pt, "sc_sim_decay2pt/F");
  tree->Branch("sc_sim_decay2e",  &sc_sim_decay2e,  "sc_sim_decay2e/F" );
  tree->Branch("sc_sim_decay2pid",  &sc_sim_decay2pid,  "sc_sim_decay2pid/I");  


  
  //converted photon info - std electron
  tree->Branch("el_con_flag", &el_con_flag, "el_con_flag/I");
  tree->Branch("el_con_numCP", &el_con_numCP, "el_con_numCP/I");
  tree->Branch("el_con_pairinvmass", &el_con_pairinvmass, "el_con_pairinvmass/F");
  tree->Branch("el_con_pairdcottheta", &el_con_pairdcottheta, "el_con_pairdcottheta/F");
  tree->Branch("el_con_pairp", &el_con_pairp, "el_con_pairp/F");
  tree->Branch("el_con_pairphi", &el_con_pairphi, "el_con_pairphi/F");
  tree->Branch("el_con_paireta", &el_con_paireta, "el_con_paireta/F");
  tree->Branch("el_con_elcottheta", &el_con_elcottheta, "el_con_elcottheta/F");
  tree->Branch("el_con_elcottheta_2", &el_con_elcottheta_2, "el_con_elcottheta_2/F");
  tree->Branch("el_con_vx", &el_con_vx, "el_con_vx/F");
  tree->Branch("el_con_vy", &el_con_vy, "el_con_vy/F");
  tree->Branch("el_con_vz", &el_con_vz, "el_con_vz/F");
  tree->Branch("el_con_vphi", &el_con_vphi, "el_con_vphi/F");  
  tree->Branch("el_con_veta", &el_con_veta, "el_con_veta/F");
  tree->Branch("el_con_vr", &el_con_vr, "el_con_vr/F");
  tree->Branch("el_con_decay1px", &el_con_decay1px, "el_con_decay1px/F");
  tree->Branch("el_con_decay1py", &el_con_decay1py, "el_con_decay1py/F");
  tree->Branch("el_con_decay1pz", &el_con_decay1pz, "el_con_decay1pz/F");
  tree->Branch("el_con_decay1pt", &el_con_decay1pt, "el_con_decay1pt/F");
  tree->Branch("el_con_decay1q", &el_con_decay1q, "el_con_decay1q/I");
  tree->Branch("el_con_decay2px", &el_con_decay2px, "el_con_decay2px/F");
  tree->Branch("el_con_decay2py", &el_con_decay2py, "el_con_decay2py/F");
  tree->Branch("el_con_decay2pz", &el_con_decay2pz, "el_con_decay2pz/F");
  tree->Branch("el_con_decay2pt", &el_con_decay2pt, "el_con_decay2pt/F");
  tree->Branch("el_con_decay2q", &el_con_decay2q, "el_con_decay2q/I");
  tree->Branch("el_con_decay1numhits",  &el_con_decay1numhits,  "el_con_decay1numhits/I");
  tree->Branch("el_con_decay1numsharedhits",  &el_con_decay1numsharedhits,  "el_con_decay1numsharedhits/I");
  tree->Branch("el_con_decay1fracsharedhits",  &el_con_decay1fracsharedhits,  "el_con_decay1fracsharedhits/F");
  tree->Branch("el_con_decay2numhits",  &el_con_decay2numhits,  "el_con_decay2numhits/I");
  tree->Branch("el_con_decay2numsharedhits",  &el_con_decay2numsharedhits,  "el_con_decay2numsharedhits/I");
  tree->Branch("el_con_decay2fracsharedhits",  &el_con_decay2fracsharedhits,  "el_con_decay2fracsharedhits/F");

  //UCSD electron1
  tree->Branch("el1_con_flag", &el1_con_flag, "el1_con_flag/I");
  tree->Branch("el1_con_numCP", &el1_con_numCP, "el1_con_numCP/I");
  tree->Branch("el1_con_pairinvmass", &el1_con_pairinvmass, "el1_con_pairinvmass/F");
  tree->Branch("el1_con_pairdcottheta", &el1_con_pairdcottheta, "el1_con_pairdcottheta/F");
  tree->Branch("el1_con_pairp", &el1_con_pairp, "el1_con_pairp/F");
  tree->Branch("el1_con_pairphi", &el1_con_pairphi, "el1_con_pairphi/F");
  tree->Branch("el1_con_paireta", &el1_con_paireta, "el1_con_paireta/F");
  tree->Branch("el1_con_el1cottheta", &el1_con_el1cottheta, "el1_con_el1cottheta/F");
  tree->Branch("el1_con_el1cottheta_2", &el1_con_el1cottheta_2, "el1_con_el1cottheta_2/F");
  tree->Branch("el1_con_vx", &el1_con_vx, "el1_con_vx/F");
  tree->Branch("el1_con_vy", &el1_con_vy, "el1_con_vy/F");
  tree->Branch("el1_con_vz", &el1_con_vz, "el1_con_vz/F");
  tree->Branch("el1_con_vphi", &el1_con_vphi, "el1_con_vphi/F");  
  tree->Branch("el1_con_veta", &el1_con_veta, "el1_con_veta/F");
  tree->Branch("el1_con_vr", &el1_con_vr, "el1_con_vr/F");
  tree->Branch("el1_con_decay1px", &el1_con_decay1px, "el1_con_decay1px/F");
  tree->Branch("el1_con_decay1py", &el1_con_decay1py, "el1_con_decay1py/F");
  tree->Branch("el1_con_decay1pz", &el1_con_decay1pz, "el1_con_decay1pz/F");
  tree->Branch("el1_con_decay1pt", &el1_con_decay1pt, "el1_con_decay1pt/F");
  tree->Branch("el1_con_decay1q", &el1_con_decay1q, "el1_con_decay1q/I");
  tree->Branch("el1_con_decay2px", &el1_con_decay2px, "el1_con_decay2px/F");
  tree->Branch("el1_con_decay2py", &el1_con_decay2py, "el1_con_decay2py/F");
  tree->Branch("el1_con_decay2pz", &el1_con_decay2pz, "el1_con_decay2pz/F");
  tree->Branch("el1_con_decay2pt", &el1_con_decay2pt, "el1_con_decay2pt/F");
  tree->Branch("el1_con_decay2q", &el1_con_decay2q, "el1_con_decay2q/I");
  tree->Branch("el1_con_decay1numhits",  &el1_con_decay1numhits,  "el1_con_decay1numhits/I");
  tree->Branch("el1_con_decay1numsharedhits",  &el1_con_decay1numsharedhits,  "el1_con_decay1numsharedhits/I");
  tree->Branch("el1_con_decay1fracsharedhits",  &el1_con_decay1fracsharedhits,  "el1_con_decay1fracsharedhits/F");
  tree->Branch("el1_con_decay2numhits",  &el1_con_decay2numhits,  "el1_con_decay2numhits/I");
  tree->Branch("el1_con_decay2numsharedhits",  &el1_con_decay2numsharedhits,  "el1_con_decay2numsharedhits/I");
  tree->Branch("el1_con_decay2fracsharedhits",  &el1_con_decay2fracsharedhits,  "el1_con_decay2fracsharedhits/F");
  
  
  //SC
  tree->Branch("sc_con_flag", &sc_con_flag, "sc_con_flag/I");
  tree->Branch("sc_con_numCP", &sc_con_numCP, "sc_con_numCP/I");
  tree->Branch("sc_con_pairinvmass", &sc_con_pairinvmass, "sc_con_pairinvmass/F");
  tree->Branch("sc_con_pairdcottheta", &sc_con_pairdcottheta, "sc_con_pairdcottheta/F");
  tree->Branch("sc_con_pairp", &sc_con_pairp, "sc_con_pairp/F");
  tree->Branch("sc_con_pairphi", &sc_con_pairphi, "sc_con_pairphi/F");
  tree->Branch("sc_con_paireta", &sc_con_paireta, "sc_con_paireta/F");
  tree->Branch("sc_con_vx", &sc_con_vx, "sc_con_vx/F");
  tree->Branch("sc_con_vy", &sc_con_vy, "sc_con_vy/F");
  tree->Branch("sc_con_vz", &sc_con_vz, "sc_con_vz/F");
  tree->Branch("sc_con_vphi", &sc_con_vphi, "sc_con_vphi/F");  
  tree->Branch("sc_con_veta", &sc_con_veta, "sc_con_veta/F");
  tree->Branch("sc_con_vr", &sc_con_vr, "sc_con_vr/F");
  tree->Branch("sc_con_decay1px", &sc_con_decay1px, "sc_con_decay1px/F");
  tree->Branch("sc_con_decay1py", &sc_con_decay1py, "sc_con_decay1py/F");
  tree->Branch("sc_con_decay1pz", &sc_con_decay1pz, "sc_con_decay1pz/F");
  tree->Branch("sc_con_decay1pt", &sc_con_decay1pt, "sc_con_decay1pt/F");
  tree->Branch("sc_con_decay1q", &sc_con_decay1q, "sc_con_decay1q/I");
  tree->Branch("sc_con_decay2px", &sc_con_decay2px, "sc_con_decay2px/F");
  tree->Branch("sc_con_decay2py", &sc_con_decay2py, "sc_con_decay2py/F");
  tree->Branch("sc_con_decay2pz", &sc_con_decay2pz, "sc_con_decay2pz/F");
  tree->Branch("sc_con_decay2pt", &sc_con_decay2pt, "sc_con_decay2pt/F");
  tree->Branch("sc_con_decay2q", &sc_con_decay2q, "sc_con_decay2q/I");
  tree->Branch("sc_con_decay1numhits",  &sc_con_decay1numhits,  "sc_con_decaynum1hits/F");
  tree->Branch("sc_con_decay2numhits",  &sc_con_decay2numhits,  "sc_con_decay2numhits/F");
    
  
  
  
}

void Conversion::endJob() {
  
  file->Write();
  file->Close();
  
}

void Conversion::analyze(const Event & event, const EventSetup& eventSetup) {

  
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

  //  Handle<SuperClusterCollection> sch1;
  //   event.getByLabel("hybridSuperClusters", sch1);
  //   const SuperClusterCollection* scb = sch1.product();

  //   Handle<SuperClusterCollection> sch2;
  //   event.getByLabel("islandSuperClusters", "islandEndcapSuperClusters", sch2);
  //   const SuperClusterCollection* sce = sch2.product();
  //   SuperClusterCollection::const_iterator itsc, itscb, itsce;
  //   SuperClusterCollection sc;
  //   sc.insert(sc.end(), scb->begin(), scb->end());
  //   sc.insert(sc.end(), sce->begin(), sce->end());


  Handle<SuperClusterCollection> sch1;
  event.getByLabel("correctedHybridSuperClusters", sch1);
  const SuperClusterCollection* scb = sch1.product();

  Handle<SuperClusterCollection> sch2;
  event.getByLabel("correctedEndcapSuperClustersWithPreshower", sch2);
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


  
  //__________________________added by PDK______________________
  //Sim Track Collection
  Handle<SimTrackContainer> simTrackCollection;
  event.getByType(simTrackCollection);
  if(!simTrackCollection.isValid()) throw cms::Exception("FatalError") << "SimTrackCollection not found!\n";
  vector<SimTrack> theSimTrks;
  theSimTrks.insert(theSimTrks.end(), 
		    simTrackCollection->begin(),
		    simTrackCollection->end());
  
  Handle<SimVertexContainer> simVertexCollection;
  event.getByType(simVertexCollection);
  if(!simVertexCollection.isValid()) throw cms::Exception("FatalError") << "SimVertexCollection not found!\n";
  vector<SimVertex> theSimVerts;
  theSimVerts.insert(theSimVerts.end(),
		     simVertexCollection->begin(),
		     simVertexCollection->end());

  InputTag blah("convertedPhotons", "", "TEST");
  Handle<ConvertedPhotonCollection> convertedPhotonCollection;
  event.getByLabel(blah, convertedPhotonCollection);
  if(!convertedPhotonCollection.isValid()) throw cms::Exception("FatalError") << "ConvertedPhotonCollection not found!\n";
  vector<ConvertedPhoton> theConPhotons;
  theConPhotons.insert(theConPhotons.end(),
                       convertedPhotonCollection->begin(),
                       convertedPhotonCollection->end());

  //____________________________________________________________


  
  logfile << endl << "Run: " << event.id().run() 
	  << " Event: " << event.id().event() << endl;
  
  //   logfile << "Size of corrected SC collection: " << sc.size() 
  // 	  << "  Size of pre-existing ConvertedPhoton Coll: " 
  // 	  << theConPhotons.size() << endl;
  //<< "   Size of pre-existing collection: " << theConPhotons2.size() << endl;
  
  
  for (itsc = sc.begin(); itsc != sc.end(); ++itsc) {
    
    int mctk_index=0;
    int mctk1_index=0;
    int mcsc_index=0;

    
    math::XYZVector scv(itsc->x(), itsc->y(), itsc->z());
    
    //save SC with Et > 5
    if (sin(scv.Theta())*itsc->energy() > 5.) {
      
      sc_e = itsc->energy(); 
      //logfile << "New SC E: " << sc_e << endl;
      sc_rawe = itsc->rawEnergy();
      sc_et = sin(scv.Theta())*itsc->energy();
      sc_eta = itsc->eta();
      sc_phi = itsc->phi(); 
      
      if (fabs(itsc->eta()) < 1.47)
        sc_type = 0; //barrel
      else
        sc_type = 1; //endcap
      
      //match SC to nearest Gen particle by dR
      double dR, dRmin = 0.1;
      HepMC::GenEvent::particle_const_iterator nearMC;
      int i=1;
      for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it,++i) { 
        
        if ((*it)->status() == 1) {      
          math::XYZVector mcv((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz());
          dR = ROOT::Math::VectorUtil::DeltaR(scv, mcv);
	  if (dR < dRmin) {
            dRmin = dR;
	    nearMC = it;
            mcsc_index = i;
          }
        }
      }//gen partcle loop

      //save Gen particle info
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
      //match standard electron to SC
      dRmin = 0.01;
      PixelMatchGsfElectronCollection::const_iterator nearElectron;	
      for(ite = electrons.begin(); ite != electrons.end(); ++ite) {
	double deta = (ite->superCluster())->eta() - itsc->eta();
	double dphi = (ite->superCluster())->phi() - itsc->phi();		
	if(fabs(dphi) > TMath::Pi() ) dphi = 2*TMath::Pi() - fabs(dphi);
        //dR = ROOT::Math::VectorUtil::DeltaR(ite->p4(), scv);
	dR = TMath::Sqrt(deta*deta+dphi*dphi);
	if (dR < dRmin) {
          dRmin = dR;
          nearElectron = ite;
        }
      }
      //logfile << "dRmin: " << dRmin << endl;
      //logfile << "El pt" << nearElectron->pt() << "El eta, phi: " << nearElectron->eta() << " , " << nearElectron->phi() << endl;
      // strore info about standard Electron if it passes matching requirement
      if (dRmin < 0.01) {
        el_pt = nearElectron->pt(); 
        el_eta = nearElectron->eta(); 
        el_e = nearElectron->energy();
        el_q = nearElectron->charge();
	el_sceta = (nearElectron->superCluster())->eta();
	el_scphi = (nearElectron->superCluster())->phi();
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
            el_xhit = hitPosition.x();
            el_yhit = hitPosition.y();
            el_zhit = hitPosition.z();
            el_rinnerhit = pos.perp(); // + pow(pos.z(),2));
            subDetector(hit, a, b);
            el_detinnerhit = a;
            el_layerinnerhit = b;
            break;
          }
          index++;
        }

        el_z0 = nearElectron->gsfTrack()->vz();
        el_d0 = nearElectron->gsfTrack()->d0();
        el_d0err = nearElectron->gsfTrack()->d0Error();
        el_tkiso = trackIsolation(nearElectron->trackMomentumAtVtx(), nearElectron->vertex(), tracks);
        el_tkpt = nearElectron->gsfTrack()->pt();
        el_tketa = nearElectron->gsfTrack()->eta(); 
        el_tkphi = nearElectron->gsfTrack()->phi();  

        double dR, dRmin = 0.1;
	double detamin = 0.1;
	HepMC::GenEvent::particle_const_iterator nearMC;
	int i=1;
        for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it, ++i) { 
          
          if ((*it)->status() == 1) {      
            math::XYZVector mcv((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz());
            dR = ROOT::Math::VectorUtil::DeltaR(nearElectron->gsfTrack()->innerMomentum(), mcv);
	    
            if (dR < dRmin) {
	      detamin = fabs(nearElectron->gsfTrack()->eta() - (*it)->momentum().eta());
              dRmin = dR;
              nearMC = it;
	      mctk_index = i;
            }
          }
        }
        
	//if MC particle is a photon, go back and look to see if theres another gen photon whose dRmin is 
	//also < 0.05, but whose deta is smaller than the deta of the closest (in dR) MC particle
	if(dRmin < 0.05 && (*nearMC)->pdg_id() == 22) {
	  i=1;
	  for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it, ++i) { 
	    
	    if ((*it)->status() == 1) {      
	      math::XYZVector mcv((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz());
	      dR = ROOT::Math::VectorUtil::DeltaR(nearElectron->gsfTrack()->innerMomentum(), mcv);
	      
	      if((*it)->pdg_id() == 22
		 && dR < 0.05
		 && fabs(nearElectron->gsfTrack()->eta() - (*it)->momentum().eta()) < detamin ) {
		//logfile << "Found a better gen photon matched to std track" << endl; 
		detamin = fabs(nearElectron->gsfTrack()->eta() - (*it)->momentum().eta());
		dRmin = dR;
		nearMC = it;
		mctk_index = i;
	      }
	    }
	  }
	} //if(dRmin < 0.1 && (*nearMC)->pdg_id == 22)
	
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
	el_q = 0;
        el_eta = 0.;
        el_phi = 0.;
	el_sceta = 0.;
	el_scphi = 0.;
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
        el_d0 = -1;
        el_d0err = -1;
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
      dRmin = 0.01;
      for(ite1 = electrons1->begin(); ite1 != electrons1->end(); ++ite1) {
	double deta = (ite1->superCluster())->eta() - itsc->eta();
	double dphi = (ite1->superCluster())->phi() - itsc->phi();		
	if(fabs(dphi) > TMath::Pi() ) dphi = 2*TMath::Pi() - fabs(dphi);
        //dR = ROOT::Math::VectorUtil::DeltaR(ite->p4(), scv);
	dR = TMath::Sqrt(deta*deta+dphi*dphi);
	if (dR < dRmin) {
          dRmin = dR;
          nearElectron1 = ite1;
        }
      }
      
      // strore info about Ele
      if (dRmin < 0.01) {
        el1_pt = nearElectron1->pt(); 
        el1_eta = nearElectron1->eta(); 
        el1_e = nearElectron1->energy();
	el1_q = nearElectron1->charge();
        el1_phi = nearElectron1->phi(); 
	el1_sceta = (nearElectron1->superCluster())->eta();
	el1_scphi = (nearElectron1->superCluster())->phi();
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
            el1_xhit = hitPosition.x();
            el1_yhit = hitPosition.y();
            el1_zhit = hitPosition.z();
            el1_rinnerhit = pos.perp(); // + pow(pos.z(),2));
            subDetector(hit, a, b);
            el1_detinnerhit = a;
            el1_layerinnerhit = b;
            break;
          }
          index++;
        }

        el1_z0 = nearElectron1->gsfTrack()->vz();
        el1_d0 = nearElectron1->gsfTrack()->d0();
        el1_d0err = nearElectron1->gsfTrack()->d0Error();
        el1_tkiso = trackIsolation(nearElectron1->trackMomentumAtVtx(), nearElectron1->vertex(), tracks);
        el1_tkpt = nearElectron1->gsfTrack()->pt();
        el1_tketa = nearElectron1->gsfTrack()->eta(); 
        el1_tkphi = nearElectron1->gsfTrack()->phi(); 

        double dR, dRmin = 0.1;
	double detamin = 0.1;
        i=1;
        HepMC::GenEvent::particle_const_iterator nearMC;
        for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it, ++i) { 
          
          if ((*it)->status() == 1) {      
            math::XYZVector mcv((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz());
            dR = ROOT::Math::VectorUtil::DeltaR(nearElectron1->gsfTrack()->innerMomentum(), mcv);
            
            if (dR < dRmin) {
              dRmin = dR;
              nearMC = it;
              mctk1_index = i;
            }
          }
        }
        

	//if MC particle is a photon, go back and look to see if theres another gen photon whose dRmin is 
	//also < 0.05, but whose deta is smaller than the deta of the closest (in dR) MC particle
	if(dRmin < 0.05 && (*nearMC)->pdg_id() == 22) {
	  i=1;
	  for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it, ++i) { 
	    
	    if ((*it)->status() == 1) {      
	      math::XYZVector mcv((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz());
	      dR = ROOT::Math::VectorUtil::DeltaR(nearElectron1->gsfTrack()->innerMomentum(), mcv);
	      
	      if((*it)->pdg_id() == 22
		 && dR < 0.05
		 && fabs(nearElectron1->gsfTrack()->eta() - (*it)->momentum().eta()) < detamin) {
		//logfile << "Found a better gen photon matched to UCSD track" << endl; 
		detamin = fabs(nearElectron1->gsfTrack()->eta() - (*it)->momentum().eta()); 
		dRmin = dR;
		nearMC = it;
		mctk1_index = i;
	      }
	    }
	  }
	} //if(dRmin < 0.1 && (*nearMC)->pdg_id == 22)


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
        el1_q = 0;
        el1_eta = 0.;
        el1_phi = 0.;
	el1_sceta = 0.;
	el1_scphi = 0.;
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
        el1_d0 = -1;
        el1_d0err = -1;
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
            
      InitSimVariables();
      
      //if there is a standard electron or a UCSD electron, fill SC info
      //now do the same for the MC particle matched closest to the SC
      if(el_dr < 0.1 || el1_dr < 0.1) {
	vector<SimTrack>::const_iterator theSimTrksIter;
	int SimIndex = 0;
	int photonTrackId = 0;
	
	for(theSimTrksIter = theSimTrks.begin();
	    theSimTrksIter!= theSimTrks.end();
	    theSimTrksIter++, SimIndex++) {
	  
	  //index of genparticle closest to the SC
	  if(theSimTrksIter->genpartIndex() == mcsc_index) { 
	
	    int i = 1;
	    HepMC::GenEvent::particle_const_iterator geniter;
	    for(geniter = myGenEvent->particles_begin(); i < mcsc_index; ++i, ++geniter) {}
	    
	    //Making sure that the Sim and Gen particle are both Photons
	    if(theSimTrksIter->type() == 22 && (*geniter)->pdg_id()==22) {
	      photonTrackId = theSimTrksIter->trackId();
	      
	      HepLorentzVector psimvect = theSimTrksIter->momentum();
	      sc_sim_gpx = psimvect.px();
	      sc_sim_gpy = psimvect.py();
	      sc_sim_gpz = psimvect.pz();
	      sc_sim_gpt = TMath::Sqrt(sc_sim_gpx*sc_sim_gpx+sc_sim_gpy*sc_sim_gpy);
	      sc_sim_ge = psimvect.e();
	    	      
	      
	      //SimVertex indicies start at 0
	      i=0;
	      vector<SimVertex>::const_iterator theSimVertsIter;
	      for(theSimVertsIter = theSimVerts.begin();
		  theSimVertsIter != theSimVerts.end();
		  theSimVertsIter++, i++) {
		
		//check to see if the parent of this vertex was the photon in question
		if(theSimVertsIter->parentIndex() == photonTrackId) {
		  		  
		  sc_sim_gvx = theSimVertsIter->position().x();
		  sc_sim_gvy = theSimVertsIter->position().y();
		  sc_sim_gvz = theSimVertsIter->position().z();
		  sc_sim_gvr = TMath::Sqrt(sc_sim_gvx*sc_sim_gvx+sc_sim_gvy*sc_sim_gvy);
		  sc_sim_gvphi = atan2(sc_sim_gvy, sc_sim_gvx);
		  //make sure that we don't get Log of 0
		  if(fabs(sc_sim_gvz) > 1E-7 && sc_sim_gvr > 1E-7) {
		    sc_sim_gveta = TMath::Log( TMath::Tan( 0.5*TMath::ATan(sc_sim_gvr/fabs(sc_sim_gvz) ) ) );
		    if(sc_sim_gvz < 0) 
		      sc_sim_gveta = -fabs(sc_sim_gvz);
		  }else sc_sim_gveta = -999.;
		  
		  
		  /*loop over the simtracks (ugh...again) to see if there are
		    any tracks which have this vertex
		  */
		  for(unsigned int k = 0; k < theSimTrks.size(); k++) {
		    if(theSimTrks.at(k).vertIndex() == i) {
		      
		      float charge = theSimTrks.at(k).charge();
		      if(el_dr < 0.1) {
			if(el_q/charge > 0 ) {
			  sc_sim_decay1px = theSimTrks.at(k).momentum().x();
			  sc_sim_decay1py = theSimTrks.at(k).momentum().y();
			  sc_sim_decay1pz = theSimTrks.at(k).momentum().z();
			  sc_sim_decay1pt = TMath::Sqrt(sc_sim_decay1px*sc_sim_decay1px+
							sc_sim_decay1py*sc_sim_decay1py);
			  sc_sim_decay1e = theSimTrks.at(k).momentum().t();
			  sc_sim_decay1pid = theSimTrks.at(k).type();
			} else {
			  sc_sim_decay2px = theSimTrks.at(k).momentum().x();
			  sc_sim_decay2py = theSimTrks.at(k).momentum().y();
			  sc_sim_decay2pz = theSimTrks.at(k).momentum().z();
			  sc_sim_decay2pt = TMath::Sqrt(sc_sim_decay1px*sc_sim_decay1px+
							sc_sim_decay1py*sc_sim_decay1py);
			  sc_sim_decay2e = theSimTrks.at(k).momentum().t();
			  sc_sim_decay2pid = theSimTrks.at(k).type();
			}
		      }  
		      
		      if(el1_dr < 0.1) {
			if(el1_q/charge > 0 ) {
			  sc_sim_decay1px = theSimTrks.at(k).momentum().x();
			  sc_sim_decay1py = theSimTrks.at(k).momentum().y();
			  sc_sim_decay1pz = theSimTrks.at(k).momentum().z();
			  sc_sim_decay1pt = TMath::Sqrt(sc_sim_decay1px*sc_sim_decay1px+
							sc_sim_decay1py*sc_sim_decay1py);
			  sc_sim_decay1e = theSimTrks.at(k).momentum().t();
			  sc_sim_decay1pid = theSimTrks.at(k).type();
			} else {
			  sc_sim_decay2px = theSimTrks.at(k).momentum().x();
			  sc_sim_decay2py = theSimTrks.at(k).momentum().y();
			  sc_sim_decay2pz = theSimTrks.at(k).momentum().z();
			  sc_sim_decay2pt = TMath::Sqrt(sc_sim_decay1px*sc_sim_decay1px+
							sc_sim_decay1py*sc_sim_decay1py);
			  sc_sim_decay2e = theSimTrks.at(k).momentum().t();
			  sc_sim_decay2pid = theSimTrks.at(k).type();
			}
		      }    
    
		    }
		  }//sim track loop 
		}//is the parent of the vertex a photon?
	      }//SimVerteces Iterator
	    }//if(theSimTrksIter->type() == 22 && (*geniter)->pdg_id()==22)  
	  }//if(theSimTrksIter->genpartIndex() == mcsc_index) {
	}//SimTracksIterator
      }//if(el_dr < 0.1 || el1_dr < 0.1) 

      
      //did we find a standard electron?
      if(el_dr < 0.1) {//existance of electron
	
	vector<SimTrack>::const_iterator theSimTrksIter;
	int SimIndex = 0;
	int photonTrackId = 0;
	for(theSimTrksIter = theSimTrks.begin();
	    theSimTrksIter!= theSimTrks.end();
	    theSimTrksIter++, SimIndex++) {
	  
	  //is the MC particle closest to the CTF track matched to the SimTrack?
	  if(theSimTrksIter->genpartIndex() == mctk_index) { 
	    
	    int i = 1;
	    HepMC::GenEvent::particle_const_iterator geniter;
	    for(geniter = myGenEvent->particles_begin(); i < mctk1_index; ++i, ++geniter) {}
	    
	    //Making sure that the Sim and Gen particle are both Photons
	    if(theSimTrksIter->type() == 22 && (*geniter)->pdg_id()==22) {
	      photonTrackId = theSimTrksIter->trackId();
	      
	      HepLorentzVector psimvect = theSimTrksIter->momentum();
	      el_sim_gpx = psimvect.px();
	      el_sim_gpy = psimvect.py();
	      el_sim_gpz = psimvect.pz();
	      el_sim_gpt = TMath::Sqrt(el_sim_gpx*el_sim_gpx+el_sim_gpy*el_sim_gpy);
	      el_sim_ge = psimvect.e();
	      
	      
	      //SimVertex indicies start at 0
	      i=0;
	      vector<SimVertex>::const_iterator theSimVertsIter;
	      for(theSimVertsIter = theSimVerts.begin();
		  theSimVertsIter != theSimVerts.end();
		  theSimVertsIter++, i++) {
		
		//check to see if the parent of this vertex was the photon in question
		if(theSimVertsIter->parentIndex() == photonTrackId) {
		  		  
		  el_sim_gvx = theSimVertsIter->position().x();
		  el_sim_gvy = theSimVertsIter->position().y();
		  el_sim_gvz = theSimVertsIter->position().z();
		  el_sim_gvr = TMath::Sqrt(el_sim_gvx*el_sim_gvx+el_sim_gvy*el_sim_gvy);
		  el_sim_gvphi = atan2(el_sim_gvy, el_sim_gvx);
		  //make sure that we don't get Log of 0
		  if(fabs(el_sim_gvz) > 1E-7 && el_sim_gvr > 1E-7) {
		    el_sim_gveta = TMath::Log( TMath::Tan( 0.5*TMath::ATan(el_sim_gvr/fabs(el_sim_gvz) ) ) );
		    if(el_sim_gvz < 0) 
		      el_sim_gveta  = -fabs(el_sim_gveta);
		  }else el_sim_gveta = -999.;

		  
		  /*loop over the simtracks (ugh...again) to see if there are
		    any tracks which have this vertex
		  */
		  for(unsigned int k = 0; k < theSimTrks.size(); k++) {
		    if(theSimTrks.at(k).vertIndex() == i) {
		      
		      float charge = theSimTrks.at(k).charge();
		      if(el1_q/charge > 0 ) {
			el_sim_decay1px = theSimTrks.at(k).momentum().x();
			el_sim_decay1py = theSimTrks.at(k).momentum().y();
			el_sim_decay1pz = theSimTrks.at(k).momentum().z();
			el_sim_decay1pt = TMath::Sqrt(el_sim_decay1px*el_sim_decay1px+
						      el_sim_decay1py*el_sim_decay1py);
			el_sim_decay1e = theSimTrks.at(k).momentum().t();
			el_sim_decay1pid = theSimTrks.at(k).type();
		      } else {
			el_sim_decay2px = theSimTrks.at(k).momentum().x();
			el_sim_decay2py = theSimTrks.at(k).momentum().y();
			el_sim_decay2pz = theSimTrks.at(k).momentum().z();
			el_sim_decay2pt = TMath::Sqrt(el_sim_decay2px*el_sim_decay2px+
						      el_sim_decay2py*el_sim_decay2py);
			el_sim_decay2e = theSimTrks.at(k).momentum().t();
			el_sim_decay2pid = theSimTrks.at(k).type();
		      }
		      
		    }
		  }//sim track loop 
		}//is the parent of the vertex a photon?
	      }//SimVerteces Iterator
	    }//if(theSimTrksIter->type() == 22 && (*geniter)->pdg_id()==22)  
	  }//if(theSimTrksIter->genpartIndex() == mctk1_index) {
	}//SimTracksIterator
      }

      
      //is there a UCSD electron?
      if(el1_dr < 0.1) {
	
	vector<SimTrack>::const_iterator theSimTrksIter;
	int SimIndex = 0;
	int photonTrackId = 0;
	for(theSimTrksIter = theSimTrks.begin();
	    theSimTrksIter!= theSimTrks.end();
	    theSimTrksIter++, SimIndex++) {
	  
	  //is the MC particle closest to the CTF track matched to the SimTrack?
	  if(theSimTrksIter->genpartIndex() == mctk1_index) { 
	    
	    int i = 1;
	    HepMC::GenEvent::particle_const_iterator geniter;
	    for(geniter = myGenEvent->particles_begin(); i < mctk1_index; ++i, ++geniter) {}
	    
	    //Making sure that the Sim and Gen particle are both Photons
	    if(theSimTrksIter->type() == 22 && (*geniter)->pdg_id()==22) {
	      photonTrackId = theSimTrksIter->trackId();
	      
	      HepLorentzVector psimvect = theSimTrksIter->momentum();
	      el1_sim_gpx = psimvect.px();
	      el1_sim_gpy = psimvect.py();
	      el1_sim_gpz = psimvect.pz();
	      el1_sim_gpt = TMath::Sqrt(el1_sim_gpx*el1_sim_gpx+el1_sim_gpy*el1_sim_gpy);
	      el1_sim_ge = psimvect.e();
	      
	      //SimVertex indicies start at 0
	      i=0;
	      vector<SimVertex>::const_iterator theSimVertsIter;
	      for(theSimVertsIter = theSimVerts.begin();
		  theSimVertsIter != theSimVerts.end();
		  theSimVertsIter++, i++) {
		
		//check to see if the parent of this vertex was the photon in question
		if(theSimVertsIter->parentIndex() == photonTrackId) {
		  el1_sim_gvx = theSimVertsIter->position().x();
		  el1_sim_gvy = theSimVertsIter->position().y();
		  el1_sim_gvz = theSimVertsIter->position().z();
		  el1_sim_gvr = TMath::Sqrt(el1_sim_gvx*el1_sim_gvx+el1_sim_gvy*el1_sim_gvy);
		  el1_sim_gvphi = atan2(el1_sim_gvy, el1_sim_gvx);
		  //make sure that we don't get Log of 0
		  if(fabs(el1_sim_gvz) > 1E-7 && el1_sim_gvr > 1E-7) {
		    el1_sim_gveta = TMath::Log( TMath::Tan( 0.5*TMath::ATan(el1_sim_gvr/fabs(el1_sim_gvz) ) ) );
		    if(el1_sim_gvz < 0) 
		      el1_sim_gveta = -fabs(el1_sim_gveta);
		  }else el1_sim_gveta = -999.;
		  
		  /*loop over the simtracks (ugh...again) to see if there are
		    any tracks which have this vertex
		  */
		  for(unsigned int k = 0; k < theSimTrks.size(); k++) {
		    if(theSimTrks.at(k).vertIndex() == i) {
		      float charge = theSimTrks.at(k).charge();
		      if(el1_q/charge > 0 ) {
			el1_sim_decay1px = theSimTrks.at(k).momentum().x();
			el1_sim_decay1py = theSimTrks.at(k).momentum().y();
			el1_sim_decay1pz = theSimTrks.at(k).momentum().z();
			el1_sim_decay1pt = TMath::Sqrt(el1_sim_decay1px*el1_sim_decay1px+
						       el1_sim_decay1py*el1_sim_decay1py);
			el1_sim_decay1e = theSimTrks.at(k).momentum().t();
			el1_sim_decay1pid = theSimTrks.at(k).type();
		      } else {
			el1_sim_decay2px = theSimTrks.at(k).momentum().x();
			el1_sim_decay2py = theSimTrks.at(k).momentum().y();
			el1_sim_decay2pz = theSimTrks.at(k).momentum().z();
			el1_sim_decay2pt = TMath::Sqrt(el1_sim_decay2px*el1_sim_decay2px+
						       el1_sim_decay2py*el1_sim_decay2py);
			el1_sim_decay2e = theSimTrks.at(k).momentum().t();
			el1_sim_decay2pid = theSimTrks.at(k).type();
		      }
		      
		    }
		  }//sim track loop 
		}//is the parent of the vertex a photon?
	      }//SimVerteces Iterator
	    }//if(theSimTrksIter->type() == 22 && (*geniter)->pdg_id()==22)  
	  }//if(theSimTrksIter->genpartIndex() == mctk1_index) {
	}//SimTracksIterator
      
      }//if(el1_dr<0.1)
      
      
      //Fill Converted photon information
      InitConvertedPhotonVariables();
      //logfile << "SC pt: " << itsc->energy() << endl;
      FillConvertedPhotonInfo(theConPhotons, *itsc,
			      &(*nearElectron), &(*nearElectron1));
      

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

void Conversion::R9_25_gsf(const Event & event, const reco::PixelMatchGsfElectron* e,
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

//
//---------------------------------------------------------------------------
//

void Conversion::InitSimVariables() {

  sc_sim_gpx = 0.0;
  sc_sim_gpy = 0.0;
  sc_sim_gpz = 0.0;
  sc_sim_gpt = 0.0;
  sc_sim_ge = 0.0; 
  sc_sim_gvx = 0.0; 
  sc_sim_gvy = 0.0;
  sc_sim_gvz = 0.0;
  sc_sim_gvr = 0.0;
  sc_sim_gvphi = 0.0;
  sc_sim_gveta = 0.0;
  sc_sim_decay1px = 0.0;
  sc_sim_decay1py = 0.0;
  sc_sim_decay1pz = 0.0;
  sc_sim_decay1pt = 0.0;  
  sc_sim_decay1e = 0.0;
  sc_sim_decay1pid = 0;
  sc_sim_decay2px = 0.0; 
  sc_sim_decay2py = 0.0;
  sc_sim_decay2pz = 0.0;
  sc_sim_decay2pt = 0.0;
  sc_sim_decay2e = 0.0;
  sc_sim_decay2pid = 0;



  el_sim_gpx = 0.0;
  el_sim_gpy = 0.0;
  el_sim_gpz = 0.0;
  el_sim_gpt = 0.0;
  el_sim_ge = 0.0; 
  el_sim_gvx = 0.0; 
  el_sim_gvy = 0.0;
  el_sim_gvz = 0.0;
  el_sim_gvr = 0.0;
  el_sim_gvphi = 0.0;
  el_sim_gveta = 0.0;
  el_sim_decay1px = 0.0;
  el_sim_decay1py = 0.0;
  el_sim_decay1pz = 0.0;
  el_sim_decay1pt = 0.0;  
  el_sim_decay1e = 0.0;
  el_sim_decay1pid = 0;
  el_sim_decay2px = 0.0; 
  el_sim_decay2py = 0.0;
  el_sim_decay2pz = 0.0;
  el_sim_decay2pt = 0.0;
  el_sim_decay2e = 0.0;
  el_sim_decay2pid = 0;

  el1_sim_gpx = 0.0;
  el1_sim_gpy = 0.0;
  el1_sim_gpz = 0.0;
  el1_sim_gpt = 0.0;
  el1_sim_ge = 0.0; 
  el1_sim_gvx = 0.0; 
  el1_sim_gvy = 0.0;
  el1_sim_gvz = 0.0;
  el1_sim_gvr = 0.0;
  el1_sim_gvphi = 0.0;
  el1_sim_gveta = 0.0;
  el1_sim_decay1px = 0.0;
  el1_sim_decay1py = 0.0;
  el1_sim_decay1pz = 0.0;
  el1_sim_decay1pt = 0.0;  
  el1_sim_decay1e = 0.0;
  el1_sim_decay1pid = 0;
  el1_sim_decay2px = 0.0; 
  el1_sim_decay2py = 0.0;
  el1_sim_decay2pz = 0.0;
  el1_sim_decay2pt = 0.0;
  el1_sim_decay2e = 0.0;
  el1_sim_decay2pid = 0;  
  
  
}

//
//---------------------------------------------------------------------------
//

void Conversion::InitConvertedPhotonVariables() {
  
  
  sc_con_flag = -1;
  sc_con_numCP = 0;
  sc_con_pairinvmass = 0.0;
  sc_con_pairdcottheta = 0.0;
  sc_con_pairp = 0.0; 
  sc_con_pairphi = 0.0;
  sc_con_paireta =0.0;
  sc_con_vx = 0.0;
  sc_con_vy = 0.0;
  sc_con_vz = 0.0;
  sc_con_vphi = 0.0;
  sc_con_veta = 0.0;
  sc_con_vr = 0.0;
  sc_con_decay1px = 0.0;
  sc_con_decay1py = 0.0;
  sc_con_decay1pz = 0.0;
  sc_con_decay1pt = 0.0;
  sc_con_decay1q  = 0;
  sc_con_decay2px = 0.0; 
  sc_con_decay2py = 0.0;
  sc_con_decay2pz = 0.0;
  sc_con_decay2pt = 0.0;
  sc_con_decay2q  = 0;
  sc_con_decay1numhits = 0;
  sc_con_decay2numhits = 0;

  
  
  el_con_flag = -1;
  el_con_numCP = 0;
  el_con_pairinvmass = 0.0;
  el_con_pairdcottheta = 0.0;
  el_con_pairp = 0.0;
  el_con_pairphi = 0.0;
  el_con_paireta = 0.0;
  el_con_elcottheta = 0.0;
  el_con_elcottheta_2 = 0.0;
  el_con_vx = 0.0;
  el_con_vy = 0.0;
  el_con_vz = 0.0;
  el_con_vphi = 0.0;
  el_con_veta = 0.0;
  el_con_vr = 0.0;
  el_con_decay1px = 0.0;
  el_con_decay1py = 0.0;
  el_con_decay1pz = 0.0;
  el_con_decay1pt = 0.0;
  el_con_decay1q  = 0;
  el_con_decay2px = 0.0;
  el_con_decay2py = 0.0;
  el_con_decay2pz = 0.0;
  el_con_decay2pt = 0.0;
  el_con_decay2q  = 0;
  el_con_decay1numhits = 0; 
  el_con_decay1numsharedhits = 0;
  el_con_decay2numhits = 0; 
  el_con_decay2numsharedhits = 0;
  el_con_decay1fracsharedhits = 0.0;
  el_con_decay2fracsharedhits = 0.0;



  el1_con_flag = -1;
  el1_con_numCP = 0;
  el1_con_pairinvmass = 0.0;
  el1_con_pairdcottheta = 0.0;
  el1_con_pairp = 0.0;
  el1_con_pairphi = 0.0;
  el1_con_paireta = 0.0 ;
  el1_con_el1cottheta = 0.0;
  el1_con_el1cottheta_2 = 0.0;
  el1_con_vx = 0.0;
  el1_con_vy = 0.0;
  el1_con_vz = 0.0;
  el1_con_vphi = 0.0;
  el1_con_veta = 0.0;
  el1_con_vr = 0.0;
  el1_con_decay1px = 0.0;
  el1_con_decay1py = 0.0;
  el1_con_decay1pz = 0.0;
  el1_con_decay1pt = 0.0;
  el1_con_decay1q  = 0;
  el1_con_decay2px = 0.0;
  el1_con_decay2py = 0.0;
  el1_con_decay2pz = 0.0;
  el1_con_decay2pt = 0.0;
  el1_con_decay2q  = 0;
  el1_con_decay1numhits = 0; 
  el1_con_decay1numsharedhits = 0;
  el1_con_decay2numhits = 0; 
  el1_con_decay2numsharedhits = 0;
  el1_con_decay1fracsharedhits = 0.0;
  el1_con_decay2fracsharedhits = 0.0;



}


//
//---------------------------------------------------------------------------
//

 
void Conversion::FillConvertedPhotonInfo(vector<ConvertedPhoton>& theConPhotons, 
                                         const SuperCluster& sc,
                                         const PixelMatchGsfElectron* El,
                                         const PixelMatchGsfElectron* El1) {

  //get the GSF tracks from the 2 electron objects
  const Track *el_tr;
  const Track *el1_tr;
  if(el_dr < 0.1) el_tr = (El->gsfTrack()).get();
  if(el1_dr < 0.1) el1_tr = (El1->gsfTrack()).get();
  

  //get ref to SC closest to given SC
  SuperClusterRef conPhotonSC;
  ConvertedPhoton matchedPhoton;
  vector<ConvertedPhoton> temp;
  for(vector<ConvertedPhoton>::const_iterator git = theConPhotons.begin();
      git != theConPhotons.end();
      git++) {
    conPhotonSC = git->superCluster();
    double dphi = conPhotonSC->phi() - sc.phi();
    
    if( fabs(dphi) > TMath::Pi() ) dphi = 2*TMath::Pi() - fabs(dphi);
    double dR = TMath::Sqrt(TMath::Power(conPhotonSC->eta() - sc.eta(), 2) +
                            TMath::Power(dphi, 2)); 
    //logfile << "dR: " << dR << endl;
    if(dR < 0.001) 
      temp.push_back(*git);
  }//convertedPhoton iterator
  

  sc_con_numCP = temp.size();
  if(el_dr < 0.1) 
    el_con_numCP = temp.size();
  if(el1_dr < 0.1) 
    el1_con_numCP = temp.size();
  
  
  
  //if no converted Photons, all flags are -1 by default 
  if(temp.size() == 0 ) {
    sc_con_flag = -1;
    el_con_flag = -1;
    el1_con_flag = -1;
  } else {
    //sc_con_flag is like a boolean - 0 means that CP exists
    // -1 means that it doesn't exist for this SC
    sc_con_flag = 0; 
  }
  
  
  //if there are more than or equal to 1 ConvertedPhoton objects
  if(temp.size() >= 1) {
    
    vector<ConvertedPhoton> sc_passedPhotons;
    vector<ConvertedPhoton> el_passedPhotons;
    vector<ConvertedPhoton> el1_passedPhotons;
    
    //case where there are no tracks in any of the CP objects = flag 0
    int numtracks = 0;
    for(vector<ConvertedPhoton>::iterator tempit = temp.begin();
	tempit != temp.end(); 
	tempit++) {
      numtracks = numtracks+(tempit->tracks()).size();
    }
    
    if(numtracks == 0) {
      if(el_dr < 0.1)
	el_con_flag = 0;
      if(el1_dr < 0.1)
	el1_con_flag = 0;
    } else {
      
      //now the case where there is atleast 1 CP with 1 track, 
      // but none of the tracks share > 50% of hits with el object - flag 1
      // if atleast 1 track shares > 50% of hits with el object - flag 2
      int el_flag2 = -999;
      int el1_flag2 = -999;
     
      for(vector<ConvertedPhoton>::iterator tempit = temp.begin();
	  tempit != temp.end(); 
	  tempit++) {
	
	int ntracks = (tempit->tracks()).size();
	// cannot have 2 tracks for flag 1, 2 - get out of loop
	if(ntracks == 2) {
	  el_flag2 = -999;
	  el1_flag2 = -999;
	  el_passedPhotons.clear();
	  el1_passedPhotons.clear();
	  break;
	}
	if(ntracks != 1) continue; //must have only 1 track
	
	if(el_flag2 == -999)
	  el_flag2 = 0;
	if(el1_flag2 == -999)
	  el1_flag2 = 0;
	
	const Track *tr = ((tempit->tracks()).at(0)).get();
	
	if(el_dr < 0.1) {
	  float el_frac = (sharedHits(*tr, *el_tr)).second;
	  el_passedPhotons.push_back(*tempit);
	  if(el_frac > 0.5) {
	    el_flag2++;
          }	
	}
	
	if(el1_dr < 0.1) {
	  float el1_frac = (sharedHits(*tr, *el1_tr)).second;
	  el1_passedPhotons.push_back(*tempit);
	  if(el1_frac > 0.5) {
	    el1_flag2++;
          }
	}
      } // for(vector<ConvertedPhoton>::iterator tempit = temp.begin();
      
      
      if( el_dr < 0.1 && el_flag2 != -999 ) {
     
	if(el_flag2 == 0) el_con_flag = 1;
	if(el_flag2 > 0) el_con_flag = 2;
      }
      if( el1_dr < 0.1 && el1_flag2 != -999 ) {
	if(el1_flag2 == 0) el1_con_flag = 1;
	if(el1_flag2 > 0) el1_con_flag = 2;
      }
          
	
     
      for(vector<ConvertedPhoton>::iterator tempit = temp.begin();
	  tempit != temp.end(); 
	  tempit++) {
      
	int ntracks = (tempit->tracks()).size();
	if(ntracks !=2) continue; //2 track requirement!!!!
	
	sc_passedPhotons.push_back(*tempit);
	
	if(el_dr < 0.1) {
	  for(int i=0; i < 2; i++) {
	    const Track *tr = ((tempit->tracks()).at(i)).get();
	    float frac = (sharedHits(*tr, *el_tr)).second;
	    
	    if(frac >=0.5) {
	      el_passedPhotons.push_back( *tempit );
	      break;
	    }
	  }
	}//if(el_dr < 0.1) {
	
	if(el1_dr < 0.1) {
	  for(int i=0; i < 2; i++) {
	    const Track *tr = ((tempit->tracks()).at(i)).get();
	    float frac = (sharedHits(*tr, *el1_tr)).second;
	    
	    if(frac >=0.5) {
	      el1_passedPhotons.push_back( *tempit );
	      break;
	    }
	  }
	}//if(el1_dr < 0.1) {
      }// for(vector<ConvertedPhoton>::iterator tempit = temp.begin();

      if( el_dr < 0.1 && el_con_flag < 0) {
	if(el_passedPhotons.size() == 0) el_con_flag = 3;
	if(el_passedPhotons.size() == 1) el_con_flag = 4;
	if(el_passedPhotons.size() > 1 ) el_con_flag = 5;
      }
      if(el1_dr < 0.1 && el1_con_flag < 0) {
	if(el1_passedPhotons.size() == 0) el1_con_flag = 3;
	if(el1_passedPhotons.size() == 1) el1_con_flag = 4;
	if(el1_passedPhotons.size() > 1 ) el1_con_flag = 5;
      }
    }

    
    /* if flag is 3, go back and loop through CP vector
       and put into the passedPhotons vector all CP objects
       with ntracks == 2
    */
    if(el_con_flag == 3) {
      for(vector<ConvertedPhoton>::iterator tempit = temp.begin();
	  tempit != temp.end(); 
	  tempit++) {
	
	if(tempit->nTracks() != 2) continue;
	el_passedPhotons.push_back(*tempit);
      }
    }
    if(el1_con_flag == 3) {
      for(vector<ConvertedPhoton>::iterator tempit = temp.begin();
	  tempit != temp.end(); 
	  tempit++) {
	
	if(tempit->nTracks() != 2) continue;
	el1_passedPhotons.push_back(*tempit);
      }
    }


    /*now go back and loop over the vector of ConPhoton objects 
      which have 2 tracks, and have atleast 1 track which 
      shares  > 50% of hits with the reco el track
    */
    
    
    if(sc_passedPhotons.size() >= 1) {
      sortConvPhotonVector(sc_passedPhotons);
      ConvertedPhoton con = sc_passedPhotons.at(0);
      sc_con_pairinvmass = con.pairInvariantMass();
      sc_con_pairdcottheta = con.pairCotThetaSeparation();
      sc_con_pairp = con.pairMomentum().mag();
      sc_con_pairphi = con.pairMomentumPhi();
      sc_con_paireta = con.pairMomentumEta();

      sc_con_vx = con.convVertexPosition().x();
      sc_con_vy = con.convVertexPosition().y();
      sc_con_vz = con.convVertexPosition().z();
      
      sc_con_vr = TMath::Sqrt(sc_con_vx*sc_con_vx+
			      sc_con_vy*sc_con_vy);
			      
      
      sc_con_vphi = atan2(sc_con_vy, sc_con_vx);
      if(fabs(sc_con_vz) > 1E-7 && sc_con_vr > 1E-7) {
	sc_con_veta = TMath::Log( TMath::Tan( 0.5*TMath::ATan(sc_con_vr/fabs(sc_con_vz) ) ) );
	if(sc_con_vz < 0) sc_con_veta = -fabs(sc_con_veta);
      } else sc_con_veta = 500;
      
      
      const Track *tr1 = ((con.tracks()).at(0)).get();
      const Track *tr2 = ((con.tracks()).at(1)).get();

      sc_con_decay1px = tr1->momentum().x();
      sc_con_decay1py = tr1->momentum().y();
      sc_con_decay1pz = tr1->momentum().z();
      sc_con_decay1pt = TMath::Sqrt(sc_con_decay1px*sc_con_decay1px+
				    sc_con_decay1py*sc_con_decay1py);
      sc_con_decay1q = tr1->charge();
      sc_con_decay1numhits = tr1->numberOfValidHits();
	    

      sc_con_decay2px = tr2->momentum().x();
      sc_con_decay2py = tr2->momentum().y();
      sc_con_decay2pz = tr2->momentum().z();
      sc_con_decay2pt = TMath::Sqrt(sc_con_decay2px*sc_con_decay2px+
				    sc_con_decay2py*sc_con_decay2py);
      sc_con_decay2q = tr2->charge();
      sc_con_decay2numhits = tr2->numberOfValidHits();
    }//if(sc_passedPhotons.size() == 1) {

    
    //deal with cases where there are more than 1 convertedPhoton objects
    if(el_passedPhotons.size() >= 1) {
      
      
      if(el_con_flag > 2 ) {
	sortConvPhotonVector(el_passedPhotons);
	ConvertedPhoton con = el_passedPhotons.at(0);

	el_con_pairdcottheta = con.pairCotThetaSeparation();
	el_con_pairinvmass = con.pairInvariantMass();
	el_con_pairp = con.pairMomentum().mag();
	el_con_pairphi = con.pairMomentumPhi();
	el_con_paireta = con.pairMomentumEta();
	
	el_con_vx = con.convVertexPosition().x();
	el_con_vy = con.convVertexPosition().y();
	el_con_vz = con.convVertexPosition().z();
	
	el_con_vr = TMath::Sqrt(el_con_vx*el_con_vx+
				el_con_vy*el_con_vy);
	
	el_con_vphi = atan2(el_con_vy, el_con_vx);
	if(fabs(el_con_vz) > 1E-7 && el_con_vr > 1E-7) {
	  el_con_veta = TMath::Log( TMath::Tan( 0.5*TMath::ATan(el_con_vr/fabs(el_con_vz) ) ) );
	  if(el_con_vz < 0) el_con_veta = -fabs(el_con_veta);
	} else el_con_veta = 500;
	
	
	const Track *tr1;
	const Track *tr2;
	
	
	if( (sharedHits( *(((con.tracks()).at(0)).get()), *el_tr)).first
	    >= (sharedHits( *(((con.tracks()).at(1)).get()), *el_tr)).first) {
	  tr1 =  ((con.tracks()).at(0)).get();
	  tr2 =  ((con.tracks()).at(1)).get();
	} else {
	  tr1 =  ((con.tracks()).at(1)).get();
	  tr2 =  ((con.tracks()).at(0)).get();
	}

	/*fill the difference in cot(theta) between the reco
	  electron and the leading (i.e. more hits) track
	*/
	float eltheta = el_tr->innerMomentum().Theta();
	float contrtheta = tr1->innerMomentum().Theta();
	
	
	if(el_passedPhotons.size() > 1) {
	  el_con_elcottheta = fabs(1/TMath::Tan(eltheta) - 1/TMath::Tan(contrtheta));
	  double secondbesttheta = (FindSecondBestTrack(el_passedPhotons, *el_tr)).innerMomentum().Theta();
	  double trtheta = el_tr->innerMomentum().Theta();
	  el_con_elcottheta_2 = fabs(1/TMath::Tan(secondbesttheta) - 1/TMath::Tan(trtheta));	 
	}
	
	el_con_decay1px = tr1->momentum().x();
	el_con_decay1py = tr1->momentum().y();
	el_con_decay1pz = tr1->momentum().z();
	el_con_decay1pt = TMath::Sqrt(el_con_decay1px*el_con_decay1px+
				      el_con_decay1py*el_con_decay1py);
	el_con_decay1q = tr1->charge();
	el_con_decay1numhits = tr1->numberOfValidHits();
	el_con_decay1numsharedhits = (sharedHits(*tr1, *el_tr)).first;
	el_con_decay1fracsharedhits = (sharedHits(*tr1, *el_tr)).second;
	
	
	el_con_decay2px = tr2->momentum().x();
	el_con_decay2py = tr2->momentum().y();
	el_con_decay2pz = tr2->momentum().z();
	el_con_decay2pt = TMath::Sqrt(el_con_decay2px*el_con_decay2px+
				      el_con_decay2py*el_con_decay2py);
	el_con_decay2q = tr2->charge();
	el_con_decay2numhits = tr2->numberOfValidHits();
	el_con_decay2numsharedhits = (sharedHits(*tr2, *el_tr)).first;
	el_con_decay2fracsharedhits = (sharedHits(*tr2, *el_tr)).second;
	
	
      } else {
	sortConvPhotonVectorByHits(el_passedPhotons, *el_tr);
	ConvertedPhoton con = el_passedPhotons.at(0);
	
	const Track *tr1 = ((con.tracks()).at(0)).get();

	/*fill the difference in cot(theta) between the reco
	  electron and the track
	*/
	float eltheta = el_tr->innerMomentum().Theta();
	float contrtheta = tr1->innerMomentum().Theta();
	el_con_elcottheta = fabs(1/TMath::Tan(eltheta) - 1/TMath::Tan(contrtheta));
	
	if(el_passedPhotons.size() > 1) {
	  double secondbesttheta = (FindSecondBestTrack(el_passedPhotons, *el_tr)).innerMomentum().Theta();
	  double trtheta = el_tr->innerMomentum().Theta();
	  el_con_elcottheta_2 = fabs(1/TMath::Tan(secondbesttheta) - 1/TMath::Tan(trtheta));	 
	}

	el_con_decay1px = tr1->momentum().x();
	el_con_decay1py = tr1->momentum().y();
	el_con_decay1pz = tr1->momentum().z();
	el_con_decay1pt = TMath::Sqrt(el_con_decay1px*el_con_decay1px+
				      el_con_decay1py*el_con_decay1py);
	el_con_decay1q = tr1->charge();
	el_con_decay1numhits = tr1->numberOfValidHits();
	el_con_decay1numsharedhits = (sharedHits(*tr1, *el_tr)).first;
	el_con_decay1fracsharedhits = (sharedHits(*tr1, *el_tr)).second;
	
      }//else
    } //if(el_passedPhotons.size() >= 1) {
      
    
    if(el1_passedPhotons.size() >= 1) {
      
      if(el1_con_flag > 2 ) {
	sortConvPhotonVector(el1_passedPhotons);
	ConvertedPhoton con = el1_passedPhotons.at(0);

	el1_con_pairdcottheta = con.pairCotThetaSeparation();
	el1_con_pairinvmass = con.pairInvariantMass();
	el1_con_pairp = con.pairMomentum().mag();
	el1_con_pairphi = con.pairMomentumPhi();
	el1_con_paireta = con.pairMomentumEta();
	
	el1_con_vx = con.convVertexPosition().x();
	el1_con_vy = con.convVertexPosition().y();
	el1_con_vz = con.convVertexPosition().z();
	
	el1_con_vr = TMath::Sqrt(el1_con_vx*el1_con_vx+
				el1_con_vy*el1_con_vy);
	
	el1_con_vphi = atan2(el1_con_vy, el1_con_vx);
	if(fabs(el1_con_vz) > 1E-7 && el1_con_vr > 1E-7) {
	  el1_con_veta = TMath::Log( TMath::Tan( 0.5*TMath::ATan(el1_con_vr/fabs(el1_con_vz) ) ) );
	  if(el1_con_vz < 0) el1_con_veta = -fabs(el1_con_veta);
	} else el1_con_veta = 500;
	
	
	const Track *tr1;
	const Track *tr2;
	
	
	if( (sharedHits( *(((con.tracks()).at(0)).get()), *el1_tr)).first
	    >= (sharedHits( *(((con.tracks()).at(1)).get()), *el1_tr)).first) {
	  tr1 =  ((con.tracks()).at(0)).get();
	  tr2 =  ((con.tracks()).at(1)).get();
	} else {
	  tr1 =  ((con.tracks()).at(1)).get();
	  tr2 =  ((con.tracks()).at(0)).get();
	}

	/*fill the difference in cot(theta) between the reco
	  electron and the leading (i.e. more hits) track
	*/
	float el1theta = el1_tr->innerMomentum().Theta();
	float contrtheta = tr1->innerMomentum().Theta();
	
	
	if(el1_passedPhotons.size() > 1) {
	  el1_con_el1cottheta = fabs(1/TMath::Tan(el1theta) - 1/TMath::Tan(contrtheta));
	  double secondbesttheta = (FindSecondBestTrack(el1_passedPhotons, *el1_tr)).innerMomentum().Theta();
	  double trtheta = el1_tr->innerMomentum().Theta();
	  el1_con_el1cottheta_2 = fabs(1/TMath::Tan(secondbesttheta) - 1/TMath::Tan(trtheta));	 
	}
	
	el1_con_decay1px = tr1->momentum().x();
	el1_con_decay1py = tr1->momentum().y();
	el1_con_decay1pz = tr1->momentum().z();
	el1_con_decay1pt = TMath::Sqrt(el1_con_decay1px*el1_con_decay1px+
				      el1_con_decay1py*el1_con_decay1py);
	el1_con_decay1q = tr1->charge();
	el1_con_decay1numhits = tr1->numberOfValidHits();
	el1_con_decay1numsharedhits = (sharedHits(*tr1, *el1_tr)).first;
	el1_con_decay1fracsharedhits = (sharedHits(*tr1, *el1_tr)).second;
	
	
	el1_con_decay2px = tr2->momentum().x();
	el1_con_decay2py = tr2->momentum().y();
	el1_con_decay2pz = tr2->momentum().z();
	el1_con_decay2pt = TMath::Sqrt(el1_con_decay2px*el1_con_decay2px+
				      el1_con_decay2py*el1_con_decay2py);
	el1_con_decay2q = tr2->charge();
	el1_con_decay2numhits = tr2->numberOfValidHits();
	el1_con_decay2numsharedhits = (sharedHits(*tr2, *el1_tr)).first;
	el1_con_decay2fracsharedhits = (sharedHits(*tr2, *el1_tr)).second;
	
	
      } else {
	sortConvPhotonVectorByHits(el1_passedPhotons, *el1_tr);
	ConvertedPhoton con = el1_passedPhotons.at(0);
	
	const Track *tr1 = ((con.tracks()).at(0)).get();

	/*fill the difference in cot(theta) between the reco
	  electron and the track
	*/
	float el1theta = el1_tr->innerMomentum().Theta();
	float contrtheta = tr1->innerMomentum().Theta();
	el1_con_el1cottheta = fabs(1/TMath::Tan(el1theta) - 1/TMath::Tan(contrtheta));
	
	if(el1_passedPhotons.size() > 1) {
	  double secondbesttheta = (FindSecondBestTrack(el1_passedPhotons, *el1_tr)).innerMomentum().Theta();
	  double trtheta = el1_tr->innerMomentum().Theta();
	  el1_con_el1cottheta_2 = fabs(1/TMath::Tan(secondbesttheta) - 1/TMath::Tan(trtheta));	 
	}

	el1_con_decay1px = tr1->momentum().x();
	el1_con_decay1py = tr1->momentum().y();
	el1_con_decay1pz = tr1->momentum().z();
	el1_con_decay1pt = TMath::Sqrt(el1_con_decay1px*el1_con_decay1px+
				      el1_con_decay1py*el1_con_decay1py);
	el1_con_decay1q = tr1->charge();
	el1_con_decay1numhits = tr1->numberOfValidHits();
	el1_con_decay1numsharedhits = (sharedHits(*tr1, *el1_tr)).first;
	el1_con_decay1fracsharedhits = (sharedHits(*tr1, *el1_tr)).second;
	
      }//else
    } //if(el1_passedPhotons.size() >= 1) {
    
    
    
  }//if more than 1 converted Photon object...if(temp.size() >= 1) {
  //logfile << "Finished ConvertedPhotonInfo" << endl;
}  
//
//------------------------------------------------------------------------------------
//


pair<unsigned int,float> Conversion::sharedHits(const Track& trackA,
                                                const Track& trackB) {
  
    
  unsigned int shared = 0;
  for(trackingRecHit_iterator tkHitA = trackA.recHitsBegin(); tkHitA !=trackA.recHitsEnd(); ++tkHitA){
    for(trackingRecHit_iterator tkHitB = trackB.recHitsBegin();
        tkHitB !=trackB.recHitsEnd(); ++tkHitB){
      if( (**tkHitA).isValid() && (**tkHitB).isValid() &&(**tkHitA).sharesInput( &(**tkHitB),TrackingRecHit::all)) {
        shared++;
        break;
      }
    }
  }
  
  float fraction = (float) shared/min(trackA.found(),trackB.found());
  return make_pair(shared,fraction);
  
}


//
//------------------------------------------------------------------------------------
//

void Conversion::sortConvPhotonVector(vector<ConvertedPhoton>& conPhotonvect) {
  
  
  /* loop through the vector using iterators. when you find the smallest delta cot theta,
     store in a different vector the corresponding ConvertedPhoton object, and then 
     erase the vector's element in the mother vector
  */

  vector<ConvertedPhoton> orderedvect;
  
  while(conPhotonvect.size() != 0) {
    
    vector<ConvertedPhoton>::iterator minCotThetaIt = conPhotonvect.begin();
    float minCotTheta = minCotThetaIt->pairCotThetaSeparation();
    for(vector<ConvertedPhoton>::iterator it = conPhotonvect.begin();
	it != conPhotonvect.end();
	it++) {
      
      if( fabs(it->pairCotThetaSeparation()) < fabs(minCotTheta) ) {
	minCotTheta = it->pairCotThetaSeparation();
	minCotThetaIt = it;
      } 
    }
    
    orderedvect.push_back(*minCotThetaIt);
    conPhotonvect.erase(minCotThetaIt);
  }
  
  conPhotonvect.clear();
  conPhotonvect = orderedvect;
  
}

//
//------------------------------------------------------------------------------------
//

void Conversion::sortConvPhotonVectorByHits(vector<ConvertedPhoton>& conPhotonvect, const reco::Track& trk) {
  
  
  /* loop through the vector using iterators. when you find the smallest delta cot theta,
     store in a different vector the corresponding ConvertedPhoton object, and then 
     erase the vector's element in the mother vector
  */

  vector<ConvertedPhoton> orderedvect;
  
  while(conPhotonvect.size() != 0) {
    
    vector<ConvertedPhoton>::iterator maxfracSharedHitsIt = conPhotonvect.begin();
    float maxfracSharedHits = sharedHits( *(((maxfracSharedHitsIt->tracks()).at(0)).get()), trk).second;
    for(vector<ConvertedPhoton>::iterator it = conPhotonvect.begin();
	it != conPhotonvect.end();
	it++) {
      
      if( sharedHits(*((it->tracks().at(0)).get()), trk).second >  maxfracSharedHits) {
	maxfracSharedHits =  sharedHits(*((it->tracks().at(0)).get()), trk).second;
	maxfracSharedHitsIt = it;
      } 
    }
    
    orderedvect.push_back(*maxfracSharedHitsIt);
    conPhotonvect.erase(maxfracSharedHitsIt);
  }
  
  conPhotonvect.clear();
  conPhotonvect = orderedvect;
  
    
}

//
//------------------------------------------------------------------------------------
//

    
    
Track Conversion::FindSecondBestTrack(const vector<ConvertedPhoton>& conPhotonvect, const Track& tr) {
  
  /* put every track into a map (except for the best conphoton object which is the best one)
     where the key is the delta cot theta between the track and 
     the electron track and the value is the track. This will get rid of duplicates, 
     and sort in order of delta cot theta. Then return the track in the second element
     in the map
  */
   
  
  map<int, Track> mp;
  
  for(vector<ConvertedPhoton>::const_iterator it = conPhotonvect.begin();
      it != conPhotonvect.end();
      it++) {
     
    if(it==conPhotonvect.begin()) continue;

    float eltheta = tr.innerMomentum().Theta();
    vector<TrackRef> conPhotontrks = it->tracks();
    for(vector<TrackRef>::iterator trkit = conPhotontrks.begin();
	trkit != conPhotontrks.end();
	trkit++) {
      float trtheta = (*trkit)->innerMomentum().theta();
      float delcottheta = fabs(1/TMath::Tan(eltheta) - 1/TMath::Tan(trtheta));
      mp[static_cast<int>(100000*delcottheta)] = *((*trkit).get());
    }
  }

  map<int, Track>::iterator mpiter = mp.begin();
   
  return mpiter->second;
}
   
   
    
    

  

	  
