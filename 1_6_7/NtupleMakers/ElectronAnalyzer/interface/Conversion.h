#ifndef Conversion_H
#define Conversion_H

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/EgammaCandidates/interface/ConvertedPhoton.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>
#include <Math/Point3D.h>
#include <vector>
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "CLHEP/HepMC/GenParticle.h"

#include <fstream>
#include <iostream>


namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

class TFile;
class TTree;

class Conversion: public edm::EDAnalyzer {
public:

  Conversion(const edm::ParameterSet& pset);

  virtual ~Conversion();

  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  void beginJob(const edm::EventSetup& eventSetup);
  void endJob();
  bool inCrack(float eta);
  int mother(HepMC::GenParticle *p); 
  void R9_25_gsf(const edm::Event & event, const reco::PixelMatchGsfElectron*,
                 float&, float&, float&, float&, float&); 
  double trackIsolation(const math::XYZVector, const math::XYZPoint,
                        const reco::TrackCollection*);
  void subDetector(TrackingRecHitRef hit, int& subdet, int& layer);
  void nHits(const reco::GsfTrackRef t, int& nPixelHits, int& nSiTkHits);
  void InitSimVariables();
  void InitConvertedPhotonVariables();
  void FillConvertedPhotonInfo(std::vector<reco::ConvertedPhoton>&, 
                               const reco::SuperCluster&,
                               const reco::PixelMatchGsfElectron*,
                               const reco::PixelMatchGsfElectron*);
  std::pair<unsigned int,float> sharedHits(const reco::Track& trackA,
					   const reco::Track& trackB);
  //method to sort the ConvertedPhoton objects 
  void sortConvPhotonVector(std::vector<reco::ConvertedPhoton>&, const reco::Track&);
  reco::Track FindSecondBestTrack(const std::vector<reco::ConvertedPhoton>&, const reco::Track&);
  
 protected:
  
 private:

  ofstream logfile;

  
  std::string baselineEleCollName;
  std::string customEleCollName;
  std::string fileName;
  
  TFile *file;
  TTree *tree;
  
  
  int run, id;
  
  //info about the Gen particle closest matched in dR to Supercluster 
  float mcsc_pt, mcsc_eta, mcsc_phi, mcsc_e, mcsc_dr;
  int mcsc_id, mcsc_mother, mcsc_crack;

  //info about the Gen particle closest matched in dR 
  //to the track belonging to the standard electron
  float mctk_pt, mctk_eta, mctk_phi, mctk_e, mctk_dr;
  int mctk_id, mctk_mother, mctk_crack;

  //info about the Gen particle closest matched in dR 
  //to the track belonging to the UCSD electron
  float mctk1_pt, mctk1_eta, mctk1_phi, mctk1_e, mctk1_dr;
  int mctk1_id, mctk1_mother, mctk1_crack;
  
  //standard electron's SC info
  float sc_e, sc_eta, sc_phi, sc_dr, sc_et;
  float sc_rawe;
  int sc_type;
  
  //UCSD electron's SC info
  float sc1_e, sc1_eta, sc1_phi, sc1_dr, sc1_et;
  float sc1_rawe;
  
  
  //standard electron information
  float el_pt, el_eta, el_phi, el_dr, el_e, el_q;
  float el_eopin, el_eopout, el_hoe, el_detain, el_dphiin;
  float el_fbrem, el_eseed, el_e3x3, el_detaout, el_dphiout;
  //sigma phi-phi, sigma eta-eta - width of the shower in phi and eta plane
  float el_spp, el_see;
  float el_e5x5, el_pout, el_z0;	
  int el_class, el_npxhits, el_nsihits;
  //innermost hit information
  int el_detinnerhit, el_layerinnerhit;
  float el_rinnerhit, el_tkpt, el_tketa, el_tkphi;

  
  //UCSD electron info
  float el1_pt, el1_eta, el1_phi, el1_dr, el1_e, el1_z0, el1_q;
  float el1_eopin, el1_eopout, el1_hoe, el1_detain, el1_dphiin;
  float el1_fbrem, el1_eseed, el1_e3x3, el1_detaout, el1_dphiout;
  //sigma phi-phi, sigma eta-eta - width of the shower in phi and eta plane
  float el1_spp, el1_see;
  float el1_e5x5, el1_pout;	
  float el1_tkiso, el_tkiso;
  float el1_xhit, el1_yhit, el1_zhit, el_xhit, el_yhit, el_zhit;
  int el1_class, el1_npxhits, el1_nsihits;
  //innermost hit information
  int el1_detinnerhit, el1_layerinnerhit;
  float el1_rinnerhit, el1_tkpt, el1_tketa, el1_tkphi;
  float el_d0, el_d0err, el1_d0, el1_d0err;
  
  
  //Sim information
  //px,py,pz,e of the sim photon closest to the ctf track for the standard electron
  float el_sim_gpx, el_sim_gpy, el_sim_gpz, el_sim_gpt, el_sim_ge; 
  //position of the vertex where the conversion happened
  float el_sim_gvx, el_sim_gvy, el_sim_gvz, el_sim_gveta, el_sim_gvphi, el_sim_gvr;
  //px,py,pz, e, pid of the decay products
  float el_sim_decay1px, el_sim_decay1py, el_sim_decay1pz, el_sim_decay1pt, el_sim_decay1e;
  int el_sim_decay1pid;
  float el_sim_decay2px, el_sim_decay2py, el_sim_decay2pz, el_sim_decay2pt, el_sim_decay2e;
  int el_sim_decay2pid;


  //for the UCSD electron
  float el1_sim_gpx, el1_sim_gpy, el1_sim_gpz, el1_sim_gpt, el1_sim_ge; 
  //position of the vertex where the conversion happened - take out tk
  float el1_sim_gvx, el1_sim_gvy, el1_sim_gvz, el1_sim_gveta, el1_sim_gvphi, el1_sim_gvr;
  //px,py,pz, e, pid of the decay products
  float el1_sim_decay1px, el1_sim_decay1py, el1_sim_decay1pz, el1_sim_decay1pt, el1_sim_decay1e;
  int el1_sim_decay1pid;
  float el1_sim_decay2px, el1_sim_decay2py, el1_sim_decay2pz, el1_sim_decay2pt, el1_sim_decay2e;
  int el1_sim_decay2pid;
  
  
  //similarily for the supercluster: 
  float sc_sim_gpx, sc_sim_gpy, sc_sim_gpz, sc_sim_gpt, sc_sim_ge; 
  float sc_sim_gvx, sc_sim_gvy, sc_sim_gvz, sc_sim_gveta, sc_sim_gvphi, sc_sim_gvr;
  float sc_sim_decay1px, sc_sim_decay1py, sc_sim_decay1pz, sc_sim_decay1pt, sc_sim_decay1e;
  int sc_sim_decay1pid;
  float sc_sim_decay2px, sc_sim_decay2py, sc_sim_decay2pz, sc_sim_decay2pt, sc_sim_decay2e; 
  int sc_sim_decay2pid;
  //_______________________________________________________________________//


  //ConvertedPhoton (convg) information
  /* what the flag means:
      flag = -1 no ConvertedPhoton (CP) object found for the SuperCluster 
           (SC) in the electron object (see comment at the end) 
  flag = 0  CP found, but no tracks in it 
  flag = 1  CP found, only one track in it, track shares < 50% of hits 
            with track in electron object 
  flag = 2  CP found, only one track in it, track shares > 50% of hits 
            with track in electron object 
  flag = 3  CP found, possibly more than one, with 2 tracks in it. 
            None of the CPs have track that shares > 50% of hits with 
            track in electron object 
  flag = 4  CP found, possibly more than one, with 2 tracks in it. 
            One and only one of the CPs have track that shares > 50% of 
            hits with track in electron object 
  flag = 5  CP found, possibly more than one, with 2 tracks in it. 
            More than one CP has track that shares > 50% of hits with 
            track in electron object 
  */
  
  
  //ConvertedPhoton information - UCSD electron
  int el_con_flag;
  //inv mass, delta cot(theta) between the 2 tracks, momentum of 2 tracks, phi , eta of 2 tracks
  float el_con_pairinvmass, el_con_pairdcottheta, el_con_pairp, el_con_pairphi, el_con_paireta;
  /*delta cot(theta) between leading track(defined as being the one whose delta cot theta between 
    it's 2 tracks is the minumum) and electron object. 
  */
  float el_con_elcottheta; 
  /* delta cot(theta) between second best track (2nd best track is defined by the track with the
     smallest delta cot(theta) between itself and the electron object, among all tracks in ConPhoton
     objects EXCEPT the leading one,
  */
  float el_con_elcottheta_2;
  //conversion vertex position
  float el_con_vx, el_con_vy, el_con_vz, el_con_vphi, el_con_veta, el_con_vr;
  //track info
  float el_con_decay1px, el_con_decay1py, el_con_decay1pz, el_con_decay1pt;
  int el_con_decay1q;
  float el_con_decay2px, el_con_decay2py, el_con_decay2pz, el_con_decay2pt;
  int el_con_decay2q;
  //num of hits in the track, number of hits shared with el object and fraction of shared hits
  int el_con_decay1numhits, el_con_decay1numsharedhits;
  int el_con_decay2numhits, el_con_decay2numsharedhits;
  float el_con_decay1fracsharedhits;
  float el_con_decay2fracsharedhits;


  //ConvertedPhoton information - UCSD electron
  int el1_con_flag;
  //cot(theta) between the 2 tracks, momentum of 2 tracks, phi , eta of pair momentum
  float el1_con_pairinvmass, el1_con_pairdcottheta, el1_con_pairp, el1_con_pairphi, el1_con_paireta;
  //delta cot(theta) between leading track in con photon object and the reco electron
  float el1_con_el1cottheta, el1_con_el1cottheta_2;
  //conversion vertex position
  float el1_con_vx, el1_con_vy, el1_con_vz, el1_con_vphi, el1_con_veta, el1_con_vr;
  //track info
  float el1_con_decay1px, el1_con_decay1py, el1_con_decay1pz, el1_con_decay1pt;
  int el1_con_decay1q;
  float el1_con_decay2px, el1_con_decay2py, el1_con_decay2pz, el1_con_decay2pt;
  int el1_con_decay2q;
  //num of hits in the track, number of hits shared with el1 object and fraction of shared hits
  int el1_con_decay1numhits, el1_con_decay1numsharedhits;
  int   el1_con_decay2numhits, el1_con_decay2numsharedhits; 
  float el1_con_decay1fracsharedhits;
  float el1_con_decay2fracsharedhits;
  



  //ConvertedPhoton information - SC
  int sc_con_flag; //flag can only be -1,0,1,2,3 or 4
  //cot(theta) between the 2 tracks, momentum of 2 tracks, phi , eta of 2 tracks
  float sc_con_pairinvmass, sc_con_pairdcottheta, sc_con_pairp, sc_con_pairphi, sc_con_paireta;
  //conversion vertex position
  float sc_con_vx, sc_con_vy, sc_con_vz, sc_con_vphi, sc_con_veta, sc_con_vr;
  //track info
  float sc_con_decay1px, sc_con_decay1py, sc_con_decay1pz, sc_con_decay1pt;
  int sc_con_decay1q;
  float sc_con_decay2px, sc_con_decay2py, sc_con_decay2pz, sc_con_decay2pt;
  int sc_con_decay2q;
  //num of hits in the track, number of hits shared with sc object and fraction of shared hits
  int sc_con_decay1numhits, sc_con_decay2numhits;

  
  HepMC::GenEvent* myGenEvent;
};
#endif

