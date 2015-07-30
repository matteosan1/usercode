#ifndef NewElectrons_H
#define NewElectrons_H

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectronFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>
#include <Math/Point3D.h>
#include <vector>
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

class TFile;
class TTree;

class NewElectrons: public edm::EDAnalyzer {
public:

  NewElectrons(const edm::ParameterSet& pset);

  virtual ~NewElectrons();

  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  void beginJob(const edm::EventSetup& eventSetup);
  void endJob();
  bool inCrack(float eta);
  int mother(HepMC::GenParticle *p); 
  void nHits(const reco::GsfTrackRef, int&, int&);
  double trackIsolation(const math::XYZVector, const math::XYZPoint,
                        const reco::TrackCollection*);
  void subDetector(TrackingRecHitRef hit, int& subdet, int& layer);  

 protected:
  
 private:
  std::string baselineEleCollName;
  std::string customEleCollName;
  std::string fileName;
  
  TFile *file;
  TTree *tree;
  int run, id;
  float mc_pt, mc_eta, mc_phi, mc_e;
  float sc_e, sc_eta, sc_phi, sc_dr, sc_et;
  float sc_rawe;
  float tk_pt, tk_eta, tk_phi, tk_dr;
  float el_pt, el_eta, el_phi, el_dr, el_e;
  float el_eopin, el_eopout, el_hoe, el_detain, el_dphiin, el_pin, el1_pin;
  float el_fbrem, el_eseed, el_e3x3, el_detaout, el_dphiout;
  float el_e5x5, el_spp, el_see, el_pout, el_z0;	
  float el1_pt, el1_eta, el1_phi, el1_dr, el1_e, el1_z0;
  float el1_eopin, el1_eopout, el1_hoe, el1_detain, el1_dphiin;
  float el1_fbrem, el1_eseed, el1_e3x3, el1_detaout, el1_dphiout;
  float el1_e5x5, el1_spp, el1_see, el1_pout;	
  float el1_tkiso, el_tkiso;
  int mc_id, sc_type, mc_mother, mc_crack;
  int tk_layer, tk_subdet, tk_nhit;
  int el_class, el1_class, el1_npxhits, el1_nsihits, el_npxhits, el_nsihits;
  float el_rinnerhit,  el1_rinnerhit;
  int el_detinnerhit, el1_detinnerhit;
  
  HepMC::GenEvent* myGenEvent;

};
#endif

