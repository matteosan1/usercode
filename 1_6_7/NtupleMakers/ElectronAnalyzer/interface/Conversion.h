#ifndef Conversion_H
#define Conversion_H

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"

#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>
#include <Math/Point3D.h>
#include <vector>
#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "CLHEP/HepMC/GenParticle.h"

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

 protected:
  
 private:
  std::string baselineEleCollName;
  std::string customEleCollName;
  std::string fileName;
  
  TFile *file;
  TTree *tree;
  int run, id;
  float mcsc_pt, mcsc_eta, mcsc_phi, mcsc_e, mcsc_dr;
  int mcsc_id, mcsc_mother, mcsc_crack;
  float mctk_pt, mctk_eta, mctk_phi, mctk_e, mctk_dr;
  int mctk_id, mctk_mother, mctk_crack;
  float mctk1_pt, mctk1_eta, mctk1_phi, mctk1_e, mctk1_dr;
  int mctk1_id, mctk1_mother, mctk1_crack;
  float sc_e, sc_eta, sc_phi, sc_dr, sc_et;
  float sc_rawe;
  float sc1_e, sc1_eta, sc1_phi, sc1_dr, sc1_et;
  float sc1_rawe;
  float el_pt, el_eta, el_phi, el_dr, el_e;
  float el_eopin, el_eopout, el_hoe, el_detain, el_dphiin;
  float el_fbrem, el_eseed, el_e3x3, el_detaout, el_dphiout;
  float el_e5x5, el_spp, el_see, el_pout, el_z0;	
  float el1_pt, el1_eta, el1_phi, el1_dr, el1_e, el1_z0;
  float el1_eopin, el1_eopout, el1_hoe, el1_detain, el1_dphiin;
  float el1_fbrem, el1_eseed, el1_e3x3, el1_detaout, el1_dphiout;
  float el1_e5x5, el1_spp, el1_see, el1_pout;	
  float el1_tkiso, el_tkiso;
  float el1_xhit, el1_yhit, el1_zhit, el_xhit, el_yhit, el_zhit;

  int el_class, el1_class, el1_npxhits, el1_nsihits, el_npxhits, el_nsihits;
  int el_detinnerhit, el1_detinnerhit, el_layerinnerhit, el1_layerinnerhit, sc_type;
  float el_rinnerhit, el1_rinnerhit, el1_tkpt, el1_tketa, el1_tkphi, el_tkpt, el_tketa, el_tkphi;

  HepMC::GenEvent* myGenEvent;
};
#endif

