#ifndef hzz_H
#define hzz_H

#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"

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

class hzz: public edm::EDAnalyzer {
public:

  hzz(const edm::ParameterSet& pset);

  virtual ~hzz();

  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  void beginJob(const edm::EventSetup& eventSetup);
  void endJob();
  bool inCrack(float eta);
  int mother(HepMC::GenParticle *p); 
  void R9_25_gsf(const edm::Event & event, const reco::PixelMatchGsfElectron*,
                 float&, float&, float&, float&, float&); 
  //void R9_25_ctf(const edm::Event & event, const reco::GlobalCtfElectron*,
  //                 float&, float&, float&, float&, float&);

protected:

private:
  std::string baselineEleCollName;
  std::string customEleCollName;
  std::string fileName;

  TFile *file;
  TTree *tree;

  int   nz; 
  float mzh, mzl;
  float mz1, mz2, mh;

  int   rnz, qnz;
  float rmh, rmz1, rmz2;
  float qmh, qmz1, qmz2;


  int nMC, run, id;
  float mc_pt[20], mc_eta[20], mc_phi[20], mc_e[20];
  float sc_e[20], sc_eta[20], sc_phi[20], sc_dr[20], sc_et[20];
  float sc_rawe[20];
  float tk_pt[20], tk_eta[20], tk_phi[20], tk_dr[20];
  float el_pt[20], el_eta[20], el_phi[20], el_dr[20], el_e[20];
  float el_eopin[20], el_eopout[20], el_hoe[20], el_detain[20], el_dphiin[20];
  float el_fbrem[20], el_eseed[20], el_e3x3[20], el_detaout[20], el_dphiout[20];
  float el_e5x5[20], el_spp[20], el_see[20], el_pout[20];	
  float el1_pt[20], el1_eta[20], el1_phi[20], el1_dr[20], el1_e[20];
  float el1_eopin[20], el1_eopout[20], el1_hoe[20], el1_detain[20], el1_dphiin[20];
  float el1_fbrem[20], el1_eseed[20], el1_e3x3[20], el1_detaout[20], el1_dphiout[20];
  float el1_e5x5[20], el1_spp[20], el1_see[20], el1_pout[20];	
  int mc_id[20], sc_type[20], mc_mother[20], mc_crack[20];
  int tk_layer[20], tk_subdet[20], tk_nhit[20];
  int el_class[20], el1_class[20];
  int el_tkhits[20], el_pixelhits[20], el1_pixelhits[20], el1_tkhits[20];
  float el_gsfpt[20], el_gsfeta[20], el_scpt[20], el_sceta[20], el_scphi[20];
  float el1_gsfpt[20], el1_gsfeta[20], el1_scpt[20], el1_sceta[20], el1_scphi[20];

  HepMC::GenEvent* myGenEvent;
};
#endif

