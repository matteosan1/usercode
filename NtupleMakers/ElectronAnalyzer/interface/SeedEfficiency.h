#ifndef SeedEfficiency_H
#define SeedEfficiency_H

#include "FWCore/Framework/interface/EDAnalyzer.h"

//#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
//#include "DataFormats/EgammaCandidates/interface/GlobalCtfElectron.h"

//#include "DataFormats/TrackReco/interface/Track.h"

#include "DataFormats/Math/interface/LorentzVector.h"
#include <Math/VectorUtil.h>
#include <Math/Point3D.h>
#include <vector>
//#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
//#include "CLHEP/HepMC/GenParticle.h"

namespace edm {
  class ParameterSet;
  class Event;
  class EventSetup;
}

class TFile;
class TTree;

class SeedEfficiency: public edm::EDAnalyzer {
public:

  SeedEfficiency(const edm::ParameterSet& pset);

  virtual ~SeedEfficiency();

  void analyze(const edm::Event & event, const edm::EventSetup& eventSetup);
  void beginJob(const edm::EventSetup& eventSetup);
  void endJob();

 protected:
  
 private:
  std::string trajectorySeedName, simHitName, fileName;
 
  TFile *file;
  TTree *tree;
  
  int nAss[3];
  int nSeed;

};
#endif

