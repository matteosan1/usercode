#ifndef ElectronRegionalPixelSeedGenerator_h
#define ElectronRegionalPixelSeedGenerator_h

//
// Class: ElectronRegionalPixelSeedGenerator

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/InputTag.h"
#include "DataFormats/Common/interface/EDProduct.h"
#include "RecoTracker/TkTrackingRegions/interface/TrackingRegionProducer.h"
#include "RecoTracker/TkTrackingRegions/interface/GlobalTrackingRegion.h"
#include "RecoTracker/TkTrackingRegions/interface/RectangularEtaPhiTrackingRegion.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/JetReco/interface/CaloJetCollection.h"
#include "DataFormats/Math/interface/Vector3D.h"
// Math
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/JetReco/interface/Jet.h"

using namespace reco;

class ElectronRegionalPixelSeedGenerator : public TrackingRegionProducer {
 public:
    
  ElectronRegionalPixelSeedGenerator(const edm::ParameterSet& conf_);
  virtual ~ElectronRegionalPixelSeedGenerator() {};
    

  std::vector<TrackingRegion* > regions(const edm::Event& e, const edm::EventSetup& es) const;
  
 private:
  edm::ParameterSet conf_;
  
  float ptmin_;
  float originRadius_;
  float halfLength_;
  float deltaEta_;
  float deltaPhi_;
  edm::InputTag candTag_, candTagEle_, vertexSrc_;
};

#endif
