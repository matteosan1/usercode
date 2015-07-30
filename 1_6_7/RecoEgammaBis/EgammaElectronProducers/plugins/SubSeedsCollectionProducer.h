#ifndef SubSeedsCollectionProducer_h
#define SubSeedsCollectionProducer_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "DataFormats/Common/interface/EDProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

//#include "RecoTracker/CkfPattern/interface/TrackerTrajectoryBuilder.h"
//#include "TrackingTools/TrajectoryCleaning/interface/TrajectoryCleaner.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"

//class TransientInitialStateEstimator;


class SubSeedsCollectionProducer : public edm::EDProducer        
{
 public:
  
  explicit SubSeedsCollectionProducer(const edm::ParameterSet& conf);
  
  virtual ~SubSeedsCollectionProducer();
  
  virtual void beginJob (edm::EventSetup const & es);
  
  virtual void produce(edm::Event& e, const edm::EventSetup& es);
  
 private:
  edm::ParameterSet conf_;
  //const SeedsFilterWithCluster*  theSeedsFilterWithCluster;
  
  edm::ESHandle<MagneticField>                theMagField;
  edm::ESHandle<GeometricSearchTracker>       theGeomSearchTracker;
};


#endif
