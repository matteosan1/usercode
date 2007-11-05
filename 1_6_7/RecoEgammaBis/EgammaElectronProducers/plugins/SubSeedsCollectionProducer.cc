#include <memory>
#include <string>
#include "RecoEgammaBis/EgammaElectronProducers/plugins/SubSeedsCollectionProducer.h"

//#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/Common/interface/OwnVector.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"

#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"


#include "RecoTracker/CkfPattern/interface/TransientInitialStateEstimator.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"

#include <TMath.h>
#include <sstream>

#include <Math/VectorUtil.h>
#include <Math/Point3D.h>


using namespace reco;
using namespace edm;
using namespace std;

SubSeedsCollectionProducer::SubSeedsCollectionProducer(edm::ParameterSet const& conf) : 
  conf_(conf)
  //,theSeedsFilterWithCluster(0):
{  
  produces<TrajectorySeedCollection>();  
  }


// Virtual destructor needed.
SubSeedsCollectionProducer::~SubSeedsCollectionProducer() {
  //delete theSeedsFilterWithCluster;  
}  

void SubSeedsCollectionProducer::beginJob (EventSetup const & es)
{
  //services
  es.get<TrackerRecoGeometryRecord>().get( theGeomSearchTracker );
  es.get<IdealMagneticFieldRecord>().get(theMagField);
  
  /*
    // get nested parameter set for the TransientInitialStateEstimator
    ParameterSet tise_params = conf_.getParameter<ParameterSet>("TransientInitialStateEstimatorParameters") ;
    theInitialState          = new TransientInitialStateEstimator( es,tise_params);
    theNavigationSchool      = new SimpleNavigationSchool(&(*theGeomSearchTracker),&(*theMagField));
    theTrajectoryCleaner     = new TrajectoryCleanerBySharedHits();          
    
    
    // set the TrajectoryBuilder
    std::string trajectoryBuilderName = conf_.getParameter<std::string>("TrajectoryBuilder");
    edm::ESHandle<TrackerTrajectoryBuilder> theTrajectoryBuilderHandle;
    es.get<CkfComponentsRecord>().get(trajectoryBuilderName,theTrajectoryBuilderHandle);
    theTrajectoryBuilder = theTrajectoryBuilderHandle.product();    
    */
}

// Functions that gets called by framework every event
void SubSeedsCollectionProducer::produce(edm::Event& e, const edm::EventSetup& es)
{        
  // Step A: Retrieve seeds and clusters from the Event
  std::string seedProducer = conf_.getParameter<std::string>("SeedProducer");
  std::string seedLabel = conf_.getParameter<std::string>("SeedLabel");
  edm::Handle<TrajectorySeedCollection> theInitialSeedColl;
  e.getByLabel(seedProducer, seedLabel, theInitialSeedColl);
  //TrajectorySeedCollection theSeedColl = *collseed;
  
  edm::Handle<SuperClusterCollection> superClustersBarrelH; 
  e.getByLabel("correctedHybridSuperClusters",superClustersBarrelH);
  
  edm::Handle<SuperClusterCollection> superClustersEndcapH; 
  e.getByLabel("islandSuperClusters", "islandEndcapSuperClusters", superClustersEndcapH);

  edm::ESHandle<TrackerGeometry> tracker;
  es.get<TrackerDigiGeometryRecord>().get(tracker);

  // Step B: Create empty output collection
  std::auto_ptr<TrajectorySeedCollection> output(new TrajectorySeedCollection);    
    
  reco::SuperClusterRefVector superClusters;
  for(int z=0; z<2; ++z) {
    
    superClusters.clear();
    if (z == 0) {
      for(reco::SuperClusterCollection::size_type i= 0; i<superClustersBarrelH->size(); ++i){
        reco::SuperClusterRef cluster(superClustersBarrelH, i);
        superClusters.push_back(cluster);
      }
      //std::cout << superClustersBarrelH->size() << std::endl;
    }
    
    if (z == 1) {
      for(reco::SuperClusterCollection::size_type i= 0; i<superClustersEndcapH->size(); ++i){
        reco::SuperClusterRef cluster(superClustersEndcapH, i);
        superClusters.push_back(cluster);
      }
      //std::cout << superClustersEndcapH->size() << std::endl;
    }
    
    //seeds selection
    for(unsigned int i=0; i< superClusters.size(); ++i) {
      reco::SuperClusterRef theClus = superClusters[i];

      double minDr = 0.3;
      std::vector<TrajectorySeed>::const_iterator seed_iter;
      for(seed_iter = theInitialSeedColl->begin(); seed_iter != theInitialSeedColl->end(); ++seed_iter) {

        GlobalPoint  gp = tracker->idToDet( DetId(seed_iter->startingState().detId()))->surface().toGlobal( seed_iter->startingState().parameters().position());
        GlobalVector gv = tracker->idToDet( DetId(seed_iter->startingState().detId()))->surface().toGlobal( seed_iter->startingState().parameters().momentum());
        
        math::XYZVector trackGlobalDir(gv.x(),gv.y(),gv.z());   
        math::XYZVector clusterGlobalPos(theClus->x() - gp.x(), theClus->y() - gp.y(), theClus->z() - gp.z());
        
        double tmpDr = ROOT::Math::VectorUtil::DeltaR(clusterGlobalPos, trackGlobalDir);
        if(tmpDr <= minDr) output->push_back(*seed_iter); 
        /*
          if(tmpDr < minDr){
          minDr = tmpDr;
          theTrack = track;
          }
        */
      }
      
    }//end loop over cluster
    
  }
  
  
           
  edm::LogVerbatim("myElectronProd") << "========== SubSeedsCollectionProducer Info ==========";
  //edm::ESHandle<TrackerGeometry> tracker;
  //es.get<TrackerDigiGeometryRecord>().get(tracker);
  edm::LogVerbatim("myElectronProd") << "number of initial seeds: " << theInitialSeedColl->size();
  edm::LogVerbatim("myElectronProd") << "number of filtered seeds: " << output->size();
  
  edm::LogVerbatim("myElectronProd") << "=================================================";
  
    // Step G: write output to file
  e.put(output);
}





