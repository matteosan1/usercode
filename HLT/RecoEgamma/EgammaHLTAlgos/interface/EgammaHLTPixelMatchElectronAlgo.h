#ifndef EgammaHLTPixelMatchElectronAlgo_H
#define EgammaHLTPixelMatchElectronAlgo_H

/** \class EgammaHLTPixelMatchElectronAlgo
 
 * Class to reconstruct electron tracks from electron pixel seeds
 *  keep track of information about the initiating supercluster
 *
 * \author Monica Vazquez Acosta (CERN)
 *
 ************************************************************/


#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
/*
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronPixelSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronPixelSeedFwd.h"
//#include "DataFormats/EgammaReco/interface/SeedSuperClusterAssociation.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/Math/interface/LorentzVector.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
*/

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoTracker/TkDetLayers/interface/GeometricSearchTracker.h"

/*
#include "TrackingTools/TrajectoryCleaning/interface/TrajectoryCleaner.h"
//#include "TrackingTools/DetLayers/interface/NavigationSetter.h"
//#include "TrackingTools/DetLayers/interface/NavigationSchool.h"


#include "TrackingTools/PatternTools/interface/TrajectoryBuilder.h"
//#include "RecoTracker/TkNavigation/interface/SimpleNavigationSchool.h"
#include "RecoTracker/MeasurementDet/interface/MeasurementTracker.h"
#include "RecoTracker/CkfPattern/interface/TransientInitialStateEstimator.h"
#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "RecoTracker/Record/interface/CkfComponentsRecord.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "TrackingTools/MaterialEffects/interface/PropagatorWithMaterial.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TrajectoryCleaning/interface/TrajectoryCleanerBySharedHits.h"
#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/PatternTools/interface/TransverseImpactPointExtrapolator.h"
#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/PatternTools/interface/TrajectoryFitter.h"

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
*/


//class TransientInitialStateEstimator;
//class NavigationSchool;


using namespace std;
using namespace edm;
using namespace reco;

class EgammaHLTPixelMatchElectronAlgo {

public:

  EgammaHLTPixelMatchElectronAlgo( );

  ~EgammaHLTPixelMatchElectronAlgo();

  void setupES(const EventSetup& setup, const ParameterSet& conf);
  void run(Event&, ElectronCollection&);

 private:

  // create electrons from tracks
  void process(edm::Handle<TrackCollection> tracksH, ElectronCollection & outEle);
  


  // input configuration
  std::string trackLabel_;
  std::string trackInstanceName_;

  //  const TrajectoryBuilder*  theCkfTrajectoryBuilder;
  //TrajectoryCleaner*               theTrajectoryCleaner;
  //TransientInitialStateEstimator*  theInitialStateEstimator;
  
  ESHandle<MagneticField>                theMagField;
  ESHandle<GeometricSearchTracker>       theGeomSearchTracker;

  //const MeasurementTracker*     theMeasurementTracker;
  //const NavigationSchool*       theNavigationSchool;

};

#endif // EgammaHLTPixelMatchElectronAlgo_H


