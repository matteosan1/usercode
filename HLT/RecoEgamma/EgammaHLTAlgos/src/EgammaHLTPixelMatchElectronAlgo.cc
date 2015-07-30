// -*- C++ -*-
//
// Package:    EgammaHLTAlgos
// Class:      EgammaHLTPixelMatchElectronAlgo.
// 
/**\class EgammaHLTPixelMatchElectronAlgo EgammaHLTAlgos/EgammaHLTPixelMatchElectronAlgo

 Description: top algorithm producing TrackCandidate and Electron objects from supercluster
              driven pixel seeded Ckf tracking for HLT
*/
//
// Original Author:  Monica Vazquez Acosta (CERN)
// $Id: EgammaHLTPixelMatchElectronAlgo.cc,v 1.10 2007/12/10 18:12:46 ghezzi Exp $
//
//
#include "RecoEgamma/EgammaHLTAlgos/interface/EgammaHLTPixelMatchElectronAlgo.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"

#include "RecoTracker/Record/interface/TrackerRecoGeometryRecord.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/SuperClusterFwd.h"
#include "DataFormats/EgammaReco/interface/ElectronPixelSeed.h"
#include "DataFormats/EgammaReco/interface/ElectronPixelSeedFwd.h"

#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"


using namespace edm;
using namespace std;
using namespace reco;
//using namespace math; // conflicts with DataFormat/Math/interface/Point3D.h!!!!

//EgammaHLTPixelMatchElectronAlgo::EgammaHLTPixelMatchElectronAlgo():  
//  theCkfTrajectoryBuilder(0), theTrajectoryCleaner(0),
//  theInitialStateEstimator(0), theNavigationSchool(0) {}

EgammaHLTPixelMatchElectronAlgo::EgammaHLTPixelMatchElectronAlgo(){}
EgammaHLTPixelMatchElectronAlgo::~EgammaHLTPixelMatchElectronAlgo() {

  // delete theInitialStateEstimator;
  //delete theNavigationSchool;
  //delete theTrajectoryCleaner; 
    
}

void EgammaHLTPixelMatchElectronAlgo::setupES(const edm::EventSetup& es, const edm::ParameterSet &conf) {

  //services
  es.get<TrackerRecoGeometryRecord>().get( theGeomSearchTracker );
  es.get<IdealMagneticFieldRecord>().get(theMagField);

  /*
  // get nested parameter set for the TransientInitialStateEstimator
  ParameterSet tise_params = conf.getParameter<ParameterSet>("TransientInitialStateEstimatorParameters") ;
  theInitialStateEstimator       = new TransientInitialStateEstimator( es,tise_params);

  edm::ESHandle<NavigationSchool> nav;
  es.get<NavigationSchoolRecord>().get("SimpleNavigationSchool", nav);
  theNavigationSchool = nav.product();
  // set the correct navigation
  NavigationSetter setter(*theNavigationSchool);

  //  theCkfTrajectoryBuilder = new CkfTrajectoryBuilder(conf,es,theMeasurementTracker);
  theTrajectoryCleaner = new TrajectoryCleanerBySharedHits();    
  std::string trajectoryBuilderName = conf.getParameter<std::string>("TrajectoryBuilder");
  edm::ESHandle<TrajectoryBuilder> theTrajectoryBuilderHandle;
  es.get<CkfComponentsRecord>().get(trajectoryBuilderName,theTrajectoryBuilderHandle);
  theCkfTrajectoryBuilder = theTrajectoryBuilderHandle.product();    
  */

  trackLabel_ = conf.getParameter<string>("TrackLabel");
  trackInstanceName_ = conf.getParameter<string>("TrackProducer");
}

void  EgammaHLTPixelMatchElectronAlgo::run(Event& e, ElectronCollection & outEle) {

  // get the input 
  edm::Handle<TrackCollection> tracksH;
  e.getByLabel(trackLabel_,trackInstanceName_,tracksH);

  process(tracksH,outEle);

  return;
}

void EgammaHLTPixelMatchElectronAlgo::process(edm::Handle<TrackCollection> tracksH, ElectronCollection & outEle) {
  const TrackCollection *tracks=tracksH.product();
  for (unsigned int i=0;i<tracks->size();++i) {
    const Track & t=(*tracks)[i];

    const TrackRef trackRef = edm::Ref<TrackCollection>(tracksH,i);
    edm::RefToBase<TrajectorySeed> seed = trackRef->extra()->seedRef();
    ElectronPixelSeedRef elseed=seed.castTo<ElectronPixelSeedRef>();
    const SuperClusterRef & scRef=elseed->superCluster();
    // Get the momentum at vertex (not at the innermost layer)
    TSCPBuilderNoMaterial tscpBuilder;
    TrajectoryStateTransform tsTransform;
    FreeTrajectoryState fts = tsTransform.innerFreeState(t,theMagField.product());
    TrajectoryStateClosestToPoint tscp = tscpBuilder(fts, Global3DPoint(0,0,0) );
    
    float scale = scRef->energy()/tscp.momentum().mag();
  
    const math::XYZTLorentzVector momentum(tscp.momentum().x()*scale,
 					   tscp.momentum().y()*scale,
 					   tscp.momentum().z()*scale,
					   scRef->energy());

    
    Electron ele(t.charge(),momentum, t.vertex() );
    ele.setSuperCluster(scRef);
    edm::Ref<TrackCollection> myRef(tracksH,i);
    ele.setTrack(myRef);
    outEle.push_back(ele);

  }  // loop over tracks
}

