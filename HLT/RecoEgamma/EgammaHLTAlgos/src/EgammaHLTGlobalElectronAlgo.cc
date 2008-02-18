#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "RecoEgamma/EgammaHLTAlgos/interface/EgammaHLTGlobalElectronAlgo.h"

#include "TrackingTools/PatternTools/interface/TSCPBuilderNoMaterial.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TrajectoryState/interface/FreeTrajectoryState.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateClosestToPoint.h"

#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "CLHEP/Units/PhysicalConstants.h"
#include <TMath.h>
#include <Math/VectorUtil.h>

#include <sstream>

using namespace edm;
using namespace std;
using namespace reco;
//using namespace math; // conflicts with DataFormat/Math/interface/Point3D.h!!!!

EgammaHLTGlobalElectronAlgo::EgammaHLTGlobalElectronAlgo():  
  theCkfTrajectoryBuilder(0), theTrajectoryCleaner(0),
  theInitialStateEstimator(0), theNavigationSchool(0) {}

EgammaHLTGlobalElectronAlgo::~EgammaHLTGlobalElectronAlgo() {
  
  delete theInitialStateEstimator;
  delete theNavigationSchool;
  delete theTrajectoryCleaner; 
  
}

void EgammaHLTGlobalElectronAlgo::setupES(const edm::EventSetup& es, const edm::ParameterSet &conf) {
  
  //services
  es.get<TrackerRecoGeometryRecord>().get( theGeomSearchTracker );
  es.get<IdealMagneticFieldRecord>().get(theMagField);
  
  // get nested parameter set for the TransientInitialStateEstimator
  ParameterSet tise_params = conf.getParameter<ParameterSet>("TransientInitialStateEstimatorParameters") ;
  theInitialStateEstimator       = new TransientInitialStateEstimator( es,tise_params);
  
  theNavigationSchool   = new SimpleNavigationSchool(&(*theGeomSearchTracker),&(*theMagField));
  
  // set the correct navigation
  NavigationSetter setter( *theNavigationSchool);
  
  //  theCkfTrajectoryBuilder = new CkfTrajectoryBuilder(conf,es,theMeasurementTracker);
  theTrajectoryCleaner = new TrajectoryCleanerBySharedHits();    
  //std::string trajectoryBuilderName = conf.getParameter<std::string>("TrajectoryBuilder");
  std::string trajectoryBuilderName =  "electronRegionalTrajectoryBuilder";

  edm::ESHandle<TrajectoryBuilder> theTrajectoryBuilderHandle;
  es.get<CkfComponentsRecord>().get(trajectoryBuilderName,theTrajectoryBuilderHandle);
  theCkfTrajectoryBuilder = theTrajectoryBuilderHandle.product();    
  
  // prendi tracce e SC

  //trackLabel_ = conf.getParameter<string>("TrackLabel");
  trackLabel_ = "electronRegionalCtfWithMaterialTracks";
  	
  //trackInstanceName_ = conf.getParameter<string>("TrackProducer");
  //scLabel_ = conf.getParameter< edm::InputTag >("SCLabel");
  //scInstanceName_ = conf.getParameter<string>("SCProducer");
  //scEndcapLabel_ = conf.getParameter<string>("SCLEndcapLabel");
  //scEndcapInstanceName_ = conf.getParameter<string>("SCLEndcapProducer");
}

void  EgammaHLTGlobalElectronAlgo::run(Event& e, ElectronCollection& outEle) {

  // get the input 
  edm::Handle<TrackCollection> tracksH;
  edm::Handle<reco::RecoEcalCandidateCollection> scH;
  //edm::Handle<SuperClusterCollection> scEndcapH;

  e.getByLabel(trackLabel_, trackInstanceName_, tracksH);
  e.getByLabel("l1IsoRecoEcalCandidate", scH);
  const reco::RecoEcalCandidateCollection* col = scH.product();
  //e.getByLabel(scEndcapLabel_, scEndcapInstanceName_, scEndcapH);

  for (reco::RecoEcalCandidateCollection::const_iterator recoecalcand = col->begin(); recoecalcand!=col->end(); 
       recoecalcand++) {

    const reco::SuperClusterRef theClus = recoecalcand->superCluster();

    process(tracksH, theClus, outEle);
  }  
  //std::cout << "Elettroni reco: " << outEle.size() << std::endl;  
  return;
}

void EgammaHLTGlobalElectronAlgo::process(edm::Handle<TrackCollection>& tracksH, const SuperClusterRef& sc, ElectronCollection & outEle) {

  reco::TrackRef theTrack = edm::Ref<reco::TrackCollection>();
  math::XYZTLorentzVector theMomentum;

  double minDr = 0.3;
  double minDeop = 5;

  for(reco::TrackCollection::size_type i=0; i<tracksH->size(); ++i){

    reco::TrackRef track(tracksH, i);
    math::XYZVector trackGlobalDir(track->momentum());
    math::XYZVector clusterGlobalDir(sc->x() - track->vx(), sc->y() - track->vy(), sc->z() - track->vz());

    double clusEt = sc->energy()*sin(clusterGlobalDir.theta());
    double clusEstimatedCurvature = clusEt/0.3/4*100;  //4 tesla (temporary solution)
    double DphiBending = sc->position().rho()/2./clusEstimatedCurvature; //ecal radius

    double tmpDr = ROOT::Math::VectorUtil::DeltaR(clusterGlobalDir, trackGlobalDir);
    if (!(tmpDr < minDr)) 
      continue;
    /*
    TrajectoryStateOnSurface innTSOS = mtsTransform_->innerStateOnSurface(*track, *(trackerHandle_.product()), theMagField.product());
    GlobalVector innMom = computeMode(innTSOS);
    
    TrajectoryStateOnSurface outTSOS = mtsTransform_->outerStateOnSurface(*track, *(trackerHandle_.product()), theMagField.product());
    if (!outTSOS.isValid())   continue;
    
    TrajectoryStateOnSurface seedTSOS = TransverseImpactPointExtrapolator(*geomPropFw_).extrapolate(outTSOS,GlobalPoint(sc->seed()->position().x(),sc->seed()->position().y(),sc->seed()->position().z()));

    if (!seedTSOS.isValid()) 
      seedTSOS=outTSOS;

    GlobalVector seedMom = computeMode(seedTSOS);
    */
    // Get the momentum at vertex (not at the innermost layer)

    TSCPBuilderNoMaterial tscpBuilder;
    TrajectoryStateTransform tsTransform;
    FreeTrajectoryState fts = tsTransform.innerFreeState(*track, theMagField.product());
    TrajectoryStateClosestToPoint tscp = tscpBuilder(fts, Global3DPoint(0,0,0));
    
    float scale = sc->energy()/tscp.momentum().mag();
    //std::cout <<  "Scale :" << scale << " = " << sc->energy() << " / " << tscp.momentum().mag() << std::endl;
   
    const math::XYZTLorentzVector momentum(tscp.momentum().x()*scale,
                                           tscp.momentum().y()*scale,
                                           tscp.momentum().z()*scale,
                                           sc->energy());

    //double eOverPin  = sc->energy()/momentum.P();
    double Deta = fabs(clusterGlobalDir.eta() - trackGlobalDir.eta());
    double dPhi = fabs(acos(cos(clusterGlobalDir.phi() - trackGlobalDir.phi())));
    float dPhi1 = fabs(dPhi - DphiBending);
    float dPhi2 = fabs(dPhi + DphiBending);
    //std::cout << "eOverPin : " << eOverPin << " =  " << sc->energy() << " / " << momentum.P() << std::endl;
    //std::cout << "Dphi: " << Dphi << " Deta: " << Deta << std::endl;

    if(!(scale<5))  
      continue;
    if(!(min(dPhi1, dPhi2) < 0.1))  
      continue;
    if(!(Deta < 0.02)) 
      continue;

    //    cout << " in matchbox, dphi, deta: " << Dphi << " , " << Deta << endl;
    //    cout << " in matchbox, E/Pin, out: " << eOverPin << " , " << eOverPout << endl;
    
    if( fabs(scale-1.) < minDeop){
      minDeop = fabs(scale-1.) ;
      theTrack = track;
      theMomentum = momentum;
    }
  }
  
  //const TrackRef trackRef = edm::Ref<TrackCollection>(tracksH,i);
  //edm::RefToBase<TrajectorySeed> seed = trackRef->extra()->seedRef();
  //ElectronPixelSeedRef elseed=seed.castTo<ElectronPixelSeedRef>();
  //const SuperClusterRef & scRef=elseed->superCluster();
  // Get the momentum at vertex (not at the innermost layer)
  /*
  TSCPBuilderNoMaterial tscpBuilder;
  TrajectoryStateTransform tsTransform;
  FreeTrajectoryState fts = tsTransform.innerFreeState(t,theMagField.product());
  TrajectoryStateClosestToPoint tscp = tscpBuilder(fts, Global3DPoint(0,0,0) );
  
  float scale = scRef->energy()/tscp.momentum().mag();
  
  const math::XYZTLorentzVector momentum(tscp.momentum().x()*scale,
                                         tscp.momentum().y()*scale,
                                         tscp.momentum().z()*scale,
                                         scRef->energy());
  */
  if(!theTrack.isNull()) {
    Electron ele(theTrack->charge(), theMomentum, theTrack->vertex());
    ele.setSuperCluster(sc);
    ele.setTrack(theTrack);
    outEle.push_back(ele);
  }
}

