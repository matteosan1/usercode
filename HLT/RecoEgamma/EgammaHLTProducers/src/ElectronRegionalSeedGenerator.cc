//
// Package:         RecoEgamma/EgammaHLTProducers
// Class:           EgammaHLTRegionalPixelSeedGeneratorProducers
//  Modified from TkSeedGeneratorFromTrk by Jeremy Werner, Princeton University, USA
// $Id: EgammaHLTRegionalPixelSeedGeneratorProducers.cc,v 1.8 2007/09/19 23:39:01 ratnik Exp $
//

#include <iostream>
#include <memory>
#include <string>

#include "RecoEgamma/EgammaHLTProducers/interface/ElectronRegionalSeedGenerator.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "RecoTracker/TkTrackingRegions/interface/RectangularEtaPhiTrackingRegion.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Math/interface/Vector3D.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"
#include "DataFormats/EgammaCandidates/interface/Electron.h"
#include "DataFormats/EgammaCandidates/interface/ElectronFwd.h"
#include "DataFormats/EgammaReco/interface/SuperCluster.h"

#include "RecoTracker/TkTrackingRegions/interface/OrderedHitsGeneratorFactory.h"
#include "RecoTracker/TkTrackingRegions/interface/OrderedHitsGenerator.h"
#include "RecoTracker/TkSeedGenerator/interface/SeedGeneratorFromRegionHits.h"

// Math
#include "Math/GenVector/VectorUtil.h"
#include "Math/GenVector/PxPyPzE4D.h"

using namespace std;
using namespace reco;

ElectronRegionalPixelSeedGenerator::ElectronRegionalPixelSeedGenerator(edm::ParameterSet const& conf):conf_(conf) {

  edm::LogInfo ("ElectronRegionalPixelSeedGenerator")<<"Enter the ElectronRegionalPixelSeedGenerator";
  
  edm::ParameterSet regionPSet = conf_.getParameter<edm::ParameterSet>("RegionPSet");
  

  ptmin_        = regionPSet.getParameter<double>("ptMin");
  originRadius_ = regionPSet.getParameter<double>("originRadius");
  halfLength_   = regionPSet.getParameter<double>("originHalfLength");
  deltaEta_     = regionPSet.getParameter<double>("deltaEtaRegion");
  deltaPhi_     = regionPSet.getParameter<double>("deltaPhiRegion");

  candTag_      = regionPSet.getParameter<edm::InputTag>("candTag");
  candTagEle_   = regionPSet.getParameter<edm::InputTag>("candTagEle");
  vertexSrc_    = regionPSet.getParameter<edm::InputTag>("vertexSrc");
}

std::vector<TrackingRegion* > ElectronRegionalPixelSeedGenerator::regions(const edm::Event& e, const edm::EventSetup& es) const {

  std::vector<TrackingRegion*> result;
  float originZ = 0.;
  double deltaZVertex;
	  
  // get Inputs 
  // Get the recoEcalCandidates
  edm::Handle<reco::RecoEcalCandidateCollection> recoecalcands;
  e.getByLabel(candTag_,recoecalcands);

  // get the primary vertex
  edm::Handle<reco::VertexCollection> h_vertices;
  e.getByLabel(vertexSrc_, h_vertices);
  const reco::VertexCollection & vertices = * h_vertices;
  if (! vertices.empty()) {
    originZ      = vertices.front().z();
    deltaZVertex = halfLength_;
  } else {
    originZ      =  0.;
    deltaZVertex = 15.;
  }
  
  //Get the HLT electrons collection if needed
  //edm::Handle<reco::ElectronCollection> electronHandle;
  //if(useZvertex_){
  // e.getByLabel(candTagEle_,electronHandle);
  //}
  
  reco::SuperClusterRef scRef;
  for (reco::RecoEcalCandidateCollection::const_iterator recoecalcand= recoecalcands->begin(); recoecalcand!=recoecalcands->end(); recoecalcand++) {
    scRef = recoecalcand->superCluster();

    //GlobalVector dirVector((recoecalcand)->px(),(recoecalcand)->py(),(recoecalcand)->pz());
    math::XYZVector clusterGlobalDir(scRef->x(), scRef->y(), scRef->z());
    
    double clusEt = scRef->energy()*sin(clusterGlobalDir.theta());
    double clusEstimatedCurvature = clusEt/0.3/4*100;  //4 tesla (temporary solution)
    double DphiBending = scRef->position().rho()/2./clusEstimatedCurvature; //ecal radius
      
    math::RhoEtaPhiVector dirVectorPlus(clusterGlobalDir);
    dirVectorPlus.SetPhi(clusterGlobalDir.phi() + DphiBending);

    math::RhoEtaPhiVector dirVectorMinus(clusterGlobalDir);
    dirVectorMinus.SetPhi(clusterGlobalDir.phi() - DphiBending);

    GlobalVector vectPlus(dirVectorPlus.x(), dirVectorPlus.y(), dirVectorPlus.z());
    GlobalVector vectMinus(dirVectorMinus.x(), dirVectorMinus.y(), dirVectorMinus.z());
    
    // FIXME
    float dphi = 0;
    if (clusEt > 15)
      dphi = 0.05;
    else
      dphi = 0.1;

    RectangularEtaPhiTrackingRegion* etaphiRegionPlus = new RectangularEtaPhiTrackingRegion(vectPlus,
                                                                                            GlobalPoint(0, 0, originZ), 
                                                                                            ptmin_,
                                                                                            originRadius_,
                                                                                            deltaZVertex,
                                                                                            deltaEta_,
                                                                                            dphi);

    result.push_back(etaphiRegionPlus);

    RectangularEtaPhiTrackingRegion* etaphiRegionMinus = new RectangularEtaPhiTrackingRegion(vectMinus,
                                                                                             GlobalPoint(0, 0, originZ), 
                                                                                             ptmin_,
                                                                                             originRadius_,
                                                                                             deltaZVertex,
                                                                                             deltaEta_,
                                                                                             dphi);

    result.push_back(etaphiRegionMinus);

  }

  return result; 
}
