/*
  logical flow of the selection:
  1) save all SC with Et > 5. GeV
  2) for each saved SC search for matched following object (within dR = 0.1):
  a) Standard electron
  b) Our electron
  c) MC particle
  d) CTF track
  e) nearest MC particle to stored CTF track
*/

#include "NtupleMakers/ElectronAnalyzer/interface/Conversion.h"

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"

#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "CLHEP/HepMC/GenParticle.h"
#include "DataFormats/EgammaReco/interface/BasicClusterShapeAssociation.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "DataFormats/EgammaCandidates/interface/PixelMatchGsfElectron.h"
#include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHit.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"
#include "DataFormats/SiStripDetId/interface/TIDDetId.h"
#include "DataFormats/SiStripDetId/interface/TECDetId.h"
#include "DataFormats/SiStripDetId/interface/TIBDetId.h"
#include "DataFormats/SiStripDetId/interface/TOBDetId.h"

#include "DataFormats/SiPixelDetId/interface/PXBDetId.h"
#include "DataFormats/SiPixelDetId/interface/PXFDetId.h"

#include "DataFormats/EgammaReco/interface/SuperCluster.h"
#include "DataFormats/TrackReco/interface/Track.h"


//__________________________added by PDK_______________________
#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
//_____________________________________________________________


#include "TFile.h"
#include "TTree.h"

using namespace std;
using namespace edm;
using namespace reco;

/// Constructor
Conversion::Conversion(const ParameterSet& pset) {
  fileName = pset.getParameter<std::string>("RootFileName");
  baselineEleCollName =  pset.getParameter<std::string>("BaselineEleCollName");
  customEleCollName   =  pset.getParameter<std::string>("CustomEleCollName");
}

Conversion::~Conversion() {}

void Conversion::beginJob(const EventSetup& eventSetup) {

  file = new TFile(fileName.c_str(), "recreate");
  tree = new TTree("event","Event data");

  tree->Branch("run", &run, "run/I");
  tree->Branch("id", &id, "id/I");

  // tree->Branch("sc_n", &nSC, "sc_n/I");
  tree->Branch("sc_e", &sc_e, "sc_e/F");
  tree->Branch("sc_rawe", &sc_rawe, "sc_rawe/F");
  tree->Branch("sc_et", &sc_et, "sc_et/F");
  tree->Branch("sc_eta", &sc_eta, "sc_eta/F");
  tree->Branch("sc_phi", &sc_phi, "sc_phi/F");
  tree->Branch("sc_dr", &sc_dr, "sc_dr/F");
  tree->Branch("sc_type", &sc_type, "sc_type/I"); // 0 barrel, 1 endcap

  tree->Branch("mcsc_pt", &mcsc_pt, "mcsc_pt/F");
  tree->Branch("mcsc_dr", &mcsc_dr, "mcsc_dr/F");
  tree->Branch("mcsc_eta", &mcsc_eta, "mcsc_eta/F");
  tree->Branch("mcsc_phi", &mcsc_phi, "mcsc_phi/F");
  tree->Branch("mcsc_id", &mcsc_id, "mcsc_id/I");
  tree->Branch("mcsc_e", &mcsc_e, "mcsc_e/F");
  tree->Branch("mcsc_mother", &mcsc_mother, "mcsc_mother/I");
  tree->Branch("mcsc_crack", &mcsc_crack, "mcsc_crack/I"); // 0 no crack, 1 crack

  tree->Branch("mctk_pt", &mctk_pt, "mctk_pt/F");
  tree->Branch("mctk_dr", &mctk_dr, "mctk_dr/F");
  tree->Branch("mctk_eta", &mctk_eta, "mctk_eta/F");
  tree->Branch("mctk_phi", &mctk_phi, "mctk_phi/F");
  tree->Branch("mctk_id", &mctk_id, "mctk_id/I");
  tree->Branch("mctk_e", &mctk_e, "mctk_e/F");
  tree->Branch("mctk_mother", &mctk_mother, "mctk_mother/I");
  tree->Branch("mctk_crack", &mctk_crack, "mctk_crack/I"); // 0 no crack, 1 crack

  tree->Branch("mctk1_pt", &mctk1_pt, "mctk1_pt/F");
  tree->Branch("mctk1_dr", &mctk1_dr, "mctk1_dr/F");
  tree->Branch("mctk1_eta", &mctk1_eta, "mctk1_eta/F");
  tree->Branch("mctk1_phi", &mctk1_phi, "mctk1_phi/F");
  tree->Branch("mctk1_id", &mctk1_id, "mctk1_id/I");
  tree->Branch("mctk1_e", &mctk1_e, "mctk1_e/F");
  tree->Branch("mctk1_mother", &mctk1_mother, "mctk1_mother/I");
  tree->Branch("mctk1_crack", &mctk1_crack, "mctk1_crack/I"); // 0 no crack, 1 crack

  tree->Branch("el_pt", &el_pt, "el_pt/F");
  tree->Branch("el_e", &el_e, "el_e/F");
  tree->Branch("el_q", &el_q, "el_q/I");
  tree->Branch("el_eta", &el_eta, "el_eta/F");
  tree->Branch("el_phi", &el_phi, "el_phi/F");
  tree->Branch("el_dr", &el_dr, "el_dr/F");
  tree->Branch("el_eopin", &el_eopin, "el_eopin/F");
  tree->Branch("el_eopout", &el_eopout, "el_eopout/F");
  tree->Branch("el_pout", &el_pout, "el_pout/F");
  tree->Branch("el_fbrem", &el_fbrem, "el_fbrem/F");
  tree->Branch("el_hoe", &el_hoe, "el_hoe/F");
  tree->Branch("el_detain", &el_detain, "el_detain/F");
  tree->Branch("el_dphiin", &el_dphiin, "el_dphiin/F");
  tree->Branch("el_detaout", &el_detaout, "el_detaout/F");
  tree->Branch("el_dphiout", &el_dphiout, "el_dphiout/F");
  tree->Branch("el_e3x3", &el_e3x3, "el_e3x3/F");
  tree->Branch("el_e5x5", &el_e5x5, "el_e5x5/F");
  tree->Branch("el_eseed", &el_eseed, "el_eseed/F");
  tree->Branch("el_spp", &el_spp, "el_spp/F");
  tree->Branch("el_see", &el_see, "el_see/F");
  tree->Branch("el_class", &el_class, "el_class/I");
  tree->Branch("el_nsihit", &el_nsihits, "el_nsihit/I");
  tree->Branch("el_npxhit", &el_npxhits, "el_npxhit/I");
  tree->Branch("el_xhit", &el_xhit, "el_xhit/F");
  tree->Branch("el_yhit", &el_yhit, "el_yhit/F");
  tree->Branch("el_zhit", &el_zhit, "el_zhit/F");
  tree->Branch("el_detinnerhit", &el_detinnerhit, "el_detinnerhit/I");
  tree->Branch("el_layerinnerhit", &el_detinnerhit, "el_layerinnerhit/I");
  tree->Branch("el_rinnerhit", &el_rinnerhit, "el_rinnerhit/F");
  tree->Branch("el_z0", &el_z0, "el_z0/F");
  tree->Branch("el_d0", &el_d0, "el1_d0/F");
  tree->Branch("el_d0err", &el_d0err, "el_d0err/F");
  tree->Branch("el_tkiso", &el_tkiso, "el_tkiso/F");
  tree->Branch("el_tkpt", &el_tkpt, "el_tkpt/F");
  tree->Branch("el_tketa", &el_tketa, "el_tketa/F");
  tree->Branch("el_tkphi", &el_tkphi, "el_tkphi/F");
  
  tree->Branch("el1_pt", &el1_pt, "el1_pt/F");
  tree->Branch("el1_e", &el1_e, "el1_e/F");
  tree->Branch("el1_q", &el1_q, "el1_q/I");
  tree->Branch("el1_eta", &el1_eta, "el1_eta/F");
  tree->Branch("el1_phi", &el1_phi, "el1_phi/F");
  tree->Branch("el1_dr", &el1_dr, "el1_dr/F");
  tree->Branch("el1_eopin", &el1_eopin, "el1_eopin/F");
  tree->Branch("el1_eopout", &el1_eopout, "el1_eopout/F");
  tree->Branch("el1_pout", &el1_pout, "el1_pout/F");
  tree->Branch("el1_fbrem", &el1_fbrem, "el1_fbrem/F");
  tree->Branch("el1_hoe", &el1_hoe, "el1_hoe/F");
  tree->Branch("el1_detain", &el1_detain, "el1_detain/F");
  tree->Branch("el1_dphiin", &el1_dphiin, "el1_dphiin/F");
  tree->Branch("el1_detaout", &el1_detaout, "el1_detaout/F");
  tree->Branch("el1_dphiout", &el1_dphiout, "el1_dphiout/F");
  tree->Branch("el1_e3x3", &el1_e3x3, "el1_e3x3/F");
  tree->Branch("el1_e5x5", &el1_e5x5, "el1_e5x5/F");
  tree->Branch("el1_eseed", &el1_eseed, "el1_eseed/F");
  tree->Branch("el1_spp", &el1_spp, "el1_spp/F");
  tree->Branch("el1_see", &el1_see, "el1_see/F");
  tree->Branch("el1_class", &el1_class, "el1_class/I");
  tree->Branch("el1_nsihit", &el1_nsihits, "el1_nsihit/I");
  tree->Branch("el1_npxhit", &el1_npxhits, "el1_npxhit/I");
  tree->Branch("el1_xhit", &el1_xhit, "el1_xhit/F");
  tree->Branch("el1_yhit", &el1_yhit, "el1_yhit/F");
  tree->Branch("el1_zhit", &el1_zhit, "el1_zhit/F");
  tree->Branch("el1_detinnerhit", &el1_detinnerhit, "el1_detinnerhit/I");
  tree->Branch("el1_layerinnerhit", &el1_detinnerhit, "el1_layerinnerhit/I");
  tree->Branch("el1_rinnerhit", &el1_rinnerhit, "el1_rinnerhit/F");
  tree->Branch("el1_z0", &el1_z0, "el1_z0/F");
  tree->Branch("el1_d0", &el1_d0, "el1_d0/F");
  tree->Branch("el1_d0err", &el1_d0err, "el1_d0err/F");
  tree->Branch("el1_tkiso", &el1_tkiso, "el1_tkiso/F");
  tree->Branch("el1_tkpt", &el1_tkpt, "el1_tkpt/F");
  tree->Branch("el1_tketa", &el1_tketa, "el1_tketa/F");
  tree->Branch("el1_tkphi", &el1_tkphi, "el1_tkphi/F");
    
  //Sim branches
  //px,py,pz and e of the photon closest to the standard electron's track
  tree->Branch("el_sim_gpx", &el_sim_gpx, "el_sim_gpx/F");
  tree->Branch("el_sim_gpy", &el_sim_gpy, "el_sim_gpy/F");
  tree->Branch("el_sim_gpz", &el_sim_gpz, "el_sim_gpz/F");
  tree->Branch("el_sim_gpt", &el_sim_gpt, "el_sim_gpt/F");
  tree->Branch("el_sim_ge", &el_sim_ge, "el_sim_ge/F");
  //position of the vertex where the conversion happened
  tree->Branch("el_sim_gvx", &el_sim_gvx, "el_sim_gvx/F");
  tree->Branch("el_sim_gvy", &el_sim_gvy, "el_sim_gvy/F");
  tree->Branch("el_sim_gvz", &el_sim_gvz, "el_sim_gvz/F");
  tree->Branch("el_sim_gvr", &el_sim_gvr, "el_sim_gvr/F");
  tree->Branch("el_sim_gveta", &el_sim_gveta, "el_sim_gveta/F");
  tree->Branch("el_sim_gvphi", &el_sim_gvphi, "el_sim_gvphi/F");
  //px,py,pz, e of the 2 decay electrons
  tree->Branch("el_sim_decay1px", &el_sim_decay1px, "el_sim_decay1px/F");
  tree->Branch("el_sim_decay1py", &el_sim_decay1py, "el_sim_decay1py/F");
  tree->Branch("el_sim_decay1pz", &el_sim_decay1pz, "el_sim_decay1pz/F");
  tree->Branch("el_sim_decay1pt", &el_sim_decay1pt, "el_sim_decay1pt/F");
  tree->Branch("el_sim_decay1e",  &el_sim_decay1e,  "el_sim_decay1e/F" );
  tree->Branch("el_sim_decay1pid",  &el_sim_decay1pid,  "el_sim_decay1pid/I" );
  tree->Branch("el_sim_decay2px", &el_sim_decay2px, "el_sim_decay2px/F");
  tree->Branch("el_sim_decay2py", &el_sim_decay2py, "el_sim_decay2py/F");
  tree->Branch("el_sim_decay2pz", &el_sim_decay2pz, "el_sim_decay2pz/F");
  tree->Branch("el_sim_decay2pt", &el_sim_decay2pt, "el_sim_decay2pt/F");
  tree->Branch("el_sim_decay2e",  &el_sim_decay2e,  "el_sim_decay2e/F" );
  tree->Branch("el_sim_decay2pid",  &el_sim_decay2pid,  "el_sim_decay2pid/I" );

  //UCSD electron
  tree->Branch("el1_sim_gpx", &el1_sim_gpx, "el1_sim_gpx/F");
  tree->Branch("el1_sim_gpy", &el1_sim_gpy, "el1_sim_gpy/F");
  tree->Branch("el1_sim_gpz", &el1_sim_gpz, "el1_sim_gpz/F");
  tree->Branch("el1_sim_gpt", &el1_sim_gpt, "el1_sim_gpt/F");
  tree->Branch("el1_sim_ge", &el1_sim_ge, "el1_sim_ge/F");
  //position of the vertex where the conversion happened
  tree->Branch("el1_sim_gvx", &el1_sim_gvx, "el1_sim_gvx/F");
  tree->Branch("el1_sim_gvy", &el1_sim_gvy, "el1_sim_gvy/F");
  tree->Branch("el1_sim_gvz", &el1_sim_gvz, "el1_sim_gvz/F");
  tree->Branch("el1_sim_gvr", &el1_sim_gvr, "el1_sim_gvr/F");
  tree->Branch("el1_sim_gveta", &el1_sim_gveta, "el1_sim_gveta/F");
  tree->Branch("el1_sim_gvphi", &el1_sim_gvphi, "el1_sim_gvphi/F");
  //px,py,pz, e of the 2 decay electrons
  tree->Branch("el1_sim_decay1px", &el1_sim_decay1px, "el1_sim_decay1px/F");
  tree->Branch("el1_sim_decay1py", &el1_sim_decay1py, "el1_sim_decay1py/F");
  tree->Branch("el1_sim_decay1pz", &el1_sim_decay1pz, "el1_sim_decay1pz/F");
  tree->Branch("el1_sim_decay1pt", &el1_sim_decay1pt, "el1_sim_decay1pt/F");
  tree->Branch("el1_sim_decay1e",  &el1_sim_decay1e,  "el1_sim_decay1e/F" );
  tree->Branch("el1_sim_decay1pid",  &el1_sim_decay1pid,  "el1_sim_decay1pid/I" );
  tree->Branch("el1_sim_decay2px", &el1_sim_decay2px, "el1_sim_decay2px/F");
  tree->Branch("el1_sim_decay2py", &el1_sim_decay2py, "el1_sim_decay2py/F");
  tree->Branch("el1_sim_decay2pz", &el1_sim_decay2pz, "el1_sim_decay2pz/F");
  tree->Branch("el1_sim_decay2pt", &el1_sim_decay2pt, "el1_sim_decay2pt/F");
  tree->Branch("el1_sim_decay2e",  &el1_sim_decay2e,  "el1_sim_decay2e/F" );
  tree->Branch("el1_sim_decay2pid",  &el1_sim_decay2pid,  "el1_sim_decay2pid/I" );

  //px,py,pz and e of the photon closest to the sc
  tree->Branch("sc_sim_gpx", &sc_sim_gpx, "sc_sim_gpx/F");
  tree->Branch("sc_sim_gpy", &sc_sim_gpy, "sc_sim_gpy/F");
  tree->Branch("sc_sim_gpz", &sc_sim_gpz, "sc_sim_gpz/F");
  tree->Branch("sc_sim_gpt", &sc_sim_gpt, "sc_sim_gpt/F");
  tree->Branch("sc_sim_ge", &sc_sim_ge, "sc_sim_ge/F");
  //position of the vertex where the conversion happened
  tree->Branch("sc_sim_gvx", &sc_sim_gvx, "sc_sim_gvx/F");
  tree->Branch("sc_sim_gvy", &sc_sim_gvy, "sc_sim_gvy/F");
  tree->Branch("sc_sim_gvz", &sc_sim_gvz, "sc_sim_gvz/F");
  tree->Branch("sc_sim_gvr", &sc_sim_gvr, "sc_sim_gvr/F");
  tree->Branch("sc_sim_gveta", &sc_sim_gveta, "sc_sim_gveta/F");
  tree->Branch("sc_sim_gvphi", &sc_sim_gvphi, "sc_sim_gvphi/F");
  //px,py,pz, e of the 2 decay electrons
  tree->Branch("sc_sim_decay1px", &sc_sim_decay1px, "sc_sim_decay1px/F");
  tree->Branch("sc_sim_decay1py", &sc_sim_decay1py, "sc_sim_decay1py/F");
  tree->Branch("sc_sim_decay1pz", &sc_sim_decay1pz, "sc_sim_decay1pz/F");
  tree->Branch("sc_sim_decay1pt", &sc_sim_decay1pt, "sc_sim_decay1pt/F");
  tree->Branch("sc_sim_decay1e",  &sc_sim_decay1e,  "sc_sim_decay1e/F" );
  tree->Branch("sc_sim_decay1pid",  &sc_sim_decay1pid,  "sc_sim_decay1pid/I" );
  tree->Branch("sc_sim_decay2px", &sc_sim_decay2px, "sc_sim_decay2px/F");
  tree->Branch("sc_sim_decay2py", &sc_sim_decay2py, "sc_sim_decay2py/F");
  tree->Branch("sc_sim_decay2pz", &sc_sim_decay2pz, "sc_sim_decay2pz/F");
  tree->Branch("sc_sim_decay2pt", &sc_sim_decay2pt, "sc_sim_decay2pt/F");
  tree->Branch("sc_sim_decay2e",  &sc_sim_decay2e,  "sc_sim_decay2e/F" );
  tree->Branch("sc_sim_decay2pid",  &sc_sim_decay2pid,  "sc_sim_decay2pid/I");  
  
}

void Conversion::endJob() {
  
  file->Write();
  file->Close();
  
}

void Conversion::analyze(const Event & event, const EventSetup& eventSetup) {

  cout << "Run: " << event.id().run() << " Event: " << event.id().event() << endl;

  // access the tracker
  edm::ESHandle<TrackerGeometry> theTrackerGeometry;
  eventSetup.get<TrackerDigiGeometryRecord>().get(theTrackerGeometry);
  const TrackerGeometry& theTracker(*theTrackerGeometry);

  Handle<HepMCProduct> evt;
  event.getByLabel("source", evt);
  myGenEvent = new HepMC::GenEvent(*(evt->GetEvent()));

  Handle<PixelMatchGsfElectronCollection> elh;
  event.getByLabel(baselineEleCollName, elh);
  const PixelMatchGsfElectronCollection* copy = elh.product();
  PixelMatchGsfElectronCollection::const_iterator ite; 

  Handle<SuperClusterCollection> sch1;
  event.getByLabel("hybridSuperClusters", sch1);
  const SuperClusterCollection* scb = sch1.product();

  Handle<SuperClusterCollection> sch2;
  event.getByLabel("islandSuperClusters", "islandEndcapSuperClusters", sch2);
  const SuperClusterCollection* sce = sch2.product();
  SuperClusterCollection::const_iterator itsc, itscb, itsce;
  SuperClusterCollection sc;
  sc.insert(sc.end(), scb->begin(), scb->end());
  sc.insert(sc.end(), sce->begin(), sce->end());

  Handle<PixelMatchGsfElectronCollection> elh1;
  event.getByLabel(customEleCollName, elh1);
  const PixelMatchGsfElectronCollection*  electrons1 = elh1.product();
  PixelMatchGsfElectronCollection::const_iterator ite1;

  Handle<TrackCollection> tkh;
  event.getByLabel("ctfWithMaterialTracks", tkh);
  const TrackCollection* tracks = tkh.product();
  TrackCollection::const_iterator itt;


  
  //__________________________added by PDK______________________
  //Sim Track Collection
  Handle<SimTrackContainer> simTrackCollection;
  event.getByType(simTrackCollection);
  if(!simTrackCollection.isValid()) throw cms::Exception("FatalError") << "SimTrackCollection not found!\n";
  vector<SimTrack> theSimTrks;
  theSimTrks.insert(theSimTrks.end(), 
		    simTrackCollection->begin(),
		    simTrackCollection->end());
  
  Handle<SimVertexContainer> simVertexCollection;
  event.getByType(simVertexCollection);
  if(!simVertexCollection.isValid()) throw cms::Exception("FatalError") << "SimVertexCollection not found!\n";
  vector<SimVertex> theSimVerts;
  theSimVerts.insert(theSimVerts.end(),
		     simVertexCollection->begin(),
		     simVertexCollection->end());

  //____________________________________________________________


  for (itsc = sc.begin(); itsc != sc.end(); ++itsc) {
    
    int mctk_index=0;
    int mctk1_index=0;
    int mcsc_index=0;

    
    math::XYZVector scv(itsc->x(), itsc->y(), itsc->z());
    
    //save SC with Et > 5
    if (sin(scv.Theta())*itsc->energy() > 5.) {
      sc_e = itsc->energy(); 
      sc_rawe = itsc->rawEnergy();
      sc_et = sin(scv.Theta())*itsc->energy();
      sc_eta = itsc->eta();
      sc_phi = itsc->phi(); 
      
      if (fabs(itsc->eta()) < 1.47)
        sc_type = 0; //barrel
      else
        sc_type = 1; //endcap
      
      //match SC to nearest Gen particle by dR
      double dR, dRmin = 0.1;
      HepMC::GenEvent::particle_const_iterator nearMC;
      int i=1;
      for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it,++i) { 
        
        if ((*it)->status() == 1) {      
          math::XYZVector mcv((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz());
          dR = ROOT::Math::VectorUtil::DeltaR(scv, mcv);
          
          if (dR < dRmin) {
            dRmin = dR;
            nearMC = it;
            mcsc_index = i;
          }
        }
      }//gen partcle loop
      
      //save Gen particle info
      if (dRmin < 0.1) {
        mcsc_dr = dRmin;
        mcsc_mother = mother(*nearMC);
        mcsc_pt = (*nearMC)->momentum().perp();
        mcsc_eta = (*nearMC)->momentum().eta();
        mcsc_phi = (*nearMC)->momentum().phi();
        mcsc_e = (*nearMC)->momentum().e();
        mcsc_id = (*nearMC)->pdg_id();
        
        // check if it is in a crack
        if (inCrack(fabs((*nearMC)->momentum().eta())))
          mcsc_crack = 1;
        else
          mcsc_crack = 0;  
      } else {
        mcsc_dr = 0.1;
        mcsc_mother = -1;
        mcsc_pt = -1;
        mcsc_eta = -1;
        mcsc_phi = -1;
        mcsc_e = -1;
        mcsc_id = 0;
        mcsc_crack = -1;
      }
      
      // remove duplicate electrons
      PixelMatchGsfElectronCollection electrons;
      PixelMatchGsfElectronCollection::const_iterator it1, it2;
      
      for(it1=copy->begin(); it1!=copy->end(); ++it1) {
        
        bool isRemoved = false;
        for(it2=copy->begin(); it2!=copy->end(); ++it2) {
          if (it1 == it2)
            continue;
          if (((*it1).superCluster().id() == (*it2).superCluster().id()) &&
              ((*it1).superCluster().index() == (*it2).superCluster().index())) {
            
            float deltaEp1 = fabs((*it1).eSuperClusterOverP() - 1.);
            float deltaEp2 = fabs((*it2).eSuperClusterOverP() - 1.);
            if (deltaEp1 > deltaEp2) {
              isRemoved = true;
              break;
            }
          }
        }
        
        if (!isRemoved)
          electrons.push_back(*it1);
      }
      
      //match standard electron to SC
      dRmin = 0.1;
      PixelMatchGsfElectronCollection::const_iterator nearElectron;	
      for(ite = electrons.begin(); ite != electrons.end(); ++ite) {
        dR = ROOT::Math::VectorUtil::DeltaR(ite->p4(), scv);
        if (dR < dRmin) {
          dRmin = dR;
          nearElectron = ite;
        }
      }
      // strore info about standard Electron if it passes matching requirement
      if (dRmin < 0.1) {
        el_pt = nearElectron->pt(); 
        el_eta = nearElectron->eta(); 
        el_e = nearElectron->energy();
        el_q = nearElectron->charge();
        el_phi = nearElectron->phi(); 
        el_dr = dRmin; 
        el_eopin = nearElectron->eSuperClusterOverP();
        el_eopout = nearElectron->eSeedClusterOverPout();
        el_hoe = nearElectron->hadronicOverEm();
        el_dphiin = nearElectron->deltaPhiSuperClusterTrackAtVtx();
        el_detain = nearElectron->deltaEtaSuperClusterTrackAtVtx();
        el_dphiout = nearElectron->deltaPhiSeedClusterTrackAtCalo();
        el_detaout = nearElectron->deltaEtaSeedClusterTrackAtCalo();
        float pin  = nearElectron->trackMomentumAtVtx().R();
        float pout = nearElectron->trackMomentumOut().R();
        el_pout = pout;
        el_fbrem = (pin-pout)/pin;
        el_class = nearElectron->classification();
        R9_25_gsf(event, &(*nearElectron), el_eseed, el_e3x3, el_e5x5, el_spp, el_see);
        int a, b;
        nHits(nearElectron->gsfTrack(), a, b);
        el_npxhits = a;
        el_nsihits = b;
        
        int index = 1;
        while(1) {
          TrackingRecHitRef hit = nearElectron->gsfTrack()->recHit(nearElectron->gsfTrack()->recHitsSize()-index);
          
          if (hit->isValid()) {
            GlobalPoint hitPosition = theTracker.idToDet(hit->geographicalId())->surface().toGlobal(hit->localPosition());
            GlobalPoint pos(hitPosition.x()-nearElectron->gsfTrack()->vx(), hitPosition.y()-nearElectron->gsfTrack()->vy(),
                            hitPosition.z()-nearElectron->gsfTrack()->vz());
            //std::cout << "Inner: " <<  HitPosition.perp() << "  " << HitPosition.z() << std::endl;
            el_xhit = hitPosition.x();
            el_yhit = hitPosition.y();
            el_zhit = hitPosition.z();
            el_rinnerhit = pos.perp(); // + pow(pos.z(),2));
            subDetector(hit, a, b);
            el_detinnerhit = a;
            el_layerinnerhit = b;
            break;
          }
          index++;
        }

        el_z0 = nearElectron->gsfTrack()->vz();
        el_d0 = nearElectron->gsfTrack()->d0();
        el_d0err = nearElectron->gsfTrack()->d0Error();
        el_tkiso = trackIsolation(nearElectron->trackMomentumAtVtx(), nearElectron->vertex(), tracks);
        el_tkpt = nearElectron->gsfTrack()->pt();
        el_tketa = nearElectron->gsfTrack()->eta(); 
        el_tkphi = nearElectron->gsfTrack()->phi();  

        double dR, dRmin = 0.1;
        HepMC::GenEvent::particle_const_iterator nearMC;
	int i=1;
        for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it, ++i) { 
          
          if ((*it)->status() == 1) {      
            math::XYZVector mcv((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz());
            dR = ROOT::Math::VectorUtil::DeltaR(nearElectron->gsfTrack()->innerMomentum(), mcv);
	    
            if (dR < dRmin) {
              dRmin = dR;
              nearMC = it;
	      mctk_index = i;
            }
          }
        }
        
        if (dRmin < 0.1) {
          mctk_dr = dRmin;
          mctk_mother = mother(*nearMC);
          mctk_pt = (*nearMC)->momentum().perp();
          mctk_eta = (*nearMC)->momentum().eta();
          mctk_phi = (*nearMC)->momentum().phi();
          mctk_e = (*nearMC)->momentum().e();
          mctk_id = (*nearMC)->pdg_id();
          
          // check if it is in a crack
          if (inCrack(fabs((*nearMC)->momentum().eta())))
            mctk_crack = 1;
          else
            mctk_crack = 0;  
        } else {
          mctk_dr = 0.1;
          mctk_mother = -1;
          mctk_pt = -1;
          mctk_eta = -1;
          mctk_phi = -1;
          mctk_e = -1;
          mctk_id = 0;
          mctk_crack = -1;
        }
      } else {
        el_pt = 0.;
        el_e = 0.;
	el_q = 0;
        el_eta = 0.;
        el_phi = 0.;
        el_dr = 0.1;
        el_eopin = 0.;
        el_eopout = 0.;
        el_pout = 0;
        el_hoe = 0.;
        el_dphiin = 0.;
        el_detain = 0.;
        el_dphiout = 0.;
        el_detaout = 0.;
        el_fbrem = 0.;
        el_eseed = 0.;
        el_e3x3 = 0.;
        el_e5x5 = 0.;
        el_spp = 0.;
        el_see = 0.;
        el_class = -1;
        el_npxhits = -1;
        el_nsihits = -1;
        el_rinnerhit = 0.;
        el_detinnerhit = -1;
        el_z0 = -1;
        el_d0 = -1;
        el_d0err = -1;
        el_tkiso = -1;
        el_tkpt = 0.;
        el_tketa = 0.; 
        el_tkphi = 0.; 
        mctk_dr = 0.1;
        mctk_mother = -1;
        mctk_pt = -1;
        mctk_eta = -1;
        mctk_phi = -1;
        mctk_e = -1;
        mctk_id = -1;
        mctk_crack = -1;
      }
      
      // new electrons collection
      //GlobalCtfElectronCollection::const_iterator nearElectron1;
      PixelMatchGsfElectronCollection::const_iterator nearElectron1;
      dRmin = 0.1;
      for(ite1 = electrons1->begin(); ite1 != electrons1->end(); ++ite1) {
        dR = ROOT::Math::VectorUtil::DeltaR(ite1->p4(), scv);
        if (dR < dRmin) {
          dRmin = dR;
          nearElectron1 = ite1;
        }
      }
      
      // strore info about Ele
      if (dRmin < 0.1) {
        el1_pt = nearElectron1->pt(); 
        el1_eta = nearElectron1->eta(); 
        el1_e = nearElectron1->energy();
	el1_q = nearElectron1->charge();
        el1_phi = nearElectron1->phi(); 
        el1_dr = dRmin; 
        el1_eopin = nearElectron1->eSuperClusterOverP();
        el1_eopout = nearElectron1->eSeedClusterOverPout();
        el1_hoe = nearElectron1->hadronicOverEm();
        el1_dphiin = nearElectron1->deltaPhiSuperClusterTrackAtVtx();
        el1_detain = nearElectron1->deltaEtaSuperClusterTrackAtVtx();
        el1_dphiout = nearElectron1->deltaPhiSeedClusterTrackAtCalo();
        el1_detaout = nearElectron1->deltaEtaSeedClusterTrackAtCalo();
        float pin  = nearElectron1->trackMomentumAtVtx().R();
        float pout = nearElectron1->trackMomentumOut().R();
        el1_pout = pout;
        el1_fbrem = (pin-pout)/pin;
        el1_class = nearElectron1->classification();
        R9_25_gsf(event, &(*nearElectron1), el1_eseed, el1_e3x3, el1_e5x5, el1_spp, el1_see);
        int a, b;

        nHits(nearElectron1->gsfTrack(), a, b);
        el1_npxhits = a;
        el1_nsihits = b; 

        int index = 1;
        while(1) {
          TrackingRecHitRef hit = nearElectron1->gsfTrack()->recHit(nearElectron1->gsfTrack()->recHitsSize()-index);
          
          if (hit->isValid()) {
            GlobalPoint hitPosition = theTracker.idToDet(hit->geographicalId())->surface().toGlobal(hit->localPosition());
            GlobalPoint pos(hitPosition.x()-nearElectron1->gsfTrack()->vx(), hitPosition.y()-nearElectron1->gsfTrack()->vy(),
                            hitPosition.z()-nearElectron1->gsfTrack()->vz());
            //std::cout << "Inner: " <<  HitPosition.perp() << "  " << HitPosition.z() << std::endl;
            el1_xhit = hitPosition.x();
            el1_yhit = hitPosition.y();
            el1_zhit = hitPosition.z();
            el1_rinnerhit = pos.perp(); // + pow(pos.z(),2));
            subDetector(hit, a, b);
            el1_detinnerhit = a;
            el1_layerinnerhit = b;
            break;
          }
          index++;
        }

        el1_z0 = nearElectron1->gsfTrack()->vz();
        el1_d0 = nearElectron1->gsfTrack()->d0();
        el1_d0err = nearElectron1->gsfTrack()->d0Error();
        el1_tkiso = trackIsolation(nearElectron1->trackMomentumAtVtx(), nearElectron1->vertex(), tracks);
        el1_tkpt = nearElectron1->gsfTrack()->pt();
        el1_tketa = nearElectron1->gsfTrack()->eta(); 
        el1_tkphi = nearElectron1->gsfTrack()->phi(); 

        double dR, dRmin = 0.1;
        i=1;
        HepMC::GenEvent::particle_const_iterator nearMC;
        for (HepMC::GenEvent::particle_const_iterator it = myGenEvent->particles_begin(); it != myGenEvent->particles_end(); ++it, ++i) { 
          
          if ((*it)->status() == 1) {      
            math::XYZVector mcv((*it)->momentum().px(), (*it)->momentum().py(), (*it)->momentum().pz());
            dR = ROOT::Math::VectorUtil::DeltaR(nearElectron1->gsfTrack()->innerMomentum(), mcv);
            
            if (dR < dRmin) {
              dRmin = dR;
              nearMC = it;
              mctk1_index = i;
            }
          }
        }
        
        if (dRmin < 0.1) {
          mctk1_dr = dRmin;
          mctk1_mother = mother(*nearMC);
          mctk1_pt = (*nearMC)->momentum().perp();
          mctk1_eta = (*nearMC)->momentum().eta();
          mctk1_phi = (*nearMC)->momentum().phi();
          mctk1_e = (*nearMC)->momentum().e();
          mctk1_id = (*nearMC)->pdg_id();
          
          // check if it is in a crack
          if (inCrack(fabs((*nearMC)->momentum().eta())))
            mctk1_crack = 1;
          else
            mctk1_crack = 0;  
        } else {
          mctk1_dr = 0.1;
          mctk1_mother = -1;
          mctk1_pt = -1;
          mctk1_eta = -1;
          mctk1_phi = -1;
          mctk1_e = -1;
          mctk1_id = 0;
          mctk1_crack = -1;
        }
      } else {
        el1_pt = 0.;
        el1_e = 0.;
        el1_q = 0;
        el1_eta = 0.;
        el1_phi = 0.;
        el1_dr = 0.1; 
        el1_eopin = 0.;
        el1_eopout = 0.;
        el1_pout = 0;
        el1_hoe = 0.;
        el1_dphiin = 0.;
        el1_detain = 0.;
        el1_dphiout = 0.;
        el1_detaout = 0.;
        el1_fbrem = 0.;
        el1_eseed = 0.;
        el1_e3x3 = 0.;
        el1_e5x5 = 0.;
        el1_spp = 0.;
        el1_see = 0.;
        el1_class = -1;
        el1_rinnerhit = 0.;
        el1_detinnerhit = -1;
        el1_npxhits = -1;
        el1_nsihits = -1;
        el1_z0 = -1;
        el1_d0 = -1;
        el1_d0err = -1;
        el1_tkiso = -1;
        el1_tkpt = 0.;
        el1_tketa = 0.; 
        el1_tkphi = 0.; 
        mctk1_dr = 0.1;
        mctk1_mother = -1;
        mctk1_pt = -1;
        mctk1_eta = -1;
        mctk1_phi = -1;
        mctk1_e = -1;
        mctk1_id = -1;
        mctk1_crack = -1;
      }
            
      InitSimVariables();
      
      //if there is a standard electron or a UCSD electron, fill SC info
      //now do the same for the MC particle matched closest to the SC
      if(el_dr < 0.1 || el1_dr < 0.1) {
	vector<SimTrack>::const_iterator theSimTrksIter;
	int SimIndex = 0;
	int photonTrackId = 0;
	
	for(theSimTrksIter = theSimTrks.begin();
	    theSimTrksIter!= theSimTrks.end();
	    theSimTrksIter++, SimIndex++) {
	  
	  //is the MC particle closest to the CTF track matched to the SimTrack?
	  if(theSimTrksIter->genpartIndex() == mcsc_index) { 
	
	    int i = 1;
	    HepMC::GenEvent::particle_const_iterator geniter;
	    for(geniter = myGenEvent->particles_begin(); i < mcsc_index; ++i, ++geniter) {}
	    
	    //Making sure that the Sim and Gen particle are both Photons
	    if(theSimTrksIter->type() == 22 && (*geniter)->pdg_id()==22) {
	      photonTrackId = theSimTrksIter->trackId();
	      
	      HepLorentzVector psimvect = theSimTrksIter->momentum();
	      sc_sim_gpx = psimvect.px();
	      sc_sim_gpy = psimvect.py();
	      sc_sim_gpz = psimvect.pz();
	      sc_sim_gpt = TMath::Sqrt(sc_sim_gpx*sc_sim_gpx+sc_sim_gpy*sc_sim_gpy);
	      sc_sim_ge = psimvect.e();
	    	      
	      
	      //SimVertex indicies start at 0
	      i=0;
	      vector<SimVertex>::const_iterator theSimVertsIter;
	      for(theSimVertsIter = theSimVerts.begin();
		  theSimVertsIter != theSimVerts.end();
		  theSimVertsIter++, i++) {
		
		//check to see if the parent of this vertex was the photon in question
		if(theSimVertsIter->parentIndex() == photonTrackId) {
		  		  
		  sc_sim_gvx = theSimVertsIter->position().x();
		  sc_sim_gvy = theSimVertsIter->position().y();
		  sc_sim_gvz = theSimVertsIter->position().z();
		  sc_sim_gvr = TMath::Sqrt(sc_sim_gvx*sc_sim_gvx+sc_sim_gvy*sc_sim_gvy);
		  sc_sim_gvphi = atan2(sc_sim_gvy, sc_sim_gvx);
		  //make sure that we don't get Log of 0
		  if(fabs(sc_sim_gvz) > 1E-7 && sc_sim_gvr > 1E-7) {
		    sc_sim_gveta = TMath::Log( TMath::Tan( 0.5*TMath::ATan(sc_sim_gvr/fabs(sc_sim_gvz) ) ) );
		    if(sc_sim_gvz < 0) 
		      sc_sim_gveta = -fabs(sc_sim_gvz);
		  }else sc_sim_gveta = -999.;
		  
		  
		  /*loop over the simtracks (ugh...again) to see if there are
		    any tracks which have this vertex
		  */
		  for(unsigned int k = 0; k < theSimTrks.size(); k++) {
		    if(theSimTrks.at(k).vertIndex() == i) {
		      
		      //if(sc_sim_decay1pid==0) {
		      float charge = theSimTrks.at(k).charge();
		      if(el_q/charge > 0) {
			sc_sim_decay1px = theSimTrks.at(k).momentum().x();
			sc_sim_decay1py = theSimTrks.at(k).momentum().y();
			sc_sim_decay1pz = theSimTrks.at(k).momentum().z();
			sc_sim_decay1pt = TMath::Sqrt(sc_sim_decay1px*sc_sim_decay1px+
						      sc_sim_decay1py*sc_sim_decay1py);
			sc_sim_decay1e = theSimTrks.at(k).momentum().t();
			sc_sim_decay1pid = theSimTrks.at(k).type();
		      } else {
			sc_sim_decay2px = theSimTrks.at(k).momentum().x();
			sc_sim_decay2py = theSimTrks.at(k).momentum().y();
			sc_sim_decay2pz = theSimTrks.at(k).momentum().z();
			sc_sim_decay2pt = TMath::Sqrt(sc_sim_decay1px*sc_sim_decay1px+
						      sc_sim_decay1py*sc_sim_decay1py);
			sc_sim_decay2e = theSimTrks.at(k).momentum().t();
			sc_sim_decay2pid = theSimTrks.at(k).type();
		      }
		      
		    }
		  }//sim track loop 
		}//is the parent of the vertex a photon?
	      }//SimVerteces Iterator
	    }//if(theSimTrksIter->type() == 22 && (*geniter)->pdg_id()==22)  
	  }//if(theSimTrksIter->genpartIndex() == mcsc_index) {
	}//SimTracksIterator
      }//if(el_dr < 0.1 || el1_dr < 0.1) 


      
      //did we find a standard electron?
      if(el_dr < 0.1) {//existance of electron
	
	vector<SimTrack>::const_iterator theSimTrksIter;
	int SimIndex = 0;
	int photonTrackId = 0;
	for(theSimTrksIter = theSimTrks.begin();
	    theSimTrksIter!= theSimTrks.end();
	    theSimTrksIter++, SimIndex++) {
	  
	  //is the MC particle closest to the CTF track matched to the SimTrack?
	  if(theSimTrksIter->genpartIndex() == mctk1_index) { 
	    
	    int i = 1;
	    HepMC::GenEvent::particle_const_iterator geniter;
	    for(geniter = myGenEvent->particles_begin(); i < mctk1_index; ++i, ++geniter) {}
	    
	    //Making sure that the Sim and Gen particle are both Photons
	    if(theSimTrksIter->type() == 22 && (*geniter)->pdg_id()==22) {
	      photonTrackId = theSimTrksIter->trackId();
	      
	      HepLorentzVector psimvect = theSimTrksIter->momentum();
	      el_sim_gpx = psimvect.px();
	      el_sim_gpy = psimvect.py();
	      el_sim_gpz = psimvect.pz();
	      el_sim_gpt = TMath::Sqrt(el_sim_gpx*el_sim_gpx+el_sim_gpy*el_sim_gpy);
	      el_sim_ge = psimvect.e();
	      
	      
	      //SimVertex indicies start at 0
	      i=0;
	      vector<SimVertex>::const_iterator theSimVertsIter;
	      for(theSimVertsIter = theSimVerts.begin();
		  theSimVertsIter != theSimVerts.end();
		  theSimVertsIter++, i++) {
		
		//check to see if the parent of this vertex was the photon in question
		if(theSimVertsIter->parentIndex() == photonTrackId) {
		  		  
		  el_sim_gvx = theSimVertsIter->position().x();
		  el_sim_gvy = theSimVertsIter->position().y();
		  el_sim_gvz = theSimVertsIter->position().z();
		  el_sim_gvr = TMath::Sqrt(el_sim_gvx*el_sim_gvx+el_sim_gvy*el_sim_gvy);
		  el_sim_gvphi = atan2(el_sim_gvy, el_sim_gvx);
		  //make sure that we don't get Log of 0
		  if(fabs(el_sim_gvz) > 1E-7 && el_sim_gvr > 1E-7) {
		    el_sim_gveta = TMath::Log( TMath::Tan( 0.5*TMath::ATan(el_sim_gvr/fabs(el_sim_gvz) ) ) );
		    if(el_sim_gvz < 0) 
		      el_sim_gveta  = -fabs(el_sim_gveta);
		  }else el_sim_gveta = -999.;

		  
		  /*loop over the simtracks (ugh...again) to see if there are
		    any tracks which have this vertex
		  */
		  for(unsigned int k = 0; k < theSimTrks.size(); k++) {
		    if(theSimTrks.at(k).vertIndex() == i) {
		      
		      float charge = theSimTrks.at(k).charge();
		      if(el1_q/charge > 0) {
			el_sim_decay1px = theSimTrks.at(k).momentum().x();
			el_sim_decay1py = theSimTrks.at(k).momentum().y();
			el_sim_decay1pz = theSimTrks.at(k).momentum().z();
			el_sim_decay1pt = TMath::Sqrt(el_sim_decay1px*el_sim_decay1px+
						      el_sim_decay1py*el_sim_decay1py);
			el_sim_decay1e = theSimTrks.at(k).momentum().t();
			el_sim_decay1pid = theSimTrks.at(k).type();
		      } else {
			el_sim_decay2px = theSimTrks.at(k).momentum().x();
			el_sim_decay2py = theSimTrks.at(k).momentum().y();
			el_sim_decay2pz = theSimTrks.at(k).momentum().z();
			el_sim_decay2pt = TMath::Sqrt(el_sim_decay2px*el_sim_decay2px+
						      el_sim_decay2py*el_sim_decay2py);
			el_sim_decay2e = theSimTrks.at(k).momentum().t();
			el_sim_decay2pid = theSimTrks.at(k).type();
		      }
		      
		    }
		  }//sim track loop 
		}//is the parent of the vertex a photon?
	      }//SimVerteces Iterator
	    }//if(theSimTrksIter->type() == 22 && (*geniter)->pdg_id()==22)  
	  }//if(theSimTrksIter->genpartIndex() == mctk1_index) {
	}//SimTracksIterator
      }

      //is there a UCSD electron?
      if(el1_dr < 0.1) {
	
	vector<SimTrack>::const_iterator theSimTrksIter;
	int SimIndex = 0;
	int photonTrackId = 0;
	for(theSimTrksIter = theSimTrks.begin();
	    theSimTrksIter!= theSimTrks.end();
	    theSimTrksIter++, SimIndex++) {
	  
	  //is the MC particle closest to the CTF track matched to the SimTrack?
	  if(theSimTrksIter->genpartIndex() == mctk1_index) { 
	    
	    int i = 1;
	    HepMC::GenEvent::particle_const_iterator geniter;
	    for(geniter = myGenEvent->particles_begin(); i < mctk1_index; ++i, ++geniter) {}
	    
	    //Making sure that the Sim and Gen particle are both Photons
	    if(theSimTrksIter->type() == 22 && (*geniter)->pdg_id()==22) {
	      photonTrackId = theSimTrksIter->trackId();
	      
	      HepLorentzVector psimvect = theSimTrksIter->momentum();
	      el1_sim_gpx = psimvect.px();
	      el1_sim_gpy = psimvect.py();
	      el1_sim_gpz = psimvect.pz();
	      el1_sim_gpt = TMath::Sqrt(el1_sim_gpx*el1_sim_gpx+el1_sim_gpy*el1_sim_gpy);
	      el1_sim_ge = psimvect.e();
	      
	      //SimVertex indicies start at 0
	      i=0;
	      vector<SimVertex>::const_iterator theSimVertsIter;
	      for(theSimVertsIter = theSimVerts.begin();
		  theSimVertsIter != theSimVerts.end();
		  theSimVertsIter++, i++) {
		
		//check to see if the parent of this vertex was the photon in question
		if(theSimVertsIter->parentIndex() == photonTrackId) {
		  el1_sim_gvx = theSimVertsIter->position().x();
		  el1_sim_gvy = theSimVertsIter->position().y();
		  el1_sim_gvz = theSimVertsIter->position().z();
		  el1_sim_gvr = TMath::Sqrt(el1_sim_gvx*el1_sim_gvx+el1_sim_gvy*el1_sim_gvy);
		  el1_sim_gvphi = atan2(el1_sim_gvy, el1_sim_gvx);
		  //make sure that we don't get Log of 0
		  if(fabs(el1_sim_gvz) > 1E-7 && el1_sim_gvr > 1E-7) {
		    el1_sim_gveta = TMath::Log( TMath::Tan( 0.5*TMath::ATan(el1_sim_gvr/fabs(el1_sim_gvz) ) ) );
		    if(el1_sim_gvz < 0) 
		      el1_sim_gveta = -fabs(el1_sim_gveta);
		  }else el1_sim_gveta = -999.;
		  
		  /*loop over the simtracks (ugh...again) to see if there are
		    any tracks which have this vertex
		  */
		  for(unsigned int k = 0; k < theSimTrks.size(); k++) {
		    if(theSimTrks.at(k).vertIndex() == i) {
		      float charge = theSimTrks.at(k).charge();
		      if(el1_q/charge > 0) {
			el1_sim_decay1px = theSimTrks.at(k).momentum().x();
			el1_sim_decay1py = theSimTrks.at(k).momentum().y();
			el1_sim_decay1pz = theSimTrks.at(k).momentum().z();
			el1_sim_decay1pt = TMath::Sqrt(el1_sim_decay1px*el1_sim_decay1px+
						       el1_sim_decay1py*el1_sim_decay1py);
			el1_sim_decay1e = theSimTrks.at(k).momentum().t();
			el1_sim_decay1pid = theSimTrks.at(k).type();
		      } else {
			el1_sim_decay2px = theSimTrks.at(k).momentum().x();
			el1_sim_decay2py = theSimTrks.at(k).momentum().y();
			el1_sim_decay2pz = theSimTrks.at(k).momentum().z();
			el1_sim_decay2pt = TMath::Sqrt(el1_sim_decay2px*el1_sim_decay2px+
						       el1_sim_decay2py*el1_sim_decay2py);
			el1_sim_decay2e = theSimTrks.at(k).momentum().t();
			el1_sim_decay2pid = theSimTrks.at(k).type();
		      }
		      
		    }
		  }//sim track loop 
		}//is the parent of the vertex a photon?
	      }//SimVerteces Iterator
	    }//if(theSimTrksIter->type() == 22 && (*geniter)->pdg_id()==22)  
	  }//if(theSimTrksIter->genpartIndex() == mctk1_index) {
	}//SimTracksIterator
      
      }//if(el1_dr<0.1)
      
      
      run = event.id().run();
      id = event.id().event();
      tree->Fill();
    }
  }
}

bool Conversion::inCrack(float eta) {

  return (eta < 0.018 ||
          (eta>0.423 && eta<0.461) ||
          (eta>0.770 && eta<0.806) ||
          (eta>1.127 && eta<1.163) ||
          (eta>1.460 && eta<1.558));
}

int Conversion::mother(HepMC::GenParticle *p) {
  
  while (p->production_vertex()) {
    HepMC::GenVertex* inVertex = p->production_vertex();
    for(std::set<HepMC::GenParticle*>::const_iterator iter = inVertex->particles_in_const_begin();
        iter != inVertex->particles_in_const_end();iter++) {
      if ((*iter)->pdg_id() != p->pdg_id()) {
        return (*iter)->pdg_id();
      } else {
        p = *iter;
        break;
      }
    }
  }
  
  return -1;
}

void Conversion::R9_25_gsf(const Event & event, const reco::PixelMatchGsfElectron* e,
			   float& eseed, float& e3x3, float& e5x5, float& spp, float& see) {
  
  reco::SuperClusterRef sclRef=e->superCluster();

  edm::Handle<reco::BasicClusterShapeAssociationCollection> bH, eH;
  event.getByLabel("hybridSuperClusters", "hybridShapeAssoc", bH);
  const reco::BasicClusterShapeAssociationCollection* barrelClShp = bH.product();
  event.getByLabel("islandBasicClusters", "islandEndcapShapeAssoc", eH);
  const reco::BasicClusterShapeAssociationCollection* endcapClShp = eH.product();

  reco::BasicClusterShapeAssociationCollection::const_iterator seedShpItr;
  DetId id = sclRef->seed()->getHitsByDetId()[0];
  if (id.subdetId() == EcalBarrel) {
    seedShpItr = barrelClShp->find(sclRef->seed());
  } else {
    seedShpItr = endcapClShp->find(sclRef->seed());
  }

  // Get the ClusterShapeRef corresponding to the BasicCluster
  const reco::ClusterShapeRef& seedShapeRef = seedShpItr->val;

  eseed = sclRef->seed()->energy();
  e3x3 = seedShapeRef->e3x3();
  e5x5 = seedShapeRef->e5x5();
  spp = sqrt(seedShapeRef->covPhiPhi());
  see = sqrt(seedShapeRef->covEtaEta());
}

void Conversion::nHits(const reco::GsfTrackRef t, int& nPixelHits, int& nSiTkHits) {

  // loop sugli hits e conta il risultato facile no ?
  nPixelHits = 0; 
  nSiTkHits = 0;

  for(size_t i = 0; i < t->recHitsSize(); ++i) {

    TrackingRecHitRef hit = t->recHit(i);
    
    if (hit->isValid()) {
      DetId detid(hit->geographicalId());       
      unsigned int subdetId = static_cast<unsigned int>(detid.subdetId());
      if ((subdetId > 2) && (subdetId < 7))
        nSiTkHits++;
      if ((subdetId == 2) || (subdetId == 1))
        nPixelHits++;

    }
  }
}

//trackRelIsolation(el->trackMomentumAtVtx(), el->vertex(), tracks, 0.3, 0.01, 0.1, 999.9, 0.5, 1.5, 7);
double Conversion::trackIsolation(const math::XYZVector momentum, 
				  const math::XYZPoint vertex,
				  const TrackCollection* tracks) {

  double dRConeMax = 0.3;
  double dRConeMin = 0.01;
  double tkVtxDMax = 0.1;
  double vtxDiffDMax = 999.9;
  double vtxDiffZMax = 0.5;
  double ptMin = 1.5;
  unsigned int nHits = 7;
  double isoResult = -10.;

  if ( tracks == 0 ) {
    return isoResult;
  }
  
  double sumPt = 0;
  
  std::vector<Track>::const_iterator iTk;
  for (iTk = tracks->begin(); iTk != tracks->end(); ++iTk){
    double dR = ROOT::Math::VectorUtil::DeltaR(momentum, iTk->momentum());
    //exclude tks in veto cone (set it to small number to 
    //exclude this track
    double dZ = fabs(vertex.z() - iTk->vz());
    double d0 = sqrt(iTk->vertex().perp2());
    double dD0 = sqrt((iTk->vertex() - vertex).perp2());
    
    if (dR < dRConeMin) 
      continue;
    
    if ( dR < dRConeMax 
         && dZ < vtxDiffZMax
         && d0 < tkVtxDMax
         && dD0 < vtxDiffDMax 
         && iTk->pt() >= ptMin
         && iTk->found() > nHits){
      sumPt += iTk->pt();
    }
  }
  
  isoResult = sumPt;

  return isoResult;
}

void Conversion::subDetector(TrackingRecHitRef hit, int& subdet, int& layer) {

  DetId detid(hit->geographicalId());       
  unsigned int subdetId = static_cast<unsigned int>(detid.subdetId());
  switch (subdetId) {
  case 1:
    {
      PXBDetId thePXBDetId(detid.rawId());
      layer = thePXBDetId.layer();
      subdet = 1;
      break;
    }
  case 2:
    {
      PXFDetId thePXFDetId(detid.rawId());
      layer = thePXFDetId.disk();
      subdet = 2;
      break;
    }
  case StripSubdetector::TIB:
    {
      TIBDetId theTIBDetId(detid.rawId());
      layer = theTIBDetId.layer();
      subdet = 3;
      break;
    }
  case StripSubdetector::TID:
    {
      TIDDetId theTIDDetId(detid.rawId());
      layer = theTIDDetId.wheel();
      subdet = 4;
      break;
    }
  case StripSubdetector::TOB:
    {
      TOBDetId theTOBDetId(detid.rawId());
      layer = theTOBDetId.layer();
      subdet = 5;
      break;
    }
  case StripSubdetector::TEC:
    {
      TECDetId theTECDetId(detid.rawId());             
      layer = theTECDetId.wheel();
      subdet = 6;
      break;
    }
  }

}



void Conversion::InitSimVariables() {

  sc_sim_gpx = 0.0;
  sc_sim_gpy = 0.0;
  sc_sim_gpz = 0.0;
  sc_sim_gpt = 0.0;
  sc_sim_ge = 0.0; 
  sc_sim_gvx = 0.0; 
  sc_sim_gvy = 0.0;
  sc_sim_gvz = 0.0;
  sc_sim_gvr = 0.0;
  sc_sim_gvphi = 0.0;
  sc_sim_gveta = 0.0;
  sc_sim_decay1px = 0.0;
  sc_sim_decay1py = 0.0;
  sc_sim_decay1pz = 0.0;
  sc_sim_decay1pt = 0.0;  
  sc_sim_decay1e = 0.0;
  sc_sim_decay1pid = 0;
  sc_sim_decay2px = 0.0; 
  sc_sim_decay2py = 0.0;
  sc_sim_decay2pz = 0.0;
  sc_sim_decay2pt = 0.0;
  sc_sim_decay2e = 0.0;
  sc_sim_decay2pid = 0;



  el_sim_gpx = 0.0;
  el_sim_gpy = 0.0;
  el_sim_gpz = 0.0;
  el_sim_gpt = 0.0;
  el_sim_ge = 0.0; 
  el_sim_gvx = 0.0; 
  el_sim_gvy = 0.0;
  el_sim_gvz = 0.0;
  el_sim_gvr = 0.0;
  el_sim_gvphi = 0.0;
  el_sim_gveta = 0.0;
  el_sim_decay1px = 0.0;
  el_sim_decay1py = 0.0;
  el_sim_decay1pz = 0.0;
  el_sim_decay1pt = 0.0;  
  el_sim_decay1e = 0.0;
  el_sim_decay1pid = 0;
  el_sim_decay2px = 0.0; 
  el_sim_decay2py = 0.0;
  el_sim_decay2pz = 0.0;
  el_sim_decay2pt = 0.0;
  el_sim_decay2e = 0.0;
  el_sim_decay2pid = 0;

  el1_sim_gpx = 0.0;
  el1_sim_gpy = 0.0;
  el1_sim_gpz = 0.0;
  el1_sim_gpt = 0.0;
  el1_sim_ge = 0.0; 
  el1_sim_gvx = 0.0; 
  el1_sim_gvy = 0.0;
  el1_sim_gvz = 0.0;
  el1_sim_gvr = 0.0;
  el1_sim_gvphi = 0.0;
  el1_sim_gveta = 0.0;
  el1_sim_decay1px = 0.0;
  el1_sim_decay1py = 0.0;
  el1_sim_decay1pz = 0.0;
  el1_sim_decay1pt = 0.0;  
  el1_sim_decay1e = 0.0;
  el1_sim_decay1pid = 0;
  el1_sim_decay2px = 0.0; 
  el1_sim_decay2py = 0.0;
  el1_sim_decay2pz = 0.0;
  el1_sim_decay2pt = 0.0;
  el1_sim_decay2e = 0.0;
  el1_sim_decay2pid = 0;  
  
  
}
