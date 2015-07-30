#ifndef RecoEcal_EgammaClusterProducers_EgammaHLTHybridClusterProducer_h_
#define RecoEcal_EgammaClusterProducers_EgammaHLTHybridClusterProducer_h_

#include <memory>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoEcal/EgammaClusterAlgos/interface/HybridClusterAlgo.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"

//


class EgammaHLTHybridClusterProducer : public edm::EDProducer 
{
  
  public:

      EgammaHLTHybridClusterProducer(const edm::ParameterSet& ps);

      ~EgammaHLTHybridClusterProducer();

      virtual void produce(edm::Event&, const edm::EventSetup&);

   private:

      int nMaxPrintout_; // max # of printouts
      int nEvt_;         // internal counter of events

      bool doIsolated_;

      std::string basicclusterCollection_;
      std::string superclusterCollection_;
      std::string hitproducer_;
      std::string hitcollection_;

      edm::InputTag l1TagIsolated_;
      edm::InputTag l1TagNonIsolated_;
      //edm::InputTag l1Tag_;
      double l1LowerThr_;
      double l1UpperThr_;
      double l1LowerThrIgnoreIsolation_;

      double regionEtaMargin_;
      double regionPhiMargin_;

      HybridClusterAlgo::DebugLevel debugL;

      HybridClusterAlgo * hybrid_p; // clustering algorithm
      PositionCalc posCalculator_; // position calculation algorithm

      bool counterExceeded() const { return ((nEvt_ > nMaxPrintout_) || (nMaxPrintout_ < 0));}
};


#endif


