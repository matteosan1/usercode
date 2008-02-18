#ifndef HLTElectronRegionalSeedFilter_h
#define HLTElectronRegionalSeedFilter_h

/** \class HLTElectronRegionalSeedFilter
 *
 *  \author Matteo Sani (UCSD)
 *
 */

#include "HLTrigger/HLTcore/interface/HLTFilter.h"

//
// class decleration
//

class HLTElectronRegionalSeedFilter : public HLTFilter {
  
   public:
      explicit HLTElectronRegionalSeedFilter(const edm::ParameterSet&);
      ~HLTElectronRegionalSeedFilter();
      virtual bool filter(edm::Event&, const edm::EventSetup&);

   private:
      edm::InputTag candTag_;     // input tag identifying product contains filtered egammas
      edm::InputTag L1IsoPixelSeedsTag_; // input tag for the pixel seed - supercluster map
      //edm::InputTag L1IsoPixelmapendcapTag_; // input tag for the pixel seed - supercluster map
      
      edm::InputTag L1NonIsoPixelSeedsTag_; // input tag for the pixel seed - supercluster map
      //edm::InputTag L1NonIsoPixelmapendcapTag_; // input tag for the pixel seed - supercluster map

      double npixelmatchcut_;     // number of pixelmatch hits
      int    ncandcut_;           // number of electrons required
      
      bool doIsolated_;

};

#endif 


