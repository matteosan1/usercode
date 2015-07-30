/**
 * MuScleFitCorrector class
 * Author M. De Mattia - 18/11/2008
 * Author S. Casasso   - 25/10/2012
 * Author E. Migliore  - 25/10/2012
 */

#ifndef MuScleFitCorrector_h
#define MuScleFitCorrector_h

#include <iostream>
#include <fstream>
#include <sstream>
#include "HtollAnalysis/interface/ScaleFunct.h"
#include "HtollAnalysis/interface/ResolFunct.h"
#include "TLorentzVector.h"
#include "TRandom3.h"


/**
 * This is used to have a common set of functions for the specialized templates to use.
 * The constructor receives the name identifying the parameters for the correction function.
 * It reads the parameters from a txt file in data/.
 *
 * It handles multiple iterations. It is also possible to use different functions in
 * different iterations.
 *
 * ATTENTION: it is important that iterations numbers in the txt file start from 0.
 */
class MuScleFitCorrector
{
 public:
  /**
   * The constructor takes a string identifying the parameters to read. It
   * parses the txt file containing the parameters, extracts the index of the
   * correction function and saves the corresponding pointer. It then fills the
   * vector of parameters.
   */
  MuScleFitCorrector(const TString& identifier)
  {
    readParameters( identifier ); 
    gRandom_ = new TRandom3();
  }

  ~MuScleFitCorrector() {}

  /// Returns a pointer to the selected function
  scaleFunctBase<double * > * getScaleFunction() { return scaleFunct_; }
  resolFunctBase<double * > * getResolFunction() { return resolFunct_; }

  /// Method to get correct pt
  double getCorrectPt( const TLorentzVector & lorentzVector , const int & chg) {
    
    // Loop on all the functions and apply them iteratively on the pt corrected by the previous function.
    double pt = lorentzVector.Pt();
    pt = ( scaleFunct_->scale( pt, lorentzVector.Eta(), lorentzVector.Phi(), chg, scaleParArray_) );
    return pt;
  }

  // Return the squared difference of relative pT resolutions data-MC
  double getSigmaPtDiffSquared(const double & pt, const double & eta){
    
    double sigmaPtData = resolFunct_->sigmaPt(pt, eta, resolDataParArray_);
    double sigmaPtMC = resolFunct_->sigmaPt(pt, eta, resolMCParArray_);
    return (sigmaPtData*sigmaPtData - sigmaPtMC*sigmaPtMC);
    
  }


  // Method to get correct pt
  double getSmearedPt( const TLorentzVector & lorentzVector , const int & chg) {
    
    // Loop on all the functions and apply them iteratively on the pt corrected by the previous function.
    double pt = lorentzVector.Pt();
    double eta = lorentzVector.Eta();
    double squaredDiff = getSigmaPtDiffSquared(pt,eta);
    double Cfact = 0.8;
    double sPar = Cfact*sqrt(squaredDiff);
    double curv = ((double)chg/pt);
    double smearedCurv = curv + fabs(curv)*(gRandom_->Gaus(0,sPar));

    return 1./((double)chg*smearedCurv);

  }
  
  // Method to apply correction to a TLorentzVector
  void applyPtCorrection( TLorentzVector& lorentzVector, const int & chg ){

    double eta = lorentzVector.Eta();
    double phi = lorentzVector.Phi();
    double m  = lorentzVector.M(); 
    double corrPt = getCorrectPt(lorentzVector, chg);
    lorentzVector.SetPtEtaPhiM(corrPt,eta,phi,m);

  }

  // Method to apply smearing to a TLorentzVector
  void applyPtSmearing(TLorentzVector& lorentzVector, const int & chg){
    double eta = lorentzVector.Eta();
    double phi = lorentzVector.Phi();
    double m  = lorentzVector.M(); 
    double smearedPt = getSmearedPt(lorentzVector, chg);
    lorentzVector.SetPtEtaPhiM(smearedPt,eta,phi,m);

  }

 protected:
  /// Convert vectors to arrays for faster random access. The first pointer is replaced, thus it is taken by reference.
  void convertToArrays();

  // Identification numbers
  int scaleFunctId_;
  int resolFunctId_;

  // Vectors of parameters
  std::vector<double> scaleParVec_;
  std::vector<double> resolDataParVec_;
  std::vector<double> resolMCParVec_;

  // We will use the array for the function calls because it is faster than the vector for random access.
  double * scaleParArray_;
  double * resolDataParArray_;
  double * resolMCParArray_;

  /// Parser of the parameters file
  void readParameters(const TString& fileName);

  // Functions
  scaleFunctBase<double * > * scaleFunct_;
  resolFunctBase<double * > * resolFunct_;

  // Pointer for TRandom3 access
  TRandom3 * gRandom_;

  // Bool for using resolution function or not (value depends from the information on the parameters txt file)
  bool useResol_;
};

#endif // MuScleFitCorrector_h
