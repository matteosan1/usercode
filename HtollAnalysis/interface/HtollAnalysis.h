#ifndef __HtollAnalysis__
#define __HtollAnalysis__

#include "PhotonAnalysis/interface/StatAnalysis.h"

// ------------------------------------------------------------------------------------
class HtollAnalysis : public StatAnalysis {
 public:
  HtollAnalysis();
  virtual ~HtollAnalysis();
  
  virtual const std::string & name() const { return name_; };
  
  // LoopAll analysis interface implementation
  virtual void Init(LoopAll&);
  virtual void Term(LoopAll&);
  
  virtual void ReducedOutputTree(LoopAll &l, TTree *);
  virtual void GetBranches(TTree *, std::set<TBranch *>& );

  virtual void ResetAnalysis();
  
  virtual void FillReductionVariables(LoopAll& l, int jentry);   
  virtual bool SelectEventsReduction(LoopAll&, int);
  
  virtual bool SkimEvents(LoopAll&, int);
  virtual bool SelectEvents(LoopAll&, int);
  virtual bool Analysis(LoopAll&, Int_t);

  bool ElectronId(LoopAll&, Int_t);
  void Tree(LoopAll& l, Int_t lept1, Int_t lept2, const TLorentzVector & Higgs,
	    Int_t cat, Int_t vbfcat, 
	    Float_t weight, Float_t pu_weight, bool isSyst, std::string name1);
  
  bool DijetPreSelection(LoopAll& l, TLorentzVector* veto_p41, TLorentzVector* veto_p42, 
      float & dijet_deta, float & dijet_mjj, float & dijet_zep, float & dijet_dphi_ll_jj, 
      float & dijet_j1pt, float & dijet_j2pt);
  
  bool doMuon;
  
 protected:
  std::string name_;
};

#endif
