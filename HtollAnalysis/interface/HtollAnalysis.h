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
  
 protected:
  std::string name_;
};

#endif
