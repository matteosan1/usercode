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
	    Float_t weight, Float_t pu_weight, bool isSyst, std::string name1, bool* passing_jets=0);
  
  bool DijetPreSelection(LoopAll& l, TLorentzVector* veto_p41, TLorentzVector* veto_p42, 
      float & dijet_deta, float & dijet_mjj, float & dijet_zep, float & dijet_dphi_ll_jj, 
      float & dijet_j1pt, float & dijet_j2pt, bool* passing_jets=0);
  
  bool doMuon;
  TMVA::Reader *tmvaReader_vbfmumu;
  Float_t tmva_vbfmumu_mjj;
  Float_t tmva_vbfmumu_zep;
  Float_t tmva_vbfmumu_deta;
  Float_t tmva_vbfmumu_dphi_ll_jj;
  Float_t tmva_vbfmumu_j1pt;
  Float_t tmva_vbfmumu_j2pt;
  Float_t tmva_vbfmumu_has2jets;
  Float_t tmva_vbfmumu_itype;
  
  
  int nCategories_;
  float massMin,massMax;
  int nDataBins;
  std::vector<int> bkgPolOrderByCat;
  std::map<int,std::string> signalLabels;
  void buildBkgModel(LoopAll& l, const std::string & postfix); 
  void FillRooContainer(LoopAll& l, int cur_type, float mass, int category, float weight);
  std::string GetSignalLabel(int id);
  void FillSignalLabelMap(LoopAll & l);
  std::vector<int> sigPointsToBook;
  
  bool doBlinding;
 protected:
  std::string name_;
};

#endif
