#include "HtollAnalysis/interface/HtollAnalysis.h"

#include <iostream>

#define HtollAnalysisDEBUG 0

using namespace std;

// ----------------------------------------------------------------------------------------------------
HtollAnalysis::HtollAnalysis()  : 
  name_("HtollAnalysis")
{}

// ----------------------------------------------------------------------------------------------------
HtollAnalysis::~HtollAnalysis() 
{}

// ----------------------------------------------------------------------------------------------------
void HtollAnalysis::Term(LoopAll& l) 
{}

// ----------------------------------------------------------------------------------------------------
void HtollAnalysis::Init(LoopAll& l) {
  /* -------------------------------------------------------------------------------------------
     Pileup Reweighting
     https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupReweighting
     ----------------------------------------------------------------------------------------------  */
  if (puHist != "" && puHist != "auto" ) {
    if(DEBUG) 
      cout << "Opening PU file"<<endl;
    TFile* puFile = TFile::Open( puHist );
    if (puFile) {
      TH1 * target = 0;
      
      if( puTarget != "" ) {
        TFile * puTargetFile = TFile::Open( puTarget ); 
        assert( puTargetFile != 0 );
        target = (TH1*)puTargetFile->Get("pileup");
        if( target == 0 ) { target = (TH1*)puTargetFile->Get("target_pu"); }
        target->Scale( 1. / target->Integral() );
      }
      
      if( puMap != "" ) {
        loadPuMap(puMap, puFile, target); 
      } else {
        loadPuWeights(0, puFile, target);
      }
      puFile->Close();
    }
    else {
      cout<<"Error opening " <<puHist<<" pileup reweighting histogram, using 1.0"<<endl; 
      weights[0].resize(50);
      for (unsigned int i=0; i<weights[0].size(); i++) weights[0][i] = 1.0;
    }
    if(DEBUG) 
      cout << "Opening PU file END"<<endl;
  } else if ( puHist == "auto" ) {
    TFile * puTargetFile = TFile::Open( puTarget ); 
    assert( puTargetFile != 0 );
    puTargetHist = (TH1*)puTargetFile->Get("pileup");
    if( puTargetHist == 0 ) { 
      puTargetHist = (TH1*)puTargetFile->Get("target_pu"); 
    }
    puTargetHist = (TH1*)puTargetHist->Clone();
    puTargetHist->SetDirectory(0);
    puTargetHist->Scale( 1. / puTargetHist->Integral() );
    puTargetFile->Close();
  }


}

// ----------------------------------------------------------------------------------------------------
bool HtollAnalysis::Analysis(LoopAll& l, Int_t jentry) {

  //apply pileup reweighting
  float weight = 1.;
  if (l.itype[l.current] != 0) {
    unsigned int n_pu = l.pu_n;
    weight = getPuWeight(l.pu_n, l.itype[l.current], &(l.sampleContainer[l.current_sample_index]), jentry == 1) * l.sampleContainer[l.current_sample_index].weight;
    //std::cout << l.sampleContainer[l.current_sample_index].weight << " " << weight/l.sampleContainer[l.current_sample_index].weight << " " << n_pu << std::endl;
  }
  
  std::vector<int> goodMuons;

  for (int i=0; i<l.mu_n; i++) {
    TLorentzVector* p4 = (TLorentzVector*)l.mu_p4->At(i);
    if (l.MuonTightID2012(i, 0))
      if (l.MuonIsolation2012(i, p4->Pt(), true))
	goodMuons.push_back(i);
  }

  if (goodMuons.size() < 2)
    return false;

  
}

// ----------------------------------------------------------------------------------------------------
void HtollAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{}

// ----------------------------------------------------------------------------------------------------
void HtollAnalysis::FillReductionVariables(LoopAll& l, int jentry) {
  if(HtollAnalysisDEBUG) 
    cout<<"myFillReduceVar START"<<endl;
  
  if(HtollAnalysisDEBUG) 
    cout<<"myFillReduceVar END"<<endl;
}

// ----------------------------------------------------------------------------------------------------
bool HtollAnalysis::SelectEventsReduction(LoopAll& l, int jentry) {

  // Two muons/electrons with pT > 25
 std::vector<int> eleIndex; 
  
  int goodEl = 0;
  for(int i =0; i<l.el_std_n; i++) {
    TLorentzVector* p4 = (TLorentzVector*)l.el_std_sc->At(i);
    if (p4->Et() > 20.)
      eleIndex.push_back(i);
  }
  
  Float_t mass = 0.;
  if (eleIndex.size() > 1) {
    for (unsigned int i=0; i<eleIndex.size()-1; i++) {
      TLorentzVector* p1 = (TLorentzVector*)l.el_std_sc->At(eleIndex[i]);
      
      for (unsigned int j=i+1; j<eleIndex.size(); j++) {
        TLorentzVector* p2 = (TLorentzVector*)l.el_std_sc->At(eleIndex[j]);
        Float_t mass_temp = ((*p1)+(*p2)).M();
        if (mass < mass_temp) {
          mass = mass_temp;
        }
      }
    }
    
    if (mass > 60.)
      return 1;
  }
    
  return 0;
  
  if(TapAnalysisDEBUG)  cout << " ****************** SelectEventsReduction " << endl;
  return true;
}

// ----------------------------------------------------------------------------------------------------
bool HtollAnalysis::SkimEvents(LoopAll& l, int jentry) {
  return true;
}

// ----------------------------------------------------------------------------------------------------
bool HtollAnalysis::SelectEvents(LoopAll& l, int jentry) {
  return true;
}

// ----------------------------------------------------------------------------------------------------
void HtollAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) 
{}

// ----------------------------------------------------------------------------------------------------
void HtollAnalysis::ResetAnalysis()
{}

bool TapAnalysis::ElectronId(LoopAll& l, Int_t eleIndex, Int_t vertexIndex, std::string type, Float_t selection) {
  bool result = false;

  if (selection == -999.)
    return true;

  if (type == "LeptonTag") {
    //l.rho = l.rho_algo1;
    //result = l.ElectronMVACuts(eleIndex, vertexIndex);

    TLorentzVector* thisel = (TLorentzVector*) l.el_std_p4->At(eleIndex);
    TLorentzVector* thissc = (TLorentzVector*) l.el_std_sc->At(eleIndex);
    float thiseta = fabs(thissc->Eta());
    float thispt = thisel->Pt();
    /*
    double Aeff=0.;
    if(thiseta<1.0)                   Aeff=0.135;
    if(thiseta>=1.0 && thiseta<1.479) Aeff=0.168;
    if(thiseta>=1.479 && thiseta<2.0) Aeff=0.068;
    if(thiseta>=2.0 && thiseta<2.2)   Aeff=0.116;
    if(thiseta>=2.2 && thiseta<2.3)   Aeff=0.162;
    if(thiseta>=2.3 && thiseta<2.4)   Aeff=0.241;
    if(thiseta>=2.4)                  Aeff=0.23;
    */
    double Aeff=0.;
    if(thiseta<1.0)                   Aeff=0.10;
    if(thiseta>=1.0 && thiseta<1.479) Aeff=0.12;
    if(thiseta>=1.479 && thiseta<2.0) Aeff=0.085;
    if(thiseta>=2.0 && thiseta<2.2)   Aeff=0.11;
    if(thiseta>=2.2 && thiseta<2.3)   Aeff=0.12;
    if(thiseta>=2.3 && thiseta<2.4)   Aeff=0.12;
    if(thiseta>=2.4)                  Aeff=0.13;    
    float thisiso=l.el_std_pfiso_charged[eleIndex]+std::max(l.el_std_pfiso_neutral[eleIndex]+l.el_std_pfiso_photon[eleIndex]-l.rho_algo1*Aeff, 0.);
    
    if(vertexIndex!=-1){
      if(fabs(l.el_std_D0Vtx[eleIndex][vertexIndex]) > 0.02) 
        return false;
      if(fabs(l.el_std_DZVtx[eleIndex][vertexIndex]) > 0.2)  
        return false;
    }
    
    if (l.el_std_hp_expin[eleIndex] > 1)
      return false;
    
    if (l.el_std_conv[eleIndex] == 0)
      return false;
        
    result = (l.el_std_mva_nontrig[eleIndex] > selection) && (thisiso/thispt<0.15);


