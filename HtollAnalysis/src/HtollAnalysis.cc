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
  if (l.typerun != l.kReduce) {
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

  if (doMuon) {
    std::vector<int> goodMuons;
    
    for (int i=0; i<l.mu_glo_n; i++) {
      TLorentzVector* p4 = (TLorentzVector*)l.mu_glo_p4->At(i);
      if (p4->Pt() > 20.) {
	if (l.MuonTightID2012(i, -1)) {
	  if (l.MuonIsolation2012(i, p4->Pt(), true)) {
	    goodMuons.push_back(i);
	  }
	}
      }
    }
    
    if (goodMuons.size() < 2)
      return false;
    
    for (unsigned int i=0; i<goodMuons.size()-1; i++) {
      TLorentzVector* p1 = (TLorentzVector*)l.mu_glo_p4->At(goodMuons[i]);
      
      for (unsigned int j=i+1; j<goodMuons.size(); j++) {
	TLorentzVector* p2 = (TLorentzVector*)l.mu_glo_p4->At(goodMuons[j]);
	TLorentzVector higgs = (*p1)+(*p2);
	Float_t mass = higgs.M();
	std::cout << higgs.X() << " " << higgs.Y() << " " << higgs.Z() << " "  << std::endl;
	std::cout << p2->X() << " " << p2->Y() << " " << p2->Z() << " "  << std::endl;
	p2->Boost(-higgs.BoostVector());
	std::cout << p2->X() << " " << p2->Y() << " " << p2->Z() << " "  << std::endl;
	std::cout << p2->Theta() << std::endl;

	l.FillHist("massMu", 0, mass, weight);
	//l.FillHist("theta2", 0, 1/4+3/2*pow(p2->Theta(), 2)+1/4*pow(p2->Theta(), 4), weight);
	l.FillHist("theta2", 0, p2->Theta(), weight);
      }
    }
  } else {
    std::vector<int> goodEles;
    
    for (int i=0; i<l.el_std_n; i++) {
      TLorentzVector* p4 = (TLorentzVector*)l.el_std_p4->At(i);
      if (p4->Pt() > 20.) 
	if (ElectronId(l, i))
	  goodEles.push_back(i);
    }
    
    if (goodEles.size() < 2)
      return false;
    
    for (unsigned int i=0; i<goodEles.size()-1; i++) {
      TLorentzVector* p1 = (TLorentzVector*)l.el_std_p4->At(goodEles[i]);
      
      for (unsigned int j=i+1; j<goodEles.size(); j++) {
	TLorentzVector* p2 = (TLorentzVector*)l.el_std_p4->At(goodEles[j]);
	Float_t mass = ((*p1)+(*p2)).M();
	l.FillHist("massEl", 0, mass, weight);
      }
    }
  }
  
  return 0;
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

  // Two muons/electrons with pT > 20
  int goodMu = 0;
  for (int i=0; i<l.mu_glo_n; i++) {
    TLorentzVector* p4 = (TLorentzVector*)l.mu_glo_p4->At(i);
    if (p4->Pt() > 20.) 
      goodMu++;
  }

  int goodEl = 0;
  for(int i =0; i<l.el_std_n; i++) {
    TLorentzVector* p4 = (TLorentzVector*)l.el_std_sc->At(i);
    if (p4->Et() > 20.)
      goodEl++;
  }

  return (goodMu > 1 || goodEl > 1); 
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

bool HtollAnalysis::ElectronId(LoopAll& l, Int_t eleIndex) {
  
  bool result = false;

  TLorentzVector* thisel = (TLorentzVector*) l.el_std_p4->At(eleIndex);
  TLorentzVector* thissc = (TLorentzVector*) l.el_std_sc->At(eleIndex);
  float thiseta = fabs(thissc->Eta());
  float thispt = thisel->Pt();

  double Aeff=0.;
  if(thiseta<1.0)                   Aeff=0.135;
  if(thiseta>=1.0 && thiseta<1.479) Aeff=0.168;
  if(thiseta>=1.479 && thiseta<2.0) Aeff=0.068;
  if(thiseta>=2.0 && thiseta<2.2)   Aeff=0.116;
  if(thiseta>=2.2 && thiseta<2.3)   Aeff=0.162;
  if(thiseta>=2.3 && thiseta<2.4)   Aeff=0.241;
  if(thiseta>=2.4)                  Aeff=0.23;
  float thisiso=l.el_std_pfiso_charged[eleIndex]+std::max(l.el_std_pfiso_neutral[eleIndex]+l.el_std_pfiso_photon[eleIndex]-l.rho_algo1*Aeff, 0.);
        
  if (l.el_std_hp_expin[eleIndex] > 1)
    return false;
  
  if (l.el_std_conv[eleIndex] == 0)
    return false;
        
  result = (l.el_std_mva_nontrig[eleIndex] > 0.9) && (thisiso/thispt<0.15);
}

