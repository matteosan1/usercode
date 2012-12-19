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
  float pu_weight = 1.;
  
  if (l.itype[l.current] != 0) {
    unsigned int n_pu = l.pu_n;
    pu_weight = getPuWeight(l.pu_n, l.itype[l.current], &(l.sampleContainer[l.current_sample_index]), jentry == 1);
    weight = pu_weight * l.sampleContainer[l.current_sample_index].weight;
    //std::cout << l.sampleContainer[l.current_sample_index].weight << " " << weight/l.sampleContainer[l.current_sample_index].weight << " " << n_pu << std::endl;
  }

  TLorentzVector higgs;
  Int_t cat = 0, vbfcat = 0;
  
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
        higgs = (*p1)+(*p2);
	Float_t mass = higgs.M();
	//std::cout << higgs.X() << " " << higgs.Y() << " " << higgs.Z() << " "  << std::endl;
	//std::cout << p2->X() << " " << p2->Y() << " " << p2->Z() << " "  << std::endl;
	//p2->Boost(-higgs.BoostVector());
	//std::cout << p2->X() << " " << p2->Y() << " " << p2->Z() << " "  << std::endl;
	//std::cout << p2->Theta() << std::endl;
	if (higgs.M() < 130. && higgs.M() > 120. && doBlinding && l.itype[l.current] == 0)
	  return false;

	l.FillHist("massMu", 0, mass, weight);
	//l.FillHist("theta2", 0, 1/4+3/2*pow(p2->Theta(), 2)+1/4*pow(p2->Theta(), 4), weight);
	l.FillHist("theta2", 0, p2->Theta(), weight);
	Tree(l, goodMuons[i], goodMuons[j], higgs, cat, vbfcat, weight, pu_weight, false, "");
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
	higgs = (*p1)+(*p2);

	if (higgs.M() < 130. && higgs.M() > 120. && doBlinding && l.itype[l.current] == 0)
	  return false;

	l.FillHist("massEl", 0, higgs.M(), weight);
	
	Tree(l, goodEles[i], goodEles[j], higgs, cat, vbfcat, weight, pu_weight, false, "");
      }
    }
  }
  
  return true;
}

// ----------------------------------------------------------------------------------------------------
void HtollAnalysis::GetBranches(TTree *t, std::set<TBranch *>& s ) 
{}

// ----------------------------------------------------------------------------------------------------
void HtollAnalysis::FillReductionVariables(LoopAll& l, int jentry) {

  if (l.itype[l.current] < 0) {
    for (int i=0; i<l.gp_n; i++) {
      if (l.gp_pdgid[i] == 25 && l.gp_status[i] == 3) {
	TLorentzVector* p4 = (TLorentzVector*)l.gp_p4->At(i);
	*(l.higgs) = *p4;
	break;
      }
    }

    if (doMuon) {
      for (int i=0; i<l.mu_glo_n; i++) {
	Int_t mc1=-1;
	Int_t mc2=-1;
	Int_t pho=-1;
	
	l.FindMCLeptons(i, mc1, mc2, pho, 13);
	// et, eta, phi
	
	if (mc1 != -1) {
	  TLorentzVector* p4 = (TLorentzVector*)l.gp_p4->At(mc1);
	  l.mc_et[i] = p4->Et();
	  l.mc_phi[i] = p4->Phi();
	  l.mc_eta[i] = p4->Eta();
	} else {
	  l.mc_et[i] = -999.;
	  l.mc_phi[i] = -999.;
	  l.mc_eta[i] = -999.;
	}
	
	if (pho != -1) {
	  TLorentzVector* p4 = (TLorentzVector*)l.gp_p4->At(pho);
	  l.fsr_et[i] = p4->Et();
	  l.fsr_phi[i] = p4->Phi();
	  l.fsr_eta[i] = p4->Eta();
	} else {
	  l.fsr_et[i] = -999.;
	  l.fsr_phi[i] = -999.;
	  l.fsr_eta[i] = -999.;
	}
      }
    } else {
      for (int i=0; i<l.el_std_n; i++) {
	Int_t mc1=-1;
	Int_t mc2=-1;
	Int_t pho=-1;
	
	l.FindMCLeptons(i, mc1, mc2, pho, 11);
	// et, eta, phi
	
	if (mc1 != -1) {
	  TLorentzVector* p4 = (TLorentzVector*)l.gp_p4->At(mc1);
	  l.mc_et[i] = p4->Et();
	  l.mc_phi[i] = p4->Phi();
	  l.mc_eta[i] = p4->Eta();
	} else {
	  l.mc_et[i] = -999.;
	  l.mc_phi[i] = -999.;
	  l.mc_eta[i] = -999.;
	}
	
	if (pho != -1) {
	  TLorentzVector* p4 = (TLorentzVector*)l.gp_p4->At(pho);
	  l.fsr_et[i] = p4->Et();
	  l.fsr_phi[i] = p4->Phi();
	  l.fsr_eta[i] = p4->Eta();
	} else {
	  l.fsr_et[i] = -999.;
	  l.fsr_phi[i] = -999.;
	  l.fsr_eta[i] = -999.;
	}
      }
    }
  } else {
    if (doMuon) {
      for (int i=0; i<l.mu_glo_n; i++) {
	l.mc_et[i] = -999.;
	l.mc_phi[i] = -999.;
	l.mc_eta[i] = -999.;
	
	l.fsr_et[i] = -999.;
	l.fsr_phi[i] = -999.;
	l.fsr_eta[i] = -999.;
      }
    } else {
      for (int i=0; i<l.el_std_n; i++) {
	l.mc_et[i] = -999.;
	l.mc_phi[i] = -999.;
	l.mc_eta[i] = -999.;
	
	l.fsr_et[i] = -999.;
	l.fsr_phi[i] = -999.;
	l.fsr_eta[i] = -999.;
      }
    }
  }
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
void HtollAnalysis::ReducedOutputTree(LoopAll &l, TTree * outputTree) {
  l.higgs = new TLorentzVector(0,0,0,0);

  l.Branch_mc_et(outputTree);
  l.Branch_mc_eta(outputTree);
  l.Branch_mc_phi(outputTree);
  l.Branch_fsr_et(outputTree);
  l.Branch_fsr_eta(outputTree);
  l.Branch_fsr_phi(outputTree);
  l.Branch_higgs(outputTree);
}

// ----------------------------------------------------------------------------------------------------
void HtollAnalysis::ResetAnalysis()
{}

void HtollAnalysis::Tree(LoopAll& l, Int_t lept1, Int_t lept2, const TLorentzVector & Higgs, Int_t cat, Int_t vbfcat, 
			 Float_t weight, Float_t pu_weight, bool isSyst, std::string name1) {
    
    l.FillTree("run", (float)l.run);
    l.FillTree("lumis", (float)l.lumis);
    l.FillTree("event", (double)l.event);
    l.FillTree("itype", (float)l.itype[l.current]);
    l.FillTree("nvtx", (float)l.vtx_std_n);
    l.FillTree("domuon", (int)doMuon);
    l.FillTree("rho", (float)l.rho_algo1);        
    l.FillTree("mass", (float)Higgs.M());
    TLorentzVector* lep1=0;
    TLorentzVector* lep2=0;

    if (doMuon) {
      lep1 = (TLorentzVector*)l.mu_glo_p4->At(lept1);
      l.FillTree("et1", (float)lep1->Et());
      l.FillTree("eta1", (float)lep1->Eta());
      l.FillTree("phi1", (float)lep1->Phi());
      l.FillTree("mutype1", (int)l.mu_glo_type[lept1]);
      l.FillTree("muchi21", (float)l.mu_glo_chi2[lept1]/l.mu_glo_dof[lept1]);
      l.FillTree("much1", (int)l.mu_glo_validChmbhits[lept1]);
      l.FillTree("munmatch1", (int)l.mu_glo_nmatches[lept1]);
      l.FillTree("mupixhit1", (int)l.mu_glo_pixelhits[lept1]);
      l.FillTree("mutklay1", (int)l.mu_tkLayers[lept1]);
      l.FillTree("mudb1", (float)l.mu_dbCorr[lept1]);
      l.FillTree("chiso1", l.mu_glo_chhadiso04[lept1]);
      l.FillTree("neiso1", l.mu_glo_nehadiso04[lept1]);
      l.FillTree("phiso1", l.mu_glo_photiso04[lept1]);
      l.FillTree("mc_et1", l.mc_et[lept1]);
      l.FillTree("mc_eta1", (float)l.mc_eta[lept1]);
      l.FillTree("mc_phi1", (float)l.mc_phi[lept1]);
      l.FillTree("fsr_et1", (float)l.fsr_et[lept1]);
      l.FillTree("fsr_eta1", (float)l.fsr_eta[lept1]);
      l.FillTree("fsr_phi1", (float)l.fsr_phi[lept1]);

      lep2 = (TLorentzVector*)l.mu_glo_p4->At(lept2);
      l.FillTree("et2", (float)lep2->Et());
      l.FillTree("eta2", (float)lep2->Eta());
      l.FillTree("phi2", (float)lep2->Phi());
      l.FillTree("mutype2", (int)l.mu_glo_type[lept2]);
      l.FillTree("muchi22", (float)l.mu_glo_chi2[lept2]/l.mu_glo_dof[lept2]);
      l.FillTree("much2", (int)l.mu_glo_validChmbhits[lept2]);
      l.FillTree("munmatch2", (int)l.mu_glo_nmatches[lept2]);
      l.FillTree("mupixhit2", (int)l.mu_glo_pixelhits[lept2]);
      l.FillTree("mutklay2", (int)l.mu_tkLayers[lept2]);
      l.FillTree("mudb2", (float)l.mu_dbCorr[lept2]);
      l.FillTree("chiso2", l.mu_glo_chhadiso04[lept2]);
      l.FillTree("neiso2", l.mu_glo_nehadiso04[lept2]);
      l.FillTree("phiso2", l.mu_glo_photiso04[lept2]);
      l.FillTree("elmisshits1", (int)99);
      l.FillTree("elconv1", (int)99);
      l.FillTree("elmva1", (float)99.);
      l.FillTree("elregr1", (float)9999.);
      l.FillTree("elregr_err1", (float)9999.);
      l.FillTree("elmisshits2", (int)99);
      l.FillTree("elconv2", (int)99);
      l.FillTree("elmva2", (float)99.);
      l.FillTree("elregr2", (float)9999.);
      l.FillTree("elregr_err2", (float)9999.);
      l.FillTree("elpt1", (float)-9999.);
      l.FillTree("elpt2", (float)-9999.);
      l.FillTree("mc_et2", (float)l.mc_et[lept2]);
      l.FillTree("mc_eta2", (float)l.mc_eta[lept2]);
      l.FillTree("mc_phi2", (float)l.mc_phi[lept2]);
      l.FillTree("fsr_et2", (float)l.fsr_et[lept2]);
      l.FillTree("fsr_eta2", (float)l.fsr_eta[lept2]);
      l.FillTree("fsr_phi2", (float)l.fsr_phi[lept2]);
      //l.FillTree("cosDphi", (float)TMath::Cos(lead_p4.Phi()-sublead_p4.Phi()));

    } else {
      lep1 = (TLorentzVector*)l.el_std_sc->At(lept1);
      l.FillTree("et1", (float)lep1->Et());
      l.FillTree("eta1", (float)lep1->Eta());
      l.FillTree("phi1", (float)lep1->Phi());
      l.FillTree("chiso1", l.el_std_pfiso_charged[lept1]);
      l.FillTree("neiso1", l.el_std_pfiso_neutral[lept1]);
      l.FillTree("phiso1", l.el_std_pfiso_photon[lept1]);
      l.FillTree("elmisshits1", (int)l.el_std_hp_expin[lept1]);
      l.FillTree("elconv1", (int)l.el_std_conv[lept1]);
      l.FillTree("elmva1", (float)l.el_std_mva_trig[lept1]);
      l.FillTree("elregr1", (float)l.el_std_regr_energy[lept1]);
      l.FillTree("elregr_err1", (float)l.el_std_regr_energyerr[lept1]);
      lep1 = (TLorentzVector*)l.el_std_p4->At(lept1);
      l.FillTree("elpt1", (float)(l.el_std_pin[lept1]*sin(lep1->Theta())));
      l.FillTree("mc_et1", (float)l.mc_et[lept1]);
      l.FillTree("mc_eta1", (float)l.mc_eta[lept1]);
      l.FillTree("mc_phi1", (float)l.mc_phi[lept1]);
      l.FillTree("fsr_et1", (float)l.fsr_et[lept1]);
      l.FillTree("fsr_eta1", (float)l.fsr_eta[lept1]);
      l.FillTree("fsr_phi1", (float)l.fsr_phi[lept1]);
      
      lep2 = (TLorentzVector*)l.el_std_sc->At(lept2);
      l.FillTree("et2", (float)lep2->Et());
      l.FillTree("eta2", (float)lep2->Eta());
      l.FillTree("phi2", (float)lep2->Phi());
      l.FillTree("chiso2", l.el_std_pfiso_charged[lept2]);
      l.FillTree("neiso2", l.el_std_pfiso_neutral[lept2]);
      l.FillTree("phiso2", l.el_std_pfiso_photon[lept2]);
      l.FillTree("elmisshits2", (int)l.el_std_hp_expin[lept2]);
      l.FillTree("elconv2", (int)l.el_std_conv[lept2]);
      l.FillTree("elmva2", (float)l.el_std_mva_trig[lept2]);
      l.FillTree("elregr2", (float)l.el_std_regr_energy[lept2]);
      l.FillTree("elregr_err2", (float)l.el_std_regr_energyerr[lept2]);
      lep2 = (TLorentzVector*)l.el_std_p4->At(lept2);
      l.FillTree("elpt2", (float)(l.el_std_pin[lept2]*sin(lep2->Theta())));
      l.FillTree("mc_et2", (float)l.mc_et[lept2]);
      l.FillTree("mc_eta2", (float)l.mc_eta[lept2]);
      l.FillTree("mc_phi2", (float)l.mc_phi[lept2]);
      l.FillTree("fsr_et2", (float)l.fsr_et[lept2]);
      l.FillTree("fsr_eta2", (float)l.fsr_eta[lept2]);
      l.FillTree("fsr_phi2", (float)l.fsr_phi[lept2]);

      l.FillTree("mutype1"  ,(int)-9999);
      l.FillTree("muchi21"  ,(float)-9999.);
      l.FillTree("much1"    ,(int)-9999);
      l.FillTree("munmatch1",(int)-9999);
      l.FillTree("mupixhit1",(int)-9999);
      l.FillTree("mutklay1" ,(int)-9999);
      l.FillTree("mudb1"    ,(float)-9999.);
      l.FillTree("mutype2"  ,(int)-9999);
      l.FillTree("muchi22"  ,(float)-9999.);
      l.FillTree("much2"    ,(int)-9999);
      l.FillTree("munmatch2",(int)-9999);
      l.FillTree("mupixhit2",(int)-9999);
      l.FillTree("mutklay2" ,(int)-9999);
      l.FillTree("mudb2"    ,(float)-9999.);

      //l.FillTree("cosDphi", (float)TMath::Cos(lead_p4.Phi()-sublead_p4.Phi()));
    }
    
    //l.FillTree("genmatch1", (float)l.pho_genmatched[diphoton_index.first]);
    //l.FillTree("genmatch2", (float)l.pho_genmatched[diphoton_index.second]);
    l.FillTree("xsec_weight", (float)l.sampleContainer[l.current_sample_index].weight);
    l.FillTree("full_weight", (float)weight);
    l.FillTree("pu_weight", (float)pu_weight);
    l.FillTree("pu_n", (float)l.pu_n);
    l.FillTree("mass", (float)Higgs.M());
    l.FillTree("vbfcat", (int)vbfcat);
    l.FillTree("MET", (float)l.met_pfmet);
    l.FillTree("MET_phi", (float)l.met_phi_pfmet);
    l.FillTree("cat", (int)cat);
    
    l.FillTree("higgs_pt", (float)l.higgs->Pt());
    l.FillTree("higgs_phi", (float)l.higgs->Phi());
    l.FillTree("higgs_eta", (float)l.higgs->Eta());
    l.FillTree("higgs_mass", (float)l.higgs->M());
    
    // FIXME 0 vtx for the moment
    TVector3* vtx = (TVector3*)l.vtx_std_xyz->At(0);
    l.FillTree("vtx_x", (float)vtx->X());
    l.FillTree("vtx_y", (float)vtx->Y());
    l.FillTree("vtx_z", (float)vtx->Z());

    //l.FillTree("issyst", (int)isSyst);
    //l.FillTree("name1", name1);
    
    //int pass_hlt = checkEventHLT(l, hltPaths);  
    //l.FillTree("pass_hlt", pass_hlt);
    float dijet_deta; 
    float dijet_mjj;
    float dijet_zep;
    float dijet_dphi_ll_jj;
    float dijet_j1pt;
    float dijet_j2pt;
    bool dijet_has2jets = DijetPreSelection(l, lep1, lep2, 
					    dijet_deta, dijet_mjj, dijet_zep, dijet_dphi_ll_jj, dijet_j1pt, dijet_j2pt);
    l.FillTree("dijet_deta",          (float)dijet_deta);
    l.FillTree("dijet_mjj",           (float)dijet_mjj);
    l.FillTree("dijet_zep",           (float)dijet_zep);
    l.FillTree("dijet_dphi_ll_jj",    (float)dijet_dphi_ll_jj);
    l.FillTree("dijet_j1pt",          (float)dijet_j1pt);
    l.FillTree("dijet_j2pt",          (float)dijet_j2pt);
    l.FillTree("dijet_has2jets",      (float)dijet_has2jets);
}

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
        
  result = (l.el_std_mva_trig[eleIndex] > 0.9) && (thisiso/thispt<0.15);
  
  return result;
}



bool HtollAnalysis::DijetPreSelection(LoopAll& l, TLorentzVector* veto_p41, TLorentzVector* veto_p42, 
				      float & dijet_deta, float & dijet_mjj, float & dijet_zep, float & dijet_dphi_ll_jj, 
				      float & dijet_j1pt, float & dijet_j2pt) {
  bool exist=false;
  
  std::pair<int, int> myjets = l.Select2HighestPtJets(*veto_p41, *veto_p42 ); // Bool_t * jetid_flags)
  
  if( myjets.first==-1 || myjets.second == -1) { // set defaults
    dijet_deta        = -99;
    dijet_mjj         = -99;
    dijet_zep         = -99;
    dijet_dphi_ll_jj  = -99;
    dijet_j1pt        = -99;
    dijet_j2pt        = -99;
    
    exist             = false;
  } else { // get jets and get values
    TLorentzVector* jet1 = (TLorentzVector*) l.jet_algoPF1_p4->At(myjets.first);
    TLorentzVector* jet2 = (TLorentzVector*) l.jet_algoPF1_p4->At(myjets.second);
    TLorentzVector jj = (*jet1) + (*jet2);
    TLorentzVector ll = *veto_p41 + *veto_p42;
    
    dijet_deta        = abs(jet1->Eta() - jet2->Eta());
    dijet_mjj         = jj.M();
    dijet_zep         = (ll.Eta() - 0.5*(jet1->Eta()+jet2->Eta()));
    dijet_dphi_ll_jj  = abs(jj.DeltaPhi(ll));
    dijet_j1pt        = jet1->Pt();
    dijet_j2pt        = jet2->Pt();
    
    exist             = true;
  }
  return exist;
}
