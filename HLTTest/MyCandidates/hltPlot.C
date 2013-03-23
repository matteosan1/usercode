#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TLorentzVector.h"
#include "TCanvas.h"
#include "TProfile.h"

#include "TROOT.h"
#include "TStyle.h"

#include <iostream> 

void drawDivision(TH1F* hlt, TH1F* offline) {
  TH1F* den = (TH1F*) hlt->Clone("hlt");
  TH1F* num = (TH1F*) offline->Clone("offline");
  num->Sumw2();
  den->Sumw2();
  num->Divide(den);
  num->GetYaxis()->SetRangeUser(0.6, 1.4);
  num->SetMarkerStyle(20);
  num->SetMarkerColor(kBlue);
  num->SetLineColor(kBlue);
  num->DrawClone("PE");
}

void hltPlot(char* filename) {

  gROOT->GetStyle("MyStyle")->SetOptStat(0);
  gROOT->GetStyle("MyStyle")->SetPalette(1);

  TH2F* diff_vs_pt[2];
  TH2F* diff_vs_eta[2];
  diff_vs_pt[0] = new TH2F("diff_vs_pt0", "diff_vs_pt0", 100, -.5, .5, 100, 0., 100);
  diff_vs_pt[1] = new TH2F("diff_vs_pt1", "diff_vs_pt1", 100, -.5, .5, 100, 0., 100);
  diff_vs_eta[0] = new TH2F("diff_vs_eta0", "diff_vs_eta0", 100, -.5, .5, 50, -2.5, 2.5);
  diff_vs_eta[1] = new TH2F("diff_vs_eta1", "diff_vs_eta1", 100, -.5, .5, 50, -2.5, 2.5);

  TH1F* diff[9][2];
  diff[0][0] = new TH1F("diff00", "diff00", 100, -.5, .5);
  diff[1][0] = new TH1F("diff10", "diff10", 100, -.5, .5);
  diff[2][0] = new TH1F("diff20", "diff20", 100, -.5, .5);
  diff[3][0] = new TH1F("diff30", "diff30", 100, -.5, .5);
  diff[4][0] = new TH1F("diff40", "diff40", 100, -.5, .5);
  diff[5][0] = new TH1F("diff50", "diff50", 100, -.5, .5);
  diff[6][0] = new TH1F("diff60", "diff60", 100, -.5, .5);
  diff[7][0] = new TH1F("diff70", "diff70", 100, -.5, .5);
  diff[8][0] = new TH1F("diff80", "diff80", 100, -.5, .5);

  diff[0][1] = new TH1F("diff01", "diff01", 100, -.5, .5);
  diff[1][1] = new TH1F("diff11", "diff11", 100, -.5, .5);
  diff[2][1] = new TH1F("diff21", "diff21", 100, -.5, .5);
  diff[3][1] = new TH1F("diff31", "diff31", 100, -.5, .5);
  diff[4][1] = new TH1F("diff41", "diff41", 100, -.5, .5);
  diff[5][1] = new TH1F("diff51", "diff51", 100, -.5, .5);
  diff[6][1] = new TH1F("diff61", "diff61", 100, -.5, .5);
  diff[7][1] = new TH1F("diff71", "diff71", 100, -.5, .5);
  diff[8][1] = new TH1F("diff81", "diff81", 100, -.5, .5);

  TH1I* hMishits[2][2];
  hMishits[0][0] = new TH1I("hMishits0",    "hMishits0", 10, 0, 10);
  hMishits[0][1] = new TH1I("hMishits1",    "hMishits1", 10, 0, 10);  
  hMishits[1][0] = new TH1I("hMishitsOff0", "hMishitsOff0", 10, 0, 10);
  hMishits[1][1] = new TH1I("hMishitsOff1", "hMishitsOff1", 10, 0, 10);

  TH1F* hSieie[2][2];
  hSieie[0][0] = new TH1F("hSieie0", "hSieie0", 100, 0.005, 0.015);
  hSieie[0][1] = new TH1F("hSieie1", "hSieie1", 100, 0.015, 0.035);  
  hSieie[1][0] = new TH1F("hSieieOff0", "hSieieOff0", 100, 0.005, 0.015);
  hSieie[1][1] = new TH1F("hSieieOff1", "hSieieOff1", 100, 0.015, 0.035);

  TH1F* hHoe[2][2];
  hHoe[0][0] = new TH1F("hHoe0", "hHoe0", 100, 0, 10);
  hHoe[0][1] = new TH1F("hHoe1", "hHoe1", 100, 0, 10);  
  hHoe[1][0] = new TH1F("hHoeOff0", "hHoeOff0", 100, 0, 10);
  hHoe[1][1] = new TH1F("hHoeOff1", "hHoeOff1", 100, 0, 10);

  TH1F* hDeta[2][2];
  hDeta[0][0] = new TH1F("hDeta0", "hDeta0", 50, 0, 0.005);
  hDeta[0][1] = new TH1F("hDeta1", "hDeta1", 50, 0, 0.005);  
  hDeta[1][0] = new TH1F("hDetaOff0", "hDetaOff0", 50, 0, 0.005);
  hDeta[1][1] = new TH1F("hDetaOff1", "hDetaOff1", 50, 0, 0.005);

  TH1F* hDphi[2][2];
  hDphi[0][0] = new TH1F("hDphi0", "hDphi0", 50, 0, 0.05);
  hDphi[0][1] = new TH1F("hDphi1", "hDphi1", 50, 0, 0.05);  
  hDphi[1][0] = new TH1F("hDphiOff0", "hDphiOff0", 50, 0, 0.05);
  hDphi[1][1] = new TH1F("hDphiOff1", "hDphiOff1", 50, 0, 0.05);

  TH1F* hTkiso[2][2];
  hTkiso[0][0] = new TH1F("hTkiso0", "hTkiso0", 50, 0, 5);
  hTkiso[0][1] = new TH1F("hTkiso1", "hTkiso1", 50, 0, 5);  
  hTkiso[1][0] = new TH1F("hTkisoOff0", "hTkisoOff0", 50, 0, 5);
  hTkiso[1][1] = new TH1F("hTkisoOff1", "hTkisoOff1", 50, 0, 5);

  TH1F* hHcal[2][2];
  hHcal[0][0] = new TH1F("hHcal0", "hHcal0", 100, 0, 4);
  hHcal[0][1] = new TH1F("hHcal1", "hHcal1", 100, 0, 4);  
  hHcal[1][0] = new TH1F("hHcalOff0", "hHcalOff0", 100, 0, 4);
  hHcal[1][1] = new TH1F("hHcalOff1", "hHcalOff1", 100, 0, 4);

  TH1F* hEt[2][2];
  hEt[0][0] = new TH1F("hEt0", "hEt0", 100, 0, 100);
  hEt[0][1] = new TH1F("hEt1", "hEt1", 100, 0, 100);  
  hEt[1][0] = new TH1F("hEtOff0", "hEtOff0", 100, 0, 100);
  hEt[1][1] = new TH1F("hEtOff1", "hEtOff1", 100, 0, 100);

  TH1F* hEta[2];
  hEta[0] = new TH1F("hEta0", "hEta0", 50, -2.5, 2.5);
  hEta[1] = new TH1F("hEtaOff0", "hEtaOff0", 50, -2.5, 2.5);

  TH1F* hIso[2][2];
  hIso[0][0] = new TH1F("hIso0", "hIso0", 100, -10, 10);
  hIso[0][1] = new TH1F("hIso1", "hIso1", 100, -10, 10);  
  hIso[1][0] = new TH1F("hIsoOff0", "hIsoOff0", 100, -10, 10);
  hIso[1][1] = new TH1F("hIsoOff1", "hIsoOff1", 100, -10, 10);

  TH2F* hcorr0 = new TH2F("hC1", "hC1", 100, -10, 10, 100, -10, 10);  
  TH2F* hcorr1 = new TH2F("hC2", "hC2", 100, -10, 10, 100, -10, 10);
  
  TProfile* effAreaReco_iso = new TProfile("profile1", "profile1", 20, 0, 40);
  TProfile* effAreaHlt_iso = new TProfile("profile2", "profile2",  20, 0, 40);
  TProfile* effAreaReco_rho = new TProfile("profile3", "profile3", 20, 0, 40);
  TProfile* effAreaHlt_rho = new TProfile("profile4", "profile4",  20, 0, 40);

  TH1F* hRho[3];
  hRho[0] = new TH1F("hRho0", "hRho0", 100, 0, 50);
  hRho[1] = new TH1F("hRho1", "hRho1", 100, 0, 50);
  hRho[2] = new TH1F("hRho2", "hRho2", 100, 0, 50);

  TFile* file = new TFile(filename);
  TTree* tree = (TTree*)file->Get("event");

  Int_t n0, n1, el_n, el_mishits[10], mishits0[8];
  Float_t eta0[8], eta1[8], pt0[8], pt1[8], phi0[8], phi1[8];
  Float_t iso0[8], iso1[8];
  Float_t deta0[8], dphi0[8], sieie0[8], hoe0[8], hcaliso0[8], tkiso0[8];
  Float_t el_deta[8], el_dphi[8], el_sieie[8], el_hoe[8], el_hcaliso[8], el_tkiso[8];
  Float_t el_et[10], el_eta[10], el_phi[10], el_iso[10];
  Int_t nvtxReco, nvtxHlt;
  Float_t rhoReco1, rhoReco2, rhoHlt;

  Int_t barrel;
  Float_t online, offline;

  tree->SetBranchAddress("rhoReco1", &rhoReco1);
  tree->SetBranchAddress("rhoReco2", &rhoReco2);
  tree->SetBranchAddress("rhoHlt",   &rhoHlt);
  tree->SetBranchAddress("nvtxReco", &nvtxReco);
  tree->SetBranchAddress("nvtxHlt",  &nvtxHlt);

  tree->SetBranchAddress("n0", &n0);
  tree->SetBranchAddress("n1", &n1);
  tree->SetBranchAddress("eta0", &eta0);
  tree->SetBranchAddress("phi0", &phi0);
  tree->SetBranchAddress("pt0", &pt0);
  tree->SetBranchAddress("eta1", &eta1);
  tree->SetBranchAddress("phi1", &phi1);
  tree->SetBranchAddress("pt1", &pt1);
  tree->SetBranchAddress("ecaliso0", &iso0);
  tree->SetBranchAddress("deta0", &deta0);
  tree->SetBranchAddress("dphi0", &dphi0);
  tree->SetBranchAddress("sieie0", &sieie0);
  tree->SetBranchAddress("hoe0", &hoe0);
  tree->SetBranchAddress("hcaliso0", &hcaliso0);
  tree->SetBranchAddress("tkiso0", &tkiso0);
  tree->SetBranchAddress("mishits0", &mishits0);

  tree->SetBranchAddress("el_n", &el_n);
  tree->SetBranchAddress("el_et", &el_et);
  tree->SetBranchAddress("el_eta", &el_eta);
  tree->SetBranchAddress("el_phi", &el_phi);
  tree->SetBranchAddress("el_iso", &el_iso);  
  tree->SetBranchAddress("el_deta", &el_deta);  
  tree->SetBranchAddress("el_dphi", &el_dphi);  
  tree->SetBranchAddress("el_sieie", &el_sieie);  
  tree->SetBranchAddress("el_hoe", &el_hoe);  
  tree->SetBranchAddress("el_hcaliso", &el_hcaliso);  
  tree->SetBranchAddress("el_tkiso", &el_tkiso);  
  tree->SetBranchAddress("el_mishits", &el_mishits);  

  TFile* output = new TFile("output_C.root", "recreate");
  TTree* outTree = new TTree("tree", "tree");
  outTree->Branch("barrel", &barrel, "barrel/I");
  outTree->Branch("online", &online, "online/F");
  outTree->Branch("offline", &offline, "offline/F");

  file->cd();
  
  Int_t entries = tree->GetEntries();
  for (int z=0; z<entries; z++) {
    tree->GetEntry(z);

    Float_t eAreaHLTEB = 0.1523;
    Float_t eAreaHLTEE = 0.1002;
    Float_t eAreaRecoEB = 0.1024;
    Float_t eAreaRecoEE = 0.0594;
    
    Float_t Zmass = 0;
    Int_t index0 = -1;
    //Int_t index1 = -1;
    for (int i=0; i<n0; i++) {
      for (int j=0; j<n1; j++) {
	TLorentzVector e1, e2;
	e1.SetPtEtaPhiM(pt0[i], eta0[i], phi0[i], 0.000511);
	e2.SetPtEtaPhiM(pt1[j], eta1[j], phi1[j], 0.000511);
	//std::cout << pt0[0] << std::endl;
	//std::cout << pt1[1] << std::endl;
	
	TLorentzVector zeta = e1+e2;
	if (zeta.M() > 70 && zeta.M()<110 && zeta.M() > Zmass) {
	  Zmass = zeta.M();
	  index0 = i;
	  //index1 = j;
	}
      }
    }    
    
    if (index0 != -1) {
      TLorentzVector hltProbe;
      hltProbe.SetPtEtaPhiM(pt0[index0], eta0[index0], phi0[index0], 0.000511);
      
      int probeIndex = -1;
      float dRmin = 0.1;
      for (int i=0; i<el_n; i++) {
	if (el_et[i] < 8.)
	  continue;
	TLorentzVector gsf1;
	gsf1.SetPtEtaPhiM(el_et[i], el_eta[i], el_phi[i], 0.000511);
		
	int dr = gsf1.DeltaR(hltProbe);
	if (dr < dRmin && ((fabs(el_et[i]-pt0[index0])/el_et[i])<0.01)) {
	  dr = dRmin;
	  probeIndex = i;
	}
      }
      TLorentzVector gsf;
      gsf.SetPtEtaPhiM(el_et[probeIndex], el_eta[probeIndex], el_phi[probeIndex], 0.000511);

      if (probeIndex != -1) {    

	hRho[0]->Fill(rhoReco1);
	hRho[1]->Fill(rhoReco2);
	hRho[2]->Fill(rhoHlt);
	
	//std::cout << nvtxReco << " " << nvtxHlt << std::endl;

    
	//iso0[index0] = iso0[index0] - rhoHlt*eAreaHLT;
	//el_iso[probeIndex] = el_iso[probeIndex] - rhoReco2*eAreaReco;
	online = iso0[index0];
	offline = el_iso[probeIndex];

	if (fabs(el_eta[probeIndex]) < 1.479) {
	  effAreaReco_iso->Fill(nvtxReco, el_iso[probeIndex]-rhoReco1*eAreaRecoEB);
	  effAreaHlt_iso->Fill(nvtxHlt, iso0[index0]-rhoHlt*eAreaHLTEB);
	  effAreaReco_rho->Fill(nvtxReco, el_iso[probeIndex]);
	  effAreaHlt_rho->Fill(nvtxHlt, iso0[index0]);
	  
	  barrel = 1;
	  hIso  [1][0]->Fill(el_iso[probeIndex]);
	  hSieie[1][0]->Fill(el_sieie[probeIndex]);
	  hDeta [1][0]->Fill(fabs(el_deta[probeIndex]));
	  hDphi [1][0]->Fill(fabs(el_dphi[probeIndex]));
	  hHoe  [1][0]->Fill(el_hoe[probeIndex]);
	  hHcal [1][0]->Fill(el_hcaliso[probeIndex]);
	  hTkiso[1][0]->Fill(el_tkiso[probeIndex]);
	  hEt   [1][0]->Fill(el_et[probeIndex]);
	  hEta  [1]->Fill(el_eta[probeIndex]);
	  hMishits[1][0]->Fill(el_mishits[probeIndex]);
	  hcorr0->Fill(el_iso[probeIndex], iso0[index0]);
	  hIso  [0][0]->Fill(iso0[index0]);
	  hSieie[0][0]->Fill(sieie0[index0]);
	  hDeta [0][0]->Fill(deta0[index0]);
	  hDphi [0][0]->Fill(dphi0[index0]);
	  hHoe  [0][0]->Fill(hoe0[index0]);
	  hHcal [0][0]->Fill(hcaliso0[index0]);
	  hTkiso[0][0]->Fill(tkiso0[index0]);
	  hEt   [0][0]->Fill(pt0[index0]);
	  hEta  [0]->Fill(eta0[index0]);
	  hMishits[0][0]->Fill(mishits0[index0]);
	  
	  	  
	  diff_vs_pt[0]->Fill((el_tkiso[probeIndex]-tkiso0[index0])/el_tkiso[probeIndex], el_et[probeIndex]);
	  diff_vs_eta[0]->Fill((el_tkiso[probeIndex]-tkiso0[index0])/el_tkiso[probeIndex], el_eta[probeIndex]);
	  diff[0][0]->Fill((el_tkiso[probeIndex]-tkiso0[index0])/el_tkiso[probeIndex]);
	  diff[1][0]->Fill((el_iso[probeIndex]-iso0[index0])/el_iso[probeIndex]);
	  diff[2][0]->Fill((el_sieie[probeIndex]-sieie0[index0])/el_sieie[probeIndex]);
	  diff[3][0]->Fill((fabs(el_dphi[probeIndex])-dphi0[index0])/fabs(el_dphi[probeIndex]));
	  diff[4][0]->Fill((fabs(el_deta[probeIndex])-deta0[index0])/fabs(el_deta[probeIndex]));
	  diff[5][0]->Fill((el_hoe[probeIndex] - (hoe0[index0]))/el_hoe[probeIndex]);
	  diff[6][0]->Fill((el_hcaliso[probeIndex]-hcaliso0[index0])/el_hcaliso[probeIndex]);
	  diff[7][0]->Fill((el_et[probeIndex]-pt0[index0])/el_et[probeIndex]);
	  diff[8][0]->Fill((el_eta[probeIndex]-eta0[index0])/el_eta[probeIndex]);
	} else {

	  barrel = 0;
	  diff_vs_pt[1]->Fill((el_tkiso[probeIndex]-tkiso0[index0])/el_tkiso[probeIndex], el_et[probeIndex]);
	  diff_vs_eta[1]->Fill((el_tkiso[probeIndex]-tkiso0[index0])/el_tkiso[probeIndex], el_eta[probeIndex]);
	  hIso[1][1]->Fill(el_iso[probeIndex]);
	  hSieie[1][1]->Fill(el_sieie[probeIndex]);
	  hDeta [1][1]->Fill(fabs(el_deta[probeIndex]));
	  hDphi [1][1]->Fill(fabs(el_dphi[probeIndex]));
	  hHoe  [1][1]->Fill(el_hoe[probeIndex]);
	  hHcal [1][1]->Fill(el_hcaliso[probeIndex]);
	  hTkiso[1][1]->Fill(el_tkiso[probeIndex]);
	  hEt   [1][1]->Fill(el_et[probeIndex]);
	  hEta  [1]->Fill(el_eta[probeIndex]);
	  hMishits[1][1]->Fill(el_mishits[probeIndex]);
	  hcorr1->Fill(el_iso[probeIndex], iso0[index0]);
	  hIso  [0][1]->Fill(iso0[index0]);
	  hSieie[0][1]->Fill(sieie0[index0]);
	  hDeta [0][1]->Fill(deta0[index0]);
	  hDphi [0][1]->Fill(dphi0[index0]);
	  hHoe  [0][1]->Fill(hoe0[index0]);
	  hHcal [0][1]->Fill(hcaliso0[index0]);
	  hTkiso[0][1]->Fill(tkiso0[index0]);
	  hEt   [0][1]->Fill(pt0[index0]);
	  hEta  [0]->Fill(eta0[index0]);
	  hMishits[0][1]->Fill(mishits0[index0]);
	  
	  diff[0][1]->Fill((el_tkiso[probeIndex]-tkiso0[index0])/el_tkiso[probeIndex]);
	  diff[1][1]->Fill((el_iso[probeIndex]-iso0[index0])/el_iso[probeIndex]);
	  diff[2][1]->Fill((el_sieie[probeIndex]-sieie0[index0])/el_sieie[probeIndex]);
	  diff[3][1]->Fill((fabs(el_dphi[probeIndex])-dphi0[index0])/fabs(el_dphi[probeIndex]));
	  diff[4][1]->Fill((fabs(el_deta[probeIndex])-deta0[index0])/fabs(el_deta[probeIndex]));
	  diff[5][1]->Fill((el_hoe[probeIndex] - (hoe0[index0]))/el_hoe[probeIndex]);
	  diff[6][1]->Fill((el_hcaliso[probeIndex]-hcaliso0[index0])/el_hcaliso[probeIndex]);
	  diff[7][1]->Fill((el_et[probeIndex]-pt0[index0])/el_et[probeIndex]);
	  diff[8][1]->Fill((el_eta[probeIndex]-eta0[index0])/el_eta[probeIndex]);
	  
	  //if (((el_iso[probeIndex]-iso0[index0])/el_iso[probeIndex]) == 0)
	  //  std::cout << el_iso[probeIndex] << " " << iso0[index0] << std::endl;
	  
	}
      }
    }
    outTree->Fill();
  }
  
  char a[100];
  /*
  TCanvas* c_1 = new TCanvas("c_1", "c_1", 800, 600);
  hMishits[1][0]->SetMarkerStyle(20);
  hMishits[1][1]->SetMarkerStyle(20);
  hMishits[0][0]->SetFillColor(kRed);
  hMishits[0][1]->SetFillColor(kRed);
  c_1->Divide(2,1);
  for(int i=0; i<2; i++) {
    c_1->cd(i+1);
    hMishits[0][i]->Draw();
    hMishits[1][i]->Draw("PESAME");
  }
  sprintf(a, "%s_%s.png", hMishits[0][0]->GetName(), filename);
  c_1->SaveAs(a);
  */
  /*
  TCanvas* c0 = new TCanvas("c0", "c0", 800, 600);
  hTkiso[1][0]->SetMarkerStyle(20);
  hTkiso[1][1]->SetMarkerStyle(20);
  hTkiso[0][0]->SetFillColor(kRed);
  hTkiso[0][1]->SetFillColor(kRed);
  c0->Divide(2,2);
  c0->GetPad(1)->SetPad(.005, .2525, .495, .995);
  c0->GetPad(2)->SetPad(.505, .2525, .995, .995);
  c0->GetPad(3)->SetPad(.005, .005, .495, .2475);
  c0->GetPad(4)->SetPad(.505, .005, .995, .2475);
  for(int i=0; i<2; i++) {
    c0->cd(i+1);
    c0->GetPad(i+1)->SetLogy(1);
    hTkiso[0][i]->Draw();
    hTkiso[1][i]->Draw("PESAME");
    c0->cd(i+1+2);
    c0->GetPad(i+1+2)->SetGridy(1);
    drawDivision(hTkiso[0][i], hTkiso[1][i]);
  }
  sprintf(a, "%s_%s.png", hTkiso[0][0]->GetName(), filename);
  c0->SaveAs(a);
  */

  //TCanvas* c1 = new TCanvas("c1", "c1", 800, 600);
  //hIso[1][0]->SetMarkerStyle(20);
  //hIso[1][1]->SetMarkerStyle(20);
  //hIso[0][0]->SetFillColor(kRed);
  //hIso[0][1]->SetFillColor(kRed);
  //c1->Divide(2,2);
  //c1->GetPad(1)->SetPad(.005, .2525, .495, .995);
  //c1->GetPad(2)->SetPad(.505, .2525, .995, .995);
  //c1->GetPad(3)->SetPad(.005, .005, .495, .2475);
  //c1->GetPad(4)->SetPad(.505, .005, .995, .2475);
  //for(int i=0; i<2; i++) {
  //  c1->cd(i+1);
  //  hIso[0][i]->Draw();
  //  hIso[1][i]->Draw("PESAME");
  //  c1->cd(i+1+2);
  //  c1->GetPad(i+1+2)->SetGridy(1);
  //  drawDivision(hIso[0][i], hIso[1][i]);
  //}
  //sprintf(a, "%s_%s.png", hIso[0][0]->GetName(), filename);
  //c1->SaveAs(a);
  //
  //
  //TCanvas* c2 = new TCanvas("c2", "c2", 800, 600);
  //hSieie[1][0]->SetMarkerStyle(20);
  //hSieie[1][1]->SetMarkerStyle(20);
  //hSieie[0][0]->SetFillColor(kRed);
  //hSieie[0][1]->SetFillColor(kRed);
  //c2->Divide(2,2);
  //c2->GetPad(1)->SetPad(.005, .2525, .495, .995);
  //c2->GetPad(2)->SetPad(.505, .2525, .995, .995);
  //c2->GetPad(3)->SetPad(.005, .005, .495, .2475);
  //c2->GetPad(4)->SetPad(.505, .005, .995, .2475);
  //for(int i=0; i<2; i++) {
  //  c2->cd(i+1);
  //  hSieie[0][i]->Draw();
  //  hSieie[0][i]->SetNdivisions(505);
  //  hSieie[1][i]->Draw("PESAME");
  //  c2->cd(i+1+2);
  //  c2->GetPad(i+1+2)->SetGridy(1);
  //  drawDivision(hSieie[0][i], hSieie[1][i]);
  //}
  //sprintf(a, "%s_%s.png", hSieie[0][0]->GetName(), filename);
  //c2->SaveAs(a);

  /*
  TCanvas* c3 = new TCanvas("c3", "c3", 800, 600);
    hDphi[1][0]->SetMarkerStyle(20);
  hDphi[1][1]->SetMarkerStyle(20);
  hDphi[0][0]->SetFillColor(kRed);
  hDphi[0][1]->SetFillColor(kRed);
  c3->Divide(2,2);
  c3->GetPad(1)->SetPad(.005, .2525, .495, .995);
  c3->GetPad(2)->SetPad(.505, .2525, .995, .995);
  c3->GetPad(3)->SetPad(.005, .005, .495, .2475);
  c3->GetPad(4)->SetPad(.505, .005, .995, .2475);
  for(int i=0; i<2; i++) {
    c3->cd(i+1);
    hDphi[0][i]->Draw();
    hDphi[0][i]->SetNdivisions(505);
    hDphi[1][i]->Draw("PESAME");
    c3->cd(i+1+2);
    c3->GetPad(i+1+2)->SetGridy(1);
    drawDivision(hDphi[0][i], hDphi[1][i]);
  }
  sprintf(a, "%s_%s.png", hDphi[0][0]->GetName(), filename);
  c3->SaveAs(a);

  TCanvas* c4 = new TCanvas("c4", "c4", 800, 600);
  hDeta[1][0]->SetMarkerStyle(20);
  hDeta[1][1]->SetMarkerStyle(20);
  hDeta[0][0]->SetFillColor(kRed);
  hDeta[0][1]->SetFillColor(kRed);
  c4->Divide(2,2);
  c4->GetPad(1)->SetPad(.005, .2525, .495, .995);
  c4->GetPad(2)->SetPad(.505, .2525, .995, .995);
  c4->GetPad(3)->SetPad(.005, .005, .495, .2475);
  c4->GetPad(4)->SetPad(.505, .005, .995, .2475);
  for(int i=0; i<2; i++) {
    c4->cd(i+1);
    hDeta[0][i]->Draw();
    hDeta[0][i]->SetNdivisions(505);
    hDeta[1][i]->Draw("PESAME");
    c4->cd(i+1+2);
    c4->GetPad(i+1+2)->SetGridy(1);
    drawDivision(hDeta[0][i], hDeta[1][i]);
  }
  sprintf(a, "%s_%s.png", hDeta[0][0]->GetName(), filename);
  c4->SaveAs(a);
  */
  /*
  TCanvas* c5 = new TCanvas("c5", "c5", 800, 600);
  hHoe[1][0]->SetMarkerStyle(20);
  hHoe[1][1]->SetMarkerStyle(20);
  hHoe[0][0]->SetFillColor(kRed);
  hHoe[0][1]->SetFillColor(kRed);
  c5->Divide(2,2);
  c5->GetPad(1)->SetPad(.005, .2525, .495, .995);
  c5->GetPad(2)->SetPad(.505, .2525, .995, .995);
  c5->GetPad(3)->SetPad(.005, .005, .495, .2475);
  c5->GetPad(4)->SetPad(.505, .005, .995, .2475);
  for(int i=0; i<2; i++) {
    c5->cd(i+1);
    c5->GetPad(i+1)->SetLogy(1);
    hHoe[0][i]->Draw();
    hHoe[0][i]->SetNdivisions(505);
    hHoe[1][i]->Draw("PESAME");
    c5->cd(i+1+2);
    c5->GetPad(i+1+2)->SetGridy(1);
    drawDivision(hHoe[0][i], hHoe[1][i]);
  }
  sprintf(a, "%s_%s.png", hHoe[0][0]->GetName(), filename);
  c5->SaveAs(a);
  
  TCanvas* c6 = new TCanvas("c6", "c6", 800, 600);
  hHcal[1][0]->SetMarkerStyle(20);
  hHcal[1][1]->SetMarkerStyle(20);
  hHcal[0][0]->SetFillColor(kRed);
  hHcal[0][1]->SetFillColor(kRed);
  c6->Divide(2,2);
  c6->GetPad(1)->SetPad(.005, .2525, .495, .995);
  c6->GetPad(2)->SetPad(.505, .2525, .995, .995);
  c6->GetPad(3)->SetPad(.005, .005, .495, .2475);
  c6->GetPad(4)->SetPad(.505, .005, .995, .2475);
  for(int i=0; i<2; i++) {
    c6->cd(i+1);
    c6->GetPad(i+1)->SetLogy(1);
    hHcal[0][i]->Draw();
    hHcal[1][i]->Draw("PESAME");
    c6->cd(i+1+2);
    c6->GetPad(i+1+2)->SetGridy(1);
    drawDivision(hHcal[0][i], hHcal[1][i]);
  }
  sprintf(a, "%s_%s.png", hHcal[0][0]->GetName(), filename);
  c6->SaveAs(a);
  */
  /*
  TCanvas* c7 = new TCanvas("c7", "c7", 800, 600);
  hEt[1][0]->SetMarkerStyle(20);
  hEt[1][1]->SetMarkerStyle(20);
  hEt[0][0]->SetFillColor(kRed);
  hEt[0][1]->SetFillColor(kRed);
  c7->Divide(2,2);
  c7->GetPad(1)->SetPad(.005, .2525, .495, .995);
  c7->GetPad(2)->SetPad(.505, .2525, .995, .995);
  c7->GetPad(3)->SetPad(.005, .005, .495, .2475);
  c7->GetPad(4)->SetPad(.505, .005, .995, .2475);
  for(int i=0; i<2; i++) {
    c7->cd(i+1);
    hEt[0][i]->Draw();
    hEt[1][i]->Draw("PESAME");
    c7->cd(i+1+2);
    c7->GetPad(i+1+2)->SetGridy(1);
    drawDivision(hEt[0][i], hEt[1][i]);
  }
  sprintf(a, "%s_%s.png", hEt[0][0]->GetName(), filename);
  c7->SaveAs(a);
  */
  /*
  TCanvas* c8 = new TCanvas("c8", "c8", 800, 600);
  c8->Divide(1,2);
  c8->GetPad(1)->SetPad(.005, .2525, .995, .995);
  c8->GetPad(2)->SetPad(.005, .005, .995, .2475);
  c8->cd(1);
  hEta[0]->Scale(hEta[1]->Integral()/hEta[0]->Integral());
  hEta[1]->SetMarkerStyle(20);
  hEta[0]->SetFillColor(kRed);
  hEta[0]->Draw();
  hEta[1]->Draw("PESAME");
  c8->cd(2);
  c8->GetPad(2)->SetGridy(1);
  drawDivision(hEta[0], hEta[1]);
  sprintf(a, "%s_%s.png", hEta[0]->GetName(), filename);
  c8->SaveAs(a);
  */
  //TCanvas* c9 = new TCanvas("c9", "c9", 800, 600);
  //c9->Divide(2,1);
  //c9->cd(1);
  //hcorr0->Draw("COLZ");
  //c9->cd(2);
  //hcorr1->Draw("COLZ");

  //TCanvas* c10 = new TCanvas("c10", "c10", 800, 600);
  //c10->Divide(2,1);
  //c10->cd(1);
  //hcorr0->ProfileX()->Draw();
  //c10->cd(2);
  //hcorr1->ProfileX()->Draw();
  
  
  //TCanvas* c_2[8];
  ////char a[100];
  //for(int i=0; i<9; i++) {
  //  //if (i != 1)
  //  //  continue;
  //  sprintf(a, "diff%d", i);
  //  c_2[i] = new TCanvas(a, a, 800, 600);
  //  c_2[i]->Divide(2,1);
  //  c_2[i]->cd(1);
  //  diff[i][0]->Draw();
  //  c_2[i]->cd(2);
  //  diff[i][1]->Draw();
  //}
  
  /*
  TCanvas* c_3[2];
  c_3[0] = new TCanvas("c_30", "c_30", 800, 600);
  c_3[0]->Divide(2,1);
  c_3[0]->cd(1);
  diff_vs_pt[0]->Draw();
  c_3[0]->cd(2);
  diff_vs_pt[1]->Draw();

  c_3[1] = new TCanvas("c_31", "c_31", 800, 600);
  c_3[1]->Divide(2,1);
  c_3[1]->cd(1);
  diff_vs_eta[0]->Draw();
  c_3[1]->cd(2);
  diff_vs_eta[1]->Draw();
  */
  /*
  TCanvas* c_4 = new TCanvas("c_4", "c_4", 600, 600);
  hRho[1]->Draw();
  hRho[1]->SetLineColor(kBlue);
  hRho[1]->SetLineStyle(2);
  hRho[0]->Draw("SAME");
  hRho[0]->SetLineColor(kBlue);
  hRho[2]->Draw("SAME");
  hRho[2]->SetLineColor(kRed);
  */
  //TCanvas* c_5 = new TCanvas("c+5", "c_5", 800, 600);
  //c_5->Divide(2, 1);
  //c_5->cd(1);
  //effAreaReco_rho->Draw();
  //effAreaReco_iso->SetLineColor(kBlue);
  //effAreaReco_rho->SetLineColor(kRed);
  //effAreaReco_iso->Draw("SAME");
  //c_5->cd(2);
  //effAreaHlt_rho->Draw();
  //effAreaHlt_iso->SetLineColor(kBlue);
  //effAreaHlt_rho->SetLineColor(kRed);
  //effAreaHlt_iso->Draw("SAME");

  TCanvas* c_5 = new TCanvas("c+5", "c_5", 800, 600);
  effAreaReco_iso->SetLineColor(kBlack);
  effAreaReco_iso->Draw();
  effAreaHlt_iso->SetLineColor(kRed);
  effAreaHlt_iso->Draw("SAME");
  effAreaReco_rho->Draw("SAME");
  effAreaReco_rho->SetMarkerStyle(20);
  effAreaHlt_rho->Draw("SAME");
  effAreaHlt_rho->SetMarkerStyle(20);
  effAreaHlt_rho->SetLineColor(kRed);
  effAreaHlt_rho->SetMarkerColor(kRed);

  output->cd();
  outTree->Write();
  output->Close();
}
