#include "TMath.h"
#include "TH1F.h"
#include "TFitter.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"

#include <iostream>

float myVar[4000000];
float myWeight[4000000];
int bins = 200;
TH1F* h, *zee;

int finalEntries = 0;

double chi2(TH1F* h1, TH1F* h2) {

  double likelihood = 0;
  for (int i=0; i<h1->GetNbinsX(); i++) {
    if (h1->GetBinContent(i) != 0) {
      //std::cout << h2->GetBinContent(i) << " " << pow((h1->GetBinContent(i) - h2->GetBinContent(i)),2) << std::endl;
      likelihood += sqrt(pow((h1->GetBinContent(i) - h2->GetBinContent(i)),2)/(h1->GetBinContent(i)));
    }
  }

  //std::cout << likelihood << std::endl;
  return likelihood;
}

double myFunc( double k, double j) {
  zee->Reset("ICESM");
  for(Int_t z=0; z<finalEntries; z++) {
    zee->Fill(myVar[z]*k+j, myWeight[z]);
  }
  //std::cout << h->Integral() << " " << zee->Integral() << std::endl;
  zee->Scale(h->Integral()/zee->Integral());
  
  return chi2(h, zee);
}

void minuitFunction(int& nDim, double* gout, double& result, double par[], int flg) {
  result = myFunc(par[0], par[1]);
}

void migrad(const char* filename = "presel_merged_v2_AB.root", const char* var = "r9", bool barrel=true, Float_t min = 0., Float_t max = 1., Double_t v1=0.85, Double_t xmin=0.80, Double_t xmax=1.0) {

  h = new TH1F("h", "h", bins, min, max);
  zee = new TH1F("zee", "zee", bins, min, max);    

  Int_t type;
  Float_t id;
  Float_t weight, et, eta;
  Float_t r9, r9_1;
  
  char a[100];
  	  
  TFile* file = new TFile(filename);
  TTree* tree = (TTree*)file->Get("opttree");
  tree->SetBranchStatus("*", 0);
  tree->SetBranchStatus("itype", 1);
  tree->SetBranchStatus("et", 1);
  tree->SetBranchStatus("eta", 1);
  tree->SetBranchStatus("weight", 1);
  tree->SetBranchStatus("id", 1);
  tree->SetBranchStatus(var, 1);
  //tree->SetBranchStatus("regr", 1);

  tree->SetBranchAddress("itype", &type);
  tree->SetBranchAddress("et", &et);
  tree->SetBranchAddress("eta", &eta);
  tree->SetBranchAddress("weight", &weight);
  tree->SetBranchAddress("id", &id);
  tree->SetBranchAddress(var, &r9);
  //tree->SetBranchAddress("regr", &r9_1);

  const Int_t entries = tree->GetEntries();

  for(Int_t z=0; z<entries; z++) {
    tree->GetEntry(z);
    if (fabs(eta) < 1.479 && !barrel)
      continue;
    
    if (fabs(eta) > 1.479 && barrel)
      continue;
    
    //std::cout << r9 << " " << r9_1 << std::endl;
    //r9 = r9/r9_1;
    //r9 = r9/et/TMath::CosH(eta);
    //std::cout << r9 << std::endl;
    if (type == 0 && id == 1)
      h->Fill(r9, 1);
    
    if (type == 41 && id == 1) {
      myVar[finalEntries] = r9;
      myWeight[finalEntries] = weight;
      finalEntries++;
    }
  }

  //h->Scale(1/h->Integral());
  
  TFitter* minimizer = new TFitter(2);

  //double p1 = -1;
  //minimizer->ExecuteCommand("SET PRINTOUT", &p1, 1);

  minimizer->SetFCN(minuitFunction);

  //0.983571 0.000143217
  //0.98659 2e-05


  //0.956291 0.000323068

  //0.891832 0.000913388
  //0.856044 0.00123527

  minimizer->SetParameter(0, "k", v1, 0.001, xmin, xmax);
  minimizer->SetParameter(1, "j", 0., 0.01, -1., 1);

  minimizer->ExecuteCommand("SIMPLEX", 0, 0);
  minimizer->ExecuteCommand("MIGRAD", 0, 0);
  
  //double e1, e2, e3, e4;
  //minimizer->GetErrors(0, e1, e2, e3, e4);
  //std::cout << e1 << " " << e2 << " " << e3 << " " << e4 << std::endl;
  //minimizer->mnmigr();
  //double p[2] = {0, 1};
  //minimizer->ExecuteCommand("CONTOUR", &p[0], (int)2);

  //minimizer->ExecuteCommand("HESSE", 0, 0)

  TCanvas* c = new TCanvas("c", "c", 1280,480);
  c->Divide(2,1);
  /*
  c->cd(3);
  //TMinuit* m = minimizer->GetMinuit();
  float p1, err1;
  //minimizer->GetParameter(0, p1, err1);
  //std::cout << p1 << " " << err1 << std::endl;
  minimizer->SetErrorDef(4); //note 4 and not 2!
  TGraph *gr2 = (TGraph*)minimizer->Contour(4,0,1);
  gr2->SetFillColor(42);
  gr2->Draw("alf");
  Get contour for parameter 0 versus parameter 2 for ERRDEF=1  
  minimizer->SetErrorDef(1);
  TGraph *gr1 = (TGraph*)minimizer->Contour(4,0,1);
  gr2->SetFillColor(36);
  gr2->Draw("lf");
  */
  double bestK = minimizer->GetParameter(0);
  double bestJ = minimizer->GetParameter(1);

  std::cout << "Values:" << bestK << " -  " << bestJ << std::endl;
  
  std::cout << myFunc(bestK, bestJ) << std::endl;
  c->cd(2);
  zee->SetFillColor(kRed);
  h->SetMarkerStyle(20);
  if (zee->GetMaximum() > h->GetMaximum()) {
    zee->DrawClone();
    h->DrawClone("PESAME");
  } else {
    h->DrawClone("PE");
    zee->DrawClone("HISTSAME");
    h->DrawClone("PESAME");
  }
  std::cout << myFunc(1, 0) << std::endl;
  c->cd(1);
  zee->SetFillColor(kRed);
  h->SetMarkerStyle(20);
  if (zee->GetMaximum() > h->GetMaximum()) {
    zee->DrawClone();
    h->DrawClone("PESAME");
  } else {
    h->DrawClone("PE");
    zee->DrawClone("HISTSAME");
    h->DrawClone("PESAME");
  }
  
  char name[100];
  sprintf(name, "%s_%d.png", var, int(barrel)); 
  c->SaveAs(name);

}

//R9
//1.00452 0.00103
//1.00861 -0.00071

// sieie 
//0.98659 2e-05
//0.99470 3e-05

//S4
//1.01894 -0.0103405
//1.04969 -0.0364221

//etawidth
//1.04302 -0.000618064
//0.903254 0.00134607

//phiwidth
//1.00002 -0.000371354
//0.999924 4.7808e-07

//effSigmaRR
//1.00101 0.870801


//sigmaE
//1.0265 -3.17466e-05
// new 1.02693 -0.00427938
//1.0326 0.0018304
// new 1.01372 0.000156943

// sieip
//1.05407 -7.83753e-08
//1.01985 9.98054e-08


void migradIso(const char* filename = "output.root", bool barrelEndcap=true, Float_t min = 0., Float_t max = 1., Double_t v1=0.85, Double_t xmin=0.80, Double_t xmax=1.0) {

  h = new TH1F("h", "h", bins, min, max);
  zee = new TH1F("zee", "zee", bins, min, max);    

  Int_t barrel;
  Float_t offline, online;
  
  	  
  TFile* file = new TFile(filename);
  TTree* tree = (TTree*)file->Get("tree");

  tree->SetBranchAddress("barrel", &barrel);
  tree->SetBranchAddress("online", &online);
  tree->SetBranchAddress("offline", &offline);

  const Int_t entries = tree->GetEntries();

  for(Int_t z=0; z<entries; z++) {
    tree->GetEntry(z);
    if (barrel == 0 && barrelEndcap)
      continue;
    
    if (barrel == 1 && !barrelEndcap)
      continue;
    
    h->Fill(online, 1);
    
    myVar[finalEntries] = offline;
    myWeight[finalEntries] = 1;
    finalEntries++;
  }

  //h->Scale(1/h->Integral());
  
  TFitter* minimizer = new TFitter(2);

  //double p1 = -1;
  //minimizer->ExecuteCommand("SET PRINTOUT", &p1, 1);

  minimizer->SetFCN(minuitFunction);

  //0.983571 0.000143217
  //0.98659 2e-05


  //0.956291 0.000323068

  //0.891832 0.000913388
  //0.856044 0.00123527

  minimizer->SetParameter(0, "k", v1, 0.001, xmin, xmax);
  minimizer->SetParameter(1, "j", 0., 0.01, -1., 1);

  minimizer->ExecuteCommand("SIMPLEX", 0, 0);
  minimizer->ExecuteCommand("MIGRAD", 0, 0);
  
  //double e1, e2, e3, e4;
  //minimizer->GetErrors(0, e1, e2, e3, e4);
  //std::cout << e1 << " " << e2 << " " << e3 << " " << e4 << std::endl;
  //minimizer->mnmigr();
  //double p[2] = {0, 1};
  //minimizer->ExecuteCommand("CONTOUR", &p[0], (int)2);

  //minimizer->ExecuteCommand("HESSE", 0, 0)

  TCanvas* c = new TCanvas("c", "c", 1280,480);
  c->Divide(2,1);
  /*
  c->cd(3);
  //TMinuit* m = minimizer->GetMinuit();
  float p1, err1;
  //minimizer->GetParameter(0, p1, err1);
  //std::cout << p1 << " " << err1 << std::endl;
  minimizer->SetErrorDef(4); //note 4 and not 2!
  TGraph *gr2 = (TGraph*)minimizer->Contour(4,0,1);
  gr2->SetFillColor(42);
  gr2->Draw("alf");
  Get contour for parameter 0 versus parameter 2 for ERRDEF=1  
  minimizer->SetErrorDef(1);
  TGraph *gr1 = (TGraph*)minimizer->Contour(4,0,1);
  gr2->SetFillColor(36);
  gr2->Draw("lf");
  */
  double bestK = minimizer->GetParameter(0);
  double bestJ = minimizer->GetParameter(1);

  std::cout << "Values:" << bestK << " /  " << bestJ << std::endl;
  
  std::cout << myFunc(bestK, bestJ) << std::endl;
  c->cd(2);
  h->SetFillColor(kRed);
  zee->SetMarkerStyle(20);
  if (h->GetMaximum() > zee->GetMaximum()) {
    h->DrawClone();
    zee->DrawClone("PESAME");
  } else {
    zee->DrawClone("PE");
    h->DrawClone("HISTSAME");
    zee->DrawClone("PESAME");
  }
  
  std::cout << myFunc(1, 0) << std::endl;
  c->cd(1);
  h->SetFillColor(kRed);
  zee->SetMarkerStyle(20);
  if (h->GetMaximum() > zee->GetMaximum()) {
    h->DrawClone();
    zee->DrawClone("PESAME");
  } else {
    zee->DrawClone("PE");
    h->DrawClone("HISTSAME");
    zee->DrawClone("PESAME");
  }
  
  char name[100];
  sprintf(name, "%s_%d.png", "iso", int(barrel)); 
  c->SaveAs(name);

}
