/*
 *  JB fake.C
 *  
 *
 *  Created by Avi Yagil on 8/15/07.
 *  Copyright 2007 __MyCompanyName__. All rights reserved.
 *
 */

{

  gROOT->ProcessLine(".L ~yagil/ana/rootTools/binerr.C");
  gROOT->ProcessLine(".L ~yagil/ana/rootTools/mk.C");   
  gSystem->Load("~yagil/data/eid/checks/MSfuncsJB_C.so");
  gROOT->ProcessLine("gStyle->SetOptStat(110011)");

TFile *feid = TFile::Open("~yagil/data/eid/hzz/bmtlk/qcd_50_80.root");
TTree *eid  = (TTree*)(gROOT->FindObject("event")); 


//Cuts:
TCut denom("sc_et>5");
TCut gel  ("el_eopin !=0 &&fabs(el_detain) <0.02&&el_hoe <0.2&&el_dphiin<0.1");
TCut gel1 ("el1_eopin!=0 &&fabs(el1_detain)<0.02&&el1_hoe<0.2&&el1_dphiin<0.1");

//TCut denom("fabs(el_detain) <0.02&&el_hoe <0.2&&el_dphiin<0.1&&el_hoe <0.2");

TCut avio1  (" ( fabs(el_eta)>1.5 ||   el_see<0.012)  && el_see <0.03");
TCut avin1  (" (fabs(el1_eta)>1.5 ||  el1_see<0.012)  && el1_see<0.03");
TCut avio2  ("fabs(el_detain) <0.012 &&el_hoe <0.1 && el_eopin>.9 &&  sqrt(el_detain**2+el_dphiin**2)<0.1");
TCut avin2  ("fabs(el1_detain)<0.012 &&el1_hoe<0.1 &&el1_eopin>.9 &&sqrt(el1_detain**2+el1_dphiin**2)<0.1");
TCut avio = avio1+avio2;
TCut avin = avin1+avin2;

TCut jimo ("cut_alajim(el_eopin, el_fbrem, el_see, el_detain, el_dphiin, el_eopout, el_hoe, el_eta) == 1");
TCut jimn ("cut_alajim(el1_eopin, el1_fbrem, el1_see, el1_detain, el1_dphiin, el1_eopout, el1_hoe, el1_eta) == 1");

TCut cmso ("cut_cmssw(0, int(el_class),  el_eopin,  el_detain,  el_dphiin,  el_hoe,  el_e3x3,  el_e5x5,  el_eopout,  el_dphiout,  el_see,  el_spp) == 1");
TCut cmsn ("cut_cmssw(0, int(el1_class), el1_eopin, el1_detain, el1_dphiin, el1_hoe, el1_e3x3, el1_e5x5, el1_eopout, el1_dphiout, el1_see, el1_spp) == 1");

TCut gel  = avio;
TCut gel1 = avin;

  //
  //Efficiencies:
  //
 TCanvas *c4 = new TCanvas("c4", " relative eff",  2900,50,800,800);
 c4->Divide(2,2);
//gROOT->ProcessLine("gStyle->SetOptStat(0)");
 TH1F *scetn0 = new TH1F("scetn0", "SC Et - num el  ", 20, 5.0, 105.);
 TH1F *scetn1 = new TH1F("scetn1", "SC Et - num el1 ", 20, 5.0, 105.);
 TH1F *scetd  = new TH1F("scetd",  "SC Et - den",      20, 5.0, 105.);
 TH1F *scetan0 = new TH1F("scetan0", "SC Eta - num el  ", 50, -2.5, 2.5);
 TH1F *scetan1 = new TH1F("scetan1", "SC Eta - num el1 ", 50, -2.5, 2.5);
 TH1F *scetad  = new TH1F("scetad",  "SC Eta - den",      50, -2.5, 2.5);
  //book Efficiency histos:
 TH1F *effptnb  = new TH1F("effptnb",  "Eff. Barrel, new",      20, 5.0, 105.);
 TH1F *effptob  = new TH1F("effptob",  "Eff. Barrel, old",      20, 5.0, 105.);
 TH1F *effetanb  = new TH1F("effetanb",  "Eff. Barrel, new",      50, -2.5, 2.5);
 TH1F *effetaob  = new TH1F("effetaob",  "Eff. Barrel, old",      50, -2.5, 2.5);
 TH1F *effptne  = new TH1F("effptne",  "Eff. EndCap, new",      20, 5.0, 105.);
 TH1F *effptoe  = new TH1F("effptoe",  "Eff. EndCap, new",      20, 5.0, 105.);
 TH1F *effetane  = new TH1F("effetane",  "Eff. EndCap, new",      50, -2.5, 2.5);
 TH1F *effetaoe  = new TH1F("effetaoe",  "Eff. EndCap, new",      50, -2.5, 2.5);

  //Et Efficiency plots:
 c4->cd(1); 
eid->Draw("sc_et>>scetd",denom);
eid->Draw("sc_et>>scetn0",denom+gel);
eid->Draw("sc_et>>scetn1",denom+gel1);
c4->cd(1); gPad->SetLogy();
scetd->Draw(); scetn0->SetLineColor(4); scetn0->Draw("sames"); scetn1->SetLineColor(2); scetn1->Draw("sames"); 
double effbo = scetn0.Integral(1,20) / scetd.Integral(1,20) ;
double effbn = scetn1.Integral(1,20) / scetd.Integral(1,20) ;
cout << "Fake Rate old = " << effbo << endl;        
cout << "Fake Rate new = " << effbn << endl;
mk(*effptnb, *scetn1, *scetd);         
mk(*effptob, *scetn0, *scetd);         
c4->cd(2); 
effptob->GetYaxis()->SetRangeUser(0.,0.25);           effptob.SetTitle("SC Et - eff"); 
effptob->SetLineColor(4); effptob->Draw("E");   effptnb->SetLineColor(2); effptnb->Draw("Esame");
c4->Update();

  //Eta Efficiency plots:
 c4->cd(3); gPad->SetLogy(); 
scetad->Draw(); scetan0->SetLineColor(4); scetan0->Draw("same"); scetan1->SetLineColor(2); scetan1->Draw("same"); 
eid->Draw("sc_eta>>scetad",denom);
eid->Draw("sc_eta>>scetan0",denom+gel);
eid->Draw("sc_eta>>scetan1",denom+gel1);
     c4->cd(3); gPad->SetLogy();scetad->SetMinimum(0.5);
scetad->Draw(); scetan0->SetLineColor(4); scetan0->Draw("sames"); scetan1->SetLineColor(2); scetan1->Draw("sames"); 
mk(*effetanb, *scetan1, *scetad);         
mk(*effetaob, *scetan0, *scetad);
 c4->cd(4); 
effetaob->GetYaxis()->SetRangeUser(0.,0.1);      effetaob.SetTitle("SC eta - eff");    
effetaob->SetLineColor(4); effetaob->Draw("E");   effetanb->SetLineColor(2); effetanb->Draw("Esame");
c4->Update();
  gROOT->ProcessLine("gStyle->SetOptStat(110011)");
c4->Print("~yagil/data/eid/checks/pics/fakerate.gif");   


}