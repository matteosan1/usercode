{

//  gROOT->ProcessLine(".L ~yagil/ana/RootStyle.cc");
  gROOT->ProcessLine(".L ~yagil/ana/rootTools/binerr.C");
  gROOT->ProcessLine(".L ~yagil/ana/rootTools/mk.C");   
  
  gSystem->Load("~yagil/data/eid/checks/MSfuncsJB_C.so");

//  gROOT->ProcessLine("Style = RootStyle()");
  gROOT->ProcessLine("gStyle->SetOptStat(110011)");

//TFile *feid1 = TFile::Open("~yagil/data/eid/checks/Zee_algoC_pixelseed_152.root"); // RED
//TFile *feid1 = TFile::Open("~yagil/data/eid/checks/newN/Zee_algoC_152_new.root"); // RED
//TFile *feid1 = TFile::Open("~yagil/data/eid/hzz/bmtlk/hzz4e_170.root");
TFile *feid1 = TFile::Open("~yagil/data/eid/hzz/xchks/h130C.root");
//TFile *feid1 = TFile::Open("~yagil/data/eid/checks/fN/zee_ootb_newntuple.root"); 
TTree *eid1  = (TTree*)(gROOT->FindObject("event")); 

//TFile *feid2 = TFile::Open("~yagil/data/eid/checks/newN/Zee_algoC_152_new.root");
//TFile *feid2 = TFile::Open("~yagil/data/eid/checks/newN/zee_algoC_152_ootb_new.root"); // BLUE
//TFile *feid2 = TFile::Open("~yagil/data/eid/hzz/bmtlk/hzz4e_170.root");
TFile *feid2 = TFile::Open("~yagil/data/eid/hzz/xchks/h130C.root");
//TFile *feid2 = TFile::Open("~yagil/data/eid/checks/fN/zee_ootb_newntuple.root"); // BLUE
TTree *eid2  = (TTree*)(gROOT->FindObject("event")); 


//Cuts:
TCut denom("sc_et>5");
TCut denom("mc_mother==23&&sc_et>5&&sc_dr<0.1");
//TCut denom("mc_mother==23&&sc_et>5&&fabs(el1_detain) <0.02&&el1_hoe <0.2&&el1_dphiin<0.1&&el1_hoe <0.2&&sc_et>5&&el1_dr<0.1");
//TCut denom("mc_mother==23&&sc_et>5&&fabs(el_detain) <0.02&&el_hoe <0.2&&el_dphiin<0.1&&el_hoe <0.2&&sc_et>5&&el_dr<0.1");


TCut gel  ("el_dr <.1 && fabs(el_detain) <0.02 && el_hoe <0.2 && fabs(el_dphiin)<0.1");
TCut gel1 ("el1_dr<.1 && fabs(el1_detain)<0.02 && el1_hoe<0.2 && fabs(el1_dphiin)<0.1");


TCut avio1  (" ( fabs(el_eta)>1.5 ||   el_see<0.012)  && el_see <0.03");
TCut avin1  (" (fabs(el1_eta)>1.5 ||  el1_see<0.012)  && el1_see<0.03");
TCut avio2  ("fabs(el_detain) <0.012 &&el_hoe <0.1 && el_eopin>.9 &&  sqrt(el_detain**2+el_dphiin**2)<0.1");
TCut avin2  ("fabs(el1_detain)<0.012 &&el1_hoe<0.1 &&el1_eopin>.9 &&sqrt(el1_detain**2+el1_dphiin**2)<0.1");
TCut avio = avio1+avio2;
TCut avin = avin1+avin2;

/*
TCut avio1  (" (fabs(el_eta)>1.5 ||   el_see<0.012)   && (fabs(el_eta)<1.5 || fabs(el_dphiin)<0.05) && fabs(el_detaout)<0.05");
TCut avin1  (" (fabs(el1_eta)>1.5 ||  el1_see<0.012)  && (fabs(el1_eta)<1.5 || fabs(el1_dphiin)<0.05 && fabs(el1_detaout)<0.05) ");
TCut avio2  ("fabs(el_detain) <0.012 &&el_hoe <0.1 && el_eopout>.9 &&  sqrt(el_detain**2+el_dphiin**2)<0.1");
TCut avin2  ("fabs(el1_detain)<0.012 &&el1_hoe<0.1 &&el1_eopout>.9 &&sqrt(el1_detain**2+el1_dphiin**2)<0.1");
TCut avio = avio1+avio2;
TCut avin = avin1+avin2;
*/

TCut jimo ("cut_alajim(el_eopin, el_fbrem, el_see, el_detain, el_dphiin, el_eopout, el_hoe, el_eta) == 1");
TCut jimn ("cut_alajim(el1_eopin, el1_fbrem, el1_see, el1_detain, el1_dphiin, el1_eopout, el1_hoe, el1_eta) == 1");

TCut cmso ("cut_cmssw(0, int(el_class),  el_eopin,  el_detain,  el_dphiin,  el_hoe,  el_e3x3,  el_e5x5,  el_eopout,  el_dphiout,  el_see,  el_spp) == 1");
TCut cmsn ("cut_cmssw(0, int(el1_class), el1_eopin, el1_detain, el1_dphiin, el1_hoe, el1_e3x3, el1_e5x5, el1_eopout, el1_dphiout, el1_see, el1_spp) == 1");

//TCut gel  = avio;
//TCut gel1 = avin;
  //
  //Efficiencies:
  //
 TCanvas *c4 = new TCanvas("c4", " relative eff",  2900,50,800,800);
  c4->Divide(1,2);
gROOT->ProcessLine("gStyle->SetOptStat(0)");
 TH1F *scetn2 = new TH1F("scetn2", "SC Et - num el2  ", 20, 5.0, 105.);
 TH1F *scetn1 = new TH1F("scetn1", "SC Et - num el1 ", 20, 5.0, 105.);
 TH1F *scetd1 = new TH1F("scetd1", "SC Et - den 1",      20, 5.0, 105.);
TH1F *scetd2 = new TH1F("scetd2", "SC Et - den 2",      20, 5.0, 105.);
 
 TH1F *scetan2 = new TH1F("scetan2", "SC Eta - num el2  ", 25, -2.5, 2.5);
 TH1F *scetan1 = new TH1F("scetan1", "SC Eta - num el1 ", 25, -2.5, 2.5);
 TH1F *scetad1 = new TH1F("scetad1",  "SC Eta - den",      25, -2.5, 2.5);
 TH1F *scetad2 = new TH1F("scetad2",  "SC Eta - den",      25, -2.5, 2.5);
 
  //book Efficiency histos:
 TH1F *effptn  = new TH1F("effptn",  "Eff(PT), new",      20, 5.0, 105.);
 TH1F *effpto  = new TH1F("effpto",  "Eff(PT), old",      20, 5.0, 105.);
 TH1F *effetan  = new TH1F("effetan",  "Eff(Eta), new",      25, -2.5, 2.5);
 TH1F *effetao  = new TH1F("effetao",  "Eff(Eta), old",      25, -2.5, 2.5);

  //Et Efficiency plots:
 c4->cd(1); 
	eid2->Draw("sc_et>>scetd2",denom);
	eid1->Draw("sc_et>>scetd1",denom);
	eid2->Draw("sc_et>>scetn2",denom+gel);
	eid1->Draw("sc_et>>scetn1",denom+gel1);
	double effbo = scetn2.Integral(1,20) / scetd2.Integral(1,20) ;
	double effbn = scetn1.Integral(1,20) / scetd1.Integral(1,20) ;
	cout << "eff - old = " << effbo << endl;        
	cout << "eff - new = " << effbn << endl;
	mk(*effptn, *scetn1, *scetd1);         
	mk(*effpto, *scetn2, *scetd2);         
    effpto->GetYaxis()->SetRangeUser(0.5,1.0);           effpto.SetTitle("Eff(PT)"); 
	effpto->SetMarkerColor(4); effpto->SetLineColor(4); effpto->Draw("E");   effptn->SetLineColor(2); effptn->SetMarkerColor(2);effptn->Draw("Esame");
	c4->Update();
	
  //Eta Efficiency plots:
 c4->cd(2); 
	eid2->Draw("sc_eta>>scetad2",denom,"q");
	eid1->Draw("sc_eta>>scetad1",denom,"q");
	eid2->Draw("sc_eta>>scetan2",denom+gel,"q");
	eid1->Draw("sc_eta>>scetan1",denom+gel1);
	mk(*effetan, *scetan1, *scetad1);         
	mk(*effetao, *scetan2, *scetad2);
	effetao->GetYaxis()->SetRangeUser(0.5,1.0);      effetao.SetTitle("Eff(Eta)");    
	effetao->SetMarkerColor(4); effetao->SetLineColor(4); effetao->Draw("E");   effetan->SetLineColor(2); effetan->SetMarkerColor(2);effetan->Draw("Esame");
	c4->Update();


/*
   TCanvas *c4 = new TCanvas("c4", " relative eff",  2900,50,800,800);
  c4->Divide(1,2);
  c4->cd(1); 
  	effptoe->GetYaxis()->SetRangeUser(0.5,1.0);       effptoe.SetTitle("Endcap, Pt");   
	effptoe->SetMarkerColor(4); effptoe->SetLineColor(4); effptoe->Draw("E");   effptne->SetLineColor(2); effptne->SetMarkerColor(2); effptne->Draw("Esame");
c4->cd(2); 
	effetaoe->GetYaxis()->SetRangeUser(0.5,1.0);       effetaoe.SetTitle("Endcap, Eta");   
	effetaoe->SetMarkerColor(4); effetaoe->SetLineColor(4); effetaoe->Draw("E");   effetane->SetLineColor(2); effetane->SetMarkerColor(2);effetane->Draw("Esame");
	c4->Update();
c4->Print("~yagil/data/eid/checks/pics/matteo eg tlk/algCpxles - 152 ootb.gif");    
*/

}