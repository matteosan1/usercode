{


gROOT->ProcessLine(".L ~yagil/ana/rootTools/binerr.C");
gROOT->ProcessLine(".L ~yagil/ana/rootTools/mk.C");   
gROOT->ProcessLine("gStyle->SetOptStat(010011)");


TFile *feid = TFile::Open("~yagil/data/eid/hzz/bmtlk/hzz4e_150.root");
TTree *eid  = (TTree*)(gROOT->FindObject("event")); 

gStyle->SetOptStat("eou");

TCut acc(" fabs(mc_eta[0])<2.5 && fabs(mc_eta[1])<2.5 && fabs(mc_eta[2])<2.5 && fabs(mc_eta[3])<2.5 ");
// simple loose ele Q-cuts
// CMS collection
TCut qel0  (" el_pt[0]>0 && el_eopout[0] >0.5 && fabs(el_detain[0]) <0.02&&el_hoe[0] <0.2 && el_dphiin[0]<0.1");
TCut qel1  (" el_pt[1]>0 && el_eopout[1] >0.5 && fabs(el_detain[1]) <0.02&&el_hoe[1] <0.2 && el_dphiin[1]<0.1");
TCut qel2  (" el_pt[2]>0 && el_eopout[2] >0.5 && fabs(el_detain[2]) <0.02&&el_hoe[2] <0.2 && el_dphiin[2]<0.1");
TCut qel3  (" el_pt[3]>0 && el_eopout[3] >0.5 && fabs(el_detain[3]) <0.02&&el_hoe[3] <0.2 && el_dphiin[3]<0.1");
TCut qel = qel0+qel1+qel2+qel3;
//SNT collection
TCut rel0  ("el1_pt[0]>0 && el1_eopout[0] >0.5 &&fabs(el1_detain[0]) <0.02&&el1_hoe[0] <0.2 && el1_dphiin[0]<0.1");
TCut rel1  ("el1_pt[1]>0 && el1_eopout[1] >0.5 &&fabs(el1_detain[1]) <0.02&&el1_hoe[1] <0.2 && el1_dphiin[1]<0.1");
TCut rel2  ("el1_pt[2]>0 && el1_eopout[2] >0.5 &&fabs(el1_detain[2]) <0.02&&el1_hoe[2] <0.2 && el1_dphiin[2]<0.1");
TCut rel3  ("el1_pt[3]>0 && el1_eopout[3] >0.5 &&fabs(el1_detain[3]) <0.02&&el1_hoe[3] <0.2 && el1_dphiin[3]<0.1");
TCut rel = rel0+rel1+rel2+rel3;


TCanvas *b = new TCanvas("b", "MC info",  50,50,800,800);
b->Divide(1,2);
TH1F* h1 = new TH1F("h1","h1",26,-6,6);
h1->SetTitle("HEPMC electron: Eta ");
TH1F* h2 = new TH1F("h2","h2",20,0.0,100.);
h2->SetTitle("HEPMC electron: PT ");
b->cd(1); eid->Draw("mc_eta>>h1");
b->cd(2); eid->Draw("mc_pt>>h2");
b->cd(1); h1->Draw();
b->cd(2); h2->Draw();
b->Print("~yagil/data/eid/hzz/pics/mckin.gif");    


TCanvas *z = new TCanvas("z", "Z/Z* HEPMC info",  50,50,800,800);
z->Divide(1,2);
TH1F* h6 = new TH1F("h6","h6",50,0,200);
h6->SetTitle("HEPMC M(Z1)");
TH1F* h7 = new TH1F("h7","h7",50,0.0,200.);
h7->SetTitle("HEPMC M(Z2) ");
z->cd(1); eid->Draw("mz1>>h6");
z->cd(2); eid->Draw("mz2>>h7");
z->cd(1); h6->Draw();
z->cd(2); h7->Draw();
z->Print("~yagil/data/eid/hzz/pics/mczz.gif");    

    


TCanvas *b1 = new TCanvas("b1", "SC matching to MC electrons",  50,50,800,800);
b1->Divide(1,2);
TH1F* h3 = new TH1F("h3","h3",50,0.0,0.21);
h3->SetTitle("DeltaR: HEPMC - Reco ECAL S.C. (all)");
TH1F* h4 = new TH1F("h4","h4",50,0.0,0.21);
h4->SetTitle("DeltaR: HEPMC - Reco ECAL S.C. (fid)");
b1->cd(1); eid->Draw("sc_dr>>h3");    
b1->cd(2); eid->Draw("sc_dr>>h4",acc);
b1->cd(1); h3->Draw();
b1->cd(2); h4->Draw();
b1->Print("~yagil/data/eid/hzz/pics/sc.gif");    



TCanvas *y = new TCanvas("y"," Reco - Mz1 Vs. Mz2 ",  50,50,800,800);
y->Divide(1,2);
TH2F* qz12 = new TH2F("qz12","qz12",50,0,200,50,0,200);
qz12->SetTitle("Reco: M(Z1) Vs. M(Z2), loose ele cuts (CMS)");
TH2F* rz12 = new TH2F("rz12","rz12",50,0,200,50,0,200);
rz12->SetTitle("Reco: M(Z1) Vs. M(Z2), loose ele cuts (S&T) ");
y->cd(1); eid->Draw("qmz1:qmz2>>qz12",qel);
y->cd(2); eid->Draw("rmz1:rmz2>>rz12",rel);
y->cd(1); qz12->Draw("box");
y->cd(2); rz12->Draw("box");
y->Print("~yagil/data/eid/hzz/pics/recozz.gif");    



TCanvas *u = new TCanvas("u", "Higgs Reco ",  50,50,800,800);
u->Divide(1,2);
TH1F* q1 = new TH1F("q1","q1",50,0,200);
q1->SetTitle("Reco: M(Z1,Z2), loose ele cuts (CMS)");
TH1F* r1 = new TH1F("r1","r1",50,0,200);
r1->SetTitle("Reco: M(Z1,Z2), loose ele cuts (S&T) ");
u->cd(1); eid->Draw("qmh>>q1",qel);
u->cd(2); eid->Draw("rmh>>r1",rel);
u->cd(1); q1->Draw();
u->cd(2); r1->Draw();
u->Print("~yagil/data/eid/hzz/pics/recomh.gif");    


TCut qpxh  ("el_pixelhits[0]>1 && el_pixelhits[1]>1 && el_pixelhits[2]>1 && el_pixelhits[3]>1 ");
TCut rpxh  ("el1_pixelhits[0]>1 && el1_pixelhits[1]>1 && el1_pixelhits[2]>1 && el1_pixelhits[3]>1 ");
TCut nrpxh  ("el1_pixelhits[0]<2 || el1_pixelhits[1]<2 || el1_pixelhits[2]<2 || el1_pixelhits[3]<2 ");

TCanvas *x = new TCanvas("xx", "H Reco Q-cuts",  50,50,800,800);
xx->Divide(1,2);
TH1F* g8 = new TH1F("g8","g8",50,0,200);
g8->SetTitle("Reco M(Z1,Z2), Q-cuts, <2 pix hits (S&T)");
TH1F* g9 = new TH1F("g9","g9",50,0.0,200.);
g9->SetTitle("Reco M(Z1,Z2), Q-cuts, >=2 pix hts (S&T) ");
xx->cd(1); eid->Draw("rmh>>g8",rel+nrpxh);
xx->cd(2); eid->Draw("rmh>>g9",rel+rpxh);
xx->cd(1); g8->Draw();
xx->cd(2); g9->Draw();
xx->Print("~yagil/data/eid/hzz/pics/mh_pxhts.gif");

TCut qgel0  ("fabs(el_detain[0]) <0.012 && el_hoe[0]<0.1 && (fabs(el_eta[0])>1.5||  sqrt(el_see[0])<0.012) && el_see[0]<0.03 && el_eopout[0]>.9 &&  sqrt(el_detain[0]**2+el_dphiin[0]**2)<0.1");
TCut qgel1  ("fabs(el_detain[1]) <0.012 && el_hoe[1]<0.1 && (fabs(el_eta[1])>1.5||  sqrt(el_see[1])<0.012) && el_see[1]<0.03 && el_eopout[1]>.9 &&  sqrt(el_detain[1]**2+el_dphiin[1]**2)<0.1");
TCut qgel2  ("fabs(el_detain[2]) <0.012 && el_hoe[2]<0.1 && (fabs(el_eta[2])>1.5||  sqrt(el_see[2])<0.012) && el_see[2]<0.03 && el_eopout[2]>.9 &&  sqrt(el_detain[2]**2+el_dphiin[2]**2)<0.1");
TCut qgel3  ("fabs(el_detain[3]) <0.012 && el_hoe[3]<0.1 && (fabs(el_eta[3])>1.5||  sqrt(el_see[3])<0.012) && el_see[3]<0.03 && el_eopout[3]>.9 &&  sqrt(el_detain[3]**2+el_dphiin[3]**2)<0.1");
TCut qgel = qgel0+qgel1+qgel2+qgel3;

TCut rgel0  ("fabs(el1_detain[0])<0.012 &&el1_hoe[0]<0.1 &&(fabs(el1_eta[0])>1.5|| sqrt(el1_see[0])<0.012) &&el1_see[0]<0.03 &&el1_eopout[0]>.9 &&sqrt(el1_detain[0]**2+el1_dphiin[0]**2)<0.1");
TCut rgel1  ("fabs(el1_detain[0])<0.012 &&el1_hoe[0]<0.1 &&(fabs(el1_eta[1])>1.5|| sqrt(el1_see[1])<0.012) &&el1_see[1]<0.03 &&el1_eopout[1]>.9 &&sqrt(el1_detain[1]**2+el1_dphiin[1]**2)<0.1");
TCut rgel2  ("fabs(el1_detain[0])<0.012 &&el1_hoe[0]<0.1 &&(fabs(el1_eta[2])>1.5|| sqrt(el1_see[2])<0.012) &&el1_see[2]<0.03 &&el1_eopout[2]>.9 &&sqrt(el1_detain[2]**2+el1_dphiin[2]**2)<0.1");
TCut rgel3  ("fabs(el1_detain[0])<0.012 &&el1_hoe[0]<0.1 &&(fabs(el1_eta[3])>1.5|| sqrt(el1_see[3])<0.012) &&el1_see[3]<0.03 &&el1_eopout[3]>.9 &&sqrt(el1_detain[3]**2+el1_dphiin[3]**2)<0.1");
TCut rgel = rgel0+rgel1+rgel2+rgel3;

TCanvas *x = new TCanvas("x", "H Reco Q-cuts",  50,50,800,800);
x->Divide(1,2);
TH1F* h8 = new TH1F("h8","h8",50,0,200);
h8->SetTitle("Reco M(Z1,Z2), Q-cuts (CMS)");
TH1F* h9 = new TH1F("h9","h9",50,0.0,200.);
h9->SetTitle("Reco M(Z1,Z2), Q-cuts (S&T) ");
x->cd(1); eid->Draw("qmh>>h8",qgel);
x->cd(2); eid->Draw("rmh>>h9",rgel);
x->cd(1); h8->Draw();
x->cd(2); h9->Draw();
x->Print("~yagil/data/eid/hzz/pics/recomhg.gif");    

}

