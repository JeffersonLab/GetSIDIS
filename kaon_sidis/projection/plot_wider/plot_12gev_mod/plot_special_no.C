#include "plot_s.C"
void plot_special_no(){
  gStyle->SetLabelSize(0.05,"Y");
  gStyle->SetTitleSize(0.05,"Y");
  //gStyle->SetTitleOffset(1.2,"Y");
  gStyle->SetLabelSize(0.05,"X");
  gStyle->SetTitleSize(0.05,"X");
  
  

  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  c1->SetFillColor(10);
  c1->cd(1);
  plot(1,2,3,1,1);
   
  //plot_feng(1,1,1);
  plot_e06010(1,1,3);
  //plot_ase(1,1);
  //plot_pretz(1,1);
  //plot_pretz(1,2);

  TLegend *l1 = new TLegend(0.1,0.6,0.45,0.9);
  Float_t x=0.,y=0.;
  TGraphErrors *g1 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g2 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g3 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g4 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g5 = new TGraphErrors(1,&x,&y,0,0);
  TGraphErrors *g6 = new TGraphErrors(1,&x,&y,0,0);
  g1->SetMarkerColor(1);
  g1->SetMarkerSize(0.7);  
  g1->SetMarkerStyle(21);
  g2->SetLineColor(3);
  g3->SetLineColor(6);
  g3->SetLineWidth(2.);
  g3->SetLineStyle(1);
  g4->SetLineColor(4);
  g4->SetLineWidth(2.5);
  g5->SetLineColor(7);
  g5->SetLineWidth(2.5);
  g6->SetMarkerColor(2);
  g6->SetMarkerStyle(20);
  g6->SetMarkerSize(0.8);
  l1->SetFillColor(10);
  l1->SetShadowColor(0);
  l1->SetLineColor(10);
  l1->AddEntry(g1,"E06010, PRL 107 072003 (2011)","P");
  // l1->AddEntry(g2,"Vogelsang and Yuan (Collins)","L");
  //l1->AddEntry(g3,"Anselmino et al.(Collins)","L");
  //l1->AddEntry(g4,"Pasquini (Collins)","L");
  //  l1->AddEntry(g5,"Pasquini (Pretzelosity)","L");
  l1->AddEntry(g6,"90 days SoLID","P");
  l1->Draw();
  
  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  c2->SetFillColor(10);
  c2->cd(1);
  plot(2,2,3,1,1);
   
  //plot_feng(2,1,1);
  plot_e06010(2,1,3);
  //plot_ase(2,1);
  //plot_pretz(2,1);
  //plot_pretz(2,2);
  l1->Draw();
}
