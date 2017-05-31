
/*C/C++ Includes{{{*/
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <iomanip>
/*}}}*/
/*ROOT Includes{{{*/
#include <TSystem.h>
#include <TString.h>
#include <TStyle.h>
#include <Riostream.h>
#include "TObjString.h"
#include <TNamed.h>
#include <TPRegexp.h>
#include <TObjArray.h>
#include <TChain.h>
#include <TMath.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TFile.h>
#include <TROOT.h>
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TCanvas.h>
#include <TDatime.h>
#include <TError.h>
#include <TVirtualFitter.h>
#include <TSQLServer.h>
#include <TSQLResult.h>
#include <TSQLRow.h>
#include <TCut.h>
#include <TMultiGraph.h>
#include <TCutG.h>
/*}}}*/

using namespace std;

//void get_rate(){
void plot(){
    gStyle->SetPaintTextFormat("5.2e");
    const double Lumi = 1.0e33/12.0;
    const double nBcm2 = 1e-33;
    const double one_day = 24*60*60.0;
    TString cut_pip = Form("(weight_hp*%f) * (x>0.08&&x<0.12)", Lumi*nBcm2*one_day);
    TString cut_pim = Form("(weight_hm*%f) * (x>0.08&&x<0.12)", Lumi*nBcm2*one_day );
    TString cut_diff= Form("(weight_hp-weight_hm)*%f * (x>0.08&&x<0.12)", Lumi*nBcm2*one_day);

    //TString cut_pip = Form("(weight_hp+1e-10) * (abs(log10(x) + 2)<0.1&&weight_hp>1e-16) * %f", Lumi*nBcm2);
    //TString cut_pim = Form("(weight_hm+1e-10) * (abs(log10(x) + 2)<0.1&&weight_hp>1e-1 ) * %f", Lumi*nBcm2);
    //TString cut_diff= Form("(weight_hp-weight_hm+1e-10) * (abs(log10(x) + 2)<0.1&&weight_hp>1e-1 ) * %f", Lumi*nBcm2);

    gStyle->SetOptStat(0);
    TChain *T0 = new TChain("T");
    for(int i=600;i<700;i++){
        for(int j=1;j<=4;j++)
            T0->Add(Form("./c12_pion_LO_free_noPt/EIC_A12_pion_10_600_%d_%d.root", j, i));
    }

    TCanvas *c1 = new TCanvas("c1","c1", 800, 600);
    TH2D *h1 = new TH2D("h1","^{12}C(e,e'#pi^{+})X, E_{e}=10GeV, E_{A}=600GeV, Count in 1-day (0.08<x<0.12); z (GeV); log(Q^{2}) (GeV^{2}))", 7, 0.15, 0.85,8, 0.0, 1.6);
    h1->GetXaxis()->CenterTitle(1);
    h1->GetYaxis()->CenterTitle(1);
    gPad->SetLogz(1);
    T0->Draw("log10(Q2):z>>h1",TCut(cut_pip),"colz");
    TH2D* h2 = (TH2D*) h1->Clone();
    T0->Draw("log10(Q2):z>>h2",TCut(cut_pip),"TEXT");
    h1->Draw("colz");
    h2->Draw("TEXT+same");
    c1->Print("c12_pip_Q2_z_A600_xbin_log.pdf");
    c1->Print("c12_pip_Q2_z_A600_xbin_log.png");
  
    TCanvas *c2 = new TCanvas("c2","c2", 800, 600);
    TH2D *h3 = new TH2D("h3","^{12}C(e,e'#pi^{-})X, E_{e}=10GeV, E_{A}=600GeV, Count in 1-day (0.08<x<0.12); z (GeV); log(Q^{2}) (GeV^{2}))", 7, 0.15, 0.85,8, 0.0, 1.6);
    h3->GetXaxis()->CenterTitle(1);
    h3->GetYaxis()->CenterTitle(1);
    gPad->SetLogz(1);
    T0->Draw("log10(Q2):z>>h3",TCut(cut_pim),"colz");
    TH2D* h4 = (TH2D*) h3->Clone();
    T0->Draw("log10(Q2):z>>h4",TCut(cut_pim),"TEXT");
    h3->Draw("colz");
    h4->Draw("TEXT+same");
    c2->Print("c12_pim_Q2_z_A600_xbin_log.pdf");
    c2->Print("c12_pim_Q2_z_A600_xbin_log.png");
   
 
    TCanvas *c3 = new TCanvas("c3","c3", 800, 600);
    TH2D *h5 = new TH2D("h5","^{12}C:e,e'#pi^{-})X, E_{e}=10GeV, E_{A}=600GeV, Count-Diff(#pi^{+})-R(#pi^{-}) in 1-day, (0.08<x<0.12); z (GeV); log(Q^{2}) (GeV^{2}))", 7, 0.15, 0.85,8, 0.0, 1.6);
    h5->GetXaxis()->CenterTitle(1);
    h5->GetYaxis()->CenterTitle(1);
    gPad->SetLogz(1);
    T0->Draw("log10(Q2):z>>h5",TCut(cut_diff),"colz");
    TH2D* h6 = (TH2D*) h5->Clone();
    T0->Draw("log10(Q2):z>>h6",TCut(cut_diff),"TEXT");
    h5->Draw("colz");
    h6->Draw("TEXT+same");
    c3->Print("c12_diff_Q2_z_A600_xbin_log.pdf");
    c3->Print("c12_diff_Q2_z_A600_xbin_log.png");
 
}

