
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
int main(){    
    gStyle->SetPaintTextFormat("5.2e");
    gStyle->SetOptStat(0);
   
    TChain *t1 = new TChain("T");
    for(int i=0;i<10;i++){
        for(int j=1;j<=4;j++)
            t1->Add(Form("./c12_pion/EIC_A12_pion_10_600_%d_%d.root", j, i));
    }

    const double Lumi = 1.0e33/12.0;/*{{{*/
    const double nBcm2 = 1e-33;
    const double one_day = 24*60*60.0;
    TString xbj_cut = "(x>0.008 && x<0.012)";
    //TString xbj_cut = "(x>0.4 && x<0.5)";
    TString cut_dif= Form("(weight_hp-weight_hm)*%f", Lumi*nBcm2*one_day);
    TString cut_sum= Form("(weight_hp+weight_hm)*%f", Lumi*nBcm2*one_day);
    TString cut_pip= Form("(weight_hp)*%f", Lumi*nBcm2*one_day);
    TString cut_pim= Form("(weight_hm)*%f", Lumi*nBcm2*one_day);
    TString cut_inc= Form("(weight_in)*%f", Lumi*nBcm2*one_day);
    /*}}}*/

    /*Histograms{{{*/
    //Q2
    TH1F *h1Q2  =new TH1F("h1Q2","h1Q2",500,0.,100.);
    //W 
    TH1F *h1W = new TH1F("h1W","h1W",500,0.,100.);
    //x 
    TH1F *h1x = new TH1F("h1x","h1x",500,0.005,0.015);
    //pt
    TH1F *h1pt = new TH1F("h1pt","h1pt",500,0.,5.);

    //z
    TH1F *h1z_inc = new TH1F("h1z_inc","h1z_inc",500,0.,1.);
    TH1F *h1z_pip = new TH1F("h1z_pip","h1z_pip",500,0.,1.);
    TH1F *h1z_pim = new TH1F("h1z_pim","h1z_pim",500,0.,1.);
    TH1F *h1z_sum = new TH1F("h1z_sum","h1z_sum",500,0.,1.);
    TH1F *h1z_dif = new TH1F("h1z_dif","h1z_dif",500,0.,1.);
    
    //Sigma
    TH1F *h1XS_pip = new TH1F("h1XS_pip","h1XS_pip (log10)",500,-12.0,6.);
    //Sigma
    TH1F *h1XS_pim = new TH1F("h1XS_pim","h1XS_pim (log10)",500,-12.0,6.);
    //Sigma
    TH1F *h1XS_inc = new TH1F("h1XS_inc","h1XS_inc (log10)",500,-12.0,6.);
    /*}}}*/
    
    /*Binning{{{*/
    TCanvas *c1 = new TCanvas("c1","c1",800,600);
    c1->Divide(4,2);

    TString histoname;
    double N_inc=0.0, N_pip=0.0,N_pim=0.0, N_dif=0.0, N_sum=0.0;
    TString  CUT_pip, CUT_pim, CUT_dif, CUT_sum, CUT_inc;

    ofstream outf("output.dat");

    double z_min = 0.0, z_max = 0.0;
    const int zbin = 8;
    const double z_cut[9] = {0.0, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};
    double z[8], Asym[8], Astat[8], mlpt_pip[8], mlpt_pim[8];
    for (Int_t i=0;i<zbin;i++){
        z_min = z_cut[i];
        z_max = z_cut[i+1];

        /*t Binning{{{*/
        h1Q2->Reset();
        h1W->Reset();
        h1pt->Reset();
        h1x->Reset();
        
        h1z_inc->Reset();
        h1z_pip->Reset();
        h1z_pim->Reset();
        h1z_dif->Reset();
        h1z_sum->Reset();
       
        h1XS_inc->Reset();
        h1XS_pip->Reset();
        h1XS_pim->Reset();

        CUT_inc = cut_inc+"*"+Form("(%s && z>%f && z<%f)", xbj_cut.Data(),z_min,z_max); 
        CUT_pip = cut_pip+"*"+Form("(%s && z>%f && z<%f)", xbj_cut.Data(),z_min,z_max); 
        CUT_pim = cut_pim+"*"+Form("(%s && z>%f && z<%f)", xbj_cut.Data(),z_min,z_max); 
        CUT_sum = cut_sum+"*"+Form("(%s && z>%f && z<%f)", xbj_cut.Data(),z_min,z_max); 
        CUT_dif = cut_dif+"*"+Form("(%s && z>%f && z<%f)", xbj_cut.Data(),z_min,z_max);
        cout<<"--- P+ Cut: "<< CUT_pip.Data()<<endl;
        cout<<"--- P- Cut: "<< CUT_pim.Data()<<endl;
        cout<<"--- P+ & P- Cut: "<< CUT_sum.Data()<<endl;
        cout<<"--- P+ - P- Cut: "<< CUT_dif.Data()<<endl;

        c1->Clear();c1->Divide(4,2);
        c1->cd(1); h1Q2->SetLineColor(1);t1->Draw("Q2>>h1Q2",TCut(CUT_sum));
        c1->cd(2); h1x->SetLineColor(1); t1->Draw("x>>h1x",TCut(CUT_sum));
        c1->cd(3); h1W->SetLineColor(1); t1->Draw("W>>h1W",TCut(CUT_sum));
        c1->cd(4); h1pt->SetLineColor(1);t1->Draw("pt>>h1pt",TCut(CUT_sum));
        
        c1->cd(5); h1z_inc->SetLineColor(1); t1->Draw("z>>h1z_inc",TCut(CUT_inc));
        c1->cd(5); h1z_pip->SetLineColor(4); t1->Draw("z>>h1z_pip",TCut(CUT_pip),"same");
        c1->cd(5); h1z_pim->SetLineColor(6); t1->Draw("z>>h1z_pim",TCut(CUT_pim),"same");
        c1->cd(5); h1z_sum->SetLineColor(2); t1->Draw("z>>h1z_sum",TCut(CUT_sum),"same");
        c1->cd(5); h1z_dif->SetLineColor(7); t1->Draw("z>>h1z_dif",TCut(CUT_dif),"same");
      
        c1->cd(6); h1XS_inc->SetLineColor(1); t1->Draw("log10(dxs_inc)>>h1XS_inc",TCut(CUT_inc));
        c1->cd(7); h1XS_pip->SetLineColor(1); t1->Draw("log10(dxs_hp)>>h1XS_pip",TCut(CUT_pip));
        c1->cd(8); h1XS_pim->SetLineColor(1); t1->Draw("log10(dxs_hm)>>h1XS_pim",TCut(CUT_pim));


        c1->Print(Form("./plot_%d.png",i));

        N_inc = h1z_inc->GetSum();
        N_pip = h1z_pip->GetSum();
        N_pim = h1z_pim->GetSum();
        N_sum = h1z_sum->GetSum();
        N_dif = h1z_dif->GetSum();

        cout<<Form("--- N_sum = %f / %f, N_dif = %f / %f", N_sum, (N_pip+N_pim), N_dif, (N_pip-N_pim))<<endl;
        if(N_sum>1.0){
            Asym[i] = N_dif/N_sum;
            Astat[i] = 1./sqrt(N_sum);
        }else{
            Asym[i] = -1.0;
            Astat[i] = 0.0;
        }

        mlpt_pip[i] = N_pip / N_inc;
        mlpt_pim[i] = N_pim / N_inc;
        
        outf<<Form("%4d %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %10.4f %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e",
                i,
                h1Q2->GetMean(),
                h1x->GetMean(),
                h1W->GetMean(),
                h1z_sum->GetMean(),
                h1pt->GetMean(),
                h1XS_inc->GetMean(),
                h1XS_pip->GetMean(),
                h1XS_pim->GetMean(),
                N_inc,
                N_pip,
                N_pim,
                mlpt_pip[i],
                mlpt_pim[i],
                Asym[i],
                Astat[i]
                )
            <<endl;
        /*}}}*/
    }
    outf.close();

    TCanvas *c2 = new TCanvas("c2","c2",800,600);
    c2->Divide(1,3);

    TH2F *hp = new TH2F("hp","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    hp->SetXTitle("z (GeV)");
    hp->SetYTitle("#pi^{+} Mulplicity");
    
    TH2F *hm = new TH2F("hm","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    hm->SetXTitle("z (GeV)");
    hm->SetYTitle("#pi^{-} Mulplicity");
    
    TH2F *hA = new TH2F("hA","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    hA->SetXTitle("z (GeV)");
    hA->SetYTitle("Charge Asymmetry");

    c2->cd(1); hp->Draw();
    TGraphErrors *gp = new TGraphErrors(zbin, z, mlpt_pip, 0, Astat);
    gp->SetMarkerStyle(20);
    gp->SetMarkerColor(2);

    c2->cd(2); hm->Draw();
    TGraphErrors *gm = new TGraphErrors(zbin, z, mlpt_pim, 0, Astat);
    gm->SetMarkerStyle(20);
    gm->SetMarkerColor(4);
 
    c2->cd(3); hA->Draw();
    TGraphErrors *gA = new TGraphErrors(zbin, z, Asym, 0, Astat);
    gA->SetMarkerStyle(20);
    gA->SetMarkerColor(1);

    c2->Print("asym_mlpt.png");
     
    /*}}}*/
}

