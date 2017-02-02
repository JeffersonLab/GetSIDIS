
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

#include "SIDIS_Lite.h" //this version doesn't include LHAPDF
using namespace std;

//void get_rate(){
int main(){    
    gStyle->SetPaintTextFormat("5.2e");
    gStyle->SetOptStat(0);
   
    SIDIS *sidis = new SIDIS("EPS09");
    sidis->SetEPS09();
  
    const int fOrder = 1; //for EPS09, 1->LO, 2->NLO
    const int fErrSet = 1;
    int fA = 2; cout<<"-- A = "; cin >> fA;
    int fZ = 1; //cout<<"-- Z = "; cin >> fZ;
    if( fA==1) fZ = 1; 
    if( fA==2) fZ = 1; 
    if( fA==12) fZ = 6; 

     TChain *t1 = new TChain("T");
    for(int i=0;i<10;i++){
        for(int j=1;j<=4;j++){
            if(fA ==12)
                t1->Add(Form("./c12_pion_LO/EIC_A12_pion_10_600_%d_%d.root", j, i));
            if(fA ==2)
                t1->Add(Form("./d2_pion_LO/EIC_A2_pion_10_100_%d_%d.root", j, i));
            if(fA ==1)
                t1->Add(Form("./prot_pion/EIC_A1_pion_10_100_%d_%d.root", j, i));
        }
    }
   
    const double Lumi = 1.0e33/fA;/*{{{*/
    const double nBcm2 = 1e-33;
    const double one_day = 24*60*60.0;
    TString xbj_cut = "(x>0.008 && x<=0.012)";
    //TString xbj_cut = "(x>0.4 && x<0.5)";
    const double xbj_center = 0.01;

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
    
    TH1F *h1u_pip = new TH1F("h1u_pip","h1u_pip",500,0.,1.);
    TH1F *h1ubar_pip = new TH1F("h1ubar_pip","h1ubar_pip",500,0.,1.);
    TH1F *h1d_pip = new TH1F("h1d_pip","h1d_pip",500,0.,1.);
    TH1F *h1dbar_pip = new TH1F("h1dbar_pip","h1dbar_pip",500,0.,1.);
    //TH1F *h1s_pip = new TH1F("h1s_pip","h1s_pip",500,0.,1.);
    //TH1F *h1sbar_pip = new TH1F("h1sbar_pip","h1sbar_pip",500,0.,1.);
    
    TH1F *h1u_pim = new TH1F("h1u_pim","h1u_pim",500,0.,1.);
    TH1F *h1ubar_pim = new TH1F("h1ubar_pim","h1ubar_pim",500,0.,1.);
    TH1F *h1d_pim = new TH1F("h1d_pim","h1d_pim",500,0.,1.);
    TH1F *h1dbar_pim = new TH1F("h1dbar_pim","h1dbar_pim",500,0.,1.);
    //TH1F *h1s_pim = new TH1F("h1s_pim","h1s_pim",500,0.,1.);
    //TH1F *h1sbar_pim = new TH1F("h1sbar_pim","h1sbar_pim",500,0.,1.);
      
    //Sigma
    TH1F *h1XS_pip = new TH1F("h1XS_pip","h1XS_pip (log10)",500,-12.0,6.);
    //Sigma
    TH1F *h1XS_pim = new TH1F("h1XS_pim","h1XS_pim (log10)",500,-12.0,6.);
    //Sigma
    TH1F *h1XS_inc = new TH1F("h1XS_inc","h1XS_inc (log10)",500,-12.0,6.);
    /*}}}*/

    /*Binning{{{*/
    //TCanvas *c1 = new TCanvas("c1","c1",800,600);
    //c1->Divide(4,2);

    TString histoname;
    double N_inc=0.0, N_pip=0.0,N_pim=0.0, N_dif=0.0, N_sum=0.0;
    TString  CUT_bin, CUT_pip, CUT_pim, CUT_dif, CUT_sum, CUT_inc;

    double z_min = 0.0, z_max = 0.0;
    const int zbin = 7;
    const double z_cut[8] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};
    double Q2log_min = 0.0, Q2log_max = 0.0;
    const int Q2bin=8;
    const double Q2_log[9] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
    double z[8], Asym[8], Astat[8], mlpt_pip[8], mlpt_pim[8], Ax[8];
    for (Int_t j=0;j<Q2bin;j++){
        Q2log_min = Q2_log[j];
        Q2log_max = Q2_log[j+1];

        TString Q2_cut = Form(" log10(Q2)>%f&&log10(Q2)<=%f", Q2log_min, Q2log_max);

        double Q2_center = pow(10., 0.5*(Q2log_max+Q2log_min));
  
        ofstream outf;
        if(fA ==12)
            outf.open(Form("c12_output_Q2_%d_lo.dat",j));
        if(fA ==2)
            outf.open(Form("d2_output_Q2_%d_lo.dat",j));
        if(fA ==1)
            outf.open(Form("prot_output_Q2_%d_lo.dat",j));

        outf<<Form("%4s %10s %10s %10s %10s %10s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s",
                "#bin","Q2","x", "W", "z", "pt", "xs_inc","xs_pip","xs_pim","N_inc","N_pip","N_pim",
            "mulp_pip","mulp_pim","Asym","Astat", "Ax", "u", "d", "ubar","dbar","s","sbar","g")<<endl;

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

        CUT_bin=Form("(%s &&%s && z>%f && z<=%f)", xbj_cut.Data(), Q2_cut.Data(),z_min,z_max); 
        CUT_inc = cut_inc+"*"+Form("(%s &&%s && z>%f && z<=%f)", xbj_cut.Data(), Q2_cut.Data(),z_min,z_max); 
        CUT_pip = cut_pip+"*"+Form("(%s &&%s && z>%f && z<=%f)", xbj_cut.Data(), Q2_cut.Data(),z_min,z_max); 
        CUT_pim = cut_pim+"*"+Form("(%s &&%s && z>%f && z<=%f)", xbj_cut.Data(), Q2_cut.Data(),z_min,z_max); 
        CUT_sum = cut_sum+"*"+Form("(%s &&%s && z>%f && z<=%f)", xbj_cut.Data(), Q2_cut.Data(),z_min,z_max); 
        CUT_dif = cut_dif+"*"+Form("(%s &&%s && z>%f && z<=%f)", xbj_cut.Data(), Q2_cut.Data(),z_min,z_max);
        cout<<"--- P+ Cut: "<< CUT_pip.Data()<<endl;
        cout<<"--- P- Cut: "<< CUT_pim.Data()<<endl;
        cout<<"--- P+ & P- Cut: "<< CUT_sum.Data()<<endl;
        cout<<"--- P+ - P- Cut: "<< CUT_dif.Data()<<endl;

        //c1->Clear();c1->Divide(4,2);
        //c1->cd(1); h1Q2->SetLineColor(1);t1->Draw("Q2>>h1Q2",TCut(CUT_pip));
        //c1->cd(2); h1x->SetLineColor(1); t1->Draw("x>>h1x",TCut(CUT_pip));
        //c1->cd(3); h1W->SetLineColor(1); t1->Draw("W>>h1W",TCut(CUT_pip));
        //c1->cd(4); h1pt->SetLineColor(1);t1->Draw("pt>>h1pt",TCut(CUT_pip));
        
        //c1->cd(5); h1z_inc->SetLineColor(1); t1->Draw("z>>h1z_inc",TCut(CUT_inc));
        //c1->cd(5); h1z_pip->SetLineColor(4); t1->Draw("z>>h1z_pip",TCut(CUT_pip),"same");
        //c1->cd(5); h1z_pim->SetLineColor(6); t1->Draw("z>>h1z_pim",TCut(CUT_pim),"same");
        //c1->cd(5); h1z_sum->SetLineColor(2); t1->Draw("z>>h1z_sum",TCut(CUT_sum),"same");
        //c1->cd(5); h1z_dif->SetLineColor(7); t1->Draw("z>>h1z_dif",TCut(CUT_dif),"same");
      
        //c1->cd(6); h1XS_inc->SetLineColor(1); t1->Draw("log10(dxs_incl)>>h1XS_inc",TCut(CUT_inc));
        //c1->cd(7); h1XS_pip->SetLineColor(1); t1->Draw("log10(dxs_hp)>>h1XS_pip",TCut(CUT_pip));
        //c1->cd(8); h1XS_pim->SetLineColor(1); t1->Draw("log10(dxs_hm)>>h1XS_pim",TCut(CUT_pim));

        t1->Project("h1Q2","Q2",TCut(CUT_pip));
        t1->Project("h1x","x",TCut(CUT_pip));
        t1->Project("h1W","W",TCut(CUT_pip));
        t1->Project("h1pt","pt",TCut(CUT_pip));
        
        t1->Project("h1z_inc","z",TCut(CUT_inc));
        t1->Project("h1z_pip","z",TCut(CUT_pip));
        t1->Project("h1z_pim","z",TCut(CUT_pim));
        t1->Project("h1z_sum","z",TCut(CUT_sum));
        t1->Project("h1z_dif","z",TCut(CUT_dif));
        
        t1->Project("h1XS_inc","log10(dxs_incl)",TCut(CUT_bin));
        t1->Project("h1XS_pip","log10(dxs_hp)",TCut(CUT_bin));
        t1->Project("h1XS_pim","log10(dxs_hm)",TCut(CUT_bin));

        t1->Project("h1u_pip","u_pdf",TCut(CUT_pip));
        t1->Project("h1d_pip","d_pdf",TCut(CUT_pip));
        //t1->Project("h1s_pip","s_pdf",TCut(CUT_pip));
        t1->Project("h1ubar_pip","ubar_pdf",TCut(CUT_pip));
        t1->Project("h1dbar_pip","dbar_pdf",TCut(CUT_pip));
        //t1->Project("h1sbar_pip","sbar_pdf",TCut(CUT_pip));
        
        t1->Project("h1u_pim","u_pdf",TCut(CUT_pim));
        t1->Project("h1d_pim","d_pdf",TCut(CUT_pim));
        //t1->Project("h1s_pim","s_pdf",TCut(CUT_pim));
        t1->Project("h1ubar_pim","ubar_pdf",TCut(CUT_pim));
        t1->Project("h1dbar_pim","dbar_pdf",TCut(CUT_pim));
        //t1->Project("h1sbar_pim","sbar_pdf",TCut(CUT_pim));
        
        if(fA==12)
            //sidis->RunEPS09(fOrder, fErrSet, fA, fZ, h1x->GetMean(), h1Q2->GetMean());
            sidis->RunEPS09(fOrder, fErrSet, fA, fZ, xbj_center, Q2_center);
        else
            //sidis->Run_CTEQPDF( h1x->GetMean(), h1Q2->GetMean());         
            sidis->Run_CTEQPDF( xbj_center, Q2_center);         

        double u = sidis->get_uA();
        double d = sidis->get_dA();
        double ubar = sidis->get_ubar();
        double dbar = sidis->get_dbar();
        double s= sidis->get_s();
        double sbar = sidis->get_sbar();
        double g= sidis->get_g();


        double u_pip = h1u_pip->GetMean();
        double d_pip = h1d_pip->GetMean();
        double ubar_pip = h1ubar_pip->GetMean();
        double dbar_pip = h1dbar_pip->GetMean();

        double u_pim = h1u_pim->GetMean();
        double d_pim = h1d_pim->GetMean();
        double ubar_pim = h1ubar_pim->GetMean();
        double dbar_pim = h1dbar_pim->GetMean();


        Ax[i] = 3./5. * (u-ubar + d-dbar) / (u+ubar + d+dbar);

        //c1->Print(Form("./plot_%d.png",i));

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
        
        outf<<Form("%4d %10.4f %10.4f %10.4f %10.4f %10.4f %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e",
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
                Astat[i],
                Ax[i],
                u,
                d,
                ubar,
                dbar,
                s,
                sbar,
                g
                )
            <<endl;
        /*}}}*/
    }
    outf.close();
    }

    //TCanvas *c2 = new TCanvas("c2","c2",800,600);
    //c2->Divide(1,3);

    //TH2F *hp = new TH2F("hp","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    //hp->SetXTitle("z (GeV)");
    //hp->SetYTitle("#pi^{+} Mulplicity");
    
    //TH2F *hm = new TH2F("hm","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    //hm->SetXTitle("z (GeV)");
    //hm->SetYTitle("#pi^{-} Mulplicity");
    
    //TH2F *hA = new TH2F("hA","", 100, 0.0, 1.0, 100, 0.0, 1.0);
    //hA->SetXTitle("z (GeV)");
    //hA->SetYTitle("Charge Asymmetry");

    //c2->cd(1); hp->Draw();
    //TGraphErrors *gp = new TGraphErrors(zbin, z, mlpt_pip, 0, Astat);
    //gp->SetMarkerStyle(20);
    //gp->SetMarkerColor(2);

    //c2->cd(2); hm->Draw();
    //TGraphErrors *gm = new TGraphErrors(zbin, z, mlpt_pim, 0, Astat);
    //gm->SetMarkerStyle(20);
    //gm->SetMarkerColor(4);
 
    //c2->cd(3); hA->Draw();
    //TGraphErrors *gA = new TGraphErrors(zbin, z, Asym, 0, Astat);
    //gA->SetMarkerStyle(20);
    //gA->SetMarkerColor(1);

    //c2->Print("asym_mlpt.png");
     
    /*}}}*/
}

