
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
#include <TH1D.h>
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

#include "../generator/SIDIS_Lite_LO.h" //this version doesn't include LHAPDF
using namespace std;

//void get_rate(){
int main(){    
    gStyle->SetPaintTextFormat("5.2e");
    gStyle->SetOptStat(0);

    int imodel = 0; TString MODEL = "";
    cout<<"-- MODEL: 1-> free, 2->EPS09 "; cin>> imodel;
    if(imodel==1) MODEL = "CTEQPDF";
    if(imodel==2) MODEL = "EPS09";

    SIDIS *sidis = new SIDIS(MODEL.Data());
    if(imodel==1) sidis->SetCTEQ(4);
    if(imodel==2) sidis->SetEPS09();

    TString CUT_RANGE="";
    cout<<"-- CUT: tight or wide? "; cin>>CUT_RANGE;

    const int fOrder = 1; //for EPS09, 1->LO, 2->NLO
    const int fErrSet = 1;
    int fA = 2; cout<<"-- A = "; cin >> fA;
    int fZ = 1; //cout<<"-- Z = "; cin >> fZ;
    if( fA==1) fZ = 1; 
    if( fA==2) fZ = 1; 
    if( fA==4) fZ = 2; 
    if( fA==12) fZ = 6; 

    const double Lumi = 1.0e33;/*{{{*/
    const double nBcm2 = 1e-33;
    const double one_day = 24*60*60.0;
    const double xbj_center = 0.01;

    TString cut_dif= Form("(weight_hp-weight_hm)*%6.1f", Lumi*nBcm2*one_day);
    TString cut_sum= Form("(weight_hp+weight_hm)*%6.1f", Lumi*nBcm2*one_day);
    TString cut_pip= Form("(weight_hp)*%6.1f", Lumi*nBcm2*one_day);
    TString cut_pim= Form("(weight_hm)*%6.1f", Lumi*nBcm2*one_day);
    TString cut_inc= Form("(weight_in)*%6.1f", Lumi*nBcm2*one_day);
    /*}}}*/

    //TCanvas *c1 = new TCanvas("c1","c1",800,600);
    //c1->Divide(4,2);

    TString histoname;
    double N_inc=0.0, N_pip=0.0,N_pim=0.0, N_dif=0.0, N_sum=0.0;
    TString  CUT_bin, CUT_unif, CUT_pip, CUT_pim, CUT_dif, CUT_sum, CUT_inc;

    double z_min = 0.0, z_max = 0.0;
    const int zbin = 7;
    const double z_cut[8] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};
    double Q2log_min = 0.0, Q2log_max = 0.0;
    const int Q2bin=8;
    const double Q2_log[9] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};
    //const int Q2bin=4;
    //const double Q2_log[5] = {0.0, 0.4, 0.8, 1.2, 1.6};
    double z[8], Asym[8], Astat[8], mlpt_pip[8], mlpt_pim[8], Ax[8];
    
    /*Binning{{{*/
    for (Int_t j=0;j<Q2bin;j++){
        /*Define{{{*/
        TString new_filename;
 //       TFile *f1 = new TFile(new_filename,"r");
  //      TTree *t1 = (TTree*) f1->Get("T");
       TChain* t1 = new TChain("T");
       if(fA ==12){
            if(imodel==1) new_filename=Form("./c12_pion_LO_free/EIC_A12_pion_10_600_skim%d_wide_free.root",j);
            if(imodel==2) new_filename=Form("./c12_pion_LO/EIC_A12_pion_10_600_skim%d_wide.root",j);
            t1->Add(new_filename.Data());
            cerr<<Form("-- Reading in file: %s", new_filename.Data())<<endl;
       }
       if(fA ==2){
            if(imodel==1) new_filename=Form("./d2_pion_LO_free/EIC_A2_pion_10_100_skim%d_wide_free.root",j);
            if(imodel==2) new_filename=Form("./d2_pion_LO/EIC_A2_pion_10_100_skim%d_wide.root",j);
            t1->Add(new_filename.Data());
            cerr<<Form("-- Reading in file: %s", new_filename.Data())<<endl;
        }

        //cerr<<Form("-- Reading in file: %s", new_filename.Data())<<endl;
        Q2log_min = Q2_log[j];
        Q2log_max = Q2_log[j+1];

        TString Q2_cut = Form(" log10(Q2)>%3.2f&&log10(Q2)<=%3.2f", Q2log_min, Q2log_max);
        double Q2_center = pow(10., 0.5*(Q2log_max+Q2log_min));
        //double Q2_center =0.5*( pow(10., Q2log_min) + pow(10.,Q2log_max));

        TString xbj_cut = "1";
        if(CUT_RANGE=="tight")
            xbj_cut = "(x>0.08 && x<=0.12)";
        else
            xbj_cut = "(x>0.00 && x<=0.3)";

        ofstream outf;/*{{{*/
        if(fA ==12){
            if(CUT_RANGE=="tight"){
                if(imodel==1)    outf.open(Form("c12_output_Q2_%d_lo_tight_free_noPT.dat",j));
                if(imodel==2)    outf.open(Form("c12_output_Q2_%d_lo_tight_noPT.dat",j));
            }else{
                if(imodel==1)     outf.open(Form("c12_output_Q2_%d_lo_wide_free_noPT.dat",j));
                if(imodel==2)     outf.open(Form("c12_output_Q2_%d_lo_wide_noPT.dat",j));
            }
        }
        if(fA ==2){
            if(CUT_RANGE=="tight"){
            if(imodel==1)     outf.open(Form("d2_output_Q2_%d_lo_tight_free_noPT.dat",j));
            if(imodel==2)     outf.open(Form("d2_output_Q2_%d_lo_tight_noPT.dat",j));
            }else{
            if(imodel==1)     outf.open(Form("d2_output_Q2_%d_lo_wide_free_noPT.dat",j));
            if(imodel==2)     outf.open(Form("d2_output_Q2_%d_lo_wide_noPT.dat",j));
            }
        }/*}}}*/

        outf<<Form("%4s %10s %10s %10s %10s %10s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %12s %4s %4s %8s",
                "#bin","Q2","x", "W", "z", "pt", "xs_inc","xs_pip","xs_pim","N_inc","N_pip","N_pim",
                "mulp_pip","mulp_pim","Asym","Astat", "Ax","u_avg","ub_avg","d_avg","db_avg","s_avg","sb_avg","g_avg", "u","ubar","d","dbar","s","sbar","g", "u_c", "ubar_c", "d_c","dbar_c","s_c","sbar_c","g_c","N_mc","A","Z","Model")<<endl;/*}}}*/
    
    /*Histograms{{{ */
    //Q2
    TH1D *h1Q2  =new TH1D("h1Q2","h1Q2",1000,0.,100.);
    //W 
    TH1D *h1W = new TH1D("h1W","h1W",1000,0.,100.);
    //x 
    TH1D *h1x = new TH1D("h1x","h1x",1000,0.005,0.15);
    //pt
    TH1D *h1pt = new TH1D("h1pt","h1pt",1000,0.,5.);

    //Fragmentation Functions
    TH1D *h1D_fav = new TH1D("h1D_fav","h1D_fav",1000,0.,3.);
    TH1D *h1D_unfav = new TH1D("h1D_unfav","h1D_unfav",1000,0.,3.);
    TH1D *h1D_s = new TH1D("h1D_s","h1D_s",1000,0.,3.);
    TH1D *h1D_g = new TH1D("h1D_g","h1D_g",1000,0.,3.);

    //z
    TH1D *h1z_inc = new TH1D("h1z_inc","h1z_inc",1000,0.,1.);
    TH1D *h1z_pip = new TH1D("h1z_pip","h1z_pip",1000,0.,1.);
    TH1D *h1z_pim = new TH1D("h1z_pim","h1z_pim",1000,0.,1.);
    TH1D *h1z_sum = new TH1D("h1z_sum","h1z_sum",1000,0.,1.);
    TH1D *h1z_dif = new TH1D("h1z_dif","h1z_dif",1000,0.,1.);

    TH1D *h1u_avg = new TH1D("h1u_avg","h1u_avg",1000,0.,1.);
    TH1D *h1ubar_avg = new TH1D("h1ubar_avg","h1ubar_avg",1000,0.,1.);
    TH1D *h1d_avg = new TH1D("h1d_avg","h1d_avg",1000,0.,1.);
    TH1D *h1dbar_avg = new TH1D("h1dbar_avg","h1dbar_avg",1000,0.,1.);
    TH1D *h1s_avg = new TH1D("h1s_avg","h1s_avg",1000,0.,1.);
    TH1D *h1sbar_avg = new TH1D("h1sbar_avg","h1sbar_avg",1000,0.,1.);
    TH1D *h1g_avg = new TH1D("h1g_avg","h1g_avg",1000,0.,1.);

    /*    //Sigma*/
    //TH1D *h1XS_pip = new TH1D("h1XS_pip","h1XS_pip (log10)",1000,-12.0,6.);
    ////Sigma
    //TH1D *h1XS_pim = new TH1D("h1XS_pim","h1XS_pim (log10)",1000,-12.0,6.);
    ////Sigma
    //TH1D *h1XS_inc = new TH1D("h1XS_inc","h1XS_inc (log10)",1000,-12.0,6.);


    //Sigma
    TH1D *h1XS_pip = new TH1D("h1XS_pip","h1XS_pip",1000,1e-12,1e6);
    //Sigma
    TH1D *h1XS_pim = new TH1D("h1XS_pim","h1XS_pim",1000,1e-12,1e6);
    //Sigma
    TH1D *h1XS_inc = new TH1D("h1XS_inc","h1XS_inc",1000,1e-12,1e6);
    /*}}}*/

        for (Int_t i=0;i<zbin;i++){
            z_min = z_cut[i];
            z_max = z_cut[i+1];

            /*z Binning{{{*/
            h1Q2->Reset();/*{{{*/
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

            h1D_fav->Reset();
            h1D_unfav->Reset();
            h1D_s->Reset();
            h1D_g->Reset();
            /*}}}*/

            CUT_bin=Form("(%s &&%s && z>%f && z<=%f)", xbj_cut.Data(), Q2_cut.Data(),z_min,z_max); 
            CUT_inc = cut_inc+"*"+Form("(%s)", CUT_bin.Data());
            CUT_pip = cut_pip+"*"+Form("(%s)", CUT_bin.Data());
            CUT_pim = cut_pim+"*"+Form("(%s)", CUT_bin.Data());
            CUT_sum = cut_sum+"*"+Form("(%s)", CUT_bin.Data()); 
            CUT_dif = cut_dif+"*"+Form("(%s)", CUT_bin.Data());
            cout<<"--- P+ Cut: "<< CUT_pip.Data()<<endl;
            cout<<"--- P- Cut: "<< CUT_pim.Data()<<endl;
            cout<<"--- P+ & P- Cut: "<< CUT_sum.Data()<<endl;
            cout<<"--- P+ - P- Cut: "<< CUT_dif.Data()<<endl;

            double N_test = (double) t1->GetEntries(TCut(CUT_pip));
            cerr<<"N_Test = "<<N_test<<endl;
            CUT_unif=Form("(%s &&%s && z>%f && z<=%f)*weight_hp/dxs_hp", xbj_cut.Data(), Q2_cut.Data(),z_min,z_max); 

            t1->Project("h1Q2","Q2",TCut(CUT_pip.Data()));
            t1->Project("h1x","x",TCut(CUT_pip.Data()));
            t1->Project("h1W","W",TCut(CUT_pip.Data()));
            t1->Project("h1pt","pt",TCut(CUT_pip.Data()));

            t1->Project("h1z_inc","z",TCut(CUT_inc));
            t1->Project("h1z_pip","z",TCut(CUT_pip));
            t1->Project("h1z_pim","z",TCut(CUT_pim));
            t1->Project("h1z_sum","z",TCut(CUT_sum));
            t1->Project("h1z_dif","z",TCut(CUT_dif));

            //t1->Project("h1XS_inc","log10(dxs_incl)",TCut(CUT_unif));
            //t1->Project("h1XS_pip","log10(dxs_hp)",TCut(CUT_unif));
            //t1->Project("h1XS_pim","log10(dxs_hm)",TCut(CUT_unif));

            t1->Project("h1XS_inc","(dxs_incl)",TCut(CUT_unif));
            t1->Project("h1XS_pip","(dxs_hp)",TCut(CUT_unif));
            t1->Project("h1XS_pim","(dxs_hm)",TCut(CUT_unif));

            t1->Project("h1u_avg","u_pdf",TCut(CUT_unif));
            t1->Project("h1d_avg","d_pdf",TCut(CUT_unif));
            t1->Project("h1s_avg","s_pdf",TCut(CUT_unif));
            t1->Project("h1g_avg","g_pdf",TCut(CUT_unif));
            t1->Project("h1ubar_avg","ubar_pdf",TCut(CUT_unif));
            t1->Project("h1dbar_avg","dbar_pdf",TCut(CUT_unif));
            t1->Project("h1sbar_avg","sbar_pdf",TCut(CUT_unif));
            
            t1->Project("h1D_fav","D_fav",TCut(CUT_unif));
            t1->Project("h1D_unfav","D_unfav",TCut(CUT_unif));
            t1->Project("h1D_s","D_s",TCut(CUT_unif));
            t1->Project("h1D_g","D_g",TCut(CUT_unif));
            
            TCanvas *c1 = new TCanvas("c1","c1", 800,800);
            c1->Divide(3,2);
            c1->cd(1); h1z_inc->Draw();//t1->Draw("z>>h1z_inc",TCut(CUT_inc));
            c1->cd(2); h1z_pip->Draw();// t1->Draw("z>>h1z_pip",TCut(CUT_pip));
            c1->cd(3); h1z_pim->Draw();// t1->Draw("z>>h1z_pim",TCut(CUT_pim));
            c1->cd(4); h1z_sum->Draw();// t1->Draw("z>>h1z_sum",TCut(CUT_sum));
            c1->cd(5); h1z_dif->Draw();// t1->Draw("z>>h1z_dif",TCut(CUT_dif));/*}}}*/

            if(imodel==1) sidis->RunCTEQPDF( xbj_center, Q2_center);         
            if(imodel==2) sidis->RunEPS09(fOrder, fErrSet, fA, fZ, xbj_center, Q2_center);
            double u_c = sidis->get_uA();
            double ubar_c = sidis->get_ubar();
            double d_c = sidis->get_dA();
            double dbar_c = sidis->get_dbar();
            double s_c= sidis->get_s();
            double sbar_c = sidis->get_sbar();
            double g_c= sidis->get_g();
            
            if(imodel==1) sidis->RunCTEQPDF( h1x->GetMean(), h1Q2->GetMean());         
            if(imodel==2) sidis->RunEPS09(fOrder, fErrSet, fA, fZ, h1x->GetMean(), h1Q2->GetMean());
            double u = sidis->get_uA();
            double ubar = sidis->get_ubar();
            double d = sidis->get_dA();
            double dbar = sidis->get_dbar();
            double s= sidis->get_s();
            double sbar = sidis->get_sbar();
            double g= sidis->get_g();

            double u_avg = h1u_avg->GetMean();
            double ubar_avg = h1ubar_avg->GetMean();
            double d_avg = h1d_avg->GetMean();
            double dbar_avg = h1dbar_avg->GetMean();
            double s_avg = h1s_avg->GetMean();
            double sbar_avg = h1sbar_avg->GetMean();
            double g_avg = h1g_avg->GetMean();

            cout<<Form("u = %e/%e,  d = %e/%e,  ubar = %e/%e, dbar = %e/%e", u, u_avg, d, d_avg, ubar, ubar_avg, dbar, dbar_avg)<<endl;

            Ax[i] = 3./5. * (u-ubar - d+dbar) / (u-ubar + d-dbar);

            //c1->Print(Form("./plot_%d.png",i));

            N_inc =(double) h1z_inc->GetSum();
            N_pip =(double) h1z_pip->GetSum();
            N_pim =(double) h1z_pim->GetSum();
            N_sum =(double) h1z_sum->GetSum();
            N_dif =(double) h1z_dif->GetSum();

            double N_inc1 = t1->GetEntries(CUT_inc);
            double N_pip1 = t1->GetEntries(CUT_pip);
            double N_pim1 = t1->GetEntries(CUT_pim);
            double N_sum1 = t1->GetEntries(CUT_sum);
            double N_dif1 = t1->GetEntries(CUT_dif);

            TString CUT_unom=Form("(%s &&%s && z>%3.2f && z<=%3.2f)", xbj_cut.Data(), Q2_cut.Data(),z_min,z_max); 
            double N_mc = (double) t1->GetEntries(TCut(CUT_unom));

            cout<<Form("--- N_sum = %f / %f, N_dif = %f / %f, N_mc = %f", N_sum, (N_pip+N_pim), N_dif, (N_pip-N_pim), N_mc)<<endl;
            cout<<Form("*** N_sum = %f / %f, N_dif = %f / %f", N_sum1, (N_pip1+N_pim1), N_dif1, (N_pip1-N_pim1))<<endl;
            if(N_sum>1.0){
                Asym[i] = N_dif/N_sum;
                Astat[i] = 1./sqrt(N_sum);
            }else{
                Asym[i] = -1.0;
                Astat[i] = 0.0;
            }

            mlpt_pip[i] = N_pip / N_inc;
            mlpt_pim[i] = N_pim / N_inc;

            outf<<Form("%4d %10.4f %10.4f %10.4f %10.4f %10.4f %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %12.4e %4d %4d %8s",
                    i,                   //0
                    h1Q2->GetMean(),     //1
                    h1x->GetMean(),      //2
                    h1W->GetMean(),      //3
                    h1z_sum->GetMean(),  //4
                    h1pt->GetMean(),     //5
                    //pow(10.,h1XS_inc->GetMean()),
                    //pow(10.,h1XS_pip->GetMean()),
                    //pow(10.,h1XS_pim->GetMean()),
                    (h1XS_inc->GetMean()), //6
                    (h1XS_pip->GetMean()), //7
                    (h1XS_pim->GetMean()), //8
                    N_inc,                 //9
                    N_pip,                 //10
                    N_pim,                 //11
                    mlpt_pip[i],           //12
                    mlpt_pim[i],           //13
                    Asym[i],               //14
                    Astat[i],              //15
                    Ax[i],                 //16
                    u_avg,                 //17
                    ubar_avg,              //18
                    d_avg,                 //19
                    dbar_avg,              //20
                    s_avg,                 //21
                    sbar_avg,              //22
                    g_avg,                 //23
                    u,                     //24
                    ubar,                  //25
                    d,                     //26
                    dbar,                  //27
                    s,                     //28
                    sbar,                  //29
                    g,                     //30
                    u_c,                   //31
                    ubar_c,                //32
                    d_c,                   //33
                    dbar_c,                //34
                    s_c,                   //35
                    sbar_c,                //36
                    g_c,                   //37
                    N_mc,                  //38
                    fA,                     //39
                    fZ,                     //40
                    MODEL.Data()           //41
                        )
                        <<endl;
            /*}}}*/
        }
        outf.close();

        //f1->Close();
    }
    /*}}}Bining*/

}

