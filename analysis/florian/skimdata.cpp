/*C/C++ Includes{{{*/
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
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
#include <TLorentzVector.h>
#include <TMath.h>
#include <TClass.h>
#include <TPaletteAxis.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <Rtypes.h>
#include <TTree.h>


using namespace std;
const double DEG = 180./3.1415926;
const double PI = 3.1415926;

int main(Int_t argc, char *argv[]){
  	int fA = 0; cerr<<"-- What nucleus: A = ? "; cin >> fA;

    TChain *T = new TChain("T");
    if(fA ==12){
          T->Add("../massproduction_CTEQfree/EIC_A12_pion_10_600_1_*.root");
    }
    else if(fA ==2) {
          T->Add("../massproduction_CTEQfree/EIC_A2_pion_10_100_1_*.root");
    }
    else { cout << "no valid nucleus" << endl; return 0;}

    //Define
    Double_t Q2, W, Wp, x, y, z, pt, nu, s, gamma, epsilon,rapidity, jacoF;
    Double_t theta_gen= 0.0 , phi_gen = 0.0, mom_gen = 0.0;
    Double_t mom_gen_ele,mom_gen_had;
    Double_t theta_gen_ele,theta_gen_had;
    Double_t phi_gen_ele,phi_gen_had;
    Double_t theta_q, theta_s,phi_h,phi_s,mom_ele,mom_had,theta_ele, theta_had,phi_ele,phi_had;
    Double_t dxs_incl,dxs_hm,dxs_hp,dilute_hp,dilute_hm;
    Double_t px_ele, py_ele,pz_ele, px_had, py_had, pz_had, E_ele,E_had;
    Double_t u_pdf, d_pdf, s_pdf, g_pdf, ubar_pdf, dbar_pdf, sbar_pdf;
    Double_t weight_hp, weight_hm, weight_in;
    ULong64_t nsim = 0;
    //For Beam Position and Vertex info
    Double_t vx_ele, vy_ele, vz_ele, vx_had, vy_had, vz_had;
    Double_t D_fav, D_unfav, D_s, D_g;


	/*Define old root file{{{*/
    T->SetBranchAddress("Q2",&Q2);
    T->SetBranchAddress("W",&W);
    T->SetBranchAddress("Wp",&Wp);
    T->SetBranchAddress("x",&x );
    T->SetBranchAddress("y",&y );
    T->SetBranchAddress("z",&z );
    T->SetBranchAddress("nu",&nu );
    T->SetBranchAddress("s",&s );
    T->SetBranchAddress("epsilon",&epsilon );
    T->SetBranchAddress("gamma",&gamma );
    T->SetBranchAddress("pt",&pt );
    T->SetBranchAddress("weight_hp",&weight_hp );
    T->SetBranchAddress("weight_hm",&weight_hm );
    T->SetBranchAddress("weight_in",&weight_in );
    T->SetBranchAddress("rapidity",&rapidity );
    T->SetBranchAddress("theta_q",&theta_q );
    T->SetBranchAddress("theta_s",&theta_s );
    T->SetBranchAddress("phi_h",&phi_h );
    T->SetBranchAddress("phi_s",&phi_s );
    T->SetBranchAddress("jacoF",&jacoF);
    T->SetBranchAddress("dxs_hm",&dxs_hm);
    T->SetBranchAddress("dxs_hp",&dxs_hp);
    T->SetBranchAddress("dxs_incl",&dxs_incl);
    T->SetBranchAddress("mom_ele",&mom_ele);
    //T->SetBranchAddress("mom_gen_ele",&mom_gen_ele);
    T->SetBranchAddress("mom_had",&mom_had);
    //T->SetBranchAddress("mom_gen_had",&mom_gen_had);
    T->SetBranchAddress("theta_ele",&theta_ele);
    //T->SetBranchAddress("theta_gen_ele",&theta_gen_ele);
    T->SetBranchAddress("theta_had",&theta_had);
    //T->SetBranchAddress("theta_gen_had",&theta_gen_had);
    T->SetBranchAddress("phi_ele",&phi_ele);
    //T->SetBranchAddress("phi_gen_ele",&phi_gen_ele);
    T->SetBranchAddress("phi_had",&phi_had);
    //T->SetBranchAddress("phi_gen_had",&phi_gen_had);
    T->SetBranchAddress("nsim",&nsim);
    T->SetBranchAddress("dilute_p",&dilute_hp );
    T->SetBranchAddress("dilute_m",&dilute_hm  );
    T->SetBranchAddress("px_ele",&px_ele);
    T->SetBranchAddress("py_ele",&py_ele);
    T->SetBranchAddress("pz_ele",&pz_ele);
    T->SetBranchAddress("E_ele",&E_ele);
    T->SetBranchAddress("px_had",&px_had);
    T->SetBranchAddress("py_had",&py_had);
    T->SetBranchAddress("pz_had",&pz_had);
    T->SetBranchAddress("E_had",&E_had);
    T->SetBranchAddress("vx_ele",&vx_ele);
    T->SetBranchAddress("vy_ele",&vy_ele);
    T->SetBranchAddress("vz_ele",&vz_ele);
    T->SetBranchAddress("vx_had",&vx_had);
    T->SetBranchAddress("vy_had",&vy_had);
    T->SetBranchAddress("vz_had",&vz_had);

    T->SetBranchAddress("u_pdf", &u_pdf);
    T->SetBranchAddress("d_pdf", &d_pdf);
    T->SetBranchAddress("s_pdf", &s_pdf);
    T->SetBranchAddress("g_pdf", &g_pdf);
    T->SetBranchAddress("ubar_pdf", &ubar_pdf);
    T->SetBranchAddress("dbar_pdf", &dbar_pdf);
    T->SetBranchAddress("sbar_pdf", &sbar_pdf);

    T->SetBranchAddress("D_g", &D_g);
    T->SetBranchAddress("D_s", &D_s);
    T->SetBranchAddress("D_fav", &D_fav);
    T->SetBranchAddress("D_unfav", &D_unfav);



    ULong64_t N_Total=T->GetEntries();
    cout << "Total number of events " << N_Total << endl;
    N_Total=10000000;
    cout << "MODIFIED NUMBER OF EVENTS TO 10M"<< endl;

//Binning definition
    //double z_min = 0.0, z_max = 1.2;
    //const int zbin = 7;
    //const double z_cut[8] = {0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 1.0};

    double z_binmin = 0.2, z_binmax = 0.8;
    const int zbin = 3;
    const double z_cut[4] = {0.2, 0.4, 0.6, 0.8};

    //double Q2log_min = 0.0, Q2log_max = 1.8;
    //const int Q2bin=8;
    //const double Q2_log[9] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6};

    double Q2binmin = 2.0, Q2binmax = 10;
    const int Q2bin=4;
    const double Q2_cut[5] = {2.0, 4.0, 6.0, 8.0, 10.0};

//xcut definition
    const double xbj_min = 0.05;
    const double xbj_max = 0.1;

    TString new_filename[zbin*Q2bin];
    TFile *outfiles[zbin*Q2bin];
    TTree *outtrees[zbin*Q2bin];
    if(fA ==12){
        for(int i=0;i<zbin;i++) {
          for(int j=0;j<Q2bin;j++) {
             new_filename[i*Q2bin+j]=Form("EIC_A12_pion_10_600_skim_Q2bin%d_zbin%d.root",j+1,i+1);
             outfiles[i*Q2bin+j] = new TFile(new_filename[i*Q2bin+j],"RECREATE");
             outtrees[i*Q2bin+j] = (TTree*) T->GetTree()->CloneTree(0);
             outtrees[i*Q2bin+j]->SetDirectory(outfiles[i*Q2bin+j]);
          }
        }
    }
    if(fA ==2){
      for(int i=0;i<zbin;i++) {
        for(int j=0;j<Q2bin;j++) {
           new_filename[i*Q2bin+j]=Form("EIC_A2_pion_10_100_skim_Q2bin%d_zbin%d.root",j+1,i+1);
           outfiles[i*Q2bin+j] = new TFile(new_filename[i*Q2bin+j],"RECREATE");
           outtrees[i*Q2bin+j] = (TTree*) T->GetTree()->CloneTree(0);
           outtrees[i*Q2bin+j]->SetDirectory(outfiles[i*Q2bin+j]);
        }
      }
    }

   for (ULong64_t i=0;i<N_Total;i++){
      T->GetEntry(i);
      if (x>=xbj_min && x<=xbj_max) {
        for(int i=0;i<zbin;i++) {
          for(int j=0;j<Q2bin;j++) {
            if (Q2>=Q2_cut[j] && Q2<Q2_cut[j+1] && z>=z_cut[i] && z<z_cut[i+1]) {
              outtrees[i*Q2bin+j]->Fill();
            }
          }
        }
      }
        if(!(i%100000))
            cerr<<Form("--- Working on evt=%d",i)<<"\r";
    }

   for (int i=0; i<zbin*Q2bin; i++) {
     outfiles[i]->Write();
     outfiles[i]->Close();
   }

}
