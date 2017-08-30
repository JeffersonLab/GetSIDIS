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

    double z_binmin = 0.2, z_binmax = 0.8;
    const int zbin = 3;
    const double z_cut[4] = {0.2, 0.4, 0.6, 0.8};
    double Q2binmin = 2.0, Q2binmax = 10;
    const int Q2bin=4;
    const double Q2_cut[5] = {2.0, 4.0, 6.0, 8.0, 10.0};

    TString histotitle;
    if (fA == 2) histotitle = "Counts per z and Q2 bin for d(e,e'#pi^{+})X; z ; Q2 (GeV^{2})";
    if (fA == 12) histotitle = "Counts per z and Q2 bin for ^{12}C(e,e'#pi^{+})X; z ; Q2 (GeV^{2})";
    TH2D *histo_counts = new TH2D("histo_counts",histotitle,zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
    histo_counts->GetXaxis()->CenterTitle(1);
    histo_counts->GetYaxis()->CenterTitle(1);
    histo_counts->SetMarkerSize(2);
    histotitle.Clear();
    if (fA == 2) histotitle = "Weighted total count rate d(e,e'#pi^{+})X; z ; Q2 (GeV^{2})";
    if (fA == 12) histotitle = "Weighted total count rate ^{12}C(e,e'#pi^{+})X; z ; Q2 (GeV^{2})";

    TH2D *pionp_cross = new TH2D("pionp_cross ",histotitle,zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
    pionp_cross->GetXaxis()->CenterTitle(1);
    pionp_cross->GetYaxis()->CenterTitle(1);
    pionp_cross->SetMarkerSize(2);
    pionp_cross->Sumw2();
    histotitle.Clear();
    if (fA == 2) histotitle = "Weighted total count rate d(e,e'#pi^{-})X; z ; Q2 (GeV^{2})";
    if (fA == 12) histotitle = "Weighted total count rate ^{12}C(e,e'#pi^{-})X; z ; Q2 (GeV^{2})";

    TH2D *pionm_cross = new TH2D("pionm_cross",histotitle,zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
    pionm_cross->GetXaxis()->CenterTitle(1);
    pionm_cross->GetYaxis()->CenterTitle(1);
    pionm_cross->SetMarkerSize(2);
    pionm_cross->Sumw2();
    histotitle.Clear();
    if (fA == 2) histotitle = "Count rate difference d(e,e'#pi^{-})X; z ; Q2 (GeV^{2})";
    if (fA == 12) histotitle = "Count rate difference ^{12}C(e,e'#pi^{-})X; z ; Q2 (GeV^{2})";

    TH2D *countratediff = new TH2D("countratediff",histotitle,zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
    countratediff->GetXaxis()->CenterTitle(1);
    countratediff->GetYaxis()->CenterTitle(1);
    countratediff->Sumw2();
    countratediff->SetMarkerSize(2);

    for(int i=0;i<zbin;i++) {
      for(int j=0;j<Q2bin;j++) {
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

        TString filename;
        if(fA ==2){
          filename= Form("EIC_A2_pion_10_100_skim_Q2bin%d_zbin%d.root",j+1,i+1);
        }
        if(fA == 12){
          filename= Form("EIC_A12_pion_10_600_skim_Q2bin%d_zbin%d.root",j+1,i+1);
        }

        TFile *input = new TFile(filename,"R");
        TTree *T = (TTree*) input->GetObjectChecked("T","TTree");
    	//Set Branches of Input files
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


        TString cut = Form("(weight_hp)*(Q2>=%f && Q2<%f && z>=%f && z<%f)",Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);
        //counts for weight_hm is the same
        double counts_in_bin = double(T->GetEntries(TCut(cut)));
        cout << "Total number of events " << counts_in_bin << " for Q2 bin " << j << " and zbin " << i << endl;
        histo_counts->SetBinContent(i+1,j+1,counts_in_bin);
        TH1D *htemp = new TH1D("htemp","Q2 dist Piplus", 80,Q2binmin,Q2binmax);
        htemp->GetXaxis()->CenterTitle(1);
        htemp->GetXaxis()->SetTitle("Q2");
        TH1D *htemp2 = new TH1D("htemp2","Q2 dist Piminus", 80,Q2binmin,Q2binmax);
        htemp2->GetXaxis()->CenterTitle(1);
        htemp2->GetXaxis()->SetTitle("Q2");

        TString pipluscut = Form("(weight_hp)*(Q2>=%f && Q2<%f && z>=%f && z<%f)",Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);
        TString piminuscut = Form("(weight_hm)*(Q2>=%f && Q2<%f && z>=%f && z<%f)",Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);
        T->Draw("Q2>>htemp",TCut(pipluscut),"");
        T->Draw("Q2>>htemp2",TCut(piminuscut),"");
      /*  for (int event=0; event<counts_in_bin; event++) {
            if (weight_hp!=0 && Q2>=Q2_cut[j] && Q2<Q2_cut[j+1] && z>=z_cut[i] && z<z_cut[i+1]) {
              htemp->Fill(Q2,weight_hp);
            }
            else if (weight_hm!=0 && Q2>=Q2_cut[j] && Q2<Q2_cut[j+1] && z>=z_cut[i] && z<z_cut[i+1]) {
                htemp2->Fill(Q2,weight_hm);
            }
       }*/
       Double_t histsum = (double) htemp->GetSum();
       Double_t histsum2 = (double) htemp2->GetSum();

       Double_t histerror = sqrt(htemp->GetSumw2()->GetSum());
       Double_t hist2error = sqrt(htemp2->GetSumw2()->GetSum());

       pionp_cross->SetBinContent(i+1,j+1,histsum);
       pionm_cross->SetBinContent(i+1,j+1,histsum2);

       pionp_cross->SetBinError(i+1,j+1,histerror);
       pionm_cross->SetBinError(i+1,j+1,hist2error);
       delete htemp;
       delete htemp2;

      }
    }

    countratediff->Add(pionp_cross,pionm_cross,1,-1);
    TFile *output;
    if (fA == 2) output= new TFile("perbinplots_d.root","RECREATE");
    if (fA == 12) output= new TFile("perbinplots_c12.root","RECREATE");
    output->cd();
    countratediff->Write();
    pionp_cross->Write();
    pionm_cross->Write();
    histo_counts->Write();

    output->Close();


}
