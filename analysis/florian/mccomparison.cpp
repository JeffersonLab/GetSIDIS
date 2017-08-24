#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
/*ROOT Includes*/
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
#include <TH2D.h>

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

#include "../../generator/SIDIS_Lite_LO.h" //this version doesn't include LHAPDF

int main()
{
  gROOT->Reset();
  gStyle->SetOptStat(0);
  //gStyle->SetTextFont(1);
  gStyle->SetPaintTextFormat("4.3f");
  //gROOT->SetStyle("BABAR");
  TChain *T0 = new TChain("T");


  const Double_t EIC_Mom_Min_e = 8.5;
  const Double_t EIC_Mom_Max_e = 10.5; //not in use
  const Double_t EIC_Mom_Min_h = 0.0;
  const Double_t EIC_Mom_Max_h = 10.0;

  const Double_t EIC_Th_Min_e = 0.0;
  const Double_t EIC_Th_Max_e  = 25.0;
  const Double_t EIC_Th_Min_h = 0.0;
  const Double_t EIC_Th_Max_h  = 180.0;

  const Double_t EIC_Ph_Min_e = 0.0;
  const Double_t EIC_Ph_Max_e  = 360.0;
  const Double_t EIC_Ph_Min_h = 0.0;
  const Double_t EIC_Ph_Max_h  = 360.0;

  const Double_t degtorad = TMath::Pi()/180.;

  int imodel = 0; TString MODEL = "";
  cout<<"-- MODEL: 1-> free, 2->EPS09 "; //cin>> imodel;
  imodel = 1;
  if(imodel==1) MODEL = "CTEQPDF";
  if(imodel==2) MODEL = "EPS09";

  SIDIS *sidis = new SIDIS(MODEL.Data());
  if(imodel==1) sidis->SetCTEQ(4);
  if(imodel==2) sidis->SetEPS09();

  const int fOrder = 1; //for EPS09, 1->LO, 2->NLO
  const int fErrSet = 1;
//  int fA = 2; cout<<"-- A = "; cin >> fA;
int fA = 12;
  int fZ = 1; //cout<<"-- Z = "; cin >> fZ;
  if( fA==1) fZ = 1;
  if( fA==2) fZ = 1;
  if( fA==4) fZ = 2;
  if( fA==12) fZ = 6;

  //T0->Add(Form("./c12_pion/EIC_A12_pion_10_600_%d_%d.root", j, i));
  //EIC:  1 -> Q2<10, pt<1
  //               2 -> Q2<10, pt>1
  //           EIC:  3 -> Q2>10, pt<1
  //                 4 -> Q2>10, pt>1
 TString outputfileending = "-";

 Double_t histo_ptmin = 0;
 Double_t histo_ptmax = 5;

 cout << "File with Q2<10 and pt<1 is processed" << endl;
 T0->Add("../../test/massproduction_CTEQ/EIC_A12_pion_10_600_1_0.root");

 outputfileending += "f1-" ;
 histo_ptmax = 1.1;

 cout << "Total events " << T0->GetEntries() << endl;

//Choice of cut for all histograms


TString xbj_cut_small = "(x>0.05 && x<=0.1)";
TString weightcut = "(weight_hp>0)";
TString weighting = "(weight_hp)";
TString allcut;
TString xandnoweighting = xbj_cut_small+"*"+weightcut;

  allcut = xbj_cut_small+"*"+weighting;
  outputfileending += "c_xsm_w1" ;


 //allcut += "*(1000000/nsim)";
 outputfileending += ".pdf" ;
 cout << "file name extension: " << outputfileending << endl;

 //**Define**
 TString savestring; //to store file name of saved canvas
  //Definitions of Tree values
 Double_t Q2, W, Wp, x, y, z, pt, nu, s, epsilon,rapidity, jacoF;
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

 T0->SetBranchAddress("Q2",&Q2);
 T0->SetBranchAddress("W",&W);
 T0->SetBranchAddress("Wp",&Wp);
 T0->SetBranchAddress("x",&x );
 T0->SetBranchAddress("y",&y );
 T0->SetBranchAddress("z",&z );
 T0->SetBranchAddress("nu",&nu );
 T0->SetBranchAddress("s",&s );
 T0->SetBranchAddress("epsilon",&epsilon );
 T0->SetBranchAddress("pt",&pt );
 T0->SetBranchAddress("weight_hp",&weight_hp );
 T0->SetBranchAddress("weight_hm",&weight_hm );
 T0->SetBranchAddress("weight_in",&weight_in );
 T0->SetBranchAddress("rapidity",&rapidity );
 T0->SetBranchAddress("theta_q",&theta_q );
 T0->SetBranchAddress("theta_s",&theta_s );
 T0->SetBranchAddress("phi_h",&phi_h );
 T0->SetBranchAddress("phi_s",&phi_s );
 T0->SetBranchAddress("jacoF",&jacoF);
 T0->SetBranchAddress("dxs_hm",&dxs_hm);
 T0->SetBranchAddress("dxs_hp",&dxs_hp);
 T0->SetBranchAddress("dxs_incl",&dxs_incl);
 T0->SetBranchAddress("mom_ele",&mom_ele);
 //T0->SetBranchAddress("mom_gen_ele",&mom_gen_ele);
 T0->SetBranchAddress("mom_had",&mom_had);
 //T0->SetBranchAddress("mom_gen_had",&mom_gen_had);
 T0->SetBranchAddress("theta_ele",&theta_ele);
 //T0->SetBranchAddress("theta_gen_ele",&theta_gen_ele);
 T0->SetBranchAddress("theta_had",&theta_had);
 //T0->SetBranchAddress("theta_gen_had",&theta_gen_had);
 T0->SetBranchAddress("phi_ele",&phi_ele);
 //T0->SetBranchAddress("phi_gen_ele",&phi_gen_ele);
 T0->SetBranchAddress("phi_had",&phi_had);
 //T0->SetBranchAddress("phi_gen_had",&phi_gen_had);
 T0->SetBranchAddress("nsim",&nsim);
 T0->SetBranchAddress("dilute_p",&dilute_hp );
 T0->SetBranchAddress("dilute_m",&dilute_hm  );
 T0->SetBranchAddress("px_ele",&px_ele);
 T0->SetBranchAddress("py_ele",&py_ele);
 T0->SetBranchAddress("pz_ele",&pz_ele);
 T0->SetBranchAddress("E_ele",&E_ele);
 T0->SetBranchAddress("px_had",&px_had);
 T0->SetBranchAddress("py_had",&py_had);
 T0->SetBranchAddress("pz_had",&pz_had);
 T0->SetBranchAddress("E_had",&E_had);
 T0->SetBranchAddress("vx_ele",&vx_ele);
 T0->SetBranchAddress("vy_ele",&vy_ele);
 T0->SetBranchAddress("vz_ele",&vz_ele);
 T0->SetBranchAddress("vx_had",&vx_had);
 T0->SetBranchAddress("vy_had",&vy_had);
 T0->SetBranchAddress("vz_had",&vz_had);

 T0->SetBranchAddress("u_pdf", &u_pdf);
 T0->SetBranchAddress("d_pdf", &d_pdf);
 T0->SetBranchAddress("s_pdf", &s_pdf);
 T0->SetBranchAddress("g_pdf", &g_pdf);
 T0->SetBranchAddress("ubar_pdf", &ubar_pdf);
 T0->SetBranchAddress("dbar_pdf", &dbar_pdf);
 T0->SetBranchAddress("sbar_pdf", &sbar_pdf);

 T0->SetBranchAddress("D_g", &D_g);
 T0->SetBranchAddress("D_s", &D_s);
 T0->SetBranchAddress("D_fav", &D_fav);
 T0->SetBranchAddress("D_unfav", &D_unfav);
  int  number_of_hminusevents = 0;
  int  number_of_hplusevents = 0;
  for (int i = 0; i<T0->GetEntries() ; i++)
  {
    T0->GetEntry(i);
    if (x>0.05 && x<=0.1){
      if (weight_hp > 0){
        number_of_hplusevents++;
      }
      if (weight_hm > 0) {
        number_of_hminusevents++;
      }

    }
  }
  cout << "For (x>0.05 && x<=0.1) and weight > 0: pos hadrons events: " << number_of_hplusevents << " and neg. hadron events: " << number_of_hminusevents  << endl;
  TH1D *h1x = new TH1D("h1x","x distribution, E_{e}=10GeV, E_{A}=600GeV; x", 100,0.03,0.1);
  h1x->GetXaxis()->CenterTitle(1);
  T0->Draw("x>>h1x",TCut(allcut),"");



  Double_t nofweightcut_nozbinQ2bincuts= (double) T0->GetEntries(TCut(xandnoweighting)) ;
  Double_t nofallcut_nozbinQ2bincuts= (double) T0->GetEntries(TCut(allcut)) ;
  cout << "Compare with TCut option weight > 0 " << nofweightcut_nozbinQ2bincuts << " and weighting " << nofallcut_nozbinQ2bincuts << endl;

 Double_t generated_electron_phase_space =(cos(EIC_Th_Min_e*degtorad) - cos(EIC_Th_Max_e*degtorad))*(EIC_Ph_Max_e*degtorad - EIC_Ph_Min_e*degtorad)*(EIC_Mom_Max_e - EIC_Mom_Min_e);
 Double_t generated_hadron_phase_space   =(cos(EIC_Th_Min_h*degtorad) - cos(EIC_Th_Max_h*degtorad))*(EIC_Ph_Max_h*degtorad - EIC_Ph_Min_h*degtorad)*(EIC_Mom_Max_h - EIC_Mom_Min_h);
 Double_t generated_phase_space=generated_electron_phase_space*generated_hadron_phase_space;

 cout << "generated PS: " << generated_phase_space << endl;
 Double_t vnratio = generated_phase_space/nsim; //Ratio PS Volume and generated Events (see Charles' Note)
 cout << "V/N ratio " << vnratio << endl;


//Further procedure for pt<1 and Q2<10, xcut -> file choice 1 and cutchoice either 2, 3, 12 or 1



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


  TH2D *c12pspart = new TH2D("c12pspart","Counts per z and Q2 bin for ^{12}C(e,e'#pi^{+})X; z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
  c12pspart->GetXaxis()->CenterTitle(1);
  c12pspart->GetYaxis()->CenterTitle(1);
  c12pspart->SetMarkerSize(2);
  TH2D *c12pspartrel = new TH2D("c12pspartrel","Relative Counts per z and Q2 bin in percent for ^{12}C(e,e'#pi^{+})X; z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
  c12pspartrel->GetXaxis()->CenterTitle(1);
  c12pspartrel->GetYaxis()->CenterTitle(1);
  c12pspartrel->SetMarkerSize(2);

  TH2D *c12pionp_cross = new TH2D("c12pionp_cross ","Integral weighted cross section with correct errors for ^{12}C(e,e'#pi^{+})X;  z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
  c12pionp_cross->GetXaxis()->CenterTitle(1);
  c12pionp_cross->GetYaxis()->CenterTitle(1);
  c12pionp_cross->SetMarkerSize(2);
  TH2D *c12pionp_crossw = new TH2D("c12pionp_crossw ","Integral weighted cross section with wrong errors for ^{12}C(e,e'#pi^{+})X;  z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
  c12pionp_crossw->GetXaxis()->CenterTitle(1);
  c12pionp_crossw->GetYaxis()->CenterTitle(1);
  c12pionp_crossw->SetMarkerSize(2);
  TH2D *c12pionm_cross = new TH2D("c12pionm_cross ","Integral weighted cross section with correct errors for ^{12}C(e,e'#pi^{-})X;  z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
  c12pionm_cross->GetXaxis()->CenterTitle(1);
  c12pionm_cross->GetYaxis()->CenterTitle(1);
  c12pionm_cross->SetMarkerSize(2);
  TH2D *c12pionm_crossw = new TH2D("c12pionm_crossw ","Integral weighted cross section with wrong errors for ^{12}C(e,e'#pi^{-})X;  z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
  c12pionm_crossw->GetXaxis()->CenterTitle(1);
  c12pionm_crossw->GetYaxis()->CenterTitle(1);
  c12pionm_crossw->SetMarkerSize(2);
  TH2D *c12nofdiff = new TH2D("c12nofdiff ","Count rate difference with correct errors for ^{12}C #pi^{+}-#pi^{-};  z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
  c12nofdiff->GetXaxis()->CenterTitle(1);
  c12nofdiff->GetYaxis()->CenterTitle(1);
  c12nofdiff->Sumw2();
  c12nofdiff->SetMarkerSize(2);
  TH2D *c12nofdiffw = new TH2D("c12nofdiffw ","Count rate difference with wrong errors for ^{12}C #pi^{+}-#pi^{-};  z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
  c12nofdiffw->GetXaxis()->CenterTitle(1);
  c12nofdiffw->GetYaxis()->CenterTitle(1);
  c12nofdiffw->Sumw2();
  c12nofdiffw->SetMarkerSize(2);

  TCanvas *canv_weighthist = new TCanvas("canv_weighthist"," ", 1200, 800);
//  canv_weighthist->Divide(zbin,Q2bin);
  canv_weighthist->Divide(2,2);
  double man_sumcross_p[zbin][Q2bin];
  double man_sumcross2_p[zbin][Q2bin];
  double man_errorintegral_p[zbin][Q2bin];
  double man_sumcross_m[zbin][Q2bin];
  double man_sumcross2_m[zbin][Q2bin];
  double man_errorintegral_m[zbin][Q2bin];
  for (Int_t j=0;j<Q2bin;j++){
 //  for (Int_t j=1;j<2;j++){//Test for one Q2 bin

     for (Int_t i=0;i<zbin;i++){
   //    for (Int_t i=2;i<3;i++){
       double weightsum_p = 0;
       man_sumcross_p[i][j] = 0;
       man_sumcross2_p[i][j] = 0;
       double weightsum_m = 0;
       man_sumcross_m[i][j] = 0;
       man_sumcross2_m[i][j] = 0;
       for (int nentry = 0; nentry<T0->GetEntries() ; nentry++)
       {
         T0->GetEntry(nentry);
         if (x>0.05 && x<=0.1 && weight_hp>0 && Q2>=Q2_cut[j] && Q2<Q2_cut[j+1] && z>=z_cut[i] && z<z_cut[i+1]){
           weightsum_p+=weight_hp;
           man_sumcross_p[i][j]+=dxs_hp;
           man_sumcross2_p[i][j]+=(dxs_hp*dxs_hp);
         }
         if (x>0.05 && x<=0.1 && weight_hm>0 && Q2>=Q2_cut[j] && Q2<Q2_cut[j+1] && z>=z_cut[i] && z<z_cut[i+1]){
           weightsum_m+=weight_hm;
           man_sumcross_m[i][j]+=dxs_hm;
           man_sumcross2_m[i][j]+=(dxs_hm*dxs_hm);
         }
       }
       man_errorintegral_p[i][j] = generated_phase_space/sqrt(nsim) * sqrt ( man_sumcross2_p[i][j]/nsim - pow (man_sumcross_p[i][j]/nsim, 2) );
       man_errorintegral_m[i][j] = generated_phase_space/sqrt(nsim) * sqrt ( man_sumcross2_m[i][j]/nsim - pow (man_sumcross_m[i][j]/nsim, 2) );

       double z_center = 0.5*(z_cut[i]+z_cut[i+1]);
       double Q2_center = pow(10., 0.5*(Q2_cut[j]+Q2_cut[j+1]));
       double xbj_center = (0.05+0.1)*0.5;

       if(imodel==1) sidis->RunCTEQPDF( xbj_center, Q2_center);
       if(imodel==2) sidis->RunEPS09(fOrder, fErrSet, fA, fZ, xbj_center, Q2_center);
     //  double u_c = sidis->get_uA();
   //    double ubar_c = sidis->get_ubar();
   //    double d_c = sidis->get_dA();
   //    double dbar_c = sidis->get_dbar();
   //    double s_c= sidis->get_s();
   //    double sbar_c = sidis->get_sbar();
   //    double g_c= sidis->get_g();

       TString cut = allcut+"*"+Form("(Q2>=%f && Q2<%f && z>=%f && z<%f)",Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);
       TString crosssectioncut = xbj_cut_small+"*"+Form("(dxs_hp)*(Q2>=%f && Q2<%f && z>=%f && z<%f)",Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);
       TString piminuscut = xbj_cut_small+"*"+Form("(weight_hm)*(Q2>=%f && Q2<%f && z>=%f && z<%f)",Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);

       double counts_in_bin = double(T0->GetEntries(TCut(cut)));
       cout << "Nof events for zbin " << i << " and Q2 bin " << j << ": "  << counts_in_bin << endl;
       c12pspart->SetBinContent(i+1,j+1,counts_in_bin);
       c12pspartrel->SetBinContent(i+1,j+1,counts_in_bin/nofallcut_nozbinQ2bincuts*100);

       TString histotitle = Form("%.1f #leq Q2 < %.1f and %.1f #leq z < %.1f",Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);
       TH1D *htemp1 = new TH1D("htemp1",histotitle, 500,0,1);
       htemp1->GetXaxis()->CenterTitle(1);
       htemp1->Sumw2();
       htemp1->GetXaxis()->SetTitle("z");
       htemp1->GetYaxis()->CenterTitle(1);
       htemp1->GetYaxis()->SetTitle("weight_hp");
       TH1D *htemp2 = new TH1D("htemp2",histotitle, 20,2,10);
       htemp2->GetXaxis()->CenterTitle(1);
       htemp2->GetXaxis()->SetTitle("Q2");
       htemp2->Sumw2();
       htemp2->GetYaxis()->CenterTitle(1);
       htemp2->GetYaxis()->SetTitle("weight_hp");
       TH1D *htemp3 = new TH1D("htemp3",histotitle, 500,0,1);
       htemp3->GetXaxis()->CenterTitle(1);
       htemp3->GetXaxis()->SetTitle("z");
       htemp3->GetYaxis()->CenterTitle(1);
       htemp3->GetYaxis()->SetTitle("dxs_hp");
       TH1D *htemp4 = new TH1D("htemp4",histotitle, 500,2,10);
       htemp4->GetXaxis()->CenterTitle(1);
       htemp4->GetXaxis()->SetTitle("Q2");
       htemp4->GetYaxis()->CenterTitle(1);
       htemp4->GetYaxis()->SetTitle("dxs_hp");
       TH1D *htemp5 = new TH1D("htemp5",histotitle, 20,2,10);
       htemp5->GetXaxis()->CenterTitle(1);
       htemp5->GetXaxis()->SetTitle("Q2");
       htemp5->Sumw2();
       htemp5->GetYaxis()->CenterTitle(1);
       htemp5->GetYaxis()->SetTitle("weight_hm");
       T0->Draw("Q2>>htemp5",TCut(piminuscut),"");

   //    canv_weighthist->cd(i+(j*zbin)+1);
       canv_weighthist->cd(1);
       T0->Draw("z>>htemp1",TCut(cut),"");
       canv_weighthist->cd(2);
       T0->Draw("Q2>>htemp2",TCut(cut),"");
       canv_weighthist->cd(3);
       T0->Draw("z>>htemp3",TCut(crosssectioncut),"");
       canv_weighthist->cd(4);
       T0->Draw("Q2>>htemp4",TCut(crosssectioncut),"");
     //  T0->Project("htemp1","z",TCut(cut));

       Double_t histsum = (double) htemp1->GetSum();
       Double_t histsum2 = (double) htemp2->GetSum();
       Double_t histsum3 = (double) htemp3->GetSum();
       Double_t histsum4 = (double) htemp4->GetSum();
       Double_t histsum5 = (double) htemp5->GetSum();
       Double_t histentries = htemp1->GetEntries();
       Double_t histentries2 = htemp2->GetEntries();
       Double_t histentries3 = htemp2->GetEntries();
       Double_t histentries4 = htemp2->GetEntries();

       Double_t hist2error = sqrt(htemp2->GetSumw2()->GetSum());
       Double_t hist5error = sqrt(htemp5->GetSumw2()->GetSum());
     //  cout << "Sum weight histo (z)" <<  histsum << " and histentries " << histentries<< endl;
     //  cout << "Sum weight histo " <<  histsum << " , my weight sum " << weightsum << endl;
     //  cout << "Sum weight histo2 (Q2) " <<  histsum2 << " and histentries2 " << histentries2<< endl;
     //  cout << "Cross section sum " << man_sumcross_p[i][j] << " , weightsum/crosssum " << weightsum/man_sumcross_p[i][j]  << endl;
     //  cout << "Sum weight histo3 (z)" <<  histsum3 << " and histentries3 " << histentries3<< endl;
     //  cout << "Sum weight histo4 (Q2) " <<  histsum4 << " and histentries4 " << histentries4<< endl;
       cout << "Pi plus results C12:" << endl;
       cout << " Histo Weighted Cross Section Integral = " << histsum2 << " +/- " <<  hist2error << endl;
       cout << " Weighted Cross Section Integral = " << weightsum_p << " +/- " << man_errorintegral_p[i][j] << endl;
       cout << "Pi minus results C12:" << endl;
       cout << " Histo Weighted Cross Section Integral = " << histsum5 << " +/- " <<  hist5error << endl;
       cout << " Weighted Cross Section Integral = " << weightsum_m << " +/- " << man_errorintegral_m[i][j] << endl;

       c12pionp_cross->SetBinContent(i+1,j+1,weightsum_p);
       c12pionp_crossw->SetBinContent(i+1,j+1,histsum2);
       c12pionm_cross->SetBinContent(i+1,j+1,weightsum_m);
       c12pionm_crossw->SetBinContent(i+1,j+1,histsum5);

       c12pionp_cross->SetBinError(i+1,j+1,man_errorintegral_p[i][j]);
       c12pionp_crossw->SetBinError(i+1,j+1,hist2error);
       c12pionm_cross->SetBinError(i+1,j+1,man_errorintegral_m[i][j]);
       c12pionm_crossw->SetBinError(i+1,j+1,hist5error);
     //  canv_weighthist->cd(i+(j*zbin)+1); //gPad->SetLogz(1);
     //  htemp1->Draw("");



     //  canv_Q2bins_ele->Update();


     }

   }


 c12nofdiff->Add(c12pionp_cross,c12pionm_cross,1,-1);
 c12nofdiffw->Add(c12pionp_crossw,c12pionm_crossw,1,-1);
 TCanvas *asymc12 = new TCanvas("asymc12","",1200,800);
 asymc12->Divide(3,2);
 asymc12->cd(1);
 c12pionp_cross->Draw("TEXTE");
 asymc12->cd(2);
 c12pionm_cross->Draw("TEXTE");
 asymc12->cd(3);
 c12nofdiff->Draw("TEXTE");
 asymc12->cd(4);
 c12pionp_crossw->Draw("TEXTE");
 asymc12->cd(5);
 c12pionm_crossw->Draw("TEXTE");
 asymc12->cd(6);
 c12nofdiffw->Draw("TEXTE");
 asymc12->SaveAs("c12difference.pdf");



 //canv_weighthist->SaveAs("weightest.pdf");
 TCanvas *ct2 = new TCanvas("ct2","pspart");
 c12pspart->Draw("TEXT");
 savestring.Clear();
 savestring = "c12_pspart" + outputfileending;
 ct2->SaveAs(savestring);

 delete T0;
 TChain *T1 = new TChain("T");
 T1->Add("../../test/massproduction_CTEQ/EIC_A2_pion_10_100_1_0.root");
 //T1->Add("../../test/EIC_A2_pion_10_100_1_0.root");
 T1->SetBranchAddress("Q2",&Q2);
 T1->SetBranchAddress("W",&W);
 T1->SetBranchAddress("Wp",&Wp);
 T1->SetBranchAddress("x",&x );
 T1->SetBranchAddress("y",&y );
 T1->SetBranchAddress("z",&z );
 T1->SetBranchAddress("nu",&nu );
 T1->SetBranchAddress("s",&s );
 T1->SetBranchAddress("epsilon",&epsilon );
 T1->SetBranchAddress("pt",&pt );
 T1->SetBranchAddress("weight_hp",&weight_hp );
 T1->SetBranchAddress("weight_hm",&weight_hm );
 T1->SetBranchAddress("weight_in",&weight_in );
 T1->SetBranchAddress("dxs_hm",&dxs_hm);
 T1->SetBranchAddress("dxs_hp",&dxs_hp);
 T1->SetBranchAddress("dxs_incl",&dxs_incl);
 Double_t d_nofallcut_nozbinQ2bincuts= (double) T1->GetEntries(TCut(allcut)) ;
 cout << "deuteron: For (x>0.05 && x<=0.1) and weight > 0: pos hadrons events: " << d_nofallcut_nozbinQ2bincuts << endl;

TH2D *d2pspart = new TH2D("d2pspart","Counts per z and Q2 bin for ^{}d(e,e'#pi^{+})X; z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
d2pspart->GetXaxis()->CenterTitle(1);
d2pspart->GetYaxis()->CenterTitle(1);
d2pspart->SetMarkerSize(2);
TH2D *d2pspartrel = new TH2D("d2pspartrel","Relative Counts per z and Q2 bin in percent for ^{}d(e,e'#pi^{+})X; z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
d2pspartrel->GetXaxis()->CenterTitle(1);
d2pspartrel->GetYaxis()->CenterTitle(1);
d2pspartrel->SetMarkerSize(2);

TH2D *d2pionp_cross = new TH2D("d2pionp_cross ","Integral weighted cross section with correct errors for ^{}d(e,e'#pi^{+})X;  z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
d2pionp_cross ->GetXaxis()->CenterTitle(1);
d2pionp_cross ->GetYaxis()->CenterTitle(1);
d2pionp_cross->SetMarkerSize(2);
TH2D *d2pionp_crossw = new TH2D("d2pionp_crossw ","Integral weighted cross section with wrong errors for ^{}d(e,e'#pi^{+})X;  z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
d2pionp_crossw->GetXaxis()->CenterTitle(1);
d2pionp_crossw->GetYaxis()->CenterTitle(1);
d2pionp_crossw->SetMarkerSize(2);
TH2D *d2pionm_cross = new TH2D("d2pionm_cross ","Integral weighted cross section with correct errors for ^{}d(e,e'#pi^{-})X;  z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
d2pionm_cross->GetXaxis()->CenterTitle(1);
d2pionm_cross->GetYaxis()->CenterTitle(1);
d2pionm_cross->SetMarkerSize(2);
TH2D *d2pionm_crossw = new TH2D("d2pionm_crossw ","Integral weighted cross section with wrong errors for ^{}d(e,e'#pi^{-})X;  z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
d2pionm_crossw ->GetXaxis()->CenterTitle(1);
d2pionm_crossw ->GetYaxis()->CenterTitle(1);
d2pionm_crossw->SetMarkerSize(2);
TH2D *d2nofdiff = new TH2D("d2nofdiff ","Count rate difference with correct errors for d #pi^{+}-#pi^{-};  z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
d2nofdiff->GetXaxis()->CenterTitle(1);
d2nofdiff->GetYaxis()->CenterTitle(1);
d2nofdiff->Sumw2();
d2nofdiff->SetMarkerSize(2);
TH2D *d2nofdiffw = new TH2D("d2nofdiffw ","Count rate difference with wrong errors for d #pi^{+}-#pi^{-};  z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
d2nofdiffw->GetXaxis()->CenterTitle(1);
d2nofdiffw->GetYaxis()->CenterTitle(1);
 d2nofdiffw->Sumw2();
 d2nofdiffw->SetMarkerSize(2);
 for (Int_t j=0;j<Q2bin;j++){
//  for (Int_t j=1;j<2;j++){//Test for one Q2 bin

    for (Int_t i=0;i<zbin;i++){
  //    for (Int_t i=2;i<3;i++){
      double weightsum_p = 0;
      man_sumcross_p[i][j] = 0;
      man_sumcross2_p[i][j] = 0;
      double weightsum_m = 0;
      man_sumcross_m[i][j] = 0;
      man_sumcross2_m[i][j] = 0;
      for (int nentry = 0; nentry<T1->GetEntries() ; nentry++)
      {
        T1->GetEntry(nentry);
        if (x>0.05 && x<=0.1 && weight_hp>0 && Q2>=Q2_cut[j] && Q2<Q2_cut[j+1] && z>=z_cut[i] && z<z_cut[i+1]){
          weightsum_p+=weight_hp;
          man_sumcross_p[i][j]+=dxs_hp;
          man_sumcross2_p[i][j]+=(dxs_hp*dxs_hp);
        }
        if (x>0.05 && x<=0.1 && weight_hm>0 && Q2>=Q2_cut[j] && Q2<Q2_cut[j+1] && z>=z_cut[i] && z<z_cut[i+1]){
          weightsum_m+=weight_hm;
          man_sumcross_m[i][j]+=dxs_hm;
          man_sumcross2_m[i][j]+=(dxs_hm*dxs_hm);
        }
      }
      man_errorintegral_p[i][j] = generated_phase_space/sqrt(nsim) * sqrt ( man_sumcross2_p[i][j]/nsim - pow (man_sumcross_p[i][j]/nsim, 2) );
      man_errorintegral_m[i][j] = generated_phase_space/sqrt(nsim) * sqrt ( man_sumcross2_m[i][j]/nsim - pow (man_sumcross_m[i][j]/nsim, 2) );

      double z_center = 0.5*(z_cut[i]+z_cut[i+1]);
      double Q2_center = pow(10., 0.5*(Q2_cut[j]+Q2_cut[j+1]));
      double xbj_center = (0.05+0.1)*0.5;

      TString cut = allcut+"*"+Form("(Q2>=%f && Q2<%f && z>=%f && z<%f)",Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);
      TString crosssectioncut = xbj_cut_small+"*"+Form("(dxs_hp)*(Q2>=%f && Q2<%f && z>=%f && z<%f)",Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);
      TString piminuscut = xbj_cut_small+"*"+Form("(weight_hm)*(Q2>=%f && Q2<%f && z>=%f && z<%f)",Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);

      double counts_in_bin = double(T1->GetEntries(TCut(cut)));
      cout << "Nof events for zbin " << i << " and Q2 bin " << j << ": "  << counts_in_bin << endl;
      d2pspart->SetBinContent(i+1,j+1,counts_in_bin);
      d2pspartrel->SetBinContent(i+1,j+1,counts_in_bin/nofallcut_nozbinQ2bincuts*100);

      TString histotitle = Form("%.1f #leq Q2 < %.1f and %.1f #leq z < %.1f",Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);
      TH1D *htemp1 = new TH1D("htemp1",histotitle, 500,0,1);
      htemp1->GetXaxis()->CenterTitle(1);
      htemp1->Sumw2();
      htemp1->GetXaxis()->SetTitle("z");
      htemp1->GetYaxis()->CenterTitle(1);
      htemp1->GetYaxis()->SetTitle("weight_hp");
      TH1D *htemp2 = new TH1D("htemp2",histotitle, 20,2,10);
      htemp2->GetXaxis()->CenterTitle(1);
      htemp2->GetXaxis()->SetTitle("Q2");
      htemp2->Sumw2();
      htemp2->GetYaxis()->CenterTitle(1);
      htemp2->GetYaxis()->SetTitle("weight_hp");
      TH1D *htemp3 = new TH1D("htemp3",histotitle, 500,0,1);
      htemp3->GetXaxis()->CenterTitle(1);
      htemp3->GetXaxis()->SetTitle("z");
      htemp3->GetYaxis()->CenterTitle(1);
      htemp3->GetYaxis()->SetTitle("dxs_hp");
      TH1D *htemp4 = new TH1D("htemp4",histotitle, 500,2,10);
      htemp4->GetXaxis()->CenterTitle(1);
      htemp4->GetXaxis()->SetTitle("Q2");
      htemp4->GetYaxis()->CenterTitle(1);
      htemp4->GetYaxis()->SetTitle("dxs_hp");
      TH1D *htemp5 = new TH1D("htemp5",histotitle, 20,2,10);
      htemp5->GetXaxis()->CenterTitle(1);
      htemp5->GetXaxis()->SetTitle("Q2");
      htemp5->Sumw2();
      htemp5->GetYaxis()->CenterTitle(1);
      htemp5->GetYaxis()->SetTitle("weight_hm");
      T1->Draw("Q2>>htemp5",TCut(piminuscut),"");

  //    canv_weighthist->cd(i+(j*zbin)+1);
      canv_weighthist->cd(1);
      T1->Draw("z>>htemp1",TCut(cut),"");
      canv_weighthist->cd(2);
      T1->Draw("Q2>>htemp2",TCut(cut),"");
      canv_weighthist->cd(3);
      T1->Draw("z>>htemp3",TCut(crosssectioncut),"");
      canv_weighthist->cd(4);
      T1->Draw("Q2>>htemp4",TCut(crosssectioncut),"");


      Double_t histsum = (double) htemp1->GetSum();
      Double_t histsum2 = (double) htemp2->GetSum();
      Double_t histsum3 = (double) htemp3->GetSum();
      Double_t histsum4 = (double) htemp4->GetSum();
      Double_t histsum5 = (double) htemp5->GetSum();
      Double_t histentries = htemp1->GetEntries();
      Double_t histentries2 = htemp2->GetEntries();
      Double_t histentries3 = htemp2->GetEntries();
      Double_t histentries4 = htemp2->GetEntries();

      Double_t hist2error = sqrt(htemp2->GetSumw2()->GetSum());
      Double_t hist5error = sqrt(htemp5->GetSumw2()->GetSum());
    //  cout << "Sum weight histo (z)" <<  histsum << " and histentries " << histentries<< endl;
    //  cout << "Sum weight histo " <<  histsum << " , my weight sum " << weightsum << endl;
    //  cout << "Sum weight histo2 (Q2) " <<  histsum2 << " and histentries2 " << histentries2<< endl;
    //  cout << "Cross section sum " << man_sumcross_p[i][j] << " , weightsum/crosssum " << weightsum/man_sumcross_p[i][j]  << endl;
    //  cout << "Sum weight histo3 (z)" <<  histsum3 << " and histentries3 " << histentries3<< endl;
    //  cout << "Sum weight histo4 (Q2) " <<  histsum4 << " and histentries4 " << histentries4<< endl;
      cout << "Pi plus results C12:" << endl;
      cout << " Histo Weighted Cross Section Integral = " << histsum2 << " +/- " <<  hist2error << endl;
      cout << " Weighted Cross Section Integral = " << weightsum_p << " +/- " << man_errorintegral_p[i][j] << endl;
      cout << "Pi minus results C12:" << endl;
      cout << " Histo Weighted Cross Section Integral = " << histsum5 << " +/- " <<  hist5error << endl;
      cout << " Weighted Cross Section Integral = " << weightsum_m << " +/- " << man_errorintegral_m[i][j] << endl;

      d2pionp_cross->SetBinContent(i+1,j+1,weightsum_p);
      d2pionp_crossw->SetBinContent(i+1,j+1,histsum2);
      d2pionm_cross->SetBinContent(i+1,j+1,weightsum_m);
      d2pionm_crossw->SetBinContent(i+1,j+1,histsum5);

      d2pionp_cross->SetBinError(i+1,j+1,man_errorintegral_p[i][j]);
      d2pionp_crossw->SetBinError(i+1,j+1,hist2error);
      d2pionm_cross->SetBinError(i+1,j+1,man_errorintegral_m[i][j]);
      d2pionm_crossw->SetBinError(i+1,j+1,hist5error);


    }

  }


  d2nofdiff->Add(d2pionp_cross,d2pionm_cross,1,-1);
  d2nofdiffw->Add(d2pionp_crossw,d2pionm_crossw,1,-1);
  TCanvas *asymd2 = new TCanvas("asymd2","",1200,800);
  asymd2->Divide(3,2);
  asymd2->cd(1);
  d2pionp_cross->Draw("TEXTE");
  d2pionp_cross->SetMarkerSize(2);
  asymd2->cd(2);
  d2pionm_cross->Draw("TEXTE");
  asymd2->cd(3);
  d2nofdiff->Draw("TEXTE");
  asymd2->cd(4);
  d2pionp_crossw->Draw("TEXTE");
  asymd2->cd(5);
  d2pionm_crossw->Draw("TEXTE");
  asymd2->cd(6);
  d2nofdiffw->Draw("TEXTE");
  asymd2->SaveAs("d2difference.pdf");

  //canv_weighthist->SaveAs("weightest.pdf");
  TCanvas *ct4 = new TCanvas("ct4","pspart");
  d2pspart->Draw("TEXT");
  savestring.Clear();
  savestring = "d2_pspart" + outputfileending;
  ct4->SaveAs(savestring);



  TH2D *asym = new TH2D("asym ","Asymmetry with correct errors;  z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
  asym->GetXaxis()->CenterTitle(1);
  asym->GetYaxis()->CenterTitle(1);
  asym->SetMarkerSize(2);
  asym->Sumw2();
  TH2D *asymw = new TH2D("asymw","Asymmetry with wrong errors;  z (GeV) ; Q2 (GeV^{2})",zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
  asymw->GetXaxis()->CenterTitle(1);
  asymw->GetYaxis()->CenterTitle(1);
  asymw->SetMarkerSize(2);
  asymw->Sumw2();
  asym->Divide(c12nofdiff,d2nofdiff,1/12.,1/2.);
  asymw->Divide(c12nofdiffw,d2nofdiffw,1/12.,1/2.);

  TCanvas *ct6 = new TCanvas("ct6","");
  ct6->Divide(2,1);
  ct6->cd(1);
  asym->Draw("TEXTE");
  ct6->cd(2);
  asymw->Draw("TEXTE");
  ct6->SaveAs("asym.pdf");
}
