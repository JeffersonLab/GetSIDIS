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
//#include <TMatrix.h>
/*}}}*/

using namespace std;
const double DEG = 180./3.1415926;
const double PI = 3.1415926;

int main(Int_t argc, char *argv[]){
	int fA = 0; cerr<<"-- What nucleus: A = ? "; cin >> fA;

    TChain *T = new TChain("T");
    for(int i=0;i<100;i++){
        for(int j=1;j<=4;j++){
            if(fA ==12)
                //T->Add(Form("./c12_pion_LO_free_noPt/EIC_A12_pion_10_600_%d_%d.root", j, i));
                T->Add(Form("./c12_pion_LO_noPt/EIC_A12_pion_10_600_%d_%d.root", j, i));
             if(fA ==2)
                //T->Add(Form("./d2_pion_LO_free_noPt/EIC_A2_pion_10_100_%d_%d.root", j, i));
                T->Add(Form("./d2_pion_LO_noPt/EIC_A2_pion_10_100_%d_%d.root", j, i));
        }
    }
    /*Define{{{*/
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
    ULong64_t nsim = 0, Nsim1 = 0, Nsim2 = 0, Nsim3 = 0,Nsim4 = 0;
    //For Beam Position and Vertex info
    Double_t vx_ele, vy_ele, vz_ele, vx_had, vy_had, vz_had;
    Double_t D_fav, D_unfav, D_s, D_g;
    /*}}}*/
 
	/*Define old root file{{{*/
    T->SetBranchAddress("Q2",&Q2);/*{{{*/
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
    /*}}}*/

    /*}}}*/

    ULong64_t N_Total=T->GetEntries();

    TString new_filename[8];
    if(fA ==12){
        for(int i=0;i<8;i++)
            //new_filename[i]=Form("./c12_pion_LO_free_noPt/EIC_A12_pion_10_600_skim%d_wide_free_noPt.root",i);
            new_filename[i]=Form("./c12_pion_LO_noPt/EIC_A12_pion_10_600_skim%d_wide_noPt_new1.root",i);
    }
    if(fA ==2){
        for(int i=0;i<8;i++)
            //new_filename[i]=Form("./d2_pion_LO_free_noPt/EIC_A2_pion_10_100_skim%d_wide_free_noPt.root",i);
            new_filename[i]=Form("./d2_pion_LO_noPt/EIC_A2_pion_10_100_skim%d_wide_noPt_new1.root",i);
    }

    /*Define new rootfile for each bin{{{*/
    TFile *file1= new TFile(new_filename[0],"RECREATE");
    TTree *t1 = new TTree("T","T");
    t1->SetDirectory(file1);
    t1->Branch("Q2",&Q2,"data/D");/*{{{*/
    t1->Branch("W",&W,"data/D");
    t1->Branch("Wp",&Wp,"data/D");
    t1->Branch("x",&x,"data/D");
    t1->Branch("y",&y,"data/D");
    t1->Branch("z",&z,"data/D");
    t1->Branch("nu",&nu,"data/D");
    t1->Branch("s",&s,"data/D");
    t1->Branch("gamma",&gamma,"data/D");
    t1->Branch("epsilon",&epsilon,"data/D");
    t1->Branch("pt",&pt,"data/D");
    t1->Branch("rapidity",&rapidity,"data/D");
    t1->Branch("theta_q",&theta_q,"data/D");
    t1->Branch("theta_s",&theta_s,"data/D");
    t1->Branch("phi_h",&phi_h,"data/D");
    t1->Branch("phi_s",&phi_s,"data/D");
    t1->Branch("jacoF",&jacoF,"jacoF/D");
    t1->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
    t1->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
    t1->Branch("dxs_incl",&dxs_incl,"dxs_incl/D");
    t1->Branch("mom_ele",&mom_ele,"mom_ele/D");
    t1->Branch("mom_gen_ele",&mom_gen_ele,"mom_gen_ele/D");
    t1->Branch("mom_had",&mom_had,"mom_had/D");
    t1->Branch("mom_gen_had",&mom_gen_had,"mom_gen_had/D");
    t1->Branch("theta_ele",&theta_ele,"theta_ele/D");
    t1->Branch("theta_gen_ele",&theta_gen_ele,"theta_gen_ele/D");
    t1->Branch("theta_had",&theta_had,"theta_had/D");
    t1->Branch("theta_gen_had",&theta_gen_had,"theta_gen_had/D");
    t1->Branch("phi_ele",&phi_ele,"phi_ele/D");
    t1->Branch("phi_gen_ele",&phi_gen_ele,"phi_gen_ele/D");
    t1->Branch("phi_had",&phi_had,"phi_had/D");
    t1->Branch("phi_gen_had",&phi_gen_had,"phi_gen_had/D");
    t1->Branch("nsim",&nsim,"nsim/l");
    t1->Branch("dilute_p",&dilute_hp,"data/D");
    t1->Branch("dilute_m",&dilute_hm ,"data/D");
    t1->Branch("px_ele",&px_ele, "px_ele/D");
    t1->Branch("py_ele",&py_ele, "py_ele/D");
    t1->Branch("pz_ele",&pz_ele, "pz_ele/D");
    t1->Branch("E_ele",&E_ele, "E_ele/D");
    t1->Branch("px_had",&px_had, "px_had/D");
    t1->Branch("py_had",&py_had, "py_had/D");
    t1->Branch("pz_had",&pz_had, "pz_had/D");
    t1->Branch("E_had",&E_had, "E_had/D");
    t1->Branch("vx_ele",&vx_ele, "vx_ele/D");
    t1->Branch("vy_ele",&vy_ele, "vy_ele/D");
    t1->Branch("vz_ele",&vz_ele, "vz_ele/D");
    t1->Branch("vx_had",&vx_had, "vx_had/D");
    t1->Branch("vy_had",&vy_had, "vy_had/D");
    t1->Branch("vz_had",&vz_had, "vz_had/D");
    
    t1->Branch("u_pdf", &u_pdf, "u_pdf/D");
    t1->Branch("d_pdf", &d_pdf, "d_pdf/D");
    t1->Branch("s_pdf", &s_pdf, "s_pdf/D");
    t1->Branch("g_pdf", &g_pdf, "g_pdf/D");
    t1->Branch("ubar_pdf", &ubar_pdf, "ubar_pdf/D");
    t1->Branch("dbar_pdf", &dbar_pdf, "dbar_pdf/D");
    t1->Branch("sbar_pdf", &sbar_pdf, "sbar_pdf/D");
    t1->Branch("weight_hp",&weight_hp,"data/D");
    t1->Branch("weight_hm",&weight_hm,"data/D");
    t1->Branch("weight_in",&weight_in,"data/D");
    
    t1->Branch("D_fav", &D_fav, "D_fav/D");
    t1->Branch("D_unfav", &D_unfav, "D_unfav/D");
    t1->Branch("D_s", &D_s, "D_s/D");
    t1->Branch("D_g", &D_g, "D_g/D");
    /*}}}*/

    TFile *file2= new TFile(new_filename[1],"RECREATE");
    TTree *t2 = new TTree("T","T");
    t2->SetDirectory(file2);
    t2->Branch("Q2",&Q2,"data/D");/*{{{*/
    t2->Branch("W",&W,"data/D");
    t2->Branch("Wp",&Wp,"data/D");
    t2->Branch("x",&x,"data/D");
    t2->Branch("y",&y,"data/D");
    t2->Branch("z",&z,"data/D");
    t2->Branch("nu",&nu,"data/D");
    t2->Branch("s",&s,"data/D");
    t2->Branch("pt",&pt,"data/D");
    t2->Branch("gamma",&gamma,"data/D");
    t2->Branch("epsilon",&epsilon,"data/D");
    t2->Branch("weight_hp",&weight_hp,"data/D");
    t2->Branch("weight_hm",&weight_hm,"data/D");
    t2->Branch("weight_in",&weight_in,"data/D");
    t2->Branch("rapidity",&rapidity,"data/D");
    t2->Branch("theta_q",&theta_q,"data/D");
    t2->Branch("theta_s",&theta_s,"data/D");
    t2->Branch("phi_h",&phi_h,"data/D");
    t2->Branch("phi_s",&phi_s,"data/D");
    t2->Branch("jacoF",&jacoF,"jacoF/D");
    t2->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
    t2->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
    t2->Branch("dxs_incl",&dxs_incl,"dxs_incl/D");
    t2->Branch("mom_ele",&mom_ele,"mom_ele/D");
    t2->Branch("mom_gen_ele",&mom_gen_ele,"mom_gen_ele/D");
    t2->Branch("mom_had",&mom_had,"mom_had/D");
    t2->Branch("mom_gen_had",&mom_gen_had,"mom_gen_had/D");
    t2->Branch("theta_ele",&theta_ele,"theta_ele/D");
    t2->Branch("theta_gen_ele",&theta_gen_ele,"theta_gen_ele/D");
    t2->Branch("theta_had",&theta_had,"theta_had/D");
    t2->Branch("theta_gen_had",&theta_gen_had,"theta_gen_had/D");
    t2->Branch("phi_ele",&phi_ele,"phi_ele/D");
    t2->Branch("phi_gen_ele",&phi_gen_ele,"phi_gen_ele/D");
    t2->Branch("phi_had",&phi_had,"phi_had/D");
    t2->Branch("phi_gen_had",&phi_gen_had,"phi_gen_had/D");
    t2->Branch("nsim",&nsim,"nsim/l");
    t2->Branch("dilute_p",&dilute_hp,"data/D");
    t2->Branch("dilute_m",&dilute_hm ,"data/D");
    t2->Branch("px_ele",&px_ele, "px_ele/D");
    t2->Branch("py_ele",&py_ele, "py_ele/D");
    t2->Branch("pz_ele",&pz_ele, "pz_ele/D");
    t2->Branch("E_ele",&E_ele, "E_ele/D");
    t2->Branch("px_had",&px_had, "px_had/D");
    t2->Branch("py_had",&py_had, "py_had/D");
    t2->Branch("pz_had",&pz_had, "pz_had/D");
    t2->Branch("E_had",&E_had, "E_had/D");
    t2->Branch("vx_ele",&vx_ele, "vx_ele/D");
    t2->Branch("vy_ele",&vy_ele, "vy_ele/D");
    t2->Branch("vz_ele",&vz_ele, "vz_ele/D");
    t2->Branch("vx_had",&vx_had, "vx_had/D");
    t2->Branch("vy_had",&vy_had, "vy_had/D");
    t2->Branch("vz_had",&vz_had, "vz_had/D");
    
    t2->Branch("u_pdf", &u_pdf, "u_pdf/D");
    t2->Branch("d_pdf", &d_pdf, "d_pdf/D");
    t2->Branch("s_pdf", &s_pdf, "s_pdf/D");
    t2->Branch("g_pdf", &g_pdf, "g_pdf/D");
    t2->Branch("ubar_pdf", &ubar_pdf, "ubar_pdf/D");
    t2->Branch("dbar_pdf", &dbar_pdf, "dbar_pdf/D");
    t2->Branch("sbar_pdf", &sbar_pdf, "sbar_pdf/D");
    
    t2->Branch("D_fav", &D_fav, "D_fav/D");
    t2->Branch("D_unfav", &D_unfav, "D_unfav/D");
    t2->Branch("D_s", &D_s, "D_s/D");
    t2->Branch("D_g", &D_g, "D_g/D");
  /*}}}*/

    TFile *file3= new TFile(new_filename[2],"RECREATE");
    TTree *t3 = new TTree("T","T");
    t3->SetDirectory(file3);
    t3->Branch("Q2",&Q2,"data/D");/*{{{*/
    t3->Branch("W",&W,"data/D");
    t3->Branch("Wp",&Wp,"data/D");
    t3->Branch("x",&x,"data/D");
    t3->Branch("y",&y,"data/D");
    t3->Branch("z",&z,"data/D");
    t3->Branch("nu",&nu,"data/D");
    t3->Branch("s",&s,"data/D");
    t3->Branch("gamma",&gamma,"data/D");
    t3->Branch("epsilon",&epsilon,"data/D");
    t3->Branch("pt",&pt,"data/D");
    t3->Branch("weight_hp",&weight_hp,"data/D");
    t3->Branch("weight_hm",&weight_hm,"data/D");
    t3->Branch("weight_in",&weight_in,"data/D");
    t3->Branch("rapidity",&rapidity,"data/D");
    t3->Branch("theta_q",&theta_q,"data/D");
    t3->Branch("theta_s",&theta_s,"data/D");
    t3->Branch("phi_h",&phi_h,"data/D");
    t3->Branch("phi_s",&phi_s,"data/D");
    t3->Branch("jacoF",&jacoF,"jacoF/D");
    t3->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
    t3->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
    t3->Branch("dxs_incl",&dxs_incl,"dxs_incl/D");
    t3->Branch("mom_ele",&mom_ele,"mom_ele/D");
    t3->Branch("mom_gen_ele",&mom_gen_ele,"mom_gen_ele/D");
    t3->Branch("mom_had",&mom_had,"mom_had/D");
    t3->Branch("mom_gen_had",&mom_gen_had,"mom_gen_had/D");
    t3->Branch("theta_ele",&theta_ele,"theta_ele/D");
    t3->Branch("theta_gen_ele",&theta_gen_ele,"theta_gen_ele/D");
    t3->Branch("theta_had",&theta_had,"theta_had/D");
    t3->Branch("theta_gen_had",&theta_gen_had,"theta_gen_had/D");
    t3->Branch("phi_ele",&phi_ele,"phi_ele/D");
    t3->Branch("phi_gen_ele",&phi_gen_ele,"phi_gen_ele/D");
    t3->Branch("phi_had",&phi_had,"phi_had/D");
    t3->Branch("phi_gen_had",&phi_gen_had,"phi_gen_had/D");
    t3->Branch("nsim",&nsim,"nsim/l");
    t3->Branch("dilute_p",&dilute_hp,"data/D");
    t3->Branch("dilute_m",&dilute_hm ,"data/D");
    t3->Branch("px_ele",&px_ele, "px_ele/D");
    t3->Branch("py_ele",&py_ele, "py_ele/D");
    t3->Branch("pz_ele",&pz_ele, "pz_ele/D");
    t3->Branch("E_ele",&E_ele, "E_ele/D");
    t3->Branch("px_had",&px_had, "px_had/D");
    t3->Branch("py_had",&py_had, "py_had/D");
    t3->Branch("pz_had",&pz_had, "pz_had/D");
    t3->Branch("E_had",&E_had, "E_had/D");
    t3->Branch("vx_ele",&vx_ele, "vx_ele/D");
    t3->Branch("vy_ele",&vy_ele, "vy_ele/D");
    t3->Branch("vz_ele",&vz_ele, "vz_ele/D");
    t3->Branch("vx_had",&vx_had, "vx_had/D");
    t3->Branch("vy_had",&vy_had, "vy_had/D");
    t3->Branch("vz_had",&vz_had, "vz_had/D");
    
    t3->Branch("u_pdf", &u_pdf, "u_pdf/D");
    t3->Branch("d_pdf", &d_pdf, "d_pdf/D");
    t3->Branch("s_pdf", &s_pdf, "s_pdf/D");
    t3->Branch("g_pdf", &g_pdf, "g_pdf/D");
    t3->Branch("ubar_pdf", &ubar_pdf, "ubar_pdf/D");
    t3->Branch("dbar_pdf", &dbar_pdf, "dbar_pdf/D");
    t3->Branch("sbar_pdf", &sbar_pdf, "sbar_pdf/D");
   
    t3->Branch("D_fav", &D_fav, "D_fav/D");
    t3->Branch("D_unfav", &D_unfav, "D_unfav/D");
    t3->Branch("D_s", &D_s, "D_s/D");
    t3->Branch("D_g", &D_g, "D_g/D");
   /*}}}*/

    TFile *file4= new TFile(new_filename[3],"RECREATE");
    TTree *t4 = new TTree("T","T");
    t4->SetDirectory(file4);
    t4->Branch("Q2",&Q2,"data/D");/*{{{*/
    t4->Branch("W",&W,"data/D");
    t4->Branch("Wp",&Wp,"data/D");
    t4->Branch("x",&x,"data/D");
    t4->Branch("y",&y,"data/D");
    t4->Branch("z",&z,"data/D");
    t4->Branch("nu",&nu,"data/D");
    t4->Branch("s",&s,"data/D");
    t4->Branch("pt",&pt,"data/D");
    t4->Branch("gamma",&gamma,"data/D");
    t4->Branch("epsilon",&epsilon,"data/D");
    t4->Branch("weight_hp",&weight_hp,"data/D");
    t4->Branch("weight_hm",&weight_hm,"data/D");
    t4->Branch("weight_in",&weight_in,"data/D");
    t4->Branch("rapidity",&rapidity,"data/D");
    t4->Branch("theta_q",&theta_q,"data/D");
    t4->Branch("theta_s",&theta_s,"data/D");
    t4->Branch("phi_h",&phi_h,"data/D");
    t4->Branch("phi_s",&phi_s,"data/D");
    t4->Branch("jacoF",&jacoF,"jacoF/D");
    t4->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
    t4->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
    t4->Branch("dxs_incl",&dxs_incl,"dxs_incl/D");
    t4->Branch("mom_ele",&mom_ele,"mom_ele/D");
    t4->Branch("mom_gen_ele",&mom_gen_ele,"mom_gen_ele/D");
    t4->Branch("mom_had",&mom_had,"mom_had/D");
    t4->Branch("mom_gen_had",&mom_gen_had,"mom_gen_had/D");
    t4->Branch("theta_ele",&theta_ele,"theta_ele/D");
    t4->Branch("theta_gen_ele",&theta_gen_ele,"theta_gen_ele/D");
    t4->Branch("theta_had",&theta_had,"theta_had/D");
    t4->Branch("theta_gen_had",&theta_gen_had,"theta_gen_had/D");
    t4->Branch("phi_ele",&phi_ele,"phi_ele/D");
    t4->Branch("phi_gen_ele",&phi_gen_ele,"phi_gen_ele/D");
    t4->Branch("phi_had",&phi_had,"phi_had/D");
    t4->Branch("phi_gen_had",&phi_gen_had,"phi_gen_had/D");
    t4->Branch("nsim",&nsim,"nsim/l");
    t4->Branch("dilute_p",&dilute_hp,"data/D");
    t4->Branch("dilute_m",&dilute_hm ,"data/D");
    t4->Branch("px_ele",&px_ele, "px_ele/D");
    t4->Branch("py_ele",&py_ele, "py_ele/D");
    t4->Branch("pz_ele",&pz_ele, "pz_ele/D");
    t4->Branch("E_ele",&E_ele, "E_ele/D");
    t4->Branch("px_had",&px_had, "px_had/D");
    t4->Branch("py_had",&py_had, "py_had/D");
    t4->Branch("pz_had",&pz_had, "pz_had/D");
    t4->Branch("E_had",&E_had, "E_had/D");
    t4->Branch("vx_ele",&vx_ele, "vx_ele/D");
    t4->Branch("vy_ele",&vy_ele, "vy_ele/D");
    t4->Branch("vz_ele",&vz_ele, "vz_ele/D");
    t4->Branch("vx_had",&vx_had, "vx_had/D");
    t4->Branch("vy_had",&vy_had, "vy_had/D");
    t4->Branch("vz_had",&vz_had, "vz_had/D");
    
    t4->Branch("u_pdf", &u_pdf, "u_pdf/D");
    t4->Branch("d_pdf", &d_pdf, "d_pdf/D");
    t4->Branch("s_pdf", &s_pdf, "s_pdf/D");
    t4->Branch("g_pdf", &g_pdf, "g_pdf/D");
    t4->Branch("ubar_pdf", &ubar_pdf, "ubar_pdf/D");
    t4->Branch("dbar_pdf", &dbar_pdf, "dbar_pdf/D");
    t4->Branch("sbar_pdf", &sbar_pdf, "sbar_pdf/D");
    t4->Branch("D_fav", &D_fav, "D_fav/D");
    t4->Branch("D_unfav", &D_unfav, "D_unfav/D");
    t4->Branch("D_s", &D_s, "D_s/D");
    t4->Branch("D_g", &D_g, "D_g/D");
  /*}}}*/

    TFile *file5= new TFile(new_filename[4],"RECREATE");
    TTree *t5 = new TTree("T","T");
    t5->SetDirectory(file5);
    t5->Branch("Q2",&Q2,"data/D");/*{{{*/
    t5->Branch("W",&W,"data/D");
    t5->Branch("Wp",&Wp,"data/D");
    t5->Branch("x",&x,"data/D");
    t5->Branch("y",&y,"data/D");
    t5->Branch("z",&z,"data/D");
    t5->Branch("nu",&nu,"data/D");
    t5->Branch("s",&s,"data/D");
    t5->Branch("gamma",&gamma,"data/D");
    t5->Branch("epsilon",&epsilon,"data/D");
    t5->Branch("pt",&pt,"data/D");
    t5->Branch("weight_hp",&weight_hp,"data/D");
    t5->Branch("weight_hm",&weight_hm,"data/D");
    t5->Branch("weight_in",&weight_in,"data/D");
    t5->Branch("rapidity",&rapidity,"data/D");
    t5->Branch("theta_q",&theta_q,"data/D");
    t5->Branch("theta_s",&theta_s,"data/D");
    t5->Branch("phi_h",&phi_h,"data/D");
    t5->Branch("phi_s",&phi_s,"data/D");
    t5->Branch("jacoF",&jacoF,"jacoF/D");
    t5->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
    t5->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
    t5->Branch("dxs_incl",&dxs_incl,"dxs_incl/D");
    t5->Branch("mom_ele",&mom_ele,"mom_ele/D");
    t5->Branch("mom_gen_ele",&mom_gen_ele,"mom_gen_ele/D");
    t5->Branch("mom_had",&mom_had,"mom_had/D");
    t5->Branch("mom_gen_had",&mom_gen_had,"mom_gen_had/D");
    t5->Branch("theta_ele",&theta_ele,"theta_ele/D");
    t5->Branch("theta_gen_ele",&theta_gen_ele,"theta_gen_ele/D");
    t5->Branch("theta_had",&theta_had,"theta_had/D");
    t5->Branch("theta_gen_had",&theta_gen_had,"theta_gen_had/D");
    t5->Branch("phi_ele",&phi_ele,"phi_ele/D");
    t5->Branch("phi_gen_ele",&phi_gen_ele,"phi_gen_ele/D");
    t5->Branch("phi_had",&phi_had,"phi_had/D");
    t5->Branch("phi_gen_had",&phi_gen_had,"phi_gen_had/D");
    t5->Branch("nsim",&nsim,"nsim/l");
    t5->Branch("dilute_p",&dilute_hp,"data/D");
    t5->Branch("dilute_m",&dilute_hm ,"data/D");
    t5->Branch("px_ele",&px_ele, "px_ele/D");
    t5->Branch("py_ele",&py_ele, "py_ele/D");
    t5->Branch("pz_ele",&pz_ele, "pz_ele/D");
    t5->Branch("E_ele",&E_ele, "E_ele/D");
    t5->Branch("px_had",&px_had, "px_had/D");
    t5->Branch("py_had",&py_had, "py_had/D");
    t5->Branch("pz_had",&pz_had, "pz_had/D");
    t5->Branch("E_had",&E_had, "E_had/D");
    t5->Branch("vx_ele",&vx_ele, "vx_ele/D");
    t5->Branch("vy_ele",&vy_ele, "vy_ele/D");
    t5->Branch("vz_ele",&vz_ele, "vz_ele/D");
    t5->Branch("vx_had",&vx_had, "vx_had/D");
    t5->Branch("vy_had",&vy_had, "vy_had/D");
    t5->Branch("vz_had",&vz_had, "vz_had/D");
    
    t5->Branch("u_pdf", &u_pdf, "u_pdf/D");
    t5->Branch("d_pdf", &d_pdf, "d_pdf/D");
    t5->Branch("s_pdf", &s_pdf, "s_pdf/D");
    t5->Branch("g_pdf", &g_pdf, "g_pdf/D");
    t5->Branch("ubar_pdf", &ubar_pdf, "ubar_pdf/D");
    t5->Branch("dbar_pdf", &dbar_pdf, "dbar_pdf/D");
    t5->Branch("sbar_pdf", &sbar_pdf, "sbar_pdf/D");
    t5->Branch("D_fav", &D_fav, "D_fav/D");
    t5->Branch("D_unfav", &D_unfav, "D_unfav/D");
    t5->Branch("D_s", &D_s, "D_s/D");
    t5->Branch("D_g", &D_g, "D_g/D");
   /*}}}*/

    TFile *file6= new TFile(new_filename[5],"RECREATE");
    TTree *t6 = new TTree("T","T");
    t6->SetDirectory(file6);
    t6->Branch("Q2",&Q2,"data/D");/*{{{*/
    t6->Branch("W",&W,"data/D");
    t6->Branch("Wp",&Wp,"data/D");
    t6->Branch("x",&x,"data/D");
    t6->Branch("y",&y,"data/D");
    t6->Branch("z",&z,"data/D");
    t6->Branch("nu",&nu,"data/D");
    t6->Branch("s",&s,"data/D");
    t6->Branch("gamma",&gamma,"data/D");
    t6->Branch("epsilon",&epsilon,"data/D");
    t6->Branch("pt",&pt,"data/D");
    t6->Branch("weight_hp",&weight_hp,"data/D");
    t6->Branch("weight_hm",&weight_hm,"data/D");
    t6->Branch("weight_in",&weight_in,"data/D");
    t6->Branch("rapidity",&rapidity,"data/D");
    t6->Branch("theta_q",&theta_q,"data/D");
    t6->Branch("theta_s",&theta_s,"data/D");
    t6->Branch("phi_h",&phi_h,"data/D");
    t6->Branch("phi_s",&phi_s,"data/D");
    t6->Branch("jacoF",&jacoF,"jacoF/D");
    t6->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
    t6->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
    t6->Branch("dxs_incl",&dxs_incl,"dxs_incl/D");
    t6->Branch("mom_ele",&mom_ele,"mom_ele/D");
    t6->Branch("mom_gen_ele",&mom_gen_ele,"mom_gen_ele/D");
    t6->Branch("mom_had",&mom_had,"mom_had/D");
    t6->Branch("mom_gen_had",&mom_gen_had,"mom_gen_had/D");
    t6->Branch("theta_ele",&theta_ele,"theta_ele/D");
    t6->Branch("theta_gen_ele",&theta_gen_ele,"theta_gen_ele/D");
    t6->Branch("theta_had",&theta_had,"theta_had/D");
    t6->Branch("theta_gen_had",&theta_gen_had,"theta_gen_had/D");
    t6->Branch("phi_ele",&phi_ele,"phi_ele/D");
    t6->Branch("phi_gen_ele",&phi_gen_ele,"phi_gen_ele/D");
    t6->Branch("phi_had",&phi_had,"phi_had/D");
    t6->Branch("phi_gen_had",&phi_gen_had,"phi_gen_had/D");
    t6->Branch("nsim",&nsim,"nsim/l");
    t6->Branch("dilute_p",&dilute_hp,"data/D");
    t6->Branch("dilute_m",&dilute_hm ,"data/D");
    t6->Branch("px_ele",&px_ele, "px_ele/D");
    t6->Branch("py_ele",&py_ele, "py_ele/D");
    t6->Branch("pz_ele",&pz_ele, "pz_ele/D");
    t6->Branch("E_ele",&E_ele, "E_ele/D");
    t6->Branch("px_had",&px_had, "px_had/D");
    t6->Branch("py_had",&py_had, "py_had/D");
    t6->Branch("pz_had",&pz_had, "pz_had/D");
    t6->Branch("E_had",&E_had, "E_had/D");
    t6->Branch("vx_ele",&vx_ele, "vx_ele/D");
    t6->Branch("vy_ele",&vy_ele, "vy_ele/D");
    t6->Branch("vz_ele",&vz_ele, "vz_ele/D");
    t6->Branch("vx_had",&vx_had, "vx_had/D");
    t6->Branch("vy_had",&vy_had, "vy_had/D");
    t6->Branch("vz_had",&vz_had, "vz_had/D");
    
    t6->Branch("u_pdf", &u_pdf, "u_pdf/D");
    t6->Branch("d_pdf", &d_pdf, "d_pdf/D");
    t6->Branch("s_pdf", &s_pdf, "s_pdf/D");
    t6->Branch("g_pdf", &g_pdf, "g_pdf/D");
    t6->Branch("ubar_pdf", &ubar_pdf, "ubar_pdf/D");
    t6->Branch("dbar_pdf", &dbar_pdf, "dbar_pdf/D");
    t6->Branch("sbar_pdf", &sbar_pdf, "sbar_pdf/D");
    t6->Branch("D_fav", &D_fav, "D_fav/D");
    t6->Branch("D_unfav", &D_unfav, "D_unfav/D");
    t6->Branch("D_s", &D_s, "D_s/D");
    t6->Branch("D_g", &D_g, "D_g/D");
  /*}}}*/

    TFile *file7= new TFile(new_filename[6],"RECREATE");
    TTree *t7 = new TTree("T","T");
    t7->SetDirectory(file7);
    t7->Branch("Q2",&Q2,"data/D");/*{{{*/
    t7->Branch("W",&W,"data/D");
    t7->Branch("Wp",&Wp,"data/D");
    t7->Branch("x",&x,"data/D");
    t7->Branch("y",&y,"data/D");
    t7->Branch("z",&z,"data/D");
    t7->Branch("nu",&nu,"data/D");
    t7->Branch("s",&s,"data/D");
    t7->Branch("gamma",&gamma,"data/D");
    t7->Branch("epsilon",&epsilon,"data/D");
    t7->Branch("pt",&pt,"data/D");
    t7->Branch("weight_hp",&weight_hp,"data/D");
   t7->Branch("weight_hm",&weight_hm,"data/D");
    t7->Branch("weight_in",&weight_in,"data/D");
    t7->Branch("rapidity",&rapidity,"data/D");
    t7->Branch("theta_q",&theta_q,"data/D");
    t7->Branch("theta_s",&theta_s,"data/D");
    t7->Branch("phi_h",&phi_h,"data/D");
    t7->Branch("phi_s",&phi_s,"data/D");
    t7->Branch("jacoF",&jacoF,"jacoF/D");
    t7->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
    t7->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
    t7->Branch("dxs_incl",&dxs_incl,"dxs_incl/D");
    t7->Branch("mom_ele",&mom_ele,"mom_ele/D");
    t7->Branch("mom_gen_ele",&mom_gen_ele,"mom_gen_ele/D");
    t7->Branch("mom_had",&mom_had,"mom_had/D");
    t7->Branch("mom_gen_had",&mom_gen_had,"mom_gen_had/D");
    t7->Branch("theta_ele",&theta_ele,"theta_ele/D");
    t7->Branch("theta_gen_ele",&theta_gen_ele,"theta_gen_ele/D");
    t7->Branch("theta_had",&theta_had,"theta_had/D");
    t7->Branch("theta_gen_had",&theta_gen_had,"theta_gen_had/D");
    t7->Branch("phi_ele",&phi_ele,"phi_ele/D");
    t7->Branch("phi_gen_ele",&phi_gen_ele,"phi_gen_ele/D");
    t7->Branch("phi_had",&phi_had,"phi_had/D");
    t7->Branch("phi_gen_had",&phi_gen_had,"phi_gen_had/D");
    t7->Branch("nsim",&nsim,"nsim/l");
    t7->Branch("dilute_p",&dilute_hp,"data/D");
    t7->Branch("dilute_m",&dilute_hm ,"data/D");
    t7->Branch("px_ele",&px_ele, "px_ele/D");
    t7->Branch("py_ele",&py_ele, "py_ele/D");
    t7->Branch("pz_ele",&pz_ele, "pz_ele/D");
    t7->Branch("E_ele",&E_ele, "E_ele/D");
    t7->Branch("px_had",&px_had, "px_had/D");
    t7->Branch("py_had",&py_had, "py_had/D");
    t7->Branch("pz_had",&pz_had, "pz_had/D");
    t7->Branch("E_had",&E_had, "E_had/D");
    t7->Branch("vx_ele",&vx_ele, "vx_ele/D");
    t7->Branch("vy_ele",&vy_ele, "vy_ele/D");
    t7->Branch("vz_ele",&vz_ele, "vz_ele/D");
    t7->Branch("vx_had",&vx_had, "vx_had/D");
    t7->Branch("vy_had",&vy_had, "vy_had/D");
    t7->Branch("vz_had",&vz_had, "vz_had/D");
    
    t7->Branch("u_pdf", &u_pdf, "u_pdf/D");
    t7->Branch("d_pdf", &d_pdf, "d_pdf/D");
    t7->Branch("s_pdf", &s_pdf, "s_pdf/D");
    t7->Branch("g_pdf", &g_pdf, "g_pdf/D");
    t7->Branch("ubar_pdf", &ubar_pdf, "ubar_pdf/D");
    t7->Branch("dbar_pdf", &dbar_pdf, "dbar_pdf/D");
    t7->Branch("sbar_pdf", &sbar_pdf, "sbar_pdf/D");
    t7->Branch("D_fav", &D_fav, "D_fav/D");
    t7->Branch("D_unfav", &D_unfav, "D_unfav/D");
    t7->Branch("D_s", &D_s, "D_s/D");
    t7->Branch("D_g", &D_g, "D_g/D");
   /*}}}*/

    TFile *file8= new TFile(new_filename[7],"RECREATE");
    TTree *t8 = new TTree("T","T");
    t8->SetDirectory(file8);
    t8->Branch("Q2",&Q2,"data/D");/*{{{*/
    t8->Branch("W",&W,"data/D");
    t8->Branch("Wp",&Wp,"data/D");
    t8->Branch("x",&x,"data/D");
    t8->Branch("y",&y,"data/D");
    t8->Branch("z",&z,"data/D");
    t8->Branch("nu",&nu,"data/D");
    t8->Branch("s",&s,"data/D");
    t8->Branch("pt",&pt,"data/D");
    t8->Branch("gamma",&gamma,"data/D");
    t8->Branch("epsilon",&epsilon,"data/D");
    t8->Branch("weight_hp",&weight_hp,"data/D");
    t8->Branch("weight_hm",&weight_hm,"data/D");
    t8->Branch("weight_in",&weight_in,"data/D");
    t8->Branch("rapidity",&rapidity,"data/D");
    t8->Branch("theta_q",&theta_q,"data/D");
    t8->Branch("theta_s",&theta_s,"data/D");
    t8->Branch("phi_h",&phi_h,"data/D");
    t8->Branch("phi_s",&phi_s,"data/D");
    t8->Branch("jacoF",&jacoF,"jacoF/D");
    t8->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
    t8->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
    t8->Branch("dxs_incl",&dxs_incl,"dxs_incl/D");
    t8->Branch("mom_ele",&mom_ele,"mom_ele/D");
    t8->Branch("mom_gen_ele",&mom_gen_ele,"mom_gen_ele/D");
    t8->Branch("mom_had",&mom_had,"mom_had/D");
    t8->Branch("mom_gen_had",&mom_gen_had,"mom_gen_had/D");
    t8->Branch("theta_ele",&theta_ele,"theta_ele/D");
    t8->Branch("theta_gen_ele",&theta_gen_ele,"theta_gen_ele/D");
    t8->Branch("theta_had",&theta_had,"theta_had/D");
    t8->Branch("theta_gen_had",&theta_gen_had,"theta_gen_had/D");
    t8->Branch("phi_ele",&phi_ele,"phi_ele/D");
    t8->Branch("phi_gen_ele",&phi_gen_ele,"phi_gen_ele/D");
    t8->Branch("phi_had",&phi_had,"phi_had/D");
    t8->Branch("phi_gen_had",&phi_gen_had,"phi_gen_had/D");
    t8->Branch("nsim",&nsim,"nsim/l");
    t8->Branch("dilute_p",&dilute_hp,"data/D");
    t8->Branch("dilute_m",&dilute_hm ,"data/D");
    t8->Branch("px_ele",&px_ele, "px_ele/D");
    t8->Branch("py_ele",&py_ele, "py_ele/D");
    t8->Branch("pz_ele",&pz_ele, "pz_ele/D");
    t8->Branch("E_ele",&E_ele, "E_ele/D");
    t8->Branch("px_had",&px_had, "px_had/D");
    t8->Branch("py_had",&py_had, "py_had/D");
    t8->Branch("pz_had",&pz_had, "pz_had/D");
    t8->Branch("E_had",&E_had, "E_had/D");
    t8->Branch("vx_ele",&vx_ele, "vx_ele/D");
    t8->Branch("vy_ele",&vy_ele, "vy_ele/D");
    t8->Branch("vz_ele",&vz_ele, "vz_ele/D");
    t8->Branch("vx_had",&vx_had, "vx_had/D");
    t8->Branch("vy_had",&vy_had, "vy_had/D");
    t8->Branch("vz_had",&vz_had, "vz_had/D");
    
    t8->Branch("u_pdf", &u_pdf, "u_pdf/D");
    t8->Branch("d_pdf", &d_pdf, "d_pdf/D");
    t8->Branch("s_pdf", &s_pdf, "s_pdf/D");
    t8->Branch("g_pdf", &g_pdf, "g_pdf/D");
    t8->Branch("ubar_pdf", &ubar_pdf, "ubar_pdf/D");
    t8->Branch("dbar_pdf", &dbar_pdf, "dbar_pdf/D");
    t8->Branch("sbar_pdf", &sbar_pdf, "sbar_pdf/D");
    t8->Branch("D_fav", &D_fav, "D_fav/D");
    t8->Branch("D_unfav", &D_unfav, "D_unfav/D");
    t8->Branch("D_s", &D_s, "D_s/D");
    t8->Branch("D_g", &D_g, "D_g/D");
    /*}}}*/
    /*}}}*/

    const double xbj_min = 0.0;
    const double xbj_max = 0.35;
    const int Q2bin=8;
    const double Q2_log[9] = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 3.0};
    //const int Q2bin=4;
    //const double Q2_log[5] = {0.0, 0.4, 0.8, 1.2, 1.6};
    double Q2_Cut[9];
    for (Int_t j=0;j<Q2bin+1;j++){
        //Q2_Cut[j] = pow(10., Q2_log[j]);
        Q2_Cut[j] = Q2_log[j];
    }

    for (ULong64_t i=0;i<N_Total;i++){
        T->GetEntry(i);
        if (x>=xbj_min&&x<xbj_max&&log10(Q2)>=Q2_Cut[0]&&log10(Q2)<Q2_Cut[1])  t1->Fill();
        if (x>=xbj_min&&x<xbj_max&&log10(Q2)>=Q2_Cut[1]&&log10(Q2)<Q2_Cut[2])  t2->Fill();
        if (x>=xbj_min&&x<xbj_max&&log10(Q2)>=Q2_Cut[2]&&log10(Q2)<Q2_Cut[3])  t3->Fill();
        if (x>=xbj_min&&x<xbj_max&&log10(Q2)>=Q2_Cut[3]&&log10(Q2)<Q2_Cut[4])  t4->Fill();
        if (x>=xbj_min&&x<xbj_max&&log10(Q2)>=Q2_Cut[4]&&log10(Q2)<Q2_Cut[5])  t5->Fill();
        if (x>=xbj_min&&x<xbj_max&&log10(Q2)>=Q2_Cut[5]&&log10(Q2)<Q2_Cut[6])  t6->Fill();
        if (x>=xbj_min&&x<xbj_max&&log10(Q2)>=Q2_Cut[6]&&log10(Q2)<Q2_Cut[7])  t7->Fill();
        if (x>=xbj_min&&x<xbj_max&&log10(Q2)>=Q2_Cut[7]&&log10(Q2)<Q2_Cut[8])  t8->Fill();

        if(!(i%10000))
            cerr<<Form("--- Working on evt=%d",i)<<"\r";
    }

    file1->Write();  file1->Close();
    file2->Write();  file2->Close();
    file3->Write();  file3->Close();
    file4->Write();  file4->Close();
    file5->Write();  file5->Close();
    file6->Write();  file6->Close();
    file7->Write();  file7->Close();
    file8->Write();  file8->Close();
}

