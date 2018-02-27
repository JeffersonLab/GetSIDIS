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
#include <TH2D.h>
#include <TH3D.h>
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

#include "acceptance/SIDIS_Acceptance.h"

using namespace std;
const double DEG = 180./3.1415926;
const double PI = 3.1415926;

int main(Int_t argc, char *argv[]){
    TString prefix = "./rootfiles/SoLID_";
    TString new_prefix = "./skim_rootfiles/skim_";

    TString target = "X"; cerr<<"-- What target (He3, NH3 (p))? "; cin >> target;
    TString particle = "X"; cerr<<"-- What particle (pion, kaon)? "; cin >> particle;
    int Ebeam = 0; cerr<<"-- What beam energy (11 or 8.8 GeV)? "; cin >> Ebeam;
    TString TA="X";
    if(target=="He3") TA = "A3";
    if(target=="NH3") TA = "A1";

    Int_t start_num  = 0;
    Int_t end_num  = 2;
    Int_t zflag = 0;
    Int_t Q2flag = 0;

    TString posfix,new_filename;
    TChain *T0 = new TChain("T","T");
    for (Int_t i=start_num; i<end_num;i++){
        posfix.Form("_%d_0_1_%d.root",int(Ebeam), i);
        new_filename = prefix + TA +"_"+ particle+posfix;
        cerr<<Form(" @@@ Adding Root File: %s", new_filename.Data())<<endl;
        T0->AddFile(new_filename);
    
        posfix.Form("_%d_0_2_%d.root",int(Ebeam), i);
        new_filename = prefix + TA +"_"+ particle+posfix;
        cerr<<Form(" @@@ Adding Root File: %s", new_filename.Data())<<endl;
        T0->AddFile(new_filename);
    }

    /*Define{{{*/
    Double_t Q2, W, Wp, x, y, z, pt, nu, s, gamma, epsilon,rapidity, jacoF;
    //Double_t theta_gen= 0.0 , phi_gen = 0.0, mom_gen = 0.0;
    //Double_t mom_gen_ele,mom_gen_had;
    //Double_t theta_gen_ele,theta_gen_had;
    //Double_t phi_gen_ele,phi_gen_had;
    Double_t theta_q, theta_s,phi_h,phi_s,mom_ele,mom_had,theta_ele, theta_had,phi_ele,phi_had;
    Double_t dxs_incl,dxs_hm,dxs_hp,dxs_hm_sidis,dxs_hp_sidis,dilute_hp,dilute_hm;
    Double_t px_ele, py_ele,pz_ele, px_had, py_had, pz_had, E_ele,E_had;
    Double_t u_pdf, d_pdf, s_pdf, g_pdf, ubar_pdf, dbar_pdf, sbar_pdf;
    Double_t weight_hp, weight_hm, weight_in;
    ULong64_t nsim = 0;// Nsim1 = 0, Nsim2 = 0, Nsim3 = 0,Nsim4 = 0;
    //For Beam Position and Vertex info
    Double_t vx_ele, vy_ele, vz_ele, vx_had, vy_had, vz_had;
    Double_t D_fav, D_unfav, D_s, D_g;
    Int_t isphy_hp, isphy_hm;
    /*}}}*/
 
	/*Define old root file{{{*/
    T0->SetBranchAddress("Q2",&Q2);
    T0->SetBranchAddress("W",&W);
    T0->SetBranchAddress("Wp",&Wp);
    T0->SetBranchAddress("x",&x );
    T0->SetBranchAddress("y",&y );
    T0->SetBranchAddress("z",&z );
    T0->SetBranchAddress("nu",&nu );
    T0->SetBranchAddress("s",&s );
    T0->SetBranchAddress("epsilon",&epsilon );
    T0->SetBranchAddress("gamma",&gamma );
    T0->SetBranchAddress("pt",&pt );
    T0->SetBranchAddress("isphy_hp",&isphy_hp );
    T0->SetBranchAddress("isphy_hm",&isphy_hm );
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
    T0->SetBranchAddress("dxs_hm_sidis",&dxs_hm_sidis);
    T0->SetBranchAddress("dxs_hp_sidis",&dxs_hp_sidis);
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
    /*}}}*/

    ULong64_t N_Total=T0->GetEntries();
    cout<<"--- Total Number of Events: "<<N_Total<<endl;

    SIDIS_Acceptance* accept = new SIDIS_Acceptance();
    accept->Init(target);

    Double_t zmin=1000,zmax=-1000;
    Double_t Q2min=1000,Q2max=-1000.;
    Double_t acc_f_ele = 0.0, acc_l_ele = 0.0;
    Double_t acc_f_hp = 0.0, acc_l_hp = 0.0;
    Double_t acc_f_hm = 0.0, acc_l_hm = 0.0;
    Double_t W_hp = 0.0, W_hm = 0.0;
    Double_t luminosity = 0.0, time = 0.0;
    Double_t rate_hp =0.0, rate_hm=0.0;

    for (Int_t i=0;i!=T0->GetEntries();i++){
        T0->GetEntry(i);
        theta_ele *= DEG; theta_had *= DEG;
        phi_ele *= DEG;   phi_had *= DEG;

        /*Get acceptance of e and pi-{{{*/
        //Make sure to use the corrected quantities for multipile scattering and eloss effects
        //Do not use the smeared quantities since we are about whether particles are in the accepntace or not, but not how good we measure
        /*Elec Acc {{{*/
        if(target=="He3"){
            acc_f_ele = accept->GetAcc("e-","forward", mom_ele, theta_ele);
            acc_l_ele = accept->GetAcc("e-","large", mom_ele, theta_ele);
        }
        else if(target=="NH3"){
            acc_f_ele = accept->GetAcc3D("e-","forward", mom_ele, theta_ele, phi_ele);
            acc_l_ele = accept->GetAcc3D("e-","large", mom_ele, theta_ele, phi_ele);
        }
        //if(mom_ele<1.0||theta_ele>14.8||theta_ele<8.0)//GeV, CLEO
        //acc_f_ele=0.0;//Farward-Angle EC Cut at 1 GeV
        //if(mom_ele<3.5||theta_ele<16.0||theta_ele>24)//GeV,CLEO
        //acc_l_ele=0.0; //Larger-Angle EC Cut at 3 GeV
        if(acc_f_ele>1.) 
            acc_f_ele=1.0; 
        if(acc_l_ele>1.) 
            acc_l_ele=1.0; 
        /*}}}*/

        /*Hadron Acc{{{*/
        if(particle=="pion"){
            if(target=="He3"){
                acc_f_hp = accept->GetAcc("pi+","forward", mom_had, theta_had);
                acc_f_hm = accept->GetAcc("pi-","forward", mom_had, theta_had);
                acc_l_hp = accept->GetAcc("pi+","large",   mom_had, theta_had);
                acc_l_hm = accept->GetAcc("pi-","large",   mom_had, theta_had);
            }
            else if(target=="NH3"){
                acc_f_hp = accept->GetAcc3D("pi+","forward", mom_had, theta_had, phi_had);
                acc_f_hm = accept->GetAcc3D("pi-","forward", mom_had, theta_had, phi_had);
                acc_l_hp = accept->GetAcc3D("pi+","large",   mom_had, theta_had, phi_had);
                acc_l_hm = accept->GetAcc3D("pi-","large",   mom_had, theta_had, phi_had);
            }
        }
        else if(particle=="kaon"){
            if(target=="He3"){
                acc_f_hp = accept->GetAcc("K+","forward", mom_had, theta_had);
                acc_f_hm = accept->GetAcc("K-","forward", mom_had, theta_had);
                acc_l_hp = accept->GetAcc("K+","large",   mom_had, theta_had);
                acc_l_hm = accept->GetAcc("K-","large",   mom_had, theta_had);
            }
            else if(target=="NH3"){
                acc_f_hp = accept->GetAcc3D("K+","forward", mom_had, theta_had, phi_had);
                acc_f_hm = accept->GetAcc3D("K-","forward", mom_had, theta_had, phi_had);
                acc_l_hp = accept->GetAcc3D("K+","large",   mom_had, theta_had, phi_had);
                acc_l_hm = accept->GetAcc3D("K-","large",   mom_had, theta_had, phi_had);
            }
        }
        else
            cerr<<"*** ERROR, I don't know the PID in skim.C, "<<particle.Data()<<endl;

        //if(theta_had>14.8||theta_had<8.0||mom_had<0.||mom_had>Ebeam){//GeV, CLEO
        //acc_f_hp=0.0; acc_f_hm=0.0;
        //}
        //if(theta_had<16.0||theta_had>24.0||mom_had<0.||mom_had>Ebeam){//GeV, CLEO
        //acc_l_hp=0.0; acc_l_hm=0.0;
        //}
        if(acc_f_hp>1.) {acc_f_hp=0.0;}
        if(acc_f_hm>1.) {acc_f_hm=0.0;}
        if(acc_l_hp>1.) {acc_l_hp=0.0;}
        if(acc_l_hm>1.) {acc_l_hm=0.0;}

        //Add the TOF Cut off here, assuming 30ps timing resolutions
        if(particle=="kaon"&&mom_had>10.0){//GeV,FIX HERE, need to conform this cut-off value
            acc_f_hp=0.0; 
        }
        if(particle=="kaon"&&mom_had>10.0){//GeV,FIX HERE, need to conform this cut-off value
            acc_l_hp=0.0; 
        }
        /*}}}*/

        /*}}}*/

        /*Get Luminosity and Beam Time{{{*/
        //1e36cm^-1*s^-1 for he3, 1e-33 is for nbar->cm
        luminosity = 1e36 * 1e-33;
        double day = 0; // days
        // nevents = dxs (nbar) * L (nucleons/cm^2/s) * T(days) 
        // 1 day = 24hr*3600s = 86400
        // nbar=10-9 barn = 10^-9*10^-28 m^2=10^-9*10^-28 *10^4 cm^2=10^-33 cm^2
        // so weight  = Lumi * nbar * 86400 *day *acc_ele*acc_had / nsim;
        if(fabs(Ebeam-11.0)<0.1)
            day  = 48; // days
        if(fabs(Ebeam-8.80)<0.1)
            day  = 21; // days
        time = day * 24 *3600; //sec
        /*}}}*/

        //Use Tianbo Liu's SIDIS XS
        if(dxs_hp>1e-33&&isphy_hp)
            weight_hp *= dxs_hp_sidis/dxs_hp/end_num;
        else
            weight_hp = 1e-33;

        if(dxs_hm>1e-33&&isphy_hm)
            weight_hm *= dxs_hm_sidis/dxs_hm/end_num;
        else
            weight_hm = 1e-33;

        if(Q2>1.0&&W>2.3&&Wp>1.6&&z>0.3&&z<0.7){
            rate_hp += weight_hp * luminosity * (acc_f_ele+acc_l_ele)*(acc_f_hp);
            rate_hm += weight_hm * luminosity * (acc_f_ele+acc_l_ele)*(acc_f_hm);
        }
    }
delete T0;

cout<<"--- Hadron+: Rate = "<<rate_hp<<endl;
cout<<"--- Hadron-: Rate = "<<rate_hm<<endl;

return 0;
}
