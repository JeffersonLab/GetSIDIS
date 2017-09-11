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
    int number_of_files = 1; cerr<<"How many files with 1M events? " ; cin >> number_of_files;
    if (number_of_files !=50 && number_of_files !=500 && number_of_files!=1 && number_of_files!=5) {
      cout << "This skim is not available. Number of input file is set to 1" << endl;
      number_of_files = 1;
    }


    const Double_t EIC_Mom_Min_e = 8.5;
    const Double_t EIC_Mom_Max_e = 10.5;
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

    Double_t generated_electron_phase_space =(cos(EIC_Th_Min_e*degtorad) - cos(EIC_Th_Max_e*degtorad))*(EIC_Ph_Max_e*degtorad - EIC_Ph_Min_e*degtorad)*(EIC_Mom_Max_e - EIC_Mom_Min_e);
    Double_t generated_hadron_phase_space   =(cos(EIC_Th_Min_h*degtorad) - cos(EIC_Th_Max_h*degtorad))*(EIC_Ph_Max_h*degtorad - EIC_Ph_Min_h*degtorad)*(EIC_Mom_Max_h - EIC_Mom_Min_h);
    Double_t generated_phase_space=generated_electron_phase_space*generated_hadron_phase_space;
    cout << "generated PS: " << generated_phase_space << endl;

    const double luminosity = 1e33; // unit 1/(cm^2 s)
    const double nbarntocm2 = 1e-33;

    double integrated_lumi1 = 1e6; // 1/fb = 10^6/nb
    double integrated_lumi2 = 1e5; // 1/10fb = 10^5/nb

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

    if (fA == 2) histotitle = "Weighted total cross section d(e,e'#pi^{+})X (MC error); z ; Q2 (GeV^{2})";
    if (fA == 12) histotitle = "Weighted total cross section ^{12}C(e,e'#pi^{+})X (MC error); z ; Q2 (GeV^{2})";
    TH2D *pionp_cross = new TH2D("pionp_cross ",histotitle,zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
    pionp_cross->GetXaxis()->CenterTitle(1);
    pionp_cross->GetYaxis()->CenterTitle(1);
    pionp_cross->SetMarkerSize(2);
    pionp_cross->Sumw2();
    histotitle.Clear();

    if (fA == 2) histotitle = "Weighted total cross section d(e,e'#pi^{-})X (MC error); z ; Q2 (GeV^{2})";
    if (fA == 12) histotitle = "Weighted total cross section ^{12}C(e,e'#pi^{-})X (MC error); z ; Q2 (GeV^{2})";
    TH2D *pionm_cross = new TH2D("pionm_cross",histotitle,zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
    pionm_cross->GetXaxis()->CenterTitle(1);
    pionm_cross->GetYaxis()->CenterTitle(1);
    pionm_cross->SetMarkerSize(2);
    pionm_cross->Sumw2();
    histotitle.Clear();

    if (fA == 2) histotitle = "Cross section difference e+d (MC error); z ; Q2 (GeV^{2})";
    if (fA == 12) histotitle = "Cross section difference e+^{12}C (MC error); z ; Q2 (GeV^{2})";
    TH2D *crosssectiondiff = new TH2D("crosssectiondiff",histotitle,zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
    crosssectiondiff->GetXaxis()->CenterTitle(1);
    crosssectiondiff->GetYaxis()->CenterTitle(1);
    crosssectiondiff->Sumw2();
    crosssectiondiff->SetMarkerSize(2);

    if (fA == 2) histotitle = "Count rate for luminosity 1*fb^{-1},  d(e,e'#pi^{+})X (statistical error); z ; Q2 (GeV^{2})";
    if (fA == 12) histotitle = "Count rate for luminosity 1*fb^{-1}, ^{12}C(e,e'#pi^{+})X (statistical error); z ; Q2 (GeV^{2})";
    TH2D *pionp_counts1 = new TH2D("pionp_counts1 ",histotitle,zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
    pionp_counts1->GetXaxis()->CenterTitle(1);
    pionp_counts1->GetYaxis()->CenterTitle(1);
    pionp_counts1->SetMarkerSize(2);
    pionp_counts1->Sumw2();
    histotitle.Clear();

    if (fA == 2) histotitle = "Count rate for luminosity 0.1*fb^{-1},  d(e,e'#pi^{+})X (statistical error); z ; Q2 (GeV^{2})";
    if (fA == 12) histotitle = "Count rate for luminosity 0.1*fb^{-1}, ^{12}C(e,e'#pi^{+})X (statistical error); z ; Q2 (GeV^{2})";
    TH2D *pionp_counts2 = new TH2D("pionp_counts2 ",histotitle,zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
    pionp_counts2->GetXaxis()->CenterTitle(1);
    pionp_counts2->GetYaxis()->CenterTitle(1);
    pionp_counts2->SetMarkerSize(2);
    pionp_counts2->Sumw2();
    histotitle.Clear();

    if (fA == 2) histotitle = "Count rate for luminosity 1*fb^{-1},  d(e,e'#pi^{-})X (statistical error); z ; Q2 (GeV^{2})";
    if (fA == 12) histotitle = "Count rate for luminosity 1*fb^{-1}, ^{12}C(e,e'#pi^{-})X (statistical error); z ; Q2 (GeV^{2})";
    TH2D *pionm_counts1 = new TH2D("pionm_counts1 ",histotitle,zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
    pionm_counts1->GetXaxis()->CenterTitle(1);
    pionm_counts1->GetYaxis()->CenterTitle(1);
    pionm_counts1->SetMarkerSize(2);
    pionm_counts1->Sumw2();
    histotitle.Clear();

    if (fA == 2) histotitle = "Count rate for luminosity 0.1*fb^{-1},  d(e,e'#pi^{-})X (statistical error); z ; Q2 (GeV^{2})";
    if (fA == 12) histotitle = "Count rate for luminosity 0.1*fb^{-1}, ^{12}C(e,e'#pi^{-})X (statistical error); z ; Q2 (GeV^{2})";
    TH2D *pionm_counts2 = new TH2D("pionm_counts2 ",histotitle,zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
    pionm_counts2->GetXaxis()->CenterTitle(1);
    pionm_counts2->GetYaxis()->CenterTitle(1);
    pionm_counts2->SetMarkerSize(2);
    pionm_counts2->Sumw2();
    histotitle.Clear();

    if (fA == 2)  histotitle = "Count rate difference for luminosity 1*fb^{-1} e+d (statistical error); z ; Q2 (GeV^{2})";
    if (fA == 12) histotitle = "Count rate difference for luminosity 1*fb^{-1} e+^{12}C (statistical) error); z ; Q2 (GeV^{2})";
    TH2D *countratediff1 = new TH2D("countratediff1",histotitle,zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
    countratediff1->GetXaxis()->CenterTitle(1);
    countratediff1->GetYaxis()->CenterTitle(1);
    countratediff1->Sumw2();
    countratediff1->SetMarkerSize(2);

    if (fA == 2)  histotitle = "Statistical relative error count rate difference (luminosity 1*fb^{-1}, e+d); z ; Q2 (GeV^{2})";
    if (fA == 12) histotitle = "Statistical relative error count rate difference (luminosity 1*fb^{-1}, e+^{12}C ; z ; Q2 (GeV^{2})";
    TH2D *relerror_countratediff1 = new TH2D("relerror_countratediff1",histotitle,zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
    relerror_countratediff1->GetXaxis()->CenterTitle(1);
    relerror_countratediff1->GetYaxis()->CenterTitle(1);
    relerror_countratediff1->Sumw2();
    relerror_countratediff1->SetMarkerSize(2);

    if (fA == 2)  histotitle = "Count rate difference for luminosity 0.1*fb^{-1} e+d (statistical error); z ; Q2 (GeV^{2})";
    if (fA == 12) histotitle = "Count rate difference for luminosity 0.1*fb^{-1} e+^{12}C (statistical) error); z ; Q2 (GeV^{2})";
    TH2D *countratediff2 = new TH2D("countratediff2",histotitle,zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
    countratediff2->GetXaxis()->CenterTitle(1);
    countratediff2->GetYaxis()->CenterTitle(1);
    countratediff2->Sumw2();
    countratediff2->SetMarkerSize(2);

    if (fA == 2)  histotitle = "Statistical relative error count rate difference (luminosity 0.1*fb^{-1}, e+d); z ; Q2 (GeV^{2})";
    if (fA == 12) histotitle = "Statistical relative error count rate difference (luminosity 0.1*fb^{-1}, e+^{12}C ; z ; Q2 (GeV^{2})";
    TH2D *relerror_countratediff2 = new TH2D("relerror_countratediff2",histotitle,zbin,z_cut[0],z_cut[zbin],Q2bin,Q2_cut[0],Q2_cut[Q2bin]);
    relerror_countratediff2->GetXaxis()->CenterTitle(1);
    relerror_countratediff2->GetYaxis()->CenterTitle(1);
    relerror_countratediff2->Sumw2();
    relerror_countratediff2->SetMarkerSize(2);

    for(int i=0;i<zbin;i++) {
      for(int j=0;j<Q2bin;j++) {
        //Define
        Double_t Q2, W, Wp, x, y, z, pt, nu, s, gamma, epsilon,rapidity, jacoF;
        Double_t mom_gen_ele,mom_gen_had;
        Double_t theta_gen_ele,theta_gen_had;
        Double_t phi_gen_ele,phi_gen_had;
        Double_t theta_q, theta_s,phi_h,phi_s,mom_ele,mom_had,theta_ele, theta_had,phi_ele,phi_had;
        Double_t dxs_incl,dxs_hm,dxs_hp,dilute_hp,dilute_hm;
        Double_t px_ele, py_ele,pz_ele, px_had, py_had, pz_had, E_ele,E_had;
        Double_t weight_hp, weight_hm, weight_in;
        ULong64_t nsim = 0;


        TString filename;
        if(fA ==2){
          if (number_of_files == 500)     filename= Form("fulldataskim/EIC_A2_pion_10_100_skim_Q2bin%d_zbin%d.root",j+1,i+1);
          else if (number_of_files == 50) filename= Form("50Mdataskim/EIC_A2_pion_10_100_skim_Q2bin%d_zbin%d.root",j+1,i+1);
          else if (number_of_files == 5) filename = "../../test/massproduction_CTEQ/EIC_A2_pion_10_100_1_total5.root";
          else     {                    //   filename = "../../test/massproduction_CTEQ/EIC_A2_pion_10_100_1_0.root";
            //  filename = "EIC_A2_pion_10_100_1_2.root";
             filename = "../massproduction_CTEQfree/EIC_A2_pion_10_100_1_2.root";
          }
        }
        if(fA == 12){
          if (number_of_files == 500)     filename= Form("fulldataskim/EIC_A12_pion_10_600_skim_Q2bin%d_zbin%d.root",j+1,i+1);
	        else if (number_of_files == 50) filename= Form("50Mdataskim/EIC_A12_pion_10_600_skim_Q2bin%d_zbin%d.root",j+1,i+1);
          else if (number_of_files == 5)  filename = "../../test/massproduction_CTEQ/EIC_A12_pion_10_600_1_total5.root";
          else                         //   filename = "../../test/massproduction_CTEQ/EIC_A12_pion_10_600_1_0.root";
           filename = "../massproduction_CTEQfree/EIC_A12_pion_10_600_1_2.root";
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

        TString pipluscut = Form("(1/%i)*(weight_hp)*(Q2>=%f && Q2<%f && z>=%f && z<%f)",number_of_files,Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);
        TString piminuscut = Form("(1/%i)*(weight_hm)*(Q2>=%f && Q2<%f && z>=%f && z<%f)",number_of_files,Q2_cut[j],Q2_cut[j+1],z_cut[i],z_cut[i+1]);
        T->Draw("Q2>>htemp",TCut(pipluscut),"");
        T->Draw("Q2>>htemp2",TCut(piminuscut),"");
      /* double man_errorintegral_p = 0;
       double man_errorintegral_m = 0;
       double man_sumcross_p = 0;
       double man_sumcross2_p = 0;
       double man_sumcross_m = 0;
       double man_sumcross2_m = 0;
       double test_m = 0;
       double test2_m = 0;
       ULong64_t nsim_real = 0;

       for (int event=0; event<counts_in_bin; event++) {
             T->GetEntry(event);
	    if (weight_hp!=0 && Q2>=Q2_cut[j] && Q2<Q2_cut[j+1] && z>=z_cut[i] && z<z_cut[i+1]) {
              htemp->Fill(Q2,weight_hp);
	      nsim_real = dxs_hp*generated_phase_space/weight_hp;
	      man_sumcross_p+= dxs_hp/pow(nsim_real,1.5);
	      man_sumcross2_p+=dxs_hp*dxs_hp/pow(nsim_real,2);
	      test_m += dxs_hp;
	      test2_m += dxs_hp*dxs_hp;
	     // cout << "Event " << event << " nsim is " << nsim << " . d*PS/weight+ " << dxs_hp*generated_phase_space/weight_hp << endl;
            }
            if (weight_hm!=0 && Q2>=Q2_cut[j] && Q2<Q2_cut[j+1] && z>=z_cut[i] && z<z_cut[i+1]) {
              htemp2->Fill(Q2,weight_hm);
	      nsim_real = dxs_hm*generated_phase_space/weight_hm;
	      man_sumcross_m+= dxs_hm/pow(nsim_real,1.5);
	      man_sumcross2_m+=dxs_hm*dxs_hm/pow(nsim_real,2);
	     // cout << "Event " << event << " nsim is " << nsim << " . d*PS/weight- " << dxs_hm*generated_phase_space/weight_hm << endl;
            }

	  //  if(!(event%100000)) cout << "Event " << event << " nsim1 is " << nsim_real << " . Test values " << test_m << " and " << test2_m << endl;
       }

       man_errorintegral_p = generated_phase_space * sqrt(man_sumcross2_p - pow (man_sumcross_p, 2) );
       man_errorintegral_m = generated_phase_space * sqrt(man_sumcross2_m - pow (man_sumcross_m, 2) );
        */
       Double_t histsum = (double) htemp->GetSum();
       Double_t histsum2 = (double) htemp2->GetSum();

       Double_t histerror = sqrt(htemp->GetSumw2()->GetSum());
       Double_t hist2error = sqrt(htemp2->GetSumw2()->GetSum());
       cout << "BIN " << i << " and " << j << " ; pi+ sum " << histsum << " , pi+ error " << histerror << " , relative error " << histerror/histsum <<  endl;
       cout << "BIN " << i << " and " << j << " ; pi- sum " << histsum2 << " , pi- error " << hist2error << " , relative error " << hist2error/histsum2 << endl;
     //  cout << "BIN " << i << " and " << j << " ; pi+ sum " << histsum << " , pi+ error " << histerror << " , relative error " << histerror/histsum << " , man error " << man_errorintegral_p<< endl;
     //  cout << "BIN " << i << " and " << j << " ; pi- sum " << histsum2 << " , pi- error " << hist2error << " , relative error " << hist2error/histsum2 << " , man error " << man_errorintegral_m<< endl;

       pionp_cross->SetBinContent(i+1,j+1,histsum);
       pionm_cross->SetBinContent(i+1,j+1,histsum2);

       pionp_cross->SetBinError(i+1,j+1,histerror);
       pionm_cross->SetBinError(i+1,j+1,hist2error);


       Double_t counts1_pion_p = histsum*luminosity*nbarntocm2*integrated_lumi1;
       Double_t counts2_pion_p = histsum*luminosity*nbarntocm2*integrated_lumi2;
       Double_t counts1_pion_m = histsum2*luminosity*nbarntocm2*integrated_lumi1;
       Double_t counts2_pion_m = histsum2*luminosity*nbarntocm2*integrated_lumi2;
       pionp_counts1->SetBinContent(i+1,j+1,counts1_pion_p);
       pionp_counts2->SetBinContent(i+1,j+1,counts2_pion_p);
       pionm_counts1->SetBinContent(i+1,j+1,counts1_pion_m);
       pionm_counts2->SetBinContent(i+1,j+1,counts2_pion_m);

       pionp_counts1->SetBinError(i+1,j+1,sqrt(counts1_pion_p));
       pionp_counts2->SetBinError(i+1,j+1,sqrt(counts2_pion_p));
       pionm_counts1->SetBinError(i+1,j+1,sqrt(counts1_pion_m));
       pionm_counts2->SetBinError(i+1,j+1,sqrt(counts2_pion_m));


       delete htemp;
       delete htemp2;

      }
    }

    crosssectiondiff->Add(pionp_cross,pionm_cross,1,-1);
    countratediff1->Add(pionp_counts1,pionm_counts1,1,-1);
    countratediff2->Add(pionp_counts2,pionm_counts2,1,-1);
    for(int i=0;i<zbin;i++) {
      for(int j=0;j<Q2bin;j++) {
        relerror_countratediff1->SetBinContent(i+1,j+1,countratediff1->GetBinError(i+1,j+1)/countratediff1->GetBinContent(i+1,j+1));
        relerror_countratediff2->SetBinContent(i+1,j+1,countratediff2->GetBinError(i+1,j+1)/countratediff2->GetBinContent(i+1,j+1));

      }
    }

    TFile *output;
    if (fA == 2) output= new TFile("perbinplots_d.root","RECREATE");
    if (fA == 12) output= new TFile("perbinplots_c12.root","RECREATE");
    output->cd();

    histo_counts->Write();
    pionp_cross->Write();
    pionm_cross->Write();
    crosssectiondiff->Write();

    pionp_counts1->Write();
    pionp_counts2->Write();
    pionm_counts1->Write();
    pionm_counts2->Write();
    countratediff1->Write();
    countratediff2->Write();

    relerror_countratediff1->Write();
    relerror_countratediff2->Write();

    output->Close();


}
