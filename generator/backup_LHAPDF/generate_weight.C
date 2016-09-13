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
  #include <TRandom3.h>
  //#include <TMatrix.h>
  /*}}}*/
using namespace std;

const double PI = 3.1415926;
const double DEG = 180.0/PI;

void generate_weight(TString inputfile){
	double momentum_ele =10.0; 

	TFile *file=new TFile(inputfile.Data(),"update");  // save Q2<10 pt<1
	TTree *T=(TTree*)file->Get("T");
	int nsim=0;
	double dxs_hp=0,dxs_hm=0;
	T->SetBranchAddress("nsim",&nsim);
	T->SetBranchAddress("dxs_hp",&dxs_hp);
	T->SetBranchAddress("dxs_hm",&dxs_hm);

	Long64_t N_entries=T->GetEntries();
	T->GetEntry(N_entries-1);          //get nsim for this rootfile
	
	double N_simulate=(double)(nsim);
	cout<<"N_simulate: "<<N_simulate<<endl;

    TString config = "EIC";

    //SoLID Acceptance
    const Double_t SoLID_Mom_Min_e = 0.5;  
    const Double_t SoLID_Mom_Max_e = 11.0;//not in use 
    const Double_t SoLID_Mom_Min_h = 0.5;  
    const Double_t SoLID_Mom_Max_h = 6.0;

    const Double_t SoLID_Th_Min_e  = 7.0;  
    const Double_t SoLID_Th_Max_e = 30.0; 
    const Double_t SoLID_Th_Min_h  =  7.0; 
    const Double_t SoLID_Th_Max_h =  30.0;

    //A rough guess but people claim EIC to be a full-acceptance device!
    const Double_t EIC_Mom_Min_e = 0.5;  
    const Double_t EIC_Mom_Max_e = 3.*10.0; //not in use 
    const Double_t EIC_Mom_Min_h = 0.0;  
    const Double_t EIC_Mom_Max_h = 50.0;

    const Double_t EIC_Th_Min_e = 0.0;  
    const Double_t EIC_Th_Max_e  = 180.0;            
    const Double_t EIC_Th_Min_h = 0.0;  
    const Double_t EIC_Th_Max_h  = 180.0;

    Double_t Mom_Max_e = 0.0, Mom_Min_e = 0.0, Mom_Max_h = 0.0,Mom_Min_h = 0.0;
    Double_t Th_Max_e = 0.0,Th_Min_e = 0.0,Th_Max_h = 0.0, Th_Min_h = 0.0;

    if(config=="SoLID" ){
        Mom_Min_e = SoLID_Mom_Min_e;  Mom_Max_e = momentum_ele; 
        Mom_Min_h = SoLID_Mom_Min_h;  Mom_Max_h = SoLID_Mom_Max_h;
        Th_Min_e = SoLID_Th_Min_e; Th_Max_e = SoLID_Th_Max_e; 
        Th_Min_h = SoLID_Th_Min_h; Th_Max_h = SoLID_Th_Max_h;
    }
    //A rough guess but people claim EIC to be a full-acceptance device!
    if(config=="EIC" ){
        Mom_Min_e = EIC_Mom_Min_e;  Mom_Max_e =  momentum_ele * 3.0; 
        Mom_Min_h = EIC_Mom_Min_h;  Mom_Max_h = EIC_Mom_Max_h;
        Th_Min_e = EIC_Th_Min_e; Th_Max_e = EIC_Th_Max_e; 
        Th_Min_h = EIC_Th_Min_h; Th_Max_h = EIC_Th_Max_h;    
    }

    double electron_phase_space=(cos(Th_Min_e/DEG) - cos(Th_Max_e/DEG))*2*PI*(Mom_Max_e - Mom_Min_e);
    double hadron_phase_space=(cos(Th_Min_h/DEG) - cos(Th_Max_h/DEG))*2*PI*(Mom_Max_h - Mom_Min_h);
    double Phase_space=electron_phase_space*hadron_phase_space;           //electron*hadron phase space eg, for electron: delta_cos_theta*delta_phi*delta_energy
	cout<<"Phase_space: "<<electron_phase_space<<"	"<<hadron_phase_space<<"	"<<Phase_space<<endl;
	double weight_hp=0;
	double weight_hm=0;

	TBranch *branch_weight_hp=T->Branch("weight_hp",&weight_hp,"weight_hp/D");
	TBranch *branch_weight_hm=T->Branch("weight_hm",&weight_hm,"weight_hm/D");

	for(Long64_t i=0;i<N_entries;i++){
		T->GetEntry(i);
		//warning: output unit is nbarn   //if calculate rate, should be translate to cm^-2     1nbarn=10^-33 cm^-2
		weight_hp=dxs_hp*Phase_space/N_simulate;   
		weight_hm=dxs_hm*Phase_space/N_simulate;
		branch_weight_hp->Fill();
		branch_weight_hm->Fill();
	}

	T->Write("",TObject::kOverwrite);
}



