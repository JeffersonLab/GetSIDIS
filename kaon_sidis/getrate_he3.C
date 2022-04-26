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
  #include <TPaletteAxis.h>
  //#include <TMatrix.h>
  /*}}}*/
using namespace std;

/*
//_________________acceptance file__________
acceptance_solid_CLEO_SIDIS_3he_negative_output.root
*/

//void generate_rato(TString particle,int hadron_flag, int energy){
int main(){
	//particle: pion    kaon
	//hadron_flag   1 for hp      2 for hm

	TString particle; cerr<<"--- Particle (pion,kaon) = "; cin >> particle;
	int hadron_flag; cerr<<"--- Hadron Type (1->plus,2->minus) = "; cin >>hadron_flag;
	int ele_energy; cerr<<"--- Energy (8->8.8GeV,11->11GeV) = "; cin >>ele_energy;

	TString charge = "A";
	if(hadron_flag==1)
		charge = "p";
	else if(hadron_flag==2)
		charge= "m";
	else
		cerr<<"********* Error"<<endl;

	gStyle->SetOptStat(1);
	TChain *T=new TChain("T");
    TString data_dir = "/work/halla/solid/yez/sidis/sidis_rate/";
	if(particle=="pion"){
		if(hadron_flag==1){
			T->Add(Form("%s/sidis_3he_pip_%d_0_1_0.root",data_dir.Data(),ele_energy));   // save Q2<10, pt<1
			T->Add(Form("%s/sidis_3he_pip_%d_0_2_0.root",data_dir.Data(),ele_energy));   // save Q2<10, pt<1
		}
		else if(hadron_flag==2){
			T->Add(Form("%s/sidis_3he_pim_%d_0_1_0.root",data_dir.Data(),ele_energy));   // save Q2<10, pt<1
			T->Add(Form("%s/sidis_3he_pim_%d_0_2_0.root",data_dir.Data(),ele_energy));   // save Q2<10, pt<1
		}
	}else if(particle=="kaon"){
		if(hadron_flag==1){
//			T->Add(Form("./sidis_3he_kp_%d_0_1_0.root",ele_energy));   // save Q2<10, pt<1
//			T->Add(Form("./sidis_3he_kp_%d_0_2_0.root",ele_energy));   // save Q2<10, pt<1
			T->Add(Form("../sidis_model/sidis_3he_kp_%d_0_1_0.root",ele_energy));   // save Q2<10, pt<1
			T->Add(Form("../sidis_model/sidis_3he_kp_%d_0_2_0.root",ele_energy));   // save Q2<10, pt<1
		}
		if(hadron_flag==2){
//			T->Add(Form("./sidis_3he_km_%d_0_1_0.root",ele_energy));   // save Q2<10, pt<1
//			T->Add(Form("./sidis_3he_km_%d_0_2_0.root",ele_energy));   // save Q2<10, pt<1
			T->Add(Form("../sidis_model/sidis_3he_km_%d_0_1_0.root",ele_energy));   // save Q2<10, pt<1
			T->Add(Form("../sidis_model/sidis_3he_km_%d_0_2_0.root",ele_energy));   // save Q2<10, pt<1
			}
	}

	//useful variables: Q2, W, Wp, x,y,z nu, s, pt, phi_h, phi_s weight_hp, weight_hm   
	//weight is xs in nbarn  weight_hp= dxs_hp*cos_theta_coverage*phi_coverage*energy_coverage/N_simulate
	const double DEG=180./3.1415926;
	double Q2,x,pt,W,Wp,z;   //kinematics
	double mom_ele,mom_had,theta_ele,theta_had;  // acceptance calculation input
	double weight_hp,weight_hm,dxs_hp,dxs_hm;
	int nsim=0;

	T->SetBranchAddress("Q2",&Q2);
	T->SetBranchAddress("x",&x);
	T->SetBranchAddress("pt",&pt);
	T->SetBranchAddress("W",&W);
	T->SetBranchAddress("Wp",&Wp);
	T->SetBranchAddress("z",&z);
	T->SetBranchAddress("mom_ele",&mom_ele);
	T->SetBranchAddress("mom_had",&mom_had);
	T->SetBranchAddress("theta_ele",&theta_ele);
	T->SetBranchAddress("theta_had",&theta_had);
//	T->SetBranchAddress("weight_hp",&weight_hp);   //unit in nbar
//	T->SetBranchAddress("weight_hm",&weight_hm);
	T->SetBranchAddress("dxs_hp",&dxs_hp);   //unit in nbar
	T->SetBranchAddress("dxs_hm",&dxs_hm);
	T->SetBranchAddress("nsim",&nsim);

	Long64_t N_entries=T->GetEntries();
	cout<<"total generated events number: "<<N_entries<<endl;

	//deal with acceptance files 
    //CLEO	
	TFile *file_negative=new TFile("/work/halla/solid/yez/dvcs/dvcs_gen/acceptance/acceptance_solid_CLEO_SIDIS_3he_negative_output.root");
	TH2F *accep_ele_forward=(TH2F*)file_negative->Get("acceptance_forwardangle");       //negative particle forward angle acceptance
	TH2F *accep_ele_large=(TH2F*)file_negative->Get("acceptance_largeangle");           //negative particle large angle acceptance
	TH2F *accep_had=(TH2F*)accep_ele_forward->Clone("acceptance_hadron");//pion is normal acceptance histogram, the same as electron acceptance
	double rate_forward=0;
	double rate_large=0;

	double beam_energy=0.0;
	if(ele_energy==11)
		beam_energy = 11.0;
	else if(ele_energy==8)
		beam_energy = 8.8;
    else{
		cerr<<"*ERROR*"<<endl;
        return -1;
	}
	T->GetEntry(N_entries-1);          //get nsim for this rootfile
	double N_simulate=(double)(nsim);
	cout<<"N_simulate: "<<N_simulate<<endl;
	double electron_phase_space=(cos(7/180.*3.1415926) - cos(30/180.*3.1415926))*2*3.14159265*(beam_energy-0.5);   // theta: 7~30 degree,  2pi phi coverage, 0.5~11 GeV Momentum coverage 	
	double hadron_phase_space=(cos(7/180.*3.1415926) - cos(30/180.*3.1415926))*2*3.14159265*(6-0.5);  //theta, 7~30 degree,  2pi phi coverage, 0.5~6 GeV Momentum coverage
	double Phase_space=electron_phase_space*hadron_phase_space;           //electron*hadron phase space eg, for electron: delta_cos_theta*delta_phi*delta_energy
	cout<<"Phase_space: "<<electron_phase_space<<"	"<<hadron_phase_space<<"	"<<Phase_space<<endl;
	
	for(Long64_t i=0;i<N_entries;i++){
		T->GetEntry(i);
		if(Q2>=1&&mom_ele>1.0 && mom_ele<11.0 && W>=2.0&&theta_had*DEG>7.5&&theta_had*DEG<24.5&&theta_ele*DEG>7.5&&theta_ele*DEG<24.5){// && Wp>=1.6 && z>0.3 && z<0.7){ //zhihong's cut
		//if(mom_ele>=1.0&&mom_ele<=10.0&&W>=2.3&&Wp>=1.6&&z>0.3&&z<0.7&&mom_had<=2.5){//any additional cuts should be added in here
			int ele_theta_bin=int(theta_ele*DEG/0.2)+1;    //0.2 degree per bin
			int ele_p_bin=int(mom_ele/0.05)+1;      //0.05 GeV per bin for mom
			double ele_forward_acceptance=accep_ele_forward->GetBinContent(ele_theta_bin,ele_p_bin);
			double ele_large_acceptance=accep_ele_large->GetBinContent(ele_theta_bin,ele_p_bin);
            if(mom_ele<1.0||theta_ele>14.8||theta_ele*DEG<8.0)//GeV, CLEO
                ele_forward_acceptance=0.0;//Farward-Angle EC Cut at 1 GeV
            if(mom_ele<3.5||theta_ele*DEG<16.0||theta_ele*DEG>24)//GeV,CLEO
                ele_large_acceptance=0.0; //Larger-Angle EC Cut at 3 GeV
            if(ele_forward_acceptance>1.) 
                ele_forward_acceptance=1.0; 
            if(ele_large_acceptance>1.) 
                ele_large_acceptance=1.0; 

			int hadron_theta_bin=int(theta_had*DEG/0.2)+1;  //0.2 degree per bin
			int hadron_p_bin=int(mom_had/0.05)+1;     //0.05 GeV per bin
			double hadron_acceptance=0;
			hadron_acceptance=accep_had->GetBinContent(hadron_theta_bin,hadron_p_bin);
            if(theta_had*DEG>14.8||theta_had*DEG<8.0||mom_had<0.||mom_had>11.)//GeV, CLEO
                hadron_acceptance=0.0;
            if(hadron_acceptance>1.) 
                hadron_acceptance=1.0; 

			double event_weight=0;   //depend on hadron charge
			weight_hp=dxs_hp*Phase_space/N_simulate; 
			//carefull here, N_simulate is different for diff. energy setting, so be aware when you chain all root files with diff. enerngy settings togeter 
			weight_hm=dxs_hm*Phase_space/N_simulate;
			if(hadron_flag==1){  //doing hp
				event_weight=weight_hp;
			}else if(hadron_flag==2){  //doing hm
				event_weight=weight_hm;
			}

			double forward_acceptance=ele_forward_acceptance*hadron_acceptance;
			double large_acceptance=ele_large_acceptance*hadron_acceptance;

			//mandatary 8 degree cut
			if(mom_ele>1.
		   // &&((theta_ele*DEG>=8.&&theta_ele*DEG<9.&&mom_ele>=4.)	
		   // ||(theta_ele*DEG>=9.&&theta_ele*DEG<10.&&mom_ele>=3.)	
		   // ||(theta_ele*DEG>=10.&&mom_ele>=2.)	
		   // )
			//&&
		    //((theta_had*DEG>=8.&&theta_had*DEG<9.&&mom_had>=4.)	
		    //||(theta_had*DEG>=9.&&theta_had*DEG<10.&&mom_had>=3.)	
		    //||(theta_had*DEG>=10.&&mom_had>=2.)	
		   // )
			)
				rate_forward+=event_weight*forward_acceptance;
			if(mom_ele>3.5
		    //&&((theta_ele*DEG<=20.&&mom_ele>=3.)	
		    //||(theta_ele*DEG>20.&&mom_ele>=2.)	
		    //)
			)
				rate_large+=event_weight*large_acceptance;

		}// if cut ends here
	}//event loop ends here

	
	/*Print&Save{{{*/
	const double Lumi = 1.0e36; // cm-2*s-1, for He3 nuclear not for nucleons
	const double KHz = 1e-3;
	const double nBcm2 = 1e-33;
	ofstream outf(Form("%s_%s_he3_rate_Q2gt1.dat", particle.Data(), charge.Data()));

	cout<<"______forward and large angle integral rate____________________"<<endl;
	cout<<"rate_forward: "<<rate_forward*Lumi*KHz*nBcm2<<endl;
	cout<<"rate_large: "<<rate_large*Lumi*KHz*nBcm2<<endl;
	cout<<"_______________________________________________________________"<<endl;
	/*}}}*/

}

