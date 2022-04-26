/*C/C++ Includes{{{*/
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
//#include "Rtypes.h"
//#include "math.h"

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
#include <TMatrixD.h>
#include <TRandom3.h>
#include <TApplication.h>
#include <Rtypes.h>
#include <TTree.h>
#include "LHAPDF/LHAPDF.h"
//#include <TMatrix.h>
/*}}}*/
#include "../solenoid_lib/Sole_inter.h"

using namespace std;
const double DEG = 180./3.1415926;
const double PI = 3.1415926;
const int N=3;

double max(double a, double b);
void getcoef(Double_t, Double_t*, Double_t*);
int makebin(int target_flag, int particle_flag, int Q2_flag);

/*int main{{{*/
int main(){
	Int_t target_flag = 3; cerr<<"--- Target (1->p, 3->he3) = "; cin >> target_flag;
	Int_t particle_flag = 2; cerr<<"--- particle (1->pip, 2->pim) = "; cin >> particle_flag;
	Int_t Q2_flag = 0;
	int err = -1000;

	for(int j=1;j<=6;j++){
		Q2_flag = j;
		err = makebin(target_flag, particle_flag, Q2_flag );
	}
   // Q2_flag = 2;//2.0 GeV ~ 3.0GeV
//	err = makebin(target_flag, particle_flag, Q2_flag );
	return err;
}
/*}}}*/

int makebin(int target_flag, int particle_flag, int Q2_flag){
	TString target = "X";
	if(target_flag==1)
		target ="p";
	else if(target_flag==2)
		target ="d2";
	else if(target_flag==3)
		target ="3he";
	else{
		cerr<<"I don't know this particle flag!"<<endl;
		return -1;
	}
	TString particle = "X";
	if(particle_flag==1)
		particle ="pip";
	else if(particle_flag==2)
		particle ="pim";
	else{
		cerr<<"I don't know this particle flag!"<<endl;
		return -1;
	}

	/*Get Q2 Bin{{{*/
	double Q2min=0., Q2max=0.;
	if (Q2_flag==1){
		Q2min = 1.; Q2max = 2.0;
	}else if (Q2_flag==2){
		Q2min = 2.0; Q2max = 3.0;
	}else if (Q2_flag==3){
		Q2min = 3.0; Q2max = 4.0;
	}else if (Q2_flag==4){
		Q2min = 4.0; Q2max = 5.0;
	}else if (Q2_flag==5){
		Q2min = 5.0; Q2max = 6.0;
	}else if (Q2_flag==6){
		Q2min = 6.0; Q2max = 8.0;
	}else if (Q2_flag==7){//Not in used!!!
		Q2min = 8.0; Q2max = 10.0;
	}
	/*}}}*/

	/*Define old root file{{{*/
	TString prefix = "../skim/skim_rootfiles_11p8/";

	TString	filename = "";	
	TString	new_filename = "";	
	TChain *T = new TChain("T","T");
	for(int z_flag=1;z_flag<=8;z_flag++){
		filename = Form("%s_skim_%s_%d_%d.root",target.Data(), particle.Data(),z_flag,Q2_flag);
		new_filename = prefix + filename;
		T->AddFile(new_filename);
	}
	cerr<<Form(" For %f<Q2<%f  Got total number of events = %e",Q2min,Q2max, (double)(T->GetEntries()))<<endl;

	Double_t Q2,W,Wp,x,y,z,pt,nu,s;
	Double_t theta_q;
	Double_t theta_s,phi_s,phi_h;
	Double_t jacoF,dxs_hp,dxs_hm;
	Double_t mom_ele,mom_had;
	Double_t theta_ele,theta_had;
	Double_t phi_ele,phi_had;
	Double_t weight,weight_hp,weight_hm;
	Double_t accp_forward, accp_large;
	Double_t dilute[2];
	Int_t nsim;

	T->SetBranchAddress("Q2",&Q2);
	T->SetBranchAddress("W",&W);
	T->SetBranchAddress("Wp",&Wp);
	T->SetBranchAddress("x",&x);
	T->SetBranchAddress("y",&y);
	T->SetBranchAddress("z",&z);
	T->SetBranchAddress("nu",&nu);
	T->SetBranchAddress("s",&s);
	T->SetBranchAddress("pt",&pt);
	T->SetBranchAddress("theta_q",&theta_q);
	T->SetBranchAddress("theta_s",&theta_s);
	T->SetBranchAddress("phi_h",&phi_h);
	T->SetBranchAddress("phi_s",&phi_s);
	T->SetBranchAddress("jacoF",&jacoF);
	T->SetBranchAddress("dxs_hm",&dxs_hm);
	T->SetBranchAddress("dxs_hp",&dxs_hp);
	T->SetBranchAddress("mom_ele",&mom_ele);
	T->SetBranchAddress("mom_had",&mom_had);
	T->SetBranchAddress("theta_ele",&theta_ele);
	T->SetBranchAddress("theta_had",&theta_had);
	T->SetBranchAddress("phi_ele",&phi_ele);
	T->SetBranchAddress("phi_had",&phi_had);
	T->SetBranchAddress("nsim",&nsim);
	T->SetBranchAddress("weight",&weight);
	T->SetBranchAddress("dilute_p",&dilute[0]);
	T->SetBranchAddress("dilute_m",&dilute[1]);
	T->SetBranchAddress("weight_hp",&weight_hp);
	T->SetBranchAddress("weight_hm",&weight_hm);
	T->SetBranchAddress("accp_forward",&accp_forward);
	T->SetBranchAddress("accp_large",&accp_large);
	/*}}}*/

	prefix = "./pretz_both_double/";
	filename = Form("%s_%s_%d.dat",target.Data(), particle.Data(),Q2_flag);
	new_filename = prefix + filename;
	ofstream outf_total(new_filename);


	// 100 days, 2e34 /cm^2/s ; the last factor is sqrt(2) for the angular separation
	//	const double dilute_factor = 0.85; //85% overall efficiency,1/3 for He3 
	double dilute_factor = 0.2;//Fix this for both pi+ (0.172) and pi-(0.267)
	if(particle_flag==1)
		dilute_factor = 0.172;//A rough guess 
	else
		dilute_factor = 0.2; //0.267;//A rough guess 

	const double polarization = 0.6/sqrt(2)  ; //60% polarization,1/2 here is for separating Asivers,Acollins,Ap..
	const double det_eff = 0.85; //85% detector efficiency for electrons and hadrons
	const double target_factor = 0.865; //neutron polarization 86.5%
	//const double Asys = 0.006; //FIX_HERE, for Collins asymmetry, 1/sqrt(N), from Min's EIC code
	const double Nsys = 4e+6; //FIX_HERE, for Collins asymmetry, 1/sqrt(N),from Min's SoLID code

	Double_t xmin,xmax,nevent_x;
	Int_t xb_max, xbin;
	const double x_step = 0.02;
	const double X_MIN = 0.0;
	const double X_MAX = 1.0;
	xbin = (int)((X_MAX-X_MIN)/x_step); //Temp, real bin size is determined by the events in each bin
	TH1F *h2 = new TH1F("h2","h2",1000,X_MIN,X_MAX);

	TString cut1,cut2;
	//"weight" is defined in skim.C
	//luminosity = 1e36, time = 48days*24hr*3600s for 11GeV and = 21days*24hr*3600s for 8.8GeV
	//weight = 1e-33*dxs_hp*ele_accp*had_accp/Nsim*luminisity*time*(solid_forward_accp+solid_large_accp) 
	if (particle_flag==1){
		cut1="weight";
	}else{
		cut1="weight";
	}
	TCut cut(cut1);

	cut2.Form("(Q2>=%f&&Q2<%f)",Q2min,Q2max);
	cut2 = cut1 + "*" + cut2;
	cut = cut2;
	h2->Reset();
	T->Project("h2","x",cut);

	//double Norm_Factor = (pow(polarization * target_factor*dilute_factor,2) * det_eff);
	//A = 1./sqrt(N*det_eff)/pol/target_factor/dilute
	//double single_bin_raw = 1.0/Norm_Factor/pow(Asys,2);
	double single_bin_raw = Nsys/det_eff;

	/*Count How many x bins{{{*/
	double x1[1000];
	int new_bin = 0;
	x1[0] = X_MIN;

	xb_max =0;
	while(xb_max<1000){
		new_bin++;
		nevent_x = 0;
		//	cerr<<Form("--- Start with x=%f, xb=%d",h2->GetBinCenter(xb_max),xb_max)<<endl;
		//while(nevent_x < h2->GetSum()/xbin&&xb_max<1000){
		while(nevent_x <=max(single_bin_raw,h2->GetSum()/4) &&xb_max<1000){
			nevent_x += h2->GetBinContent(xb_max++);
		}
		xmax = h2->GetBinCenter(xb_max);

		x1[new_bin] = xmax;
		cerr<<Form("    #%d bin: x=%f, xb=%d, N=%e ", new_bin, xmax, xb_max, nevent_x)<<endl;
	}

	xbin = new_bin;
	if (xbin<1) xbin = 1; //Counts as one bin if the statistic is really low
	/*}}}*/

	outf_total << Q2_flag << "\t" << xbin << endl;
	cerr    << Form("--- Total x-bin in #%d Q2 bin is %d with N=%e .vs. %e" ,Q2_flag,xbin,single_bin_raw, h2->GetSum()/100) << endl;
	TString histoname;
	xmin = X_MIN;
	/*Binning in X{{{*/
	for (Int_t i=0;i<xbin;i++){
		xmin = x1[i];
		xmax = x1[i+1];

		/*Histograms{{{*/
		//Q2
		TH1F *h1Q2  =new TH1F("h1Q2","h1Q2",1000,1.,8.);
		//x 
		TH1F *h1x = new TH1F("h1x","h1x",1000,X_MIN,X_MAX);
		// z
		TH1F *h1z = new TH1F("h1z","h1z",1000,0.3,0.7);
		// y 
		TH1F *h1y = new TH1F("h1y","h1y",1000,0.,1.);
		// pt
		TH1F *h1pt = new TH1F("h1pt","h1pt",1000.,0.0,1.6);
		// dilution factor
		TH1F *h1di = new TH1F("h1di","h1di",1000.,0.,10.);

		TH1F *h1phi_h = new TH1F("h1phi_h","h1phi_h",360,0.,360.);
		TH1F *h2phi_h = new TH1F("h2phi_h","h2phi_h",360,0.,360.);
		/*}}}*/

		/*Fill and Save{{{*/
		h1di->Reset();
		h1Q2->Reset();
		h1x->Reset();
		h1y->Reset();
		h1z->Reset();
		h1pt->Reset();
		h1phi_h->Reset();
		h1phi_h->Reset();

		cut2.Form("(Q2>=%f&&Q2<%f&&x>=%f&&x<%f)",Q2min,Q2max,xmin,xmax);
		cut2 = cut1 + "*" + cut2;
		cut = cut2;

		T->Project("h1Q2","Q2",cut);
		T->Project("h1x","x",cut);
		T->Project("h1z","z",cut);
		T->Project("h1y","y",cut);
		T->Project("h1pt","pt",cut);

		//In Collider, dilute = dxs_proton /dxs_neutron
		if (particle_flag==1){
			T->Project("h1di","dilute_p",cut);
		}else{
			T->Project("h1di","dilute_m",cut);
		}

		Double_t dilution = 0.0; 
		Double_t temp_dilute = h1di->GetMean();

		cut2.Form("(x>=%f&&x<%f&&phi_h>0.)",xmin,xmax);
		cut2 = cut1 + "*" + cut2;
		cut = cut2;
		T->Project("h1phi_h","phi_h*180./3.1415926",cut);
		cut2.Form("(x>=%f&&x<%f&&phi_h<0.)",xmin,xmax);
		cut2 = cut1 + "*" + cut2;
		cut = cut2;
		T->Project("h2phi_h","phi_h*180./3.1415926+360.",cut);

		histoname.Form("h_%d",i);
		TH1F *hhh = (TH1F*)h1phi_h->Clone(histoname);
		hhh->Add(h1phi_h);
		hhh->Add(h2phi_h);
		//hhh->SetDirectory(file);
		//hhh->Write();

        double phi_sum = 0.0;
		double phi_bin[360],phi_weight[360];
	    double coverage = 0.0;
		double angle[4];
	    double coef[3];
		angle[0]= 0.0; angle[1]= 0.0; angle[2]= 0.0; angle[3]= 0.0;
		coef[0] =0.0; coef[1] =0.0; coef[2] =0.0;
		phi_sum =0;
	
		for (Int_t k=0;k<360;k++){
			phi_bin[k] = hhh->GetBinContent(k+1);
			phi_sum += phi_bin[k];
		}
		if(phi_sum<1){
			coverage = 0.0;
			coef[0] =0.0; coef[1] =0.0; coef[2] =0.0;
		}
		else{
			sole_inter sole;
			sole.init(360,phi_bin);
			sole.connect();
			sole.get_new(phi_bin);
			sole.get_limit(angle);
			for (Int_t k=0;k!=360;k++){
				phi_weight[k] = sole.get_accep(k+0.5);
			//	cerr<<Form("phi_weight[%d] = %f", k, phi_weight[k])<<endl;
			}
			coverage = (angle[1]-angle[0]+angle[3]-angle[2])*PI/180.;
			getcoef(coverage,coef,phi_weight);

			cerr<< Form(" --- coef: 1 = %f, 2=%f, 3=%f ", coef[0],coef[1],coef[2])<<endl;
		}

		//A = 1./sqrt(N_out) = 1./(pol*target_factor*dilut*sqrt(N_raw*det_eff))	
		double N_out = 0;
		if (target_flag ==1){
			dilution = dilute_factor;	
			//A = 1./sqrt(N_out) = 1./(pol*target_factor*dilut*sqrt(N_raw*det_eff))	
			N_out = h1Q2->GetSum()*(pow(polarization * target_factor * dilution,2) * det_eff);
		}else if (target_flag == 2 ){
			// Need to take into account the proton polarization
			//dilution = dilute_factor/(temp_dilute*temp_dilute*4.+(1+temp_dilute)*(1+temp_dilute));	
			//dilution = dilute_factor;	
			dilution = 1./(1.+temp_dilute);
			//A = 1./sqrt(N_out) = 1./(pol*target_factor*dilut*sqrt(N_raw*det_eff))	
			N_out = h1Q2->GetSum()*(pow(polarization * target_factor * dilution,2) * det_eff);
		}else if (target_flag == 3){
			//dilution = dilute_factor/((1+2.*temp_dilute)*(1+2.*temp_dilute)) ; 
			dilution = 1./(1.+2.*temp_dilute);
			//A = 1./sqrt(N_out) = 1./(pol*target_factor*dilut*sqrt(N_raw*det_eff))	
			N_out = h1Q2->GetSum()*(pow(polarization * target_factor * dilution,2) * det_eff);
		}
		cerr<<Form("---- Working on Q2=(%d)%f, x=(%d)%f, N=%e",Q2_flag,h1Q2->GetMean(),i,h1x->GetMean(),N_out)<<endl;
		if(N_out>0.){
			outf_total << Q2_flag <<" \t" << i << " \t"
				<< h1z->GetMean() << " \t" 
				<< h1Q2->GetMean() << " \t" 
				<< h1pt->GetMean() << " \t"
				<< h1x->GetMean() << " \t" 
				<< h1y->GetMean() << " \t" 
				<< 1/sqrt(N_out)<< " \t"
				<< h1x->GetSum()<< " \t"
				<< coverage <<" \t"
				<< coef[0] <<" \t"
				<< coef[1] <<" \t"
				<< coef[2] <<" \t"
				<< endl;
		}
		else{
			outf_total << Q2_flag <<" \t" << i << " \t"
				<< h1z->GetMean() << " \t" 
				<< h1Q2->GetMean() << " \t" 
				<< h1pt->GetMean() << " \t"
				<< h1x->GetMean() << " \t" 
				<< h1y->GetMean() << " \t" 
				<< 0 << " \t"
				<< 0 << " \t"
				<< coverage <<" \t"
				<< coef[0] <<" \t"
				<< coef[1] <<" \t"
				<< coef[2] <<" \t"
				<< endl;
		}
	/*}}}*/
	
		h1z->Delete();
		h1x->Delete();
		h1y->Delete();
		h1Q2->Delete();
		h1pt->Delete();
		h1di->Delete();
		h1phi_h->Delete();
		h2phi_h->Delete();	
        hhh->Delete();		
	}
	/*}}}*/

	h2->Delete();
	outf_total.close();
	return 0;
}

/*getcoef(Double_t angle, Double_t* coef, Double_t* acpt){{{*/
void getcoef(Double_t angle, Double_t* coef, Double_t* acpt){
	//x->phi y->phih 3 terms,w/ acpt
	Double_t a_matrix[N*N]={};
	Double_t x,y;
	for(Int_t i=0;i<180;i++){
		x = PI*(i+0.5)/180.; //phi

		for(Int_t j=0;j<360;j++){
			y = PI*(j+0.5)/180.; //phi_h
			a_matrix[0] += (sin(x)*sin(x))*PI/180.*PI/180.*acpt[j];
			a_matrix[1] += (sin(x)*sin(2*y-x))*PI/180.*PI/180.*acpt[j];
			a_matrix[2] += (sin(x)*sin(2*y+x))*PI/180.*PI/180.*acpt[j];

			a_matrix[N+1] += (sin(2*y-x)*sin(2*y-x))*PI/180.*PI/180.*acpt[j];
			a_matrix[N+2] += (sin(2*y-x)*sin(2*y+x))*PI/180.*PI/180.*acpt[j];

			a_matrix[2*N+2] += (sin(2*y+x)*sin(2*y+x))*PI/180.*PI/180.*acpt[j];
		}
	}
	a_matrix[N] = a_matrix[1];
	a_matrix[2*N] = a_matrix[2];
	a_matrix[2*N+1] = a_matrix[N+2];
	//       cout << "\t Initialize a_matrix done!" << endl;

	TMatrixD a(N,N,a_matrix);
	TMatrixD asav = a;
	a.Invert();
	//   cout << "\t Initialize a invert done!" << endl;

	Double_t qx[N]={};
	for(Int_t i=0;i<180;i++){
		x = PI*(i+0.5)/180.; //phi
		for(Int_t j=0;j<360;j++){
			y = PI*(j+0.5)/180.; //phi_h
			for(Int_t k=0;k<N;k++){
				qx[k] += (pow((a(k,0)*sin(x)+a(k,1)*sin(2*y-x)+a(k,2)*sin(2*y+x)),2))*PI/180.*PI/180.*acpt[j];
			}
		}
	}
	for(Int_t i=0;i<N;i++){
		qx[i] *= PI*angle;
		coef[i] = sqrt(qx[i]/2);
	}
	return;
}
/*}}}*/

double max(double a, double b){
	if(a>b)
		return a;
	else
		return b;
}
