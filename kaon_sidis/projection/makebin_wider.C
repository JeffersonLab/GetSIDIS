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
#include <TRandom3.h>
#include <TApplication.h>
#include <Rtypes.h>
#include <TTree.h>
//#include <TMatrix.h>
/*}}}*/

using namespace std;
const double DEG = 180./3.1415926;
const double PI = 3.1415926;

double max(double a, double b);
int makebin(int Ebeam, int target_flag, int particle_flag, int z_flag, int Q2_flag);

/*int main{{{*/
int main(){
    Int_t target_flag = 3; cerr<<"--- Target (1->p, 3->he3) = "; cin >> target_flag;
    Int_t particle_flag = 2; cerr<<"--- particle (1->pion, 2->kaon) = "; cin >> particle_flag;
    Int_t Ebeam = 0; cerr<<"--- particle (11 or 8) = "; cin >> Ebeam;

    Int_t z_flag = 0;
    Int_t Q2_flag = 0;
    int err = -1000;
    for(int i=1;i<=6;i++){
        for(int j=1;j<=4;j++){
            z_flag = i;
            Q2_flag = j;
            err = makebin(Ebeam, target_flag, particle_flag, z_flag, Q2_flag );
        }
    }
    return err;
}
/*}}}*/

int makebin(int Ebeam, int target_flag, int particle_flag, int z_flag, int Q2_flag){
	TString target = "X";/*{{{*/
	if(target_flag==1)
		target ="NH3";
	else if(target_flag==2)
		target ="D2";
	else if(target_flag==3)
		target ="He3";
	else{
		cerr<<"I don't know this particle flag!"<<endl;
		return -1;
	}
	TString particle = "X";
	if(particle_flag==1)
		particle ="pion";
	else if(particle_flag==2)
		particle ="kaon";
	else{
		cerr<<"I don't know this particle flag!"<<endl;
		return -1;
	}/*}}}*/
	
	/*Get Z and Q2 Bin{{{*/
	double zmin =0., zmax = 0., Q2min=0., Q2max=0.;
	if (z_flag==1){
		zmin = 0.3; zmax = 0.35;
	}else if (z_flag==2){
		zmin = 0.35; zmax = 0.4;
	}else if (z_flag==3){
		zmin = 0.4; zmax = 0.45;
	}else if (z_flag==4){
		zmin = 0.45; zmax = 0.5;
	}else if (z_flag==5){
		zmin = 0.5; zmax = 0.6;
	}else if (z_flag==6){
		zmin = 0.6; zmax = 0.7;
	}
	if (Q2_flag==1){
		Q2min = 1.; Q2max = 2.0;
	}else if (Q2_flag==2){
		Q2min = 2.0; Q2max = 3.5;
	}else if (Q2_flag==3){
		Q2min = 3.5; Q2max = 5.5;
	}else if (Q2_flag==4){
		Q2min = 5.5; Q2max = 10.;
	}
	/*}}}*/
	
	/*Define old root file{{{*/
	TString prefix = "../skim_rootfiles/";
	TString	filename = Form("skim_%s_%s_E%d_z%d_Qsq%d.root",target.Data(), particle.Data(),Ebeam,z_flag,Q2_flag);

	TString new_filename = prefix + filename;
	TChain *T = new TChain("T","T");
	T->AddFile(new_filename);
	cerr<<Form(" For %3.2f<z<%3.2f, %2.0f<Q2<%2.0f  Got total number of events = %e",zmin,zmax,Q2min,Q2max, (double)(T->GetEntries()))<<endl;

	Double_t theta_gen, phi_gen,mom_gen;
	Double_t Q2,W,Wp,x,y,z,pt,nu,s;
	Double_t theta_q,phi_q;
	Double_t theta_s,phi_s,phi_h;
	Double_t jacoF,dxs_hp,dxs_hm,dxs_hp_sidis,dxs_hm_sidis;
	Double_t mom_ele,mom_had;
	Double_t theta_ele,theta_had;
	Double_t phi_ele,phi_had;
	Double_t mom_pro,energy_pro;
	Double_t mom_ini_ele,energy_ini_ele,weight_hp,weight_hm, W_hp, W_hm;
	Double_t acc_f_ele, acc_l_ele, acc_f_hp, acc_l_hp, acc_f_hm,acc_l_hm;
    Double_t luminosity, time;
	Double_t dilute[2];
	Int_t nsim, isphy_hp, isphy_hm;

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
    T->SetBranchAddress("dxs_hm_sidis",&dxs_hm_sidis);
    T->SetBranchAddress("dxs_hp_sidis",&dxs_hp_sidis);
	T->SetBranchAddress("mom_ele",&mom_ele);
	T->SetBranchAddress("mom_had",&mom_had);
	T->SetBranchAddress("theta_ele",&theta_ele);
	T->SetBranchAddress("theta_had",&theta_had);
	T->SetBranchAddress("phi_ele",&phi_ele);
	T->SetBranchAddress("phi_had",&phi_had);
	T->SetBranchAddress("nsim",&nsim);
	T->SetBranchAddress("dilute_hp",&dilute[0]);
	T->SetBranchAddress("dilute_hm",&dilute[1]);
	T->SetBranchAddress("weight_hp",&weight_hp);
	T->SetBranchAddress("weight_hm",&weight_hm);
	T->SetBranchAddress("acc_f_ele",&acc_f_ele);
	T->SetBranchAddress("acc_l_ele",&acc_l_ele);
	T->SetBranchAddress("acc_f_hp",&acc_f_hp);
	T->SetBranchAddress("acc_l_hp",&acc_l_hp);
    T->SetBranchAddress("acc_f_hm",&acc_f_hm);
    T->SetBranchAddress("acc_l_hm",&acc_l_hm);
    T->SetBranchAddress("W_hp",&W_hp);
    T->SetBranchAddress("W_hm",&W_hm);
	T->SetBranchAddress("isphy_hp",&isphy_hp);
	T->SetBranchAddress("isphy_hm",&isphy_hm);
	T->SetBranchAddress("time",&time);
	T->SetBranchAddress("luminosity",&luminosity);
	/*}}}*/

	// 100 days, 2e34 /cm^2/s ; the last factor is sqrt(2) for the angular separation
//	const double dilute_factor = 0.85; //85% overall efficiency,1/3 for He3 

	double dilute_factor = 0.2;//Fix this for both pi+ (0.172) and pi-(0.267)
	if(particle_flag==1)
		dilute_factor = 0.172;//A rough guess 
	else
		dilute_factor = 0.267;//A rough guess 

	const double polarization = 0.6/sqrt(2)  ; //60% polarization,1/2 here is for separating Asivers,Acollins,Ap..
	const double det_eff_e = 0.85; //85% detector efficiency for electrons 
	const double det_eff_h = 0.85; //85% detector efficiency for hadrons   
	const double target_factor = 0.865; //neutron polarization 86.5%
	const double Asys = 0.01; //FIX_HERE, for Collins asymmetry, 1/sqrt(N), from Min's EIC code
	const double Nsys = 1e+4; //FIX_HERE, for Collins asymmetry, 1/sqrt(N),from Min's SoLID code

	Double_t ptmin,ptmax,xmin,xmax,nevent_pt,nevent_x;
	Int_t xb_min,xb_max,ptb_min,ptb_max;
	Int_t xbin, ptbin;

	const double pt_step = 0.4;
	const double PT_MIN  =0.0;
	const double PT_MAX = 1.6;
	ptbin = (int)((PT_MAX-PT_MIN)/pt_step); 

	const double x_step = 0.02;
	const double X_MIN = 0.0;
	const double X_MAX = 1.0;
	xbin = (int)((X_MAX-X_MIN)/x_step); //Temp, real bin size is determined by the events in each bin

	prefix = "./databases_wider/";
	filename = Form("%s_%s_E%d_%d_%d.dat",target.Data(), particle.Data(),Ebeam, z_flag,Q2_flag);
	new_filename = prefix + filename;
	ofstream outf_total(new_filename);
	outf_total << ptbin << endl;

    filename = Form("%s_%s_E%d_%d_%d.dat",target.Data(), particle.Data(),Ebeam, z_flag,Q2_flag);
    prefix = "./databases_p_wider/";
    new_filename = prefix + filename;
    ofstream outfile_p(new_filename);
    outfile_p << ptbin<<endl;

    prefix = "./databases_m_wider/";
    new_filename = prefix + filename;
    ofstream outfile_m(new_filename);
    outfile_m << ptbin<<endl;

	ptmin = PT_MIN; ptb_min = 1;ptb_max = 1;
	TString histoname;
	//loop through x and pt bin
	for (Int_t i=0;i<ptbin;i++){
		ptmin = PT_MIN + i* pt_step;
		if(ptmin>(PT_MAX-pt_step)) i=ptbin; //Stop
		ptmax = PT_MIN + (i+1)* pt_step;
		if(ptmax>PT_MAX) continue; //Stop

		prefix = "./out_rootfiles_wider/";
        filename = Form("%s_%s_E%d_%d_%d_%d.root",target.Data(), particle.Data(),Ebeam, z_flag,Q2_flag,i);
		new_filename = prefix + filename;
        /*Define new rootfile for each bin{{{*/
        TFile *file = new TFile(new_filename,"RECREATE");
        TTree *t1 = new TTree("T","T");
		t1->SetDirectory(file);

		t1->Branch("Q2",&Q2,"data/D");
		t1->Branch("W",&W,"data/D");
		t1->Branch("Wp",&Wp,"data/D");
		t1->Branch("x",&x,"data/D");
		t1->Branch("y",&y,"data/D");
		t1->Branch("z",&z,"data/D");
		t1->Branch("nu",&nu,"data/D");
		t1->Branch("s",&s,"data/D");
		t1->Branch("pt",&pt,"data/D");
		t1->Branch("theta_q",&theta_q,"data/D");
		t1->Branch("theta_s",&theta_s,"data/D");
		t1->Branch("phi_h",&phi_h,"data/D");
		t1->Branch("phi_s",&phi_s,"data/D");
		t1->Branch("jacoF",&jacoF,"jacoF/D");
		t1->Branch("dxs_hm",&dxs_hm,"dxs_hm/D");
		t1->Branch("dxs_hp",&dxs_hp,"dxs_hp/D");
		t1->Branch("mom_ele",&mom_ele,"mom_ele/D");
		t1->Branch("mom_had",&mom_had,"mom_had/D");
		t1->Branch("theta_ele",&theta_ele,"theta_ele/D");
		t1->Branch("theta_had",&theta_had,"theta_had/D");
		t1->Branch("phi_ele",&phi_ele,"phi_ele/D");
		t1->Branch("phi_had",&phi_had,"phi_had/D");
		t1->Branch("nsim",&nsim,"nsim/I");
		t1->Branch("dilute_hp",&dilute[0],"data/D");
		t1->Branch("dilute_hm",&dilute[1],"data/D");
		t1->Branch("weight_hp",&weight_hp,"data/D");
        t1->Branch("weight_hm",&weight_hm,"data/D");
        t1->Branch("acc_f_ele",&acc_f_ele,"data/D");
        t1->Branch("acc_l_ele",&acc_l_ele,"data/D");
        t1->Branch("acc_f_hp",&acc_f_hp,"data/D");
        t1->Branch("acc_l_hp",&acc_l_hp,"data/D");
        t1->Branch("acc_f_hm",&acc_f_hm,"data/D");
        t1->Branch("acc_l_hm",&acc_l_hm,"data/D");
        t1->Branch("W_hp",&W_hp,"data/D");
        t1->Branch("W_hm",&W_hm,"data/D");
        t1->Branch("isphy_hp",&isphy_hp,"data/D");
        t1->Branch("isphy_hm",&isphy_hm,"data/D");
        t1->Branch("time",&time,"data/D");
        t1->Branch("luminosity",&luminosity,"data/D");

		for(int k=0;k<T->GetEntries();k++){
			if(pt>=ptmin && pt<ptmax 
					&& Q2>=Q2min && Q2<Q2max
					&& z>=zmin && z<zmax 	
			  )
				t1->Fill();
		}
		t1->Write();
        /*}}}*/

		/*Histograms{{{*/
		//pt
		TH1F *h1p  =new TH1F("h1p","h1p",1000,PT_MIN,PT_MAX);
		TH1F *h1m  =new TH1F("h1m","h1m",1000,PT_MIN,PT_MAX);
		//x
		TH1F *h2p = new TH1F("h2p","h2p",1000,X_MIN,X_MAX);
		TH1F *h2m = new TH1F("h2m","h2m",1000,X_MIN,X_MAX);
		/*}}}*/

		TString cut1p,cut2p;
		TString cut1m,cut2m;
		//"weight" is defined in skim.C
		//luminosity = 1e36, time = 48days*24hr*3600s for 11GeV and = 21days*24hr*3600s for 8.8GeV
		//weight = 1e-33*dxs_hp*ele_accp*had_accp/Nsim*luminisity*time*(solid_forward_accp+solid_large_accp) 
        
		/*X Binning{{{*/
        cut1p="W_hp";
        cut1m="W_hm";

        TString cut2=Form("(pt>=%f&&pt<%f&&(isphy_hp>0))",ptmin,ptmax);
        cut2p = cut1p+ "*" + cut2;
		cut2m = cut1m+ "*" + cut2;

		h1p->Reset();
		T->Project("h1p","pt",TCut(cut2p));
		h2p->Reset();
		T->Project("h2p","x",TCut(cut2p));

        h1m->Reset();
		T->Project("h1m","pt",TCut(cut2m));
		h2m->Reset();
		T->Project("h2m","x",TCut(cut2m));

        h1p->SetDirectory(file);
		h1m->SetDirectory(file);

        h2p->SetDirectory(file);
        h2m->SetDirectory(file);

		double Norm_Factor = (pow(polarization * target_factor*dilute_factor,2) * det_eff_e * det_eff_h );
		//A = 1./sqrt(N*det_eff_e*det_eff_h)/pol/target_factor/dilute
		//double single_bin_raw = 1.0/Norm_Factor/pow(Asys,2);
		double single_bin_raw = Nsys/det_eff_e/det_eff_h;
		Double_t total_eve_hp = h2p->GetSum();
		Double_t total_eve_hm = h2m->GetSum();
        Double_t total_eve = max(total_eve_hp, total_eve_hm);

		/*Count How many x bins{{{*/
		int xbin_fix = (int)((X_MAX-X_MIN)/x_step);
		double x1[1000];
		int new_bin = 0;
		x1[0] = X_MIN;

		xb_max =0;
		while(xb_max<1000){
			new_bin++;
			nevent_x = 0;

			//	cerr<<Form("--- Start with x=%f, xb=%d",h2p->GetBinCenter(xb_max),xb_max)<<endl;
			//while(nevent_x < total_eve)/xbin_fix&&xb_max<1000){
			while(nevent_x <=max(single_bin_raw,total_eve/8) &&xb_max<1000){
				nevent_x += h2p->GetBinContent(xb_max++);
			}
			xmax = h2p->GetBinCenter(xb_max);

			x1[new_bin] = xmax;
			cerr<<Form("    #%d bin: x=%f, xb=%d, N=%e ", new_bin, xmax, xb_max, nevent_x)<<endl;
		}

		xbin = new_bin;
		/*}}}*/

	//	xbin = total_eve/single_bin;
	//	if (xbin<1) xbin = 1; //Counts as one bin if the statistic is really low

		outfile_p << i << "\t" << xbin << endl;
		outfile_m << i << "\t" << xbin << endl;
		outf_total << i << "\t" << xbin << endl;
		cerr    << Form("--- Total x-bin in #%d pt bin is %d with N=%e (set) .vs. %e(h+) %e(h-)" ,i,xbin,single_bin_raw, total_eve_hp, total_eve_hm) << endl;
		xmin = X_MIN;
		for (Int_t j=0;j<xbin;j++){
			xmin = x1[j];
			xmax = x1[j+1];

			//			cerr<<Form("---- Working on pt=%f, x=%f",h1z->GetMean(),(xmax+xmin)/2.)<<endl;

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
			TH1F *h1pt = new TH1F("h1pt","h1pt",1000.,PT_MIN,PT_MAX);
			// dilution factor
			TH1F *h1di = new TH1F("h1di","h1di",1000.,0.,10.);

			TH1F *h1phi_hm = new TH1F("h1phi_hm","h1phi_hm",360,0.,360.);
			TH1F *h1phi_hp = new TH1F("h1phi_hp","h1phi_hp",360,0.,360.);
            TH1F *h2phi_hm = new TH1F("h2phi_hm","h2phi_hm",360,0.,360.);
            TH1F *h2phi_hp = new TH1F("h2phi_hp","h2phi_hp",360,0.,360.);
			/*}}}*/

			/*Fill and Save for hadron+{{{*/
			h1di->Reset();
			h1Q2->Reset();
			h1x->Reset();
			h1y->Reset();
			h1z->Reset();
			h1pt->Reset();
			h1phi_hp->Reset();
			h2phi_hp->Reset();

			cut2.Form("(pt>=%f&&pt<%f&&x>=%f&&x<%f)",ptmin,ptmax,xmin,xmax);
			cut2 = cut1p + "*" + cut2;
			TCut cut_p = TCut(cut2);

			T->Project("h1Q2","Q2",cut_p);
			T->Project("h1x","x",cut_p);
			T->Project("h1z","z",cut_p);
			T->Project("h1y","y",cut_p);
			T->Project("h1pt","pt",cut_p);

			//In Collider, dilute = dxs_proton /dxs_neutron
            T->Project("h1di","dilute_hp",cut_p);

			Double_t dilution = 0.0; 
			Double_t temp_dilute = h1di->GetMean();
			//A = 1./sqrt(N_out) = 1./(pol*target_factor*dilut*sqrt(N_raw*det_eff_e*det_eff*h))	
			double N_out = 0;
			double N_raw = 0;
			if (target_flag ==1){
				dilution = dilute_factor;	
			    //A = 1./sqrt(N_out) = 1./(pol*target_factor*dilut*sqrt(N_raw*det_eff_e*det_eff*h))	
				N_out = h1Q2->GetSum()*(pow(polarization * target_factor * dilution,2) * det_eff_e * det_eff_h );
				N_raw = h1Q2->GetSum()*(det_eff_e * det_eff_h );
				outfile_p << i << " \t" << j << " \t" 
					<< h1z->GetMean() << " \t" 
					<< h1Q2->GetMean() << " \t" 
					<< h1pt->GetMean() << " \t"
					<< h1x->GetMean() << " \t" 
					<< h1y->GetMean() << " \t" 
					<< dilution << " \t" 
					<< N_out << " \t" 
					<< N_raw
					<< endl;
				outf_total <<"had+  " << i << " \t" << j << " \t" 
					<< h1z->GetMean() << " \t" 
					<< h1Q2->GetMean() << " \t" 
					<< h1pt->GetMean() << " \t"
					<< h1x->GetMean() << " \t" 
					<< h1y->GetMean() << " \t" 
					<< dilution << " \t" 
					<< N_out << " \t" 
					<< N_raw
					<< endl;
			}else if (target_flag == 2 ){
				// Need to take into account the proton polarization
				//dilution = dilute_factor/(temp_dilute*temp_dilute*4.+(1+temp_dilute)*(1+temp_dilute));	
				//dilution = dilute_factor;	
				dilution = 1./(1.+temp_dilute);
			    //A = 1./sqrt(N_out) = 1./(pol*target_factor*dilut*sqrt(N_raw*det_eff_e*det_eff*h))	
				N_out = h1Q2->GetSum()*(pow(polarization * target_factor * dilution,2) * det_eff_e * det_eff_h );
				N_raw = h1Q2->GetSum()*(det_eff_e * det_eff_h );
				outfile_p << i << " \t" << j << " \t" 
					<< h1z->GetMean() << " \t" 
					<< h1Q2->GetMean() << " \t" 
					<< h1pt->GetMean() << " \t" 
					<< h1x->GetMean() << " \t" 
					<< h1y->GetMean() << " \t" 
					<< dilution << " \t" 
					<< N_out << " \t" 
					<< N_raw
					<< endl;
				outf_total <<"had+  " << i << " \t" << j << " \t" 
					<< h1z->GetMean() << " \t" 
					<< h1Q2->GetMean() << " \t" 
					<< h1pt->GetMean() << " \t" 
					<< h1x->GetMean() << " \t" 
					<< h1y->GetMean() << " \t" 
					<< dilution << " \t" 
					<< N_out << " \t" 
					<< N_raw
					<< endl;
			}else if (target_flag == 3){
				//dilution = dilute_factor/((1+2.*temp_dilute)*(1+2.*temp_dilute)) ; 
				dilution = 1./(1.+2.*temp_dilute);
			    //A = 1./sqrt(N_out) = 1./(pol*target_factor*dilut*sqrt(N_raw*det_eff_e*det_eff*h))	
				N_out = h1Q2->GetSum()*(pow(polarization * target_factor * dilution,2) * det_eff_e * det_eff_h );
				N_raw = h1Q2->GetSum()*(det_eff_e * det_eff_h );
				outfile_p << i << " \t" << j << " \t"
					<< h1z->GetMean() << " \t" 
					<< h1Q2->GetMean() << " \t" 
					<< h1pt->GetMean() << " \t" 
					<< h1x->GetMean() << " \t" 
					<< h1y->GetMean() << " \t" 
					<< dilution << " \t" 
					<< N_out << " \t" 
					<< N_raw
					<< endl;
				outf_total <<"had+  " << i << " \t" << j << " \t"
					<< h1z->GetMean() << " \t" 
					<< h1Q2->GetMean() << " \t" 
					<< h1pt->GetMean() << " \t" 
					<< h1x->GetMean() << " \t" 
					<< h1y->GetMean() << " \t" 
					<< dilution << " \t" 
					<< N_out << " \t" 
					<< N_raw
					<< endl;
				cerr<<Form("---- Working on pt=(%d)%f, x=(%d)%f, N=%e",i,h1pt->GetMean(),j,h1x->GetMean(),N_out)<<endl;
			}

            histoname.Form("hp_%d_%d",i,j);
			TH1F *hp = (TH1F*)h1phi_hp->Clone(histoname);

			cut2.Form("(pt>=%f&&pt<%f&&x>=%f&&x<%f&&phi_h>0.)",ptmin,ptmax,xmin,xmax);
			cut2 = cut1p + "*" + cut2;
			cut_p = cut2;
			T->Project("h1phi_hp","phi_h*180./3.1415926",cut_p);
			cut2.Form("(pt>=%f&&pt<%f&&x>=%f&&x<%f&&phi_h<0.)",ptmin,ptmax,xmin,xmax);
			cut2 = cut1p + "*" + cut2;
			cut_p = cut2;
			T->Project("h2phi_hp","phi_h*180./3.1415926+360.",cut_p);

			hp->Add(h1phi_hp);
			hp->Add(h2phi_hp);
			hp->SetDirectory(file);
			hp->Write();
            /*}}}*/
			
            /*Fill and Save for hadron-{{{*/
			h1di->Reset();
			h1Q2->Reset();
			h1x->Reset();
			h1y->Reset();
			h1z->Reset();
			h1pt->Reset();
			h1phi_hm->Reset();
			h2phi_hm->Reset();

			cut2.Form("(pt>=%f&&pt<%f&&x>=%f&&x<%f)",ptmin,ptmax,xmin,xmax);
			cut2 = cut1m + "*" + cut2;
			TCut cut_m = TCut(cut2);

			T->Project("h1Q2","Q2",cut_m);
			T->Project("h1x","x",cut_m);
			T->Project("h1z","z",cut_m);
			T->Project("h1y","y",cut_m);
			T->Project("h1pt","pt",cut_m);


			//In Collider, dilute = dxs_proton /dxs_neutron
            T->Project("h1di","dilute_hm",cut_m);

			temp_dilute = h1di->GetMean();
			//A = 1./sqrt(N_out) = 1./(pol*target_factor*dilut*sqrt(N_raw*det_eff_e*det_eff*h))	
			N_out = 0;
			N_raw = 0;
			if (target_flag ==1){
				dilution = dilute_factor;	
			    //A = 1./sqrt(N_out) = 1./(pol*target_factor*dilut*sqrt(N_raw*det_eff_e*det_eff*h))	
				N_out = h1Q2->GetSum()*(pow(polarization * target_factor * dilution,2) * det_eff_e * det_eff_h );
				N_raw = h1Q2->GetSum()*(det_eff_e * det_eff_h );
				outfile_m << i << " \t" << j << " \t" 
					<< h1z->GetMean() << " \t" 
					<< h1Q2->GetMean() << " \t" 
					<< h1pt->GetMean() << " \t"
					<< h1x->GetMean() << " \t" 
					<< h1y->GetMean() << " \t" 
					<< dilution << " \t" 
					<< N_out << " \t" 
					<< N_raw
					<< endl;
				outf_total <<"had-  " << i << " \t" << j << " \t" 
					<< h1z->GetMean() << " \t" 
					<< h1Q2->GetMean() << " \t" 
					<< h1pt->GetMean() << " \t"
					<< h1x->GetMean() << " \t" 
					<< h1y->GetMean() << " \t" 
					<< dilution << " \t" 
					<< N_out << " \t" 
					<< N_raw
					<< endl;
			}else if (target_flag == 2 ){
				// Need to take into account the proton polarization
				//dilution = dilute_factor/(temp_dilute*temp_dilute*4.+(1+temp_dilute)*(1+temp_dilute));	
				//dilution = dilute_factor;	
				dilution = 1./(1.+temp_dilute);
			    //A = 1./sqrt(N_out) = 1./(pol*target_factor*dilut*sqrt(N_raw*det_eff_e*det_eff*h))	
				N_out = h1Q2->GetSum()*(pow(polarization * target_factor * dilution,2) * det_eff_e * det_eff_h );
				N_raw = h1Q2->GetSum()*(det_eff_e * det_eff_h );
				outfile_m << i << " \t" << j << " \t" 
					<< h1z->GetMean() << " \t" 
					<< h1Q2->GetMean() << " \t" 
					<< h1pt->GetMean() << " \t" 
					<< h1x->GetMean() << " \t" 
					<< h1y->GetMean() << " \t" 
					<< dilution << " \t" 
					<< N_out << " \t" 
					<< N_raw
					<< endl;
				outf_total <<"had-  " << i << " \t" << j << " \t" 
					<< h1z->GetMean() << " \t" 
					<< h1Q2->GetMean() << " \t" 
					<< h1pt->GetMean() << " \t" 
					<< h1x->GetMean() << " \t" 
					<< h1y->GetMean() << " \t" 
					<< dilution << " \t" 
					<< N_out << " \t" 
					<< N_raw
					<< endl;
			}else if (target_flag == 3){
				//dilution = dilute_factor/((1+2.*temp_dilute)*(1+2.*temp_dilute)) ; 
				dilution = 1./(1.+2.*temp_dilute);
			    //A = 1./sqrt(N_out) = 1./(pol*target_factor*dilut*sqrt(N_raw*det_eff_e*det_eff*h))	
				N_out = h1Q2->GetSum()*(pow(polarization * target_factor * dilution,2) * det_eff_e * det_eff_h );
				N_raw = h1Q2->GetSum()*(det_eff_e * det_eff_h );
				outfile_m << i << " \t" << j << " \t"
					<< h1z->GetMean() << " \t" 
					<< h1Q2->GetMean() << " \t" 
					<< h1pt->GetMean() << " \t" 
					<< h1x->GetMean() << " \t" 
					<< h1y->GetMean() << " \t" 
					<< dilution << " \t" 
					<< N_out << " \t" 
					<< N_raw
					<< endl;
				outf_total <<"had-  "<< i << " \t" << j << " \t"
					<< h1z->GetMean() << " \t" 
					<< h1Q2->GetMean() << " \t" 
					<< h1pt->GetMean() << " \t" 
					<< h1x->GetMean() << " \t" 
					<< h1y->GetMean() << " \t" 
					<< dilution << " \t" 
					<< N_out << " \t" 
					<< N_raw
					<< endl;
				cerr<<Form("---- Working on pt=(%d)%f, x=(%d)%f, N=%e",i,h1pt->GetMean(),j,h1x->GetMean(),N_out)<<endl;
			}

			histoname.Form("hm_%d_%d",i,j);
			TH1F *hm = (TH1F*)h1phi_hm->Clone(histoname);

			cut2.Form("(pt>=%f&&pt<%f&&x>=%f&&x<%f&&phi_h>0.)",ptmin,ptmax,xmin,xmax);
			cut2 = cut1m + "*" + cut2;
			cut_m = cut2;
			T->Project("h1phi_hm","phi_h*180./3.1415926",cut_m);
			cut2.Form("(pt>=%f&&pt<%f&&x>=%f&&x<%f&&phi_h<0.)",ptmin,ptmax,xmin,xmax);
			cut2 = cut1m + "*" + cut2;
			cut_m = cut2;
			T->Project("h2phi_hm","phi_h*180./3.1415926+360.",cut_m);

			hm->Add(h1phi_hm);
			hm->Add(h2phi_hm);
			hm->SetDirectory(file);
			hm->Write();
            /*}}}*/

			h1z->Delete();
			h1x->Delete();
			h1y->Delete();
			h1Q2->Delete();
			h1pt->Delete();
			h1di->Delete();
			h1phi_hp->Delete();
			h2phi_hp->Delete();			
			h1phi_hm->Delete();
			h2phi_hm->Delete();			
		}
		/*}}}*/

		h1p->Delete();
		h1m->Delete();
		h2p->Delete();
		h2m->Delete();

		file->Write();
		file->Close();
	}
    outfile_p.close();
    outfile_m.close();
    outf_total.close();
	return 0;
}

double max(double a, double b){
 if(a>b)
	 return a;
 else
	 return b;
}
