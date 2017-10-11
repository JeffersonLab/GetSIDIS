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
#include <TObjString.h>
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
#include <TRandom3.h>
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

#include "../../generator/SIDIS_Lite_LO.h" //this version doesn't include LHAPDF

using namespace std;

//Input values to the function are
//1. the mass of the nucleus under consideration, 2 for deuterium and 12 for Carbon
//2. the number of the rootfile for smearing (usually from 0-500)
int main(Int_t argc, char *argv[]){


   if (argc<2) { cout<< "Not enough arguments given to the function. Program ended" << endl; return 0; }
   else if (argc>3) { cout << "Too many arguments for the function. Program ended" << endl; return 0;}
   else {
        try
        {
           const string str = argv[1] ;
           std::size_t pos ;
           // http://en.cppreference.com/w/cpp/string/basic_string/stol
           const int value =  std::stoi( str, &pos, 10 ) ; // decimal integer literal
           if( pos == str.size() ) cout << "integer value of first input is " << value << '\n' ;
           else { cout << "partial conversion to int; value: " << value << " (non-numeric character found at position " << pos << ")\n" ; return -1; }
        }
        catch( const std::invalid_argument& ) { cout << "invalid characters: no conversion could be performed\n" ; return -1; }
        catch( const std::out_of_range& ) {  cout << "integer value is out of the range of int\n" ; return -1; }
        try
        {
           const string str = argv[2] ;
           std::size_t pos ;
           // http://en.cppreference.com/w/cpp/string/basic_string/stol
           const int value =  std::stoi( str, &pos, 10 ) ; // decimal integer literal
           if( pos == str.size() ) cout << "integer value of second input is " << value << '\n' ;
           else { cout << "partial conversion to int; value: " << value << " (non-numeric character found at position " << pos << ")\n" ; return -1; }
        }
        catch( const std::invalid_argument& ) { cout << "invalid characters: no conversion could be performed\n" ; return -1; }
        catch( const std::out_of_range& ) {  cout << "integer value is out of the range of int\n" ; return -1; }

    }

  	int fA = 0;
    fA = atoi(argv[1]); //nucleus number
    int fZ = 1;
    double momentum_ele = 10;
    double momentum_ion = 0;
    int particle_flag = 1; //fixed for pions
    double  ion_mass = 0;

    //mass
    const Double_t mass_p = 0.93827;//GeV
    const Double_t mass_n = 0.939566;//GeV
    const Double_t mass_u = 0.931494;//GeV

 //Definition of smearing
    double ele_mom_reso = 0.01;
    double had_mom_reso = 0.02;
    double ele_the_reso = 0.002;
    double had_the_reso = ele_the_reso;

    TFile *inputfile;

    int runnumber = atoi(argv[2]); //number of input root file
    if(fA ==12){
            cout << "Input nucleus is carbon" << endl;
          //  TString filename = Form("../../test/massproduction_CTEQ/EIC_A12_pion_10_600_1_%i.root",runnumber);
            TString filename = Form("../massproduction_CTEQfree/EIC_A12_pion_10_600_1_%i.root",runnumber);
            inputfile = new TFile(filename,"R");
            momentum_ele = 10; //GeV
            momentum_ion = 600; //GeV
            fZ = 6;
            ion_mass = fZ*mass_p+(fA-fZ)*mass_n;
            cout << "Input nucleus is carbon" << endl;
            if (inputfile->IsZombie()) { cout << "Input file not found. Exit program." << endl; return 0;}
    }
    else if(fA ==2) {
            cout << "Input nucleus is deuterium" << endl;
      //     TString filename = Form("../../test/massproduction_CTEQ/EIC_A2_pion_10_100_1_%i.root",runnumber);
            TString filename = Form("../massproduction_CTEQfree/EIC_A2_pion_10_100_1_%i.root",runnumber);
            inputfile = new TFile(filename,"R");
            momentum_ele = 10; //GeV
            momentum_ion = 100; //GeV
            fZ = 1;
            ion_mass = fZ*mass_p+(fA-fZ)*mass_n;
            if (inputfile->IsZombie()) { cout << "Input file not found. Exit program." << endl; return 0;}
    }
    else { cout << "no valid nucleus: " << fA << endl; return 0;}

    cout << "Resolution values for smearing: p_ele " << ele_mom_reso << " , p_had " << had_mom_reso <<
    " , theta " << ele_the_reso << " [rad]" << endl;

    cout << " A is " << fA << " , Z is " << fZ << " , ionmass is " << ion_mass << endl;
    TTree *T = (TTree*)inputfile->GetObjectChecked("T", "TTree");

    if (T->IsZombie()) { cout << "Tree not found. Exit program" << endl; return 0;}
    //Define
    Double_t Q2, W, Wp, x, y, z, pt, nu, s, gamma, epsilon,rapidity;
    Double_t theta_q, theta_s,phi_h,phi_s,mom_ele,mom_had,theta_ele, theta_had,phi_ele,phi_had;
    Double_t dxs_incl,dxs_hm,dxs_hp;
    Double_t px_ele, py_ele,pz_ele, px_had, py_had, pz_had, E_ele,E_had;
    Double_t weight_hp, weight_hm, weight_in;
    ULong64_t nsim = 0;

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
    T->SetBranchAddress("rapidity",&rapidity );
    T->SetBranchAddress("weight_hp",&weight_hp );
    T->SetBranchAddress("weight_hm",&weight_hm );
    T->SetBranchAddress("weight_in",&weight_in );
    T->SetBranchAddress("theta_q",&theta_q );
    T->SetBranchAddress("theta_s",&theta_s );
    T->SetBranchAddress("phi_h",&phi_h );
    T->SetBranchAddress("phi_s",&phi_s );
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
    T->SetBranchAddress("px_ele",&px_ele);
    T->SetBranchAddress("py_ele",&py_ele);
    T->SetBranchAddress("pz_ele",&pz_ele);
    T->SetBranchAddress("E_ele",&E_ele);
    T->SetBranchAddress("px_had",&px_had);
    T->SetBranchAddress("py_had",&py_had);
    T->SetBranchAddress("pz_had",&pz_had);
    T->SetBranchAddress("E_had",&E_had);

//Definitions of new branches and variables
    Double_t Q2_smeared, W_smeared, Wp_smeared, x_smeared, y_smeared, z_smeared, pt_smeared;
    Double_t nu_smeared, s_smeared, gamma_smeared, epsilon_smeared, rapidity_smeared;
    Double_t theta_q_smeared, theta_s_smeared, phi_h_smeared, phi_s_smeared;
    Double_t px_ele_smeared, py_ele_smeared, pz_ele_smeared, E_ele_smeared;
    Double_t px_had_smeared, py_had_smeared, pz_had_smeared, E_had_smeared;
    Double_t mom_ele_smeared,mom_had_smeared,theta_ele_smeared,theta_had_smeared,phi_ele_smeared,phi_had_smeared;

    TBranch *branch_Q2_smeared=T->Branch("Q2_smeared",&Q2_smeared,"Q2_smeared/D");
    TBranch *branch_W_smeared=T->Branch("W_smeared",&W_smeared,"W_smeared/D");
    TBranch *branch_Wp_smeared=T->Branch("Wp_smeared",&Wp_smeared,"Wp_smeared/D");
    TBranch *branch_x_smeared=T->Branch("x_smeared",&x_smeared,"x_smeared/D");
    TBranch *branch_y_smeared=T->Branch("y_smeared",&y_smeared,"y_smeared/D");
    TBranch *branch_z_smeared=T->Branch("z_smeared",&z_smeared,"z_smeared/D");
    TBranch *branch_pt_smeared=T->Branch("pt_smeared",&pt_smeared,"pt_smeared/D");
    TBranch *branch_nu_smeared=T->Branch("nu_smeared",&nu_smeared,"nu_smeared/D");
    TBranch *branch_s_smeared=T->Branch("s_smeared",&s_smeared,"s_smeared/D");
    TBranch *branch_gamma_smeared=T->Branch("gamma_smeared",&gamma_smeared,"gamma_smeared/D");
    TBranch *branch_epsilon_smeared=T->Branch("epsilon_smeared",&epsilon_smeared,"epsilon_smeared/D");
    TBranch *branch_rapidity_smeared=T->Branch("rapidity_smeared",&rapidity_smeared,"rapidity_smeared/D");
    TBranch *branch_theta_q_smeared=T->Branch("theta_q_smeared",&theta_q_smeared,"theta_q_smeared/D");
    TBranch *branch_theta_s_smeared=T->Branch("theta_s_smeared",&theta_s_smeared,"theta_s_smeared/D");
    TBranch *branch_phi_h_smeared=T->Branch("phi_h_smeared",&phi_h_smeared,"phi_h_smeared/D");
    TBranch *branch_phi_s_smeared=T->Branch("phi_s_smeared",&phi_s_smeared,"phi_s_smeared/D");
    TBranch *branch_px_ele_smeared=T->Branch("px_ele_smeared",&px_ele_smeared,"px_ele_smeared/D");
    TBranch *branch_py_ele_smeared=T->Branch("py_ele_smeared",&py_ele_smeared,"py_ele_smeared/D");
    TBranch *branch_pz_ele_smeared=T->Branch("pz_ele_smeared",&pz_ele_smeared,"pz_ele_smeared/D");
    TBranch *branch_E_ele_smeared=T->Branch("E_ele_smeared",&E_ele_smeared,"E_ele_smeared/D");
    TBranch *branch_px_had_smeared=T->Branch("px_had_smeared",&px_had_smeared,"px_had_smeared/D");
    TBranch *branch_py_had_smeared=T->Branch("py_had_smeared",&py_had_smeared,"py_had_smeared/D");
    TBranch *branch_pz_had_smeared=T->Branch("pz_had_smeared",&pz_had_smeared,"pz_had_smeared/D");
    TBranch *branch_E_had_smeared=T->Branch("E_had_smeared",&E_had_smeared,"E_had_smeared/D");
    TBranch *branch_mom_ele_smeared=T->Branch("mom_ele_smeared",&mom_ele_smeared,"mom_ele_smeared/D");
    TBranch *branch_mom_had_smeared=T->Branch("mom_had_smeared",&mom_had_smeared,"mom_had_smeared/D");
    TBranch *branch_theta_ele_smeared=T->Branch("theta_ele_smeared",&theta_ele_smeared,"theta_ele_smeared/D");
    TBranch *branch_theta_had_smeared=T->Branch("theta_had_smeared",&theta_had_smeared,"theta_had_smeared/D");
    TBranch *branch_phi_ele_smeared=T->Branch("phi_ele_smeared",&phi_ele_smeared,"phi_ele_smeared/D");
    TBranch *branch_phi_had_smeared=T->Branch("phi_had_smeared",&phi_had_smeared,"phi_had_smeared/D");



    ULong64_t N_Total=T->GetEntries();
    cout << "Total number of events " << N_Total << endl;
    N_Total=10000;
  //  cout << "MODIFIED NUMBER OF EVENTS TO 10M"<< endl;

    gRandom = new TRandom3();
    gRandom->SetSeed(0);// uses time for seed
    cout <<"SEED number for Smearing is " << gRandom->GetSeed()<<endl;


    TString new_filename = inputfile->GetName();
    new_filename.Replace(new_filename.Sizeof()-6,6,"_smeared.root");
    TFile *outfile = new TFile(new_filename,"RECREATE");
    TTree *outtree = (TTree*) T->CloneTree(0);
    outtree->SetDirectory(outfile);
    cout << new_filename << endl;

    TH1F *diffhad = new TH1F("diffhad"," ; #Delta p_{h}/p ; counts",1000,-10*had_mom_reso,10*had_mom_reso);
    TH1F *diffele = new TH1F("diffele"," ; #Delta p_{e}/p ; counts",1000,-10*ele_mom_reso,10*ele_mom_reso);
    TH1F *diffhad_the = new TH1F("diffhad_the"," ; #Delta theta_{h}; counts",1000,-10*had_the_reso,10*had_the_reso);
    TH1F *diffele_the = new TH1F("diffele_the"," ; #Delta theta_{e}; counts",1000,-10*ele_the_reso,10*ele_the_reso);

    TH1F *diffQ2= new TH1F("diffQ2"," ; Q^{2}_{MC}-Q^{2}_{smeared} [GeV^{2}]; counts ",1000,-1,1);
    TH1F *diffz = new TH1F("diffz"," ; z_{MC}-z_{smeared} ; counts ",1000,-1,1);
    TH1F *diffx = new TH1F("diffx"," ; x_{MC}-x_{smeared} ; counts ",1000,-0.05,0.05);
    TH1F *diffpt = new TH1F("diffpt"," ; p_{t,MC}-p_{t,smeared} [GeV]; counts ",1000,-1,1);
    TH2F *zvsQ2 = new TH2F("zvsQ2"," ; Q^{2}_{MC} ; z_{MC} ",100,0,11,100,0,1);
    TH2F *diffdQ2vQ2= new TH2F("diffdQ2vQ2"," ; Q^{2}_{MC}-Q^{2}_{smeared} [GeV^{2}]; Q^{2}_{MC} [GeV^{2}]",100,-1,1,100,1,11);
    TH2F *diffdzvz = new TH2F("diffdzvz"," ; z_{MC}-z_{smeared} ; z_{MC}",100,-1,1,100,0,1);
    TH2F *diffdxvx = new TH2F("diffdxvx"," ; x_{MC}-x_{smeared} ; x_{MC}",100,-0.05,0.05,100,0.0,0.2);
    TH2F *diffdptvpt = new TH2F("diffdptvpt"," ; p_{t,MC}-p_{t,smeared} [GeV]; p_{t,MC} [GeV]",100,-1,1,100,0,1);
    for (ULong64_t i=0;i<N_Total;i++){
       T->GetEntry(i);
       mom_had_smeared = gRandom->Gaus(1,had_mom_reso)*mom_had;
       mom_ele_smeared = gRandom->Gaus(1,ele_mom_reso)*mom_ele;
       theta_had_smeared = gRandom->Gaus(0,had_the_reso)+theta_had;
       theta_ele_smeared = gRandom->Gaus(0,ele_the_reso)+theta_ele;
       phi_ele_smeared = phi_ele;
       phi_had_smeared = phi_had;


       diffhad->Fill((mom_had_smeared - mom_had)/mom_had);
       diffele->Fill((mom_ele_smeared - mom_ele)/mom_ele);
       diffhad_the->Fill(theta_had_smeared - theta_had);
       diffele_the->Fill(theta_ele_smeared - theta_ele);

       sidis->SetKin(momentum_ele, momentum_ion,
                 mom_ele_smeared, theta_ele_smeared*DEG, phi_ele*DEG,
                 mom_had_smeared, theta_had_smeared*DEG, phi_had*DEG,
                 ion_mass, fA, fZ, particle_flag);

       Q2_smeared = sidis->fQ2;
       W_smeared = sidis->fW;
       Wp_smeared = sidis->fWp;
       x_smeared = sidis->fXb;
       y_smeared = sidis->fY;
       z_smeared = sidis->fZ_h;
       pt_smeared = sidis->fPt;
       nu_smeared = sidis->fNu;
       s_smeared = sidis->fS;
       gamma_smeared = sidis->fGamma;
       epsilon_smeared = sidis->fEpsilon;
       rapidity_smeared = sidis->fRapidity;
       theta_q_smeared = sidis->fTheta_q;
       theta_s_smeared = sidis->fTheta_s;
       phi_h_smeared = sidis->fPhi_h;
       phi_s_smeared = sidis->fPhi_s;
       px_ele_smeared = sidis->fPx_ele;
       py_ele_smeared = sidis->fPy_ele;
       pz_ele_smeared = sidis->fPz_ele;
       E_ele_smeared = sidis->fE_ele;
       px_had_smeared = sidis->fPx_had;
       py_had_smeared = sidis->fPy_had;
       pz_had_smeared = sidis->fPz_had;
       E_had_smeared = sidis->fE_had;


      // cout << "Event " << i << " x_new = " << x_smeared << " , Q2_new = " << Q2_smeared << " , z_new = " << z_smeared << " , x old = " << x << " , Q2 old = " << Q2 << endl;
       diffQ2->Fill(Q2-Q2_smeared);
       diffz->Fill(z-z_smeared);
       diffx->Fill(x-x_smeared);
       diffpt->Fill(pt-pt_smeared);
       zvsQ2->Fill(Q2,z,weight_hp);
       diffdQ2vQ2->Fill(Q2-Q2_smeared, Q2, weight_hp);
       diffdzvz->Fill(z-z_smeared,z, weight_hp);
       diffdxvx->Fill(x-x_smeared,x,weight_hp);
       diffdptvpt->Fill(pt-pt_smeared,pt, weight_hp);
       outtree->Fill();
       if(!(i%100000))
             cerr<<Form("--- Working on evt=%d",i)<<"\r";
    }

    outfile->cd();
  //  diffhad->Write();
  //  diffele->Write();
    outfile->Write();
    outfile->Close();

  /*  sidis->SetKin(momentum_ele, momentum_ion,
            mom_gen_ele, theta_gen_ele, phi_gen_ele,
            mom_gen_had, theta_gen_had, phi_gen_had,
            ion_mass, A, Z, particle_flag);

    mom_ele = sidis->fMom_ele; theta_ele = sidis->fTheta_ele; phi_ele = sidis->fPhi_ele;
    mom_had = sidis->fMom_had; theta_had = sidis->fTheta_had; phi_had = sidis->fPhi_had;
    theta_q=sidis->fTheta_q;  theta_s=sidis->fTheta_s; phi_s= sidis->fPhi_s; phi_h = sidis->fPhi_h;
    px_ele = sidis->fPx_ele;	py_ele = sidis->fPy_ele;	pz_ele = sidis->fPz_ele; E_ele = sidis->fE_ele;
    px_had = sidis->fPx_had;	py_had = sidis->fPy_had;	pz_had = sidis->fPz_had; E_had = sidis->fE_had;

    x=sidis->fXb; y=sidis->fY; z=sidis->fZ_h; Q2=sidis->fQ2; W=sidis->fW; Wp=sidis->fWp;
    s=sidis->fS; nu=sidis->fNu; pt=sidis->fPt; gamma=sidis->fGamma; epsilon=sidis->fEpsilon;*/
}
