#include <iostream> 
#include <fstream>
#include <cmath> 
#include <math.h> 
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TF1.h>
#include <TLorentzVector.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TMinuit.h>
#include <TPaveText.h>
#include <TText.h>
#include <TSystem.h>
#include <TArc.h>

using namespace std;

void acceptance(string input_filename)
{
gROOT->Reset();
gStyle->SetPalette(1);
gStyle->SetOptStat(0);
// gStyle->SetOptStat(11111);

const double DEG=180./3.1415926;

bool Is_PVDIS=false,Is_SIDIS_He3=false,Is_SIDIS_NH3=false,Is_DDVCS_JPsi_LH2=false,Is_DDVCS_PVDIS_LH2=false,Is_JPsi=false;
if (input_filename.find("SIDIS_He3",0) != string::npos) Is_SIDIS_He3=true;
else if (input_filename.find("SIDIS_NH3",0) != string::npos) Is_SIDIS_NH3=true;
else if (input_filename.find("DDVCS_JPsi_LH2",0) != string::npos) Is_DDVCS_JPsi_LH2=true;
else if (input_filename.find("JPsi",0) != string::npos) Is_JPsi=true;
else if (input_filename.find("DDVCS_PVDIS_LH2",0) != string::npos) Is_DDVCS_PVDIS_LH2=true;
else if (input_filename.find("PVDIS",0) != string::npos) Is_PVDIS=true;
else {cout << "not PVDIS or SIDIS or JPsi " << endl; return;}

///radius cut standard as particles can travel 30cm in Z before leave calorimeter
///SIDIS largeangle EC outerradius is at 140cm
double rout_cut_FA,rin_cut_FA,rout_cut_LA,rin_cut_LA;
if (Is_PVDIS){
  rout_cut_FA = 250;  //target at 10,ec front at 320 back at 370 with angle 36
  rin_cut_FA = 110;  //as the detector edge
  cout << " PVDIS rcut " << rin_cut_FA << " " << rout_cut_FA <<  endl;
}
else if (Is_SIDIS_He3 || Is_SIDIS_NH3){
  rout_cut_FA=220;     //target at -350,ec front at 415 with angle 15
//   rin_cut_FA=98;   //cut at the actual edge 
  rin_cut_FA=105;   //cut at the actual edge
//   rout_cut_FA=250;     //target at -350,ec front at 415 with angle 15  
//   rin_cut_FA=80;   //cut at the actual edge
  rout_cut_LA=127;   //target at -350,ec front at -65 with angle 24
  rin_cut_LA=83;   //cut at the actual edge
  cout << " SIDIS_He3 and SIDIS_NH3 rcut " << rin_cut_FA << " " << rout_cut_FA << " " << rin_cut_LA << " " << rout_cut_LA <<  endl;
}
else if (Is_JPsi){
    rout_cut_FA=220;     //target at -350,ec front at 405 with angle 15
//     rin_cut_FA=98;   	//cut at the actual edge
    rin_cut_FA=105;   	 //cut at the actual edge
    
//   rout_cut_LA=127.1;   //target at -360,ec front at -65 with angle 23.3    
//   rout_cut_LA=126.3;   //target at -350,ec front at -65 with angle 23.9
//   rout_cut_LA=125.9;   //target at -340,ec front at -65 with angle 24.6
//   rout_cut_LA=125.3;   //target at -330,ec front at -65 with angle 25.3    
//   rout_cut_LA=124.9;   //target at -320,ec front at -65 with angle 26.1
//   rout_cut_LA=124.5;   //target at -315,ec front at -65 with angle 26.5
  rout_cut_LA=135;   //target at -315,ec front at -65 with angle 28.4
//   rout_cut_LA=124.3;   //target at -310,ec front at -65 with angle 26.9
//      rout_cut_LA=124.0;   //target at -305,ec front at -65 with angle 27.3
//      rout_cut_LA=123.9;   //target at -300,ec front at -65 with angle 27.8
//      rout_cut_LA=123.8;   //target at -295,ec front at -65 with angle 28.3
//   rout_cut_LA=123.2;   //target at -290,ec front at -65 with angle 28.7
//      rout_cut_LA=122.6;   //target at -280,ec front at -65 with angle 29.7
//      rout_cut_LA=121.7;   //target at -270,ec front at -65 with angle 30.7
//   rout_cut_LA=120.1;   //target at -250,ec front at -65 with angle 33.0
  rin_cut_LA=83;   //cut at the actual edge  
  cout << " JPsi rcut " << rin_cut_FA << " " << rout_cut_FA << " " << rin_cut_LA << " " << rout_cut_LA <<  endl;
}
else if(Is_DDVCS_JPsi_LH2 || Is_DDVCS_PVDIS_LH2){
  rout_cut_FA=1000;     //target at -350,ec front at 415 with angle 15
  rin_cut_FA=0;   //cut at the actual edge
  rout_cut_LA=1000;   //target at -350,ec front at -65 with angle 24
  rin_cut_LA=0;   //cut at the actual edge
  cout << " DDVCS rcut " << rin_cut_FA << " " << rout_cut_FA << " " << rin_cut_LA << " " << rout_cut_LA <<  endl;
}
else {cout << "not PVDIS or SIDIS or JPsi or DDVCS" << endl; return;}

double theta_max,theta_min;
double vz_max,vz_min;
if (Is_SIDIS_He3){
    theta_max=40;
    theta_min=0;    
    vz_max=-320;    
    vz_min=-380;
}
else if(Is_SIDIS_NH3) {
    theta_max=50;
    theta_min=0; 
    vz_max=-345;    
    vz_min=-355;    
}
else if(Is_PVDIS) {
    theta_max=50;
    theta_min=10;  
    vz_max=40;    
    vz_min=-20;    
}
else if (Is_JPsi){
    theta_max=40;
    theta_min=0;    
    vz_max=-305;    
    vz_min=-325;
}
else if(Is_DDVCS_JPsi_LH2) {
    theta_max=40;
    theta_min=0;    
    vz_max=-305;    
    vz_min=-325;    
}
else if(Is_DDVCS_PVDIS_LH2) {
    theta_max=50;
    theta_min=10;    
    vz_max=40;    
    vz_min=-20;        
}
else {cout << "not PVDIS or SIDIS or JPsi " << endl; return;}

char the_filename[200];
sprintf(the_filename, "%s",input_filename.substr(0,input_filename.rfind(".")).c_str());

char output_filename[200];
sprintf(output_filename, "%s_output.root",the_filename);
TFile *outputfile=new TFile(output_filename, "recreate");

// int binfactor_theta=5;
// int binfactor_p=20;
// int binfactor_phi=1;

int binfactor_theta=2;
int binfactor_p=100;
int binfactor_phi=1;
int binfactor_vz=2;

TH2F *hacceptance_ThetaP_forwardangle,*hacceptance_ThetaP_largeangle,*hacceptance_ThetaP_overall;
TH2F *hacceptance_ThetaPhi_forwardangle,*hacceptance_ThetaPhi_largeangle;
TH2F *hacceptance_PhiP_forwardangle,*hacceptance_PhiP_largeangle;
TH3F *hacceptance_ThetaPhiP_forwardangle,*hacceptance_ThetaPhiP_largeangle;

TH1F *hgen_P=new TH1F("gen_P","gen_P",int(11*binfactor_p),0,11);
TH1F *hgen_Theta=new TH1F("gen_Theta","gen_Theta",int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max);
TH1F *hgen_Phi=new TH1F("gen_Phi","gen_Phi",int(360*binfactor_phi),-180,180);
TH2F *hgen_ThetaP=new TH2F("gen_ThetaP","gen_ThetaP",int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max,int(11*binfactor_p),0,11);
TH2F *hgen_ThetaPhi=new TH2F("gen_ThetaPhi","gen_ThetaPhi",int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max,int(360*binfactor_phi),-180,180);     
TH2F *hgen_PhiP=new TH2F("gen_PhiP","gen_PhiP",int(360*binfactor_phi),-180,180,int(11*binfactor_p),0,11);
TH3F *hgen_ThetaPhiP=new TH3F("gen_ThetaPhiP","gen_ThetaPhiP",int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max,180,-180,180,110,0,11);   

TH2F *hgen_ThetaVz=new TH2F("gen_ThetaVz","gen_ThetaVz",int((vz_max-vz_min)*binfactor_vz),vz_min,vz_max,int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max);   
TH2F *hgen_ThetaVr=new TH2F("gen_ThetaVr","gen_ThetaVr",20,0,1,int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max);

const int m=2;

TH1F *hflux_P[m],*hflux_Theta[m],*hflux_Phi[m];
TH2F *hflux_ThetaP[m],*hflux_ThetaPhi[m],*hflux_PhiP[m];
TH3F *hflux_ThetaPhiP[m];
TH2F *hflux_ThetaVz[m],*hflux_ThetaVr[m];
TH1F *hacceptance_P[m],*hacceptance_Theta[m],*hacceptance_Phi[m];;
TH2F *hacceptance_ThetaP[m],*hacceptance_ThetaPhi[m],*hacceptance_PhiP[m];
TH3F *hacceptance_ThetaPhiP[m];
TH2F *hacceptance_ThetaVz[m],*hacceptance_ThetaVr[m];

char *title[m]={"FA","LA"};
for(int i=0;i<m;i++){
   char hstname[100];
   
   sprintf(hstname,"flux_P_%i",i);
   hflux_P[i]=new TH1F(hstname,hstname,int(11*binfactor_p),0,11);
   sprintf(hstname,"flux_Theta_%i",i);
   hflux_Theta[i]=new TH1F(hstname,hstname,int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max); 
   sprintf(hstname,"flux_Phi_%i",i);
   hflux_Phi[i]=new TH1F(hstname,hstname,int(360*binfactor_phi),-180,180); 
   
   sprintf(hstname,"flux_ThetaP_%i",i);
   hflux_ThetaP[i]=new TH2F(hstname,hstname,int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max,int(11*binfactor_p),0,11);        
   hflux_ThetaP[i]->SetTitle(Form("particles detected by %s;vertex #theta (deg);vertex P (GeV)",title[i]));   
   sprintf(hstname,"flux_ThetaPhi_%i",i);   
   hflux_ThetaPhi[i]=new TH2F(hstname,hstname,int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max,int(360*binfactor_phi),-180,180); 
   hflux_ThetaPhi[i]->SetTitle(Form("particles detected by %s;vertex #theta (deg);vertex Phi (GeV)",title[i]));   
   sprintf(hstname,"flux_PhiP_%i",i);   
   hflux_PhiP[i]=new TH2F(hstname,hstname,int(360*binfactor_phi),-180,180,int(11*binfactor_p),0,11);
   hflux_PhiP[i]->SetTitle(Form("particles detected by %s;vertex #phi (deg);vertex P (GeV)",title[i]));  
   sprintf(hstname,"flux_ThetaPhiP_%i",i);   
   hflux_ThetaPhiP[i]=new TH3F(hstname,hstname,int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max,180,-180,180,110,0,11);   
   hflux_ThetaPhiP[i]->SetTitle(Form("particles detected by %s;vertex #theta (deg);vertex #phi (deg);vertex P (GeV)",title[i]));   
    
   sprintf(hstname,"flux_ThetaVz_%i",i);
   hflux_ThetaVz[i]=new TH2F(hstname,hstname,int((vz_max-vz_min)*binfactor_vz),vz_min,vz_max,int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max);   
   hflux_ThetaVz[i]->SetTitle("flux;vertex Z (cm);vertex #theta (deg)");      
   sprintf(hstname,"flux_ThetaVr_%i",i);
   hflux_ThetaVr[i]=new TH2F(hstname,hstname,20,0,1,int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max);
   hflux_ThetaVr[i]->SetTitle("flux;vertex R (cm);vertex #theta (deg)");            
  
   sprintf(hstname,"acceptance_P_%i",i);
   hacceptance_P[i]=new TH1F(hstname,hstname,int(11*binfactor_p),0,11);
   hacceptance_P[i]->SetTitle(Form("acceptance by %s;vertex P (GeV);",title[i]));
   sprintf(hstname,"acceptance_Theta_%i",i);
   hacceptance_Theta[i]=new TH1F(hstname,hstname,int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max);    
   hacceptance_Theta[i]->SetTitle(Form("acceptance by %s;vertex #theta (deg);",title[i]));   
   sprintf(hstname,"acceptance_Phi_%i",i);
   hacceptance_Phi[i]=new TH1F(hstname,hstname,int(360*binfactor_phi),-180,180);    
   hacceptance_Phi[i]->SetTitle(Form("acceptance by %s;vertex #phi (deg));",title[i]));   
   
   sprintf(hstname,"acceptance_ThetaP_%i",i);
   hacceptance_ThetaP[i]=new TH2F(hstname,hstname,int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max,int(11*binfactor_p),0,11);     
   hacceptance_ThetaP[i]->SetTitle(Form("acceptance by %s;vertex #theta (deg);vertex P (GeV)",title[i]));
   sprintf(hstname,"acceptance_ThetaPhi_%i",i);
   hacceptance_ThetaPhi[i]=new TH2F(hstname,hstname,int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max,int(360*binfactor_phi),-180,180);     
   hacceptance_ThetaPhi[i]->SetTitle(Form("acceptance by %s;vertex #theta (deg);vertex #phi (deg)",title[i]));
   sprintf(hstname,"acceptance_PhiP_%i",i);
   hacceptance_PhiP[i]=new TH2F(hstname,hstname,int(360*binfactor_phi),-180,180,int(11*binfactor_p),0,11);
   hacceptance_PhiP[i]->SetTitle(Form("acceptance by %s;vertex #phi (deg);vertex P (GeV)",title[i]));      
   sprintf(hstname,"acceptance_ThetaPhiP_%i",i);   
   hacceptance_ThetaPhiP[i]=new TH3F(hstname,hstname,int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max,180,-180,180,110,0,11);
   hacceptance_ThetaPhiP[i]->SetTitle(Form("acceptance by %s;vertex #theta (deg);vertex #phi (deg);vertex P (GeV)",title[i]));

   sprintf(hstname,"acceptance_ThetaVz_%i",i);
   hacceptance_ThetaVz[i]=new TH2F(hstname,hstname,int((vz_max-vz_min)*binfactor_vz),vz_min,vz_max,int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max);   
   hacceptance_ThetaVz[i]->SetTitle("acceptance;vertex Z (cm);vertex #theta (deg)");      
   sprintf(hstname,"acceptance_ThetaVr_%i",i);
   hacceptance_ThetaVr[i]=new TH2F(hstname,hstname,20,0,1,int((theta_max-theta_min)*binfactor_theta),theta_min,theta_max);
   hacceptance_ThetaVr[i]->SetTitle("acceptance;vertex R (cm);vertex #theta (deg)");            
}

TH2F *hhit_PR[2][2];
hhit_PR[0][0]=new TH2F("hit_PR_0_0","all electrons reach FAEC and GEM;R(cm);Mom(GeV)",160,80,240,1100,0,11);
hhit_PR[0][1]=new TH2F("hit_PR_0_1","electrons (Q2>1) reach FAEC and GEM;R(cm);Mom(GeV)",160,80,240,1100,0,11);
hhit_PR[1][0]=new TH2F("hit_PR_1_0","all electrons reach LAEC and GEM;R(cm);Mom(GeV)",80,70,150,1100,0,11);
hhit_PR[1][1]=new TH2F("hit_PR_1_1","electrons (Q2>1) reach LAEC and GEM;R(cm);Mom(GeV)",80,70,150,1100,0,11);  

TH2F *hhit_rz=new TH2F("hit_rz","hit_rz",1000,-400,600,300,0,300);  

TH2F *hlinearity_GEM34=new TH2F("linearity_GEM34","linearity_GEM34",1000,theta_min,theta_max,2200,0,11);

TH1F *hmissingGEM_forwardangle=new TH1F("missingGEM_forwardangle","missingGEM_forwardangle",6,0.5,6.5);
TH1F *hmissingGEM_largeangle=new TH1F("missingGEM_largeangle","missingGEM_largeangle",6,0.5,6.5);

const int n=15;

TH2F *hhit_rMom[n];
TH2F *hhit_phidiffMom[n],*hhit_thetadiffMom[n],*hhit_anglediffMom[n];
TH2F *hmom_phidiffMom[n],*hmom_thetadiffMom[n],*hmom_anglediffMom[n];

for(int i=0;i<n;i++){
   char hstname[100];

   sprintf(hstname,"hit_rMom_%i",i);
   hhit_rMom[i]=new TH2F(hstname,hstname,300,0,300,int(11*binfactor_p),0,11);  
   
//    sprintf(hstname,"hit_phidiffMom_%i",i);
//    hhit_phidiffMom[i]=new TH2F(hstname,hstname,7200,-360,360,int(11*binfactor_p),0,11);  
//    sprintf(hstname,"hit_thetadiffMom_%i",i);
//    hhit_thetadiffMom[i]=new TH2F(hstname,hstname,3600,-180,180,int(11*binfactor_p),0,11);  
//    sprintf(hstname,"hit_anglediffMom_%i",i);
//    hhit_anglediffMom[i]=new TH2F(hstname,hstname,1800,0,180,int(11*binfactor_p),0,11);     
//    
//    sprintf(hstname,"mom_phidiffMom_%i",i);
//    hmom_phidiffMom[i]=new TH2F(hstname,hstname,7200,-360,360,int(11*binfactor_p),0,11);  
//    sprintf(hstname,"mom_thetadiffMom_%i",i);
//    hmom_thetadiffMom[i]=new TH2F(hstname,hstname,3600,-180,180,int(11*binfactor_p),0,11);  
//    sprintf(hstname,"mom_anglediffMom_%i",i);
//    hmom_anglediffMom[i]=new TH2F(hstname,hstname,1800,0,180,int(11*binfactor_p),0,11);        
   sprintf(hstname,"hit_phidiffMom_%i",i);
   hhit_phidiffMom[i]=new TH2F(hstname,hstname,7200,-10,10,int(11*binfactor_p),0,11);  
   sprintf(hstname,"hit_thetadiffMom_%i",i);
   hhit_thetadiffMom[i]=new TH2F(hstname,hstname,3600,-10,10,int(11*binfactor_p),0,11);  
   sprintf(hstname,"hit_anglediffMom_%i",i);
   hhit_anglediffMom[i]=new TH2F(hstname,hstname,1800,0,20,int(11*binfactor_p),0,11);     
   
   sprintf(hstname,"mom_phidiffMom_%i",i);
   hmom_phidiffMom[i]=new TH2F(hstname,hstname,100,-360,360,int(11*binfactor_p),0,11);  
   sprintf(hstname,"mom_thetadiffMom_%i",i);
   hmom_thetadiffMom[i]=new TH2F(hstname,hstname,100,-5,5,int(11*binfactor_p),0,11);  
   sprintf(hstname,"mom_anglediffMom_%i",i);
   hmom_anglediffMom[i]=new TH2F(hstname,hstname,100,0,20,int(11*binfactor_p),0,11);        
   
}

TH2F *hhit_xy[n][2],*hhit_PhiR[n][2],*hhit_PhiZ[n][2],*hPloss[n][2];
for(int i=0;i<n;i++){
  for(int j=0;j<2;j++){
   char hstname[100];
   sprintf(hstname,"hit_xy_%i_%i",i,j);
   hhit_xy[i][j]=new TH2F(hstname,hstname,600,-300,300,600,-300,300);        
   sprintf(hstname,"hit_PhiR_%i_%i",i,j);
   hhit_PhiR[i][j]=new TH2F(hstname,hstname,int(360*binfactor_phi),-180,180,300,0,300);       
   sprintf(hstname,"hit_PhiZ_%i_%i",i,j);
   hhit_PhiZ[i][j]=new TH2F(hstname,hstname,int(360*binfactor_phi),-180,180,400,-300,700);          
   sprintf(hstname,"Ploss_%i_%i",i,j);   
   hPloss[i][j]=new TH2F(hstname,hstname,110,0,1.1,110,0,11);
}}
  
const int Nplate=22;
TH2F *hbaffleplate[Nplate],*hbaffleplate_observer[Nplate];
for (int i=0;i<Nplate;i++){
  hbaffleplate[i]=new TH2F(Form("baffleplate_%i", i),Form("baffleplate_%i", i),300,-150,150,300,-150,150);
  hbaffleplate_observer[i]=new TH2F(Form("hbaffleplate_observer_%i", i),Form("hbaffleplate_observer_%i", i),300,-150,150,300,-150,150);  
}

  TFile *file=new TFile(input_filename.c_str());
    if (file->IsZombie()) {
       cout << "Error opening file" << input_filename << endl;
       exit(-1);
    }
    else cout << "open file " << input_filename << endl;

TTree *tree_header = (TTree*) file->Get("header");
vector <int> *evn=0,*evn_type=0;
vector <double> *beamPol=0;
vector <double> *var1=0,*var2=0,*var3=0,*var4=0,*var5=0,*var6=0,*var7=0,*var8=0;
tree_header->SetBranchAddress("evn",&evn);
tree_header->SetBranchAddress("evn_type",&evn_type);
tree_header->SetBranchAddress("beamPol",&beamPol);
tree_header->SetBranchAddress("var1",&var1);
tree_header->SetBranchAddress("var2",&var2);
tree_header->SetBranchAddress("var3",&var3);
tree_header->SetBranchAddress("var4",&var4);
tree_header->SetBranchAddress("var5",&var5);
tree_header->SetBranchAddress("var6",&var6);
tree_header->SetBranchAddress("var7",&var7);
tree_header->SetBranchAddress("var8",&var8);

TTree *tree_generated = (TTree*) file->Get("generated");
vector <int> *gen_pid=0;
vector <double> *gen_px=0,*gen_py=0,*gen_pz=0,*gen_vx=0,*gen_vy=0,*gen_vz=0;
tree_generated->SetBranchAddress("pid",&gen_pid);
tree_generated->SetBranchAddress("px",&gen_px);
tree_generated->SetBranchAddress("py",&gen_py);
tree_generated->SetBranchAddress("pz",&gen_pz);
tree_generated->SetBranchAddress("vx",&gen_vx);
tree_generated->SetBranchAddress("vy",&gen_vy);
tree_generated->SetBranchAddress("vz",&gen_vz);

TTree *tree_flux = (TTree*) file->Get("flux");
vector<int> *flux_id=0,*flux_hitn=0;
// vector<double> *flux_id=0;
// vector<int> *flux_hitn=0;
vector<int> *flux_pid=0,*flux_mpid=0,*flux_tid=0,*flux_mtid=0,*flux_otid=0;
vector<double> *flux_trackE=0,*flux_totEdep=0,*flux_avg_x=0,*flux_avg_y=0,*flux_avg_z=0,*flux_avg_lx=0,*flux_avg_ly=0,*flux_avg_lz=0,*flux_px=0,*flux_py=0,*flux_pz=0,*flux_vx=0,*flux_vy=0,*flux_vz=0,*flux_mvx=0,*flux_mvy=0,*flux_mvz=0,*flux_avg_t=0;
tree_flux->SetBranchAddress("hitn",&flux_hitn);
tree_flux->SetBranchAddress("id",&flux_id);
tree_flux->SetBranchAddress("pid",&flux_pid);
tree_flux->SetBranchAddress("mpid",&flux_mpid);
tree_flux->SetBranchAddress("tid",&flux_tid);
tree_flux->SetBranchAddress("mtid",&flux_mtid);
tree_flux->SetBranchAddress("otid",&flux_otid);
tree_flux->SetBranchAddress("trackE",&flux_trackE);
tree_flux->SetBranchAddress("totEdep",&flux_totEdep);
tree_flux->SetBranchAddress("avg_x",&flux_avg_x);
tree_flux->SetBranchAddress("avg_y",&flux_avg_y);
tree_flux->SetBranchAddress("avg_z",&flux_avg_z);
tree_flux->SetBranchAddress("avg_lx",&flux_avg_lx);
tree_flux->SetBranchAddress("avg_ly",&flux_avg_ly);
tree_flux->SetBranchAddress("avg_lz",&flux_avg_lz);
tree_flux->SetBranchAddress("px",&flux_px);
tree_flux->SetBranchAddress("py",&flux_py);
tree_flux->SetBranchAddress("pz",&flux_pz);
tree_flux->SetBranchAddress("vx",&flux_vx);
tree_flux->SetBranchAddress("vy",&flux_vy);
tree_flux->SetBranchAddress("vz",&flux_vz);
tree_flux->SetBranchAddress("mvx",&flux_mvx);
tree_flux->SetBranchAddress("mvy",&flux_mvy);
tree_flux->SetBranchAddress("mvz",&flux_mvz);
tree_flux->SetBranchAddress("avg_t",&flux_avg_t);

// TTree *tree_solid_lgc = (TTree*) file->Get("solid_lgc");
// vector<int> *solid_lgc_hitn=0;
// vector<int> *solid_lgc_sector=0,*solid_lgc_pmt=0,*solid_lgc_pixel=0,*solid_lgc_nphe=0;
// vector<double> *solid_lgc_avg_t=0;
// tree_solid_lgc->SetBranchAddress("hitn",&solid_lgc_hitn);
// tree_solid_lgc->SetBranchAddress("sector",&solid_lgc_sector);
// tree_solid_lgc->SetBranchAddress("pmt",&solid_lgc_pmt);
// tree_solid_lgc->SetBranchAddress("pixel",&solid_lgc_pixel);
// tree_solid_lgc->SetBranchAddress("nphe",&solid_lgc_nphe);
// tree_solid_lgc->SetBranchAddress("avg_t",&solid_lgc_avg_t);
// 
// TTree *tree_solid_hgc = (TTree*) file->Get("solid_hgc");
// vector<int> *solid_hgc_id=0,*solid_hgc_hitn=0;
// vector<int> *solid_hgc_pid=0,*solid_hgc_mpid=0,*solid_hgc_tid=0,*solid_hgc_mtid=0,*solid_hgc_otid=0;
// vector<double> *solid_hgc_trackE=0,*solid_hgc_totEdep=0,*solid_hgc_avg_x=0,*solid_hgc_avg_y=0,*solid_hgc_avg_z=0,*solid_hgc_avg_lx=0,*solid_hgc_avg_ly=0,*solid_hgc_avg_lz=0,*solid_hgc_px=0,*solid_hgc_py=0,*solid_hgc_pz=0,*solid_hgc_vx=0,*solid_hgc_vy=0,*solid_hgc_vz=0,*solid_hgc_mvx=0,*solid_hgc_mvy=0,*solid_hgc_mvz=0,*solid_hgc_avg_t=0;
// tree_solid_hgc->SetBranchAddress("hitn",&solid_hgc_hitn);
// tree_solid_hgc->SetBranchAddress("id",&solid_hgc_id);
// tree_solid_hgc->SetBranchAddress("pid",&solid_hgc_pid);
// tree_solid_hgc->SetBranchAddress("mpid",&solid_hgc_mpid);
// tree_solid_hgc->SetBranchAddress("tid",&solid_hgc_tid);
// tree_solid_hgc->SetBranchAddress("mtid",&solid_hgc_mtid);
// tree_solid_hgc->SetBranchAddress("otid",&solid_hgc_otid);
// tree_solid_hgc->SetBranchAddress("trackE",&solid_hgc_trackE);
// tree_solid_hgc->SetBranchAddress("totEdep",&solid_hgc_totEdep);
// tree_solid_hgc->SetBranchAddress("avg_x",&solid_hgc_avg_x);
// tree_solid_hgc->SetBranchAddress("avg_y",&solid_hgc_avg_y);
// tree_solid_hgc->SetBranchAddress("avg_z",&solid_hgc_avg_z);
// tree_solid_hgc->SetBranchAddress("avg_lx",&solid_hgc_avg_lx);
// tree_solid_hgc->SetBranchAddress("avg_ly",&solid_hgc_avg_ly);
// tree_solid_hgc->SetBranchAddress("avg_lz",&solid_hgc_avg_lz);
// tree_solid_hgc->SetBranchAddress("px",&solid_hgc_px);
// tree_solid_hgc->SetBranchAddress("py",&solid_hgc_py);
// tree_solid_hgc->SetBranchAddress("pz",&solid_hgc_pz);
// tree_solid_hgc->SetBranchAddress("vx",&solid_hgc_vx);
// tree_solid_hgc->SetBranchAddress("vy",&solid_hgc_vy);
// tree_solid_hgc->SetBranchAddress("vz",&solid_hgc_vz);
// tree_solid_hgc->SetBranchAddress("mvx",&solid_hgc_mvx);
// tree_solid_hgc->SetBranchAddress("mvy",&solid_hgc_mvy);
// tree_solid_hgc->SetBranchAddress("mvz",&solid_hgc_mvz);
// tree_solid_hgc->SetBranchAddress("avg_t",&solid_hgc_avg_t);

// cout << tree_flux->GetEntries() << " " << tree_header->GetEntries() << " " << tree_generated->GetEntries() << endl;

int nevent = (int)tree_generated->GetEntries();
int nselected = 0;
cout << "nevent " << nevent << endl;

int counter_decay=0;
int counter_loss=0;

double rate=1;

for (Int_t i=0;i<nevent;i++) {
  cout << i << "\r";
//   cout << i << "\n";

  tree_header->GetEntry(i);
//   rate=var8->at(0);
  
  tree_generated->GetEntry(i);  
  
      TVector3 coor_gen(0,0,0);
      TVector3 mom_gen(0,0,0);
      
      TLorentzVector beam(0,0,11000,11000);      
      TLorentzVector mom4_gen(0,0,0,0);
      double Q2=0;
      
  int pid_gen=0;
  double theta_gen=0,phi_gen=0,p_gen=0,px_gen=0,py_gen=0,pz_gen=0,vx_gen=0,vy_gen=0,vz_gen=0;         
  for (int j=0;j<gen_pid->size();j++) {
//       cout << gen_pid->at(j) << " " << gen_px->at(j) << " " << gen_py->at(j) << " " << gen_pz->at(j) << " " << gen_vx->at(j) << " " << gen_vy->at(j) << " " << gen_vz->at(j) << endl; 
      pid_gen=gen_pid->at(j);
      px_gen=gen_px->at(j);    	//in MeV
      py_gen=gen_py->at(j);	//in MeV
      pz_gen=gen_pz->at(j);     //in MeV
      vx_gen=gen_vx->at(j);    	//in mm
      vy_gen=gen_vy->at(j);	//in mm
      vz_gen=gen_vz->at(j);     //in mm

      coor_gen.SetXYZ(vx_gen,vy_gen,vz_gen);
      mom_gen.SetXYZ(px_gen,py_gen,pz_gen);
//       p_gen=sqrt(px_gen*px_gen+py_gen*py_gen+pz_gen*pz_gen);
//       theta_gen=acos(pz_gen/p_gen);  	//in rad
//       phi_gen=atan2(py_gen,px_gen);     //in rad
      p_gen=mom_gen.Mag();
      theta_gen=mom_gen.Theta()*DEG;
      phi_gen=mom_gen.Phi()*DEG;     //in deg

      if (pid_gen==11) {
	mom4_gen.SetXYZT(px_gen,py_gen,pz_gen,sqrt(0.511*0.511+p_gen*p_gen));                
	Q2=-(beam-mom4_gen).Mag2()/1e6;
      }
      
      hgen_P->Fill(p_gen/1e3);      
      hgen_Theta->Fill(theta_gen);            
      hgen_Phi->Fill(phi_gen);      
      hgen_ThetaP->Fill(theta_gen,p_gen/1e3);      
      hgen_ThetaPhi->Fill(theta_gen,phi_gen);            
      hgen_PhiP->Fill(phi_gen,p_gen/1e3);
      hgen_ThetaPhiP->Fill(theta_gen,phi_gen,p_gen/1e3);   
      
      hgen_ThetaVz->Fill(vz_gen/1e1,theta_gen,rate);
      hgen_ThetaVr->Fill(sqrt(vx_gen*vx_gen+vy_gen*vy_gen)/1e1,theta_gen,rate);      
  }  
  
    tree_flux->GetEntry(i);    

    //check all hits, begin ============================================================    
      TVector3 coor_acc(0,0,0);
      TVector3 mom_acc(0,0,0);
      TVector3 coor_GEM3(0,0,0);      
      TVector3 mom_GEM3(0,0,0);      
      TVector3 coor_GEM4(0,0,0);      
      TVector3 mom_GEM4(0,0,0);     
      
      double R_hit_EC=0;
      
    bool Is_out=false;
    bool Is_decay=false;
    int acc[n]={0};
    for (Int_t j=0;j<flux_hitn->size();j++) {
//           cout << j << " !!! " << flux_id->at(j) << " " << flux_pid->at(j) << " " << flux_mpid->at(j) << " " << flux_tid->at(j) << " " << flux_mtid->at(j) << " " << flux_trackE->at(j) << " " << flux_totEdep->at(j) << " " << flux_avg_x->at(j) << " " << flux_avg_y->at(j) << " " << flux_avg_z->at(j) << " " << flux_avg_lx->at(j) << " " << flux_avg_ly->at(j) << " " << flux_avg_lz->at(j) << " " << flux_px->at(j) << " " << flux_py->at(j) << " " << flux_pz->at(j) << " " << flux_vx->at(j) << " " << flux_vy->at(j) << " " << flux_vz->at(j) << " " << flux_mvx->at(j) << " " << flux_mvy->at(j) << " " << flux_mvz->at(j) << " " << flux_avg_t->at(j) << endl;           
    
    // remove decay product
    if (flux_tid->at(j)!= 1) {
//     if (flux_pid->at(j)!= pid_gen) {
//     if (abs(flux_pid->at(j)) != 13 && abs(flux_pid->at(j)) != 211 ) {      
      Is_decay=true; 
      counter_decay++;            
// 	cout << "pid " << pid_gen << " change to " << flux_pid->at(j) << ", "; 
      continue;
    }
    
    double p_flux=sqrt(flux_px->at(j)*flux_px->at(j)+flux_py->at(j)*flux_py->at(j)+flux_pz->at(j)*flux_pz->at(j));
    if (p_flux/p_gen<0.8) {
      counter_loss++;            
// 	cout << "mom " << p_gen << " change to " << p_flux << " at " << flux_id->at(j) << endl;       
//       continue;
    }
    
    int hit_id=-1;
    if (flux_id->at(j)==1110000) hit_id=0;
    else if (flux_id->at(j)==1210000) hit_id=1;
    else if (flux_id->at(j)==1310000) hit_id=2;
    else if (flux_id->at(j)==1410000) hit_id=3;
    else if (flux_id->at(j)==1510000) hit_id=4;
    else if (flux_id->at(j)==1610000) hit_id=5;    
    else if (flux_id->at(j)==2110000) hit_id=6;
    else if (flux_id->at(j)==2210000) hit_id=7;
    else if (flux_id->at(j)==5110000) hit_id=8;
    else if (flux_id->at(j)==5210000) hit_id=9;	
    else if (flux_id->at(j)==3110000) hit_id=10;
    else if (flux_id->at(j)==3210000) hit_id=11;
    else if (flux_id->at(j)==4110000) hit_id=12;	        
    else if (flux_id->at(j)==6110000) hit_id=13;
    else if (flux_id->at(j)==6100000) hit_id=14;	    
//     else if (flux_id->at(j)==6110000) hit_id=8;
//     else if (flux_id->at(j)==6120000) hit_id=8;    
//     else if (flux_id->at(j)==6130000) hit_id=8;        
//     else if (flux_id->at(j)==6140000) hit_id=8;            
//     else if (flux_id->at(j)==6210000) hit_id=9;
//     else if (flux_id->at(j)==6100000) hit_id=9;    
    else continue;  //skip other detector for now

//     if (hit_id==9) cout << "flux_id->at(j) " << flux_avg_z->at(j) << endl;

//     int detector_ID=flux_id->at(j)/1000000;
//     int subdetector_ID=(flux_id->at(j)%1000000)/100000;
//     int subsubdetector_ID=((flux_id->at(j)%1000000)%100000)/10000;
    int detector_ID=int(flux_id->at(j))/1000000;
    int subdetector_ID=(int(flux_id->at(j))%1000000)/100000;
    int subsubdetector_ID=((int(flux_id->at(j))%1000000)%100000)/10000;
//     cout << detector_ID << " " << subdetector_ID << " "  << subsubdetector_ID << endl;  
        
    coor_acc.SetXYZ(flux_avg_x->at(j),flux_avg_y->at(j),flux_avg_z->at(j));
    mom_acc.SetXYZ(flux_px->at(j),flux_py->at(j),flux_pz->at(j));   
    
    if (hit_id==2){ //hit on GEM3
      coor_GEM3.SetXYZ(flux_avg_x->at(j),flux_avg_y->at(j),flux_avg_z->at(j));    
      mom_GEM3.SetXYZ(flux_px->at(j),flux_py->at(j),flux_pz->at(j));   
    }
    if (hit_id==3){ //hit on GEM4
      coor_GEM4.SetXYZ(flux_avg_x->at(j),flux_avg_y->at(j),flux_avg_z->at(j));    
      mom_GEM4.SetXYZ(flux_px->at(j),flux_py->at(j),flux_pz->at(j));   
    }
    
     double hit_r=sqrt(pow(flux_avg_x->at(j),2)+pow(flux_avg_y->at(j),2));
     double hit_y=flux_avg_y->at(j),hit_x=flux_avg_x->at(j),hit_z=flux_avg_z->at(j);
     
//      double hit_phi=atan2(flux_avg_y->at(j)-vy_gen,flux_avg_x->at(j)-vx_gen);
// //      double hit_theta=atan((hit_r/1e1-sqrt(vx_gen*vx_gen+vy_gen*vy_gen))/(320-vz_gen))*DEG;        
//      double hit_theta=atan(sqrt(pow(flux_avg_x->at(j)-vx_gen,2)+pow(flux_avg_y->at(j)-vy_gen,2))/(flux_avg_z->at(j)-vz_gen))*DEG;

     double hit_phi=(coor_acc-coor_gen).Phi()*DEG;
     double hit_theta=(coor_acc-coor_gen).Theta()*DEG;     
//      cout << hit_phi << " " << hit_theta << endl;
     
     hhit_rz->Fill(flux_avg_z->at(j)/1e1,hit_r/1e1);
     
      hhit_phidiffMom[hit_id]->Fill((hit_phi-phi_gen),p_gen/1e3);
      hhit_thetadiffMom[hit_id]->Fill((hit_theta-theta_gen),p_gen/1e3);
      hhit_anglediffMom[hit_id]->Fill((coor_acc-coor_gen).Angle(mom_gen)*DEG,p_gen/1e3);
      
      hmom_phidiffMom[hit_id]->Fill((mom_acc.Phi()-mom_gen.Phi())*DEG,p_gen/1e3);
      hmom_thetadiffMom[hit_id]->Fill((mom_acc.Theta()-mom_gen.Theta())*DEG,p_gen/1e3);
      hmom_anglediffMom[hit_id]->Fill(mom_acc.Angle(mom_gen)*DEG,p_gen/1e3);      
      
      hhit_rMom[hit_id]->Fill(hit_r/1e1,p_gen/1e3);    
      
     hhit_xy[hit_id][0]->Fill(hit_x/1e1,hit_y/1e1);    
     hhit_PhiR[hit_id][0]->Fill(hit_phi,hit_r/1e1);
     hhit_PhiZ[hit_id][0]->Fill(hit_phi,hit_z/1e1);
     hPloss[hit_id][0]->Fill(p_flux/p_gen,p_gen/1e3);      
     
    //check hit on baffleplate
//     if (detector_ID==0 && subsubdetector_ID==0) hbaffleplate[subdetector_ID-1]->Fill(flux_avg_x->at(j)/1e1,flux_avg_y->at(j)/1e1);
//     if (detector_ID==0 && subsubdetector_ID==1) hbaffleplate_observer[subdetector_ID-1]->Fill(flux_avg_x->at(j)/1e1,flux_avg_y->at(j)/1e1);
    if (flux_id->at(j)<100) hbaffleplate_observer[flux_id->at(j)-1]->Fill(hit_x/1e1,hit_y/1e1);    
    
    //check hit on GEM
    if (detector_ID==1) {
      // some low mom tracks spiral and travel back in field and go through one plane twice or more,   flux bank average these steps to get hit position which is wrong and can be outside of the virtual plane.
//       for PVDIS
//  my @Rin = (48,59,65,105,109);
//  my @Rout = (122,143,143,230,237);	      
//       for SIDIS
//        my @Rin = (36,21,25,32,42,55);
//        my @Rout = (87,98,112,135,100,123);
      double Rin[6]={0,0,0,0,0,0},Rout[6]={0,0,0,0,0,0};
      if (Is_PVDIS || Is_DDVCS_PVDIS_LH2){
	Rin[0]=48;Rin[1]=59;Rin[2]=65;Rin[3]=105;Rin[4]=109;Rin[5]=0;     
        Rout[0]=122;Rout[1]=143;Rout[2]=143;Rout[3]=230;Rout[4]=237;Rout[5]=300;
      }
      else {
	Rin[0]=36;Rin[1]=21;Rin[2]=25;Rin[3]=32;Rin[4]=42;Rin[5]=55;     
        Rout[0]=87;Rout[1]=98;Rout[2]=112;Rout[3]=135;Rout[4]=100;Rout[5]=123;
      }

      if (hit_r/1e1<Rin[subdetector_ID-1] || Rout[subdetector_ID-1]<hit_r/1e1) {
	Is_out=true;
// 	cout << flux_id->at(j) << endl; 	
	continue;	
      }  
                      
    }       

    //check hit on EC    
    if (detector_ID==3) {  
      //cut out some part of viritual plane to consider "good" part of EC
      if ((detector_ID==3 && subdetector_ID==1) && (hit_r/1e1 < rin_cut_FA || rout_cut_FA < hit_r/1e1)) continue;        
      if ((detector_ID==3 && subdetector_ID==2) && (hit_r/1e1 < rin_cut_LA || rout_cut_LA < hit_r/1e1)) continue;
//       if ((detector_ID==3 && subdetector_ID==1) && (hit_r/1e1 < rin_cut_FA)) continue;        
//       if ((detector_ID==3 && subdetector_ID==2) && (hit_r/1e1 < rin_cut_LA)) continue;

//       if(detector_ID==3 && subdetector_ID==1){	
// 	if (105<=hit_r/1e1 && hit_r/1e1<120)     {if(p_gen/1e3<5) continue;}
// 	else if(120<=hit_r/1e1 && hit_r/1e1<135) {if(p_gen/1e3<4) continue;}
// 	else if(135<=hit_r/1e1 && hit_r/1e1<149) {if(p_gen/1e3<3) continue;}
// 	else if(149<=hit_r/1e1 && hit_r/1e1<163) {if(p_gen/1e3<3) continue;}
// 	else if(163<=hit_r/1e1 && hit_r/1e1<177) {if(p_gen/1e3<1.5) continue;}
// 	else if(177<=hit_r/1e1 && hit_r/1e1<191) {if(p_gen/1e3<1.5) continue;}
// 	else if(191<=hit_r/1e1 && hit_r/1e1<235) {if(p_gen/1e3<1) continue;}
// 	else {cout << "hit_r " << hit_r << " out of range"  << endl;}	
//       }      
      
      R_hit_EC=hit_r/1e1;
    }
    
    if (Is_PVDIS){
      acc[hit_id]=1;      
    }
    else if (Is_SIDIS_He3){
//       if(theta_gen<8) continue;  //only acceptance 8 deg above      
      acc[hit_id]=1;
    }
    else if (Is_SIDIS_NH3){     
      // only acceptance area not in sheet of flame
      
      //wrong ptarget field
//       if ((detector_ID==3 && subdetector_ID==1) && ((-74<hit_phi && hit_phi<-38 && hit_r/1e1<195)||(-92<hit_phi && hit_phi<-88 && hit_r/1e1<120)||(50<hit_phi && hit_phi<80 && hit_r/1e1<195))) continue;                
//       if ((detector_ID==3 && subdetector_ID==2) && ((-85<hit_phi && hit_phi<-60)||(65<hit_phi && hit_phi<85))) continue;                
//       
//       if ((detector_ID==1 && subdetector_ID==1) && ((-93<hit_phi && hit_phi<-82)||(80<hit_phi && hit_phi<92))) continue;          
//       if ((detector_ID==1 && subdetector_ID==2) && ((-93<hit_phi && hit_phi<-82)||(80<hit_phi && hit_phi<92))) continue;          
//       if ((detector_ID==1 && subdetector_ID==3) && ((-92<hit_phi && hit_phi<-75)||(77<hit_phi && hit_phi<92))) continue;          
//       if ((detector_ID==1 && subdetector_ID==4) && ((-92<hit_phi && hit_phi<-62)||(66<hit_phi && hit_phi<92))) continue;          
//       if ((detector_ID==1 && subdetector_ID==5) && ((-92<hit_phi && hit_phi<-65)||(66<hit_phi && hit_phi<92))) continue;          
//       if ((detector_ID==1 && subdetector_ID==6) && ((-92<hit_phi && hit_phi<-58)||(66<hit_phi && hit_phi<95))) continue;    

      //oxford ptarget field
      //cut gem      
      if ((detector_ID==1 && subdetector_ID==1) && ((-95<hit_phi && hit_phi<-75)||(75<hit_phi && hit_phi<95))) continue;          
      if ((detector_ID==1 && subdetector_ID==2) && ((-95<hit_phi && hit_phi<-75)||(75<hit_phi && hit_phi<95))) continue;          
      if ((detector_ID==1 && subdetector_ID==3) && ((-95<hit_phi && hit_phi<-75)||(75<hit_phi && hit_phi<95))) continue;          
      if ((detector_ID==1 && subdetector_ID==4) && ((-95<hit_phi && hit_phi<-75)||(75<hit_phi && hit_phi<95))) continue;          
      if ((detector_ID==1 && subdetector_ID==5) && ((-85<hit_phi && hit_phi<-55)||(70<hit_phi && hit_phi<90))) continue;          
      if ((detector_ID==1 && subdetector_ID==6) && ((-85<hit_phi && hit_phi<-55)||(70<hit_phi && hit_phi<90))) continue;  
      //cut LGC
      if ((detector_ID==2 && subdetector_ID==1) && ((-85<hit_phi && hit_phi<-55)||(70<hit_phi && hit_phi<90))) continue;
      //cut HGC
      if ((detector_ID==2 && subdetector_ID==1) && ((-85<hit_phi && hit_phi<-40 && hit_r/1e1<180)||(65<hit_phi && hit_phi<90 && hit_r/1e1<180))) continue;                
      //cut FASPD
      if ((detector_ID==5 && subdetector_ID==1) && ((-75<hit_phi && hit_phi<-40)||(70<hit_phi && hit_phi<95))) continue; 
      //cut LASPD
      if ((detector_ID==5 && subdetector_ID==2) && ((-90<hit_phi && hit_phi<-60)||(60<hit_phi && hit_phi<90))) continue;       
      //cut FAMRPC
      if ((detector_ID==4 && subdetector_ID==1) && ((-75<hit_phi && hit_phi<-40)||(70<hit_phi && hit_phi<95))) continue; 
      //cut EC    
      if ((detector_ID==3 && subdetector_ID==1) && ((-75<hit_phi && hit_phi<-40 && hit_r/1e1<195)||(70<hit_phi && hit_phi<90 && hit_r/1e1<195))) continue;                
      if ((detector_ID==3 && subdetector_ID==2) && ((-90<hit_phi && hit_phi<-60)||(60<hit_phi && hit_phi<95))) continue;
      
// double cut1[12]={-95,-95,-95,-95,-85,-85,-85,-85,-75,-90,-75,-90};
// double cut2[12]={-75,-75,-75,-75,-55,-55,-55,-40,-40,-60,-40,-60};
// double cut3[12]={ 75, 75, 75, 75, 70, 70, 70, 65, 70, 60, 70, 60};
// double cut4[12]={ 95, 95, 95, 95, 90, 90, 90, 90, 90, 95, 90, 95};
      
      acc[hit_id]=1; 
    }
    else if (Is_JPsi){
      acc[hit_id]=1;      
    }
    else if (Is_DDVCS_JPsi_LH2){
      acc[hit_id]=1;      
    }
    else if (Is_DDVCS_PVDIS_LH2){
      acc[hit_id]=1;

    }    
       
       ///DIRC cut
//        double x=flux_avg_x->at(j)/1e1,y=flux_avg_y->at(j)/1e1;
//         if (detector_ID==3 && subdetector_ID==1) {
// // 	  if ((fabs(x)<98 && fabs(y) <98) || fabs(x)>181 || fabs(y)>181) continue;  // 4 side
// //           if ( (y>-98 && y<196+sqrt(3)*x && y<196-sqrt(3)*x) || (y<-181 || y>362+sqrt(3)*x || y>362-sqrt(3)*x) ) continue;  // 3 side, y=196+-sqrt(3)*x,y=362+-sqrt(3)*x
// 	  if (fabs(y)<98 || fabs(y)>181) continue;  // 2 side, up and down
// 	}
       
//         if (vx_gen != flux_vx->at(j)/1e1) cout << "vx " << vx_gen << " " << flux_vx->at(j)/1e1 <<endl;
//         if (vy_gen != flux_vy->at(j)/1e1) cout << "vy " << vy_gen << " " << flux_vy->at(j)/1e1 <<endl;
//         if (vz_gen != flux_vz->at(j)/1e1) cout << "vz " << vz_gen << " " << flux_vz->at(j)/1e1 <<endl;		  
    }        
    //check all hits, end ============================================================

    //analysis hit result, begin ============================================================    
//     if(Is_out){
//     for (Int_t j=0;j<flux_hitn->size();j++) {
//       cout << j << " !!! " << flux_id->at(j) << " " << flux_pid->at(j) << " " << flux_mpid->at(j) << " " << flux_tid->at(j) << " " << flux_mtid->at(j) << " " << flux_trackE->at(j) << " " << flux_totEdep->at(j) << " (" << flux_avg_x->at(j) << " " << flux_avg_y->at(j) << " " << flux_avg_z->at(j) << ") " << flux_avg_lx->at(j) << " " << flux_avg_ly->at(j) << " " << flux_avg_lz->at(j) << " " << flux_px->at(j) << " " << flux_py->at(j) << " " << flux_pz->at(j) << " " << flux_vx->at(j) << " " << flux_vy->at(j) << " " << flux_vz->at(j) << " " << flux_mvx->at(j) << " " << flux_mvy->at(j) << " " << flux_mvz->at(j) << " " << flux_avg_t->at(j) << endl;           
//     }
//     for (int j=0;j<gen_pid->size();j++) {
//       cout << gen_pid->at(j) << " " << gen_px->at(j) << " " << gen_py->at(j) << " " << gen_pz->at(j) << " " << gen_vx->at(j) << " " << gen_vy->at(j) << " " << gen_vz->at(j) << endl; 
//     }    
//     }

    if (mom_GEM3.Mag()>0 && mom_GEM4.Mag()>0){ // hit on both GEM3 and GEM4
//       hlinearity_GEM34->Fill(theta_gen,p_gen/1e3,(coor_GEM4-coor_GEM3).Angle(mom_GEM3)*DEG);      
      hlinearity_GEM34->SetBinContent(int((theta_gen-theta_min)/((theta_max-theta_min)/1000.)),p_gen/1e3/(11/2200.),(coor_GEM4-coor_GEM3).Angle(mom_GEM3)*DEG);            
    }
    
    int pattern_id=-1;       
    if (Is_PVDIS){    
//       if (acc[6]==1) pattern_id=0;  //hit on FAEC
      if (acc[0]==1&&acc[1]==1&&acc[2]==1&&acc[3]==1&&acc[4]==1&&acc[10]==1) pattern_id=0;
      //hit on FAEC and GEM      
    }
    else if (Is_SIDIS_He3){
//       if (acc[6]==1) pattern_id=0;  //hit on FAEC
//       if (acc[7]==1) pattern_id=1;	  //hit on LAEC
      if (acc[1]==1&&acc[2]==1&&acc[3]==1&&acc[4]==1&&acc[5]==1&&acc[10]==1) pattern_id=0; //hit on FAEC and GEM             
      if (acc[0]==1&&acc[1]==1&&acc[2]==1&&acc[3]==1&&acc[11]==1) pattern_id=1; //hit on LAEC and GEM  

      if (acc[10]==1) {
	for (Int_t k=0;k<6;k++) {      
	    if (acc[k]==1) {
	      hmissingGEM_forwardangle->Fill(k+1);
	    }
	}      
      }
      if (acc[11]==1) {
	for (Int_t k=0;k<6;k++) {      
	    if (acc[k]==1) {
	      hmissingGEM_largeangle->Fill(k+1);
	    }
	}      
      }
      
    }
    else if (Is_SIDIS_NH3){
//       if (acc[6]==1) pattern_id=0;
//       if (acc[7]==1) pattern_id=1;	        
      if (abs(pid_gen)==11){      
      if (acc[1]==1&&acc[2]==1&&acc[3]==1&&acc[4]==1&&acc[5]==1&&acc[6]==1&&acc[8]==1&&acc[10]==1) pattern_id=0; //hit on GEM, LGC, FASPD, FAEC for electron 
      if (acc[0]==1&&acc[1]==1&&acc[2]==1&&acc[3]==1&&acc[9]==1&&acc[11]==1) pattern_id=1; //hit on GEM, LASPD, LAEC for electron
      }
      
      if (abs(pid_gen)==211){      
      if (acc[1]==1&&acc[2]==1&&acc[3]==1&&acc[4]==1&&acc[5]==1&&acc[7]==1&&acc[8]==1) pattern_id=0; //hit on GEM,HGC,FASPD for pion 
      if (acc[0]==1&&acc[1]==1&&acc[2]==1&&acc[3]==1&&acc[9]==1) pattern_id=1; //hit on GEM, LASPD for pion
      }
      
      if (abs(pid_gen)==321){      
      if (acc[1]==1&&acc[2]==1&&acc[3]==1&&acc[4]==1&&acc[5]==1&&acc[8]==1) pattern_id=0; //hit on GEM, FASPD for kaon 
      if (acc[0]==1&&acc[1]==1&&acc[2]==1&&acc[3]==1&&acc[9]==1) pattern_id=1; //hit on GEM, LASPD for kaon
      }
      
      if (abs(pid_gen)==22){      
      if (acc[10]==1) pattern_id=0; //hit on FAEC for photon
      if (acc[11]==1) pattern_id=1; //hit on LAEC for photon      
      }      
      
    }
    else if (Is_JPsi){
//       if (acc[6]==1) pattern_id=0;  //hit on FAEC
//       if (acc[7]==1) pattern_id=1;	  //hit on LAEC      
      
      if (acc[1]==1&&acc[2]==1&&acc[3]==1&&acc[4]==1&&acc[5]==1&&acc[10]==1) pattern_id=0; //hit on FAEC and GEM
      if (acc[0]==1&&acc[1]==1&&acc[2]==1&&acc[3]==1&&acc[11]==1) pattern_id=1; //hit on LAEC and GEM
    }
    else if (Is_DDVCS_JPsi_LH2){
//       if (acc[8]==1) pattern_id=0;  //hit on forward angle muon det behind endcap
//       if (acc[9]==1) pattern_id=1;	  //hit on large angle muon det outside of endcap donut
      if (acc[1]==1&&acc[2]==1&&acc[3]==1&&acc[4]==1&&acc[5]==1&&acc[13]==1) pattern_id=0;  //hit on forward angle muon det behind endcap and all gem
      if (acc[0]==1&&acc[1]==1&&acc[2]==1&&acc[3]==1&&acc[14]==1) pattern_id=1;	  //hit on large angle muon det outside of endcap donut and all gem      
    }    
    else if (Is_DDVCS_PVDIS_LH2){
//       if (acc[8]==1) pattern_id=0;  //hit on FAEC      
      if (acc[0]==1&&acc[1]==1&&acc[2]==1&&acc[3]==1&&acc[4]==1&&acc[13]==1) pattern_id=0; //hit on forward angle muon det behind endcap and all gem
    }        
    

    int counter_hit[2]={0,0};
    if(pattern_id != -1){   
      counter_hit[pattern_id]++;      
      
      hflux_P[pattern_id]->Fill(p_gen/1e3);
      hflux_Theta[pattern_id]->Fill(theta_gen);            
      hflux_Phi[pattern_id]->Fill(phi_gen);	      
      hflux_ThetaP[pattern_id]->Fill(theta_gen,p_gen/1e3);
      hflux_ThetaPhi[pattern_id]->Fill(theta_gen,phi_gen);            
      hflux_PhiP[pattern_id]->Fill(phi_gen,p_gen/1e3);	
      hflux_ThetaPhiP[pattern_id]->Fill(theta_gen,phi_gen,p_gen/1e3);
      
      hflux_ThetaVz[pattern_id]->Fill(vz_gen/1e1,theta_gen,rate);
      hflux_ThetaVr[pattern_id]->Fill(sqrt(vx_gen*vx_gen+vy_gen*vy_gen)/1e1,theta_gen,rate);
    }
      
    if (counter_hit[0]>1 || counter_hit[1]>1) {
      cout << endl;
      cout << "more than 1 hit? " << counter_hit[0] << " " << counter_hit[1] << endl;
    }
    
    if(pattern_id != -1){
//       hhit_PR[pattern_id][0]->Fill(R_hit_EC,p_gen/1e3);
//       if(Q2>1.7) hhit_PR[pattern_id][1]->Fill(R_hit_EC,p_gen/1e3);
      hhit_PR[pattern_id][0]->Fill(R_hit_EC,p_gen/1e3,rate);
      if(Q2>1.7) hhit_PR[pattern_id][1]->Fill(R_hit_EC,p_gen/1e3,rate);      
    }
    
    
//     for (Int_t j=0;j<flux_hitn->size();j++) {
//     int detector_ID=flux_id->at(j)/1000000;
//     int subdetector_ID=(flux_id->at(j)%1000000)/100000;
//     int subsubdetector_ID=((flux_id->at(j)%1000000)%100000)/1000;
// 
//       double hit_r=sqrt(pow(flux_avg_x->at(j),2)+pow(flux_avg_y->at(j),2));
//       double hit_y=flux_avg_y->at(j),hit_x=flux_avg_x->at(j);  
//       double hit_phi=atan2(flux_avg_y->at(j),flux_avg_x->at(j))*DEG;
//       double hit_theta=atan((hit_r/1e1-sqrt(vx_gen*vx_gen+vy_gen*vy_gen))/(320-vz_gen))*DEG;    
//            
//       if ((pattern_id==0 && (flux_id->at(j)==1410000 || flux_id->at(j)==1510000 || flux_id->at(j)==1610000 || flux_id->at(j)==3110000)) || (pattern_id==1 && (flux_id->at(j)==1110000 || flux_id->at(j)==1210000 || flux_id->at(j)==1310000 || flux_id->at(j)==1410000 || flux_id->at(j)==3210000))) {	
// 	//check hit on GEM
// 	if (detector_ID==1) {
//     //       if (flux_pid->at(j)==11) {
// 	  hhit_xy_gem[subdetector_ID-1][1]->Fill(hit_x/1e1,hit_y/1e1);    
// 	  hhit_PhiR_gem[subdetector_ID-1][1]->Fill(hit_phi,hit_r/1e1);
//     //       }
// 	}       
// 
// 	if (detector_ID==3) {
// 	hhit_xy_ec[subdetector_ID-1][1]->Fill(hit_x/1e1,hit_y/1e1);    
// 	hhit_PhiR_ec[subdetector_ID-1][1]->Fill(hit_phi,hit_r/1e1);
// 	}
//       }
//     }    
    
//     for(int k=0;k<n;k++){
//       if (counter[pattern_id] > 1) cout << counter[pattern_id] << " " << pattern_id << endl;
//       counter[pattern_id]=0;
//     }
    //analysis hit result, end ============================================================       
}
file->Close();

cout << "counter_decay " << counter_decay << endl;
cout << "counter_loss " << counter_loss << endl;

if (Is_PVDIS){ 
//   double planeZ[6]={40,68,96,124,152,180};
//   double Rin[6]={2.11, 12.86, 23.61, 34.36, 45.10, 55.85};
//   double Rout[6]={39.60, 59.94, 80.28, 100.63, 120.97, 141.31};
double Rin[11]={5.00,  7.88, 13.53, 19.19, 24.85, 30.50, 36.16, 41.82, 47.47, 53.13, 58.79 };
double Rout[11]={ 38.16, 47.96, 57.77, 67.57, 77.37, 87.18, 96.98, 106.78, 116.58, 126.39, 136.19 };
TCanvas *c_baffleplate = new TCanvas("baffleplate","baffleplate",1350,900);
c_baffleplate->Divide(4,3);
// for (int i=0;i<11;i++){
//   c_baffleplate->cd(i+1);      
//   hbaffleplate_observer[i]->SetMarkerColor(kRed);  
//   hbaffleplate_observer[i]->Draw("box");    
// //   hbaffleplate[i]->Draw("same");  
// //   hbaffleplate_observer[i]->Add(hbaffleplate[i],-1);
//   TArc *c_in=new TArc(0,0,Rin[i]);
//   c_in->SetLineColor(kBlack);
//   c_in->SetFillStyle(0);  
//   c_in->Draw();
//   TArc *c_out=new TArc(0,0,Rout[i]);
//   c_out->SetLineColor(kBlack);
//   c_out->SetFillStyle(0);    
//   c_out->Draw();
// }
for (int i=11;i<22;i++){
  c_baffleplate->cd(i-11+1);      
  hbaffleplate_observer[i]->SetMarkerColor(kBlue);  
//   hbaffleplate_observer[i]->Draw("box same");    
    hbaffleplate_observer[i]->Draw();    
}
c_baffleplate->SaveAs(Form("%s_%s",the_filename,"baffleplate.png"));
}

for(int i=0;i<m;i++) {
  hacceptance_P[i]->Divide(hflux_P[i],hgen_P);  
  hacceptance_P[i]->SetMinimum(0);  
  hacceptance_P[i]->SetMaximum(1);   
  hacceptance_Theta[i]->Divide(hflux_Theta[i],hgen_Theta);    
  hacceptance_Theta[i]->SetMinimum(0);  
  hacceptance_Theta[i]->SetMaximum(1);    
  hacceptance_Phi[i]->Divide(hflux_Phi[i],hgen_Phi);    
  hacceptance_Phi[i]->SetMinimum(0);  
  hacceptance_Phi[i]->SetMaximum(1);    
  hacceptance_ThetaP[i]->Divide(hflux_ThetaP[i],hgen_ThetaP);  
  hacceptance_ThetaP[i]->SetMinimum(0);  
  hacceptance_ThetaP[i]->SetMaximum(1);  
  hacceptance_ThetaPhi[i]->Divide(hflux_ThetaPhi[i],hgen_ThetaPhi);  
  hacceptance_ThetaPhi[i]->SetMinimum(0);  
  hacceptance_ThetaPhi[i]->SetMaximum(1);  
  hacceptance_PhiP[i]->Divide(hflux_PhiP[i],hgen_PhiP);  
  hacceptance_PhiP[i]->SetMinimum(0);  
  hacceptance_PhiP[i]->SetMaximum(1);    
  hacceptance_ThetaPhiP[i]->Divide(hflux_ThetaPhiP[i],hgen_ThetaPhiP);  
  hacceptance_ThetaPhiP[i]->SetMinimum(0);  
  hacceptance_ThetaPhiP[i]->SetMaximum(1); 
  hacceptance_ThetaVz[i]->Divide(hflux_ThetaVz[i],hgen_ThetaVz); 
  hacceptance_ThetaVz[i]->SetMinimum(0);  
  hacceptance_ThetaVz[i]->SetMaximum(1);      
  hacceptance_ThetaVr[i]->Divide(hflux_ThetaVr[i],hgen_ThetaVr);  
  hacceptance_ThetaVr[i]->SetMinimum(0);  
  hacceptance_ThetaVr[i]->SetMaximum(1);      
}

TCanvas *c_acc = new TCanvas("acc","acc",1800,1000);
c_acc->Divide(4,2);
c_acc->cd(1);
gPad->SetLogy(1);
hacceptance_ThetaP[0]->Draw("colz");
hacceptance_ThetaP[0]->SetAxisRange(0.1,11,"Y");
c_acc->cd(2);
hacceptance_ThetaPhi[0]->Draw("colz");
c_acc->cd(3);
gPad->SetLogy(1);
hacceptance_PhiP[0]->Draw("colz");
hacceptance_PhiP[0]->SetAxisRange(0.1,11,"Y");
c_acc->cd(4);
hacceptance_ThetaPhiP[0]->Draw("box");
c_acc->cd(5);
gPad->SetLogy(1);
hacceptance_ThetaP[1]->Draw("colz");
hacceptance_ThetaP[1]->SetAxisRange(0.1,11,"Y");
c_acc->cd(6);
hacceptance_ThetaPhi[1]->Draw("colz");
c_acc->cd(7);
gPad->SetLogy(1);
hacceptance_PhiP[1]->Draw("colz");
hacceptance_PhiP[1]->SetAxisRange(0.1,11,"Y");
c_acc->cd(8);
hacceptance_ThetaPhiP[1]->Draw("box");
// c_acc->cd(7);
// hacceptance_P[0]->Draw();
// c_acc->cd(8);
// hacceptance_Theta[0]->Draw();
// c_acc->cd(9);
// hacceptance_Phi[0]->Draw();
// c_acc->cd(10);
// hacceptance_P[1]->Draw();
// c_acc->cd(11);
// hacceptance_Theta[1]->Draw();
// c_acc->cd(12);
// hacceptance_Phi[1]->Draw();
c_acc->SaveAs(Form("%s_%s",the_filename,"acc.png"));

TCanvas *c_acc_vertex = new TCanvas("acc_vertex","acc_vertex",1600,800);
c_acc_vertex->Divide(2,1);
c_acc_vertex->cd(1);
hacceptance_ThetaVz[0]->Draw("colz");
hacceptance_ThetaVz[0]->SetTitle(Form("FA sum %i",int(hacceptance_ThetaVz[0]->Integral())));
c_acc_vertex->cd(2);
hacceptance_ThetaVz[1]->Draw("colz");
hacceptance_ThetaVz[1]->SetTitle(Form("LA sum %i",int(hacceptance_ThetaVz[1]->Integral())));
// c_acc_vertex->cd(3);
// hacceptance_ThetaVr[0]->Draw("colz");
// c_acc_vertex->cd(4);
// hacceptance_ThetaVr[1]->Draw("colz");
c_acc_vertex->SaveAs(Form("%s_%s",the_filename,"acc_vertex.png"));

TCanvas *c_gen = new TCanvas("gen","gen",1800,600);
c_gen->Divide(3,3);
c_gen->cd(1);
hgen_P->Draw();
c_gen->cd(2);
hgen_Theta->Draw();
c_gen->cd(3);
hgen_Phi->Draw();
c_gen->cd(4);
hgen_ThetaP->Draw("colz");
c_gen->cd(5);
hgen_ThetaPhi->Draw("colz");
c_gen->cd(6);
hgen_PhiP->Draw("colz");
c_gen->cd(7);
hgen_ThetaPhiP->Draw("box");

TCanvas *c_gen_vertex = new TCanvas("gen_vertex","gen_vertex",1200,600);
c_gen_vertex->Divide(2,1);
c_gen_vertex->cd(1);
hgen_ThetaVz->Draw("colz");
c_gen_vertex->cd(2);
hgen_ThetaVr->Draw("colz");

// TCanvas *c_flux = new TCanvas("flux","flux",1200,600);
// c_flux->Divide(3,2);
// c_flux->cd(1);
// hflux_ThetaP[0]->Draw("colz");
// c_flux->cd(2);
// hflux_ThetaVz[0]->Draw("colz");
// c_flux->cd(3);
// hflux_ThetaVr[0]->Draw("colz");
// c_flux->cd(4);
// hflux_ThetaP[1]->Draw("colz");
// c_flux->cd(5);
// hflux_ThetaVz[1]->Draw("colz");
// c_flux->cd(6);
// hflux_ThetaVr[1]->Draw("colz");
// c_flux->SaveAs(Form("%s_%s",the_filename,"flux.png"));


TCanvas *c_hit_PR = new TCanvas("hit_PR","hit_PR",1800,800);
c_hit_PR->Divide(2,1);
c_hit_PR->cd(1);
hhit_PR[0][0]->Draw("colz");
c_hit_PR->cd(2);
hhit_PR[0][1]->Draw("colz");
// c_hit_PR->cd(3);
// hhit_PR[1][0]->Draw("colz");
// c_hit_PR->cd(4);
// hhit_PR[1][1]->Draw("colz");

TCanvas *c_hit_rz = new TCanvas("hit_rz","hit_rz",1800,800);
hhit_rz->Draw("colz");

TCanvas *c_hit_rMom = new TCanvas("hit_rMom","hit_rMom",1800,800);
c_hit_rMom->Divide(n,1);
for(int k=0;k<n;k++){
c_hit_rMom->cd(k+1);
gPad->SetLogz(1);
gPad->SetLogy(1);
hhit_rMom[k]->Draw("colz");
}
// c_hit_rMom->SaveAs(Form("%s_%s",the_filename,"hit_rMom.png"));

TCanvas *c_hit_diffMom = new TCanvas("hit_diffMom","hit_diffMom",1800,800);
c_hit_diffMom->Divide(n,3);
for(int k=0;k<n;k++){
c_hit_diffMom->cd(k+1);
gPad->SetLogz(1);
hhit_phidiffMom[k]->Draw("colz");
}
for(int k=0;k<n;k++){
c_hit_diffMom->cd(n+k+1);
gPad->SetLogz(1);
hhit_thetadiffMom[k]->Draw("colz");
}
for(int k=0;k<n;k++){
c_hit_diffMom->cd(n*2+k+1);
gPad->SetLogz(1);
hhit_anglediffMom[k]->Draw("colz");
}
// c_hit_thetadiffMom->SaveAs(Form("%s_%s",the_filename,"hit_thetadiffMom.png"));

TCanvas *c_mom_diffMom = new TCanvas("mom_diffMom","mom_diffMom",1800,800);
c_mom_diffMom->Divide(n,3);
for(int k=0;k<n;k++){
c_mom_diffMom->cd(k+1);
gPad->SetLogz(1);
hmom_phidiffMom[k]->Draw("colz");
}
for(int k=0;k<n;k++){
c_mom_diffMom->cd(n+k+1);
gPad->SetLogz(1);
hmom_thetadiffMom[k]->Draw("colz");
}
for(int k=0;k<n;k++){
c_mom_diffMom->cd(n*2+k+1);
gPad->SetLogz(1);
hmom_anglediffMom[k]->Draw("colz");
}

TCanvas *c_hit = new TCanvas("hit","hit",1800,800);
c_hit->Divide(n,4);
for(int j=0;j<1;j++){
for(int i=0;i<n;i++){
c_hit->cd(j*6+0*n+i+1);
gPad->SetLogz(1);
hhit_xy[i][j]->Draw("colz");
c_hit->cd(j*6+1*n+i+1);
gPad->SetLogz(1);
hhit_PhiR[i][j]->Draw("colz");
c_hit->cd(j*6+2*n+i+1);
gPad->SetLogz(1);
hhit_PhiZ[i][j]->Draw("colz");
c_hit->cd(j*6+3*n+i+1);
gPad->SetLogz(1);
hPloss[i][j]->Draw("colz");
}}

// TCanvas *c_acceptance_all_gem = new TCanvas("acceptance_gem","acceptance_gem",1800,800);
// c_acceptance_all_gem->Divide(2,3);
// for(int k=0;k<6;k++){
// c_acceptance_all_gem->cd(k+1);
// hacceptance_ThetaP[k]->Draw("colz");
// char hsttitle[80];
// sprintf(hsttitle,"acceptance at GEM plane %i;theta (deg);vertex P (GeV)",k+1);
// hacceptance_ThetaP[k]->SetTitle(hsttitle);
// }
// c_acceptance_all_gem->SaveAs("acceptance_gem.png");

// TCanvas *c_acceptance_all_gem_1D = new TCanvas("acceptance_gem_1D","acceptance_gem_1D",1800,800);
// c_acceptance_all_gem_1D->Divide(2,6);
// for(int k=0;k<6;k++){
// c_acceptance_all_gem_1D->cd(k+1);
// hacceptance_P[k]->Draw();
// c_acceptance_all_gem_1D->cd(6+k+1);
// hacceptance_theta[k]->Draw();
// }

hacceptance_ThetaP_forwardangle=(TH2F*) hacceptance_ThetaP[0]->Clone();
hacceptance_ThetaP_forwardangle->SetName("acceptance_ThetaP_forwardangle");
hacceptance_ThetaP_largeangle=(TH2F*) hacceptance_ThetaP[1]->Clone();
hacceptance_ThetaP_largeangle->SetName("acceptance_ThetaP_largeangle");

hacceptance_ThetaP_overall=(TH2F*) hacceptance_ThetaP_forwardangle->Clone();
hacceptance_ThetaP_overall->Add(hacceptance_ThetaP_largeangle);
hacceptance_ThetaP_overall->SetNameTitle("acceptance_ThetaP_overall","acceptance;vertex #theta (deg);vertex P (GeV)");
hacceptance_ThetaP_overall->SetMinimum(0);  
hacceptance_ThetaP_overall->SetMaximum(1);  

hacceptance_ThetaPhi_forwardangle=(TH2F*) hacceptance_ThetaPhi[0]->Clone();
hacceptance_ThetaPhi_forwardangle->SetName("acceptance_ThetaPhi_forwardangle");
hacceptance_ThetaPhi_largeangle=(TH2F*) hacceptance_ThetaPhi[1]->Clone();
hacceptance_ThetaPhi_largeangle->SetName("acceptance_ThetaPhi_largeangle");

hacceptance_PhiP_forwardangle=(TH2F*) hacceptance_PhiP[0]->Clone();
hacceptance_PhiP_forwardangle->SetName("acceptance_PhiP_forwardangle");
hacceptance_PhiP_largeangle=(TH2F*) hacceptance_PhiP[1]->Clone();
hacceptance_PhiP_largeangle->SetName("acceptance_PhiP_largeangle");

hacceptance_ThetaPhiP_forwardangle=(TH3F*) hacceptance_ThetaPhiP[0]->Clone();
hacceptance_ThetaPhiP_forwardangle->SetNameTitle("acceptance_ThetaPhiP_forwardangle","acceptance by FA;vertex #theta (deg);vertex #phi (deg);vertex P (GeV)");
hacceptance_ThetaPhiP_largeangle=(TH3F*) hacceptance_ThetaPhiP[1]->Clone();
hacceptance_ThetaPhiP_largeangle->SetNameTitle("acceptance_ThetaPhiP_largeangle","acceptance by LA;vertex #theta (deg);vertex #phi (deg);vertex P (GeV)");

// gStyle->SetOptStat(0);

if (Is_PVDIS || Is_DDVCS_PVDIS_LH2){ 
TCanvas *c_acceptance_2D = new TCanvas("acceptance_2D","acceptance_2D",800,600);
// c_acceptance_all->Divide(1,3);
// c_acceptance_all->cd(1);
// gPad->SetLogy(1);  
hacceptance_ThetaP_forwardangle->Draw("colz");
c_acceptance_2D->SaveAs(Form("%s_%s",the_filename,"acceptance_2D.png"));
}
else{
TCanvas *c_acceptance_2D = new TCanvas("acceptance_2D","acceptance_2D",500,900);
c_acceptance_2D->Divide(1,3);
c_acceptance_2D->cd(1);
gPad->SetLogy();
gPad->SetGrid();
hacceptance_ThetaP_forwardangle->Draw("colz");
c_acceptance_2D->cd(2);
gPad->SetLogy();
gPad->SetGrid();
hacceptance_ThetaP_largeangle->Draw("colz");
c_acceptance_2D->cd(3);
gPad->SetLogy();
gPad->SetGrid();
hacceptance_ThetaP_overall->Draw("colz");
c_acceptance_2D->SaveAs(Form("%s_%s",the_filename,"acceptance_2D.png"));

TCanvas *c_acceptance_2D_sep = new TCanvas("acceptance_2D_sep","acceptance_2D_sep",1200,500);
c_acceptance_2D_sep->Divide(2,1);
c_acceptance_2D_sep->cd(1);
gPad->SetLogy();
gPad->SetGrid();
hacceptance_ThetaP_forwardangle->Draw("colz");
c_acceptance_2D_sep->cd(2);
gPad->SetLogy();
gPad->SetGrid();
hacceptance_ThetaP_largeangle->Draw("colz");
c_acceptance_2D_sep->SaveAs(Form("%s_%s",the_filename,"acceptance_2D_sep.png"));

TCanvas *c_acceptance_2D_both = new TCanvas("acceptance_2D_both","acceptance_2D_both",1200,1000);
c_acceptance_2D_both->Divide(1,1);
c_acceptance_2D_both->cd(1);
gPad->SetLogy();
gPad->SetGrid();
hacceptance_ThetaP_overall->Draw("colz");
c_acceptance_2D_both->SaveAs(Form("%s_%s",the_filename,"acceptance_2D_both.png"));

TCanvas *c_acceptance_3D = new TCanvas("acceptance_3D","acceptance_3D",1200,800);
c_acceptance_3D->Divide(2,1);
c_acceptance_3D->cd(1);
hacceptance_ThetaPhiP_forwardangle->Draw("box");
c_acceptance_3D->cd(2);
hacceptance_ThetaPhiP_largeangle->Draw("box");
c_acceptance_3D->SaveAs(Form("%s_%s",the_filename,"acceptance_3D.png"));
}

TCanvas *c_linearity_GEM34 = new TCanvas("linearity_GEM34","linearity_GEM34",800,600);
hlinearity_GEM34->Draw("colz");

TCanvas *c_missingGEM = new TCanvas("missingGEM","missingGEM",1700,800);
c_missingGEM->Divide(2,1);
c_missingGEM->cd(1);
hmissingGEM_forwardangle->Draw();
c_missingGEM->cd(2);
hmissingGEM_largeangle->Draw();


cout << "checking holes in acceptance plots" << endl;
int NbinsX,NbinsY,NbinsZ;
cout << "checking hacceptance_ThetaP_forwardangle" << ",";
NbinsX=hacceptance_ThetaP_forwardangle->GetXaxis()->GetNbins();
NbinsY=hacceptance_ThetaP_forwardangle->GetYaxis()->GetNbins();
for(int j=1;j<NbinsY;j++){
  int count=0;
  bool Is_acc=false;
  for(int i=1;i<NbinsX;i++){
    if (hacceptance_ThetaP_forwardangle->GetBinContent(i,j)>0 && !Is_acc) {Is_acc=true;count++;}
    if (hacceptance_ThetaP_forwardangle->GetBinContent(i,j)==0 && Is_acc) {Is_acc=false;count++;} 
  }
  if (count !=2 && count !=0) cout << j << " " << count << ",";
}
cout << endl;
cout << "checking hacceptance_ThetaP_largeangle" << ",";
NbinsX=hacceptance_ThetaP_largeangle->GetXaxis()->GetNbins();
NbinsY=hacceptance_ThetaP_largeangle->GetYaxis()->GetNbins();
for(int j=1;j<NbinsY;j++){
  int count=0;
  bool Is_acc=false;
  for(int i=1;i<NbinsX;i++){
    if (hacceptance_ThetaP_largeangle->GetBinContent(i,j)>0 && !Is_acc) {Is_acc=true;count++;}
    if (hacceptance_ThetaP_largeangle->GetBinContent(i,j)==0 && Is_acc) {Is_acc=false;count++;} 
  }
  if (count !=2) cout << j << " " << count << ",";
}
cout << endl;
// cout << "checking hacceptance_ThetaPhiP_forwardangle" << ",";
// NbinsX=hacceptance_ThetaPhiP_forwardangle->GetXaxis()->GetNbins();
// NbinsY=hacceptance_ThetaPhiP_forwardangle->GetYaxis()->GetNbins();
// NbinsZ=hacceptance_ThetaPhiP_forwardangle->GetZaxis()->GetNbins();
// for(int k=1;k<NbinsZ;k++){
// for(int j=1;j<NbinsY;j++){
//   int count=0;
//   bool Is_acc=false;
//   for(int i=1;i<NbinsX;i++){
//     if (hacceptance_ThetaPhiP_forwardangle->GetBinContent(i,j,k)>0 && !Is_acc) {Is_acc=true;count++;}
//     if (hacceptance_ThetaPhiP_forwardangle->GetBinContent(i,j,k)==0 && Is_acc) {Is_acc=false;count++;} 
//   }
//   if (count !=2 && count !=0) cout << j << " " << k << " " << count << ",";
// }}
// cout << endl;
// cout << "checking hacceptance_ThetaPhiP_largeangle" << ",";
// NbinsX=hacceptance_ThetaPhiP_largeangle->GetXaxis()->GetNbins();
// NbinsY=hacceptance_ThetaPhiP_largeangle->GetYaxis()->GetNbins();
// NbinsZ=hacceptance_ThetaPhiP_largeangle->GetZaxis()->GetNbins();
// for(int k=1;k<NbinsZ;k++){
// for(int j=1;j<NbinsY;j++){
//   int count=0;
//   bool Is_acc=false;
//   for(int i=1;i<NbinsX;i++){
//     if (hacceptance_ThetaPhiP_largeangle->GetBinContent(i,j,k)>0 && !Is_acc) {Is_acc=true;count++;}
//     if (hacceptance_ThetaPhiP_largeangle->GetBinContent(i,j,k)==0 && Is_acc) {Is_acc=false;count++;} 
//   }
//   if (count !=2 && count !=0) cout << j << " " << k << " " << count << ",";
// }}
// cout << endl;

outputfile->Write();
outputfile->Flush();

char output_filename_final[200];
sprintf(output_filename_final, "%s_output_final.root",the_filename);
TFile *outputfile_final=new TFile(output_filename_final, "recreate");

hacceptance_ThetaP_forwardangle->SetDirectory(outputfile_final);
hacceptance_ThetaP_largeangle->SetDirectory(outputfile_final);
hacceptance_ThetaP_overall->SetDirectory(outputfile_final);
hacceptance_ThetaPhiP_forwardangle->SetDirectory(outputfile_final);
hacceptance_ThetaPhiP_largeangle->SetDirectory(outputfile_final);

outputfile_final->Write();
outputfile_final->Flush();
}
