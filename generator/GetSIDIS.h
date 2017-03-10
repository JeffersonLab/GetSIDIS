#define NOCASE 1

/*C/C++ Includes{{{*/
#include <stdio.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <iomanip>
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
/*}}}*/

using namespace std;

const double mm = 0.001;//default is cm
const double cm = 1.0;//default is cm
const double m = 1000.0;//default is cm
//char *LHAPDF_path = getenv("LHAPDF");

//mass
const Double_t mass_p = 0.93827;//GeV
const Double_t mass_n = 0.939566;//GeV
const Double_t mass_u = 0.931494;//GeV

////////////////////////
//SoLID Acceptance
////////////////////////
const Double_t SoLID_Mom_Min_e = 0.5;  
const Double_t SoLID_Mom_Max_e = 11.0;//not in use 
const Double_t SoLID_Mom_Min_h = 0.5;  
const Double_t SoLID_Mom_Max_h = 6.0;

const Double_t SoLID_Th_Min_e  = 7.0;  
const Double_t SoLID_Th_Max_e = 30.0; 
const Double_t SoLID_Th_Min_h  =  7.0; 
const Double_t SoLID_Th_Max_h =  30.0;

const Double_t SoLID_Ph_Min_e  = 0.0;  
const Double_t SoLID_Ph_Max_e = 360.0; 
const Double_t SoLID_Ph_Min_h  =  0.0; 
const Double_t SoLID_Ph_Max_h =  360.0;

//SoLID BeamSize Info
const Double_t SoLID_BeamSizeX_ele = 0.5 * cm;
const Double_t SoLID_BeamSizeY_ele = 0.5 * cm;
//const Double_t SoLID_BeamSizeX_ion = 0.5 * cm;
//const Double_t SoLID_BeamSizeY_ion = 0.5 * cm;
const Double_t SoLID_Target_Center = -350.0 * cm;
const Double_t SoLID_Target_Length = 40.0 *cm;
////////////////////////

////////////////////////
//EIC Acceptance
////////////////////////
//A rough guess but people claim EIC to be a full-acceptance device!
const Double_t EIC_Mom_Min_e = 0.5;  
const Double_t EIC_Mom_Max_e = 3.*10.0; //not in use 
const Double_t EIC_Mom_Min_h = 0.0;  
//const Double_t EIC_Mom_Max_h = 50.0;
const Double_t EIC_Mom_Max_h = 10.0;

const Double_t EIC_Th_Min_e = 0.0;  
//const Double_t EIC_Th_Max_e  = 180.0;            
const Double_t EIC_Th_Max_e  = 140.0;            
const Double_t EIC_Th_Min_h = 0.0;  
const Double_t EIC_Th_Max_h  = 180.0;

const Double_t EIC_Ph_Min_e = 0.0;  
const Double_t EIC_Ph_Max_e  = 360.0;            
const Double_t EIC_Ph_Min_h = 0.0;  
const Double_t EIC_Ph_Max_h  = 360.0;

//EIC BeamSize Info, to be updated
const Double_t EIC_BeamSizeX_ele = 0. * cm;
const Double_t EIC_BeamSizeY_ele = 0. * cm;
//const Double_t EIC_BeamSizeX_ion = 0. * cm;
//const Double_t EIC_BeamSizeY_ion = 0. * cm;
const Double_t EIC_Vertex_Center = 0.0 * cm;
const Double_t EIC_Vertex_Length = 0.0 *cm;


////////////////////////
//SPECT Acceptance
////////////////////////
const Double_t SPECT_Mom_Min_e = 2.0;  
const Double_t SPECT_Mom_Max_e = 6.0;
const Double_t SPECT_Mom_Min_h = 2.0;  
const Double_t SPECT_Mom_Max_h = 6.0;

const Double_t SPECT_Th_Min_e  = 5.0;  
const Double_t SPECT_Th_Max_e = 40.0; 
const Double_t SPECT_Th_Min_h  = 5.0; 
const Double_t SPECT_Th_Max_h =  40.0;

const Double_t SPECT_Ph_Min_e  =-5.0;  
const Double_t SPECT_Ph_Max_e = 5.0; 
const Double_t SPECT_Ph_Min_h  = -5.0; 
const Double_t SPECT_Ph_Max_h =  5.0;


//SPECT BeamSize Info
const Double_t SPECT_BeamSizeX_ele = 0.5 * cm;
const Double_t SPECT_BeamSizeY_ele = 0.5 * cm;
//const Double_t SPECT_BeamSizeX_ion = 0.5 * cm;
//const Double_t SPECT_BeamSizeY_ion = 0.5 * cm;
const Double_t SPECT_Target_Center = -0.0 * cm;
const Double_t SPECT_Target_Length = 10.0 *cm;
 

//////////////////////////////////
//For the input parameters
//////////////////////////////////
Int_t A=1;
Int_t Z=1;
Int_t particle_flag=1;//particle_flag = 1 for pion+, -1 for pion- and 2 for kaon+, -2 for kaon-
Double_t momentum_ele=10.0;//GeV
Double_t momentum_ion=100.0;//GeV, set ion_mom=0 for fix targets
Int_t FileNo = 0;//FileNo is the file number of output, used for batch
Int_t number_of_events = 1000;//events is number of event in each file
TString config = "SoLID";//config is 'EIC' or 'SoLID' which needs ion_mom=0

TString model = "EPS09";
//TString model = "LHAPDF";
Double_t cdxs_max =1000.;

Bool_t bLUND=kFALSE;
Bool_t bXSMode=kFALSE;
TString Output_FileName;

Bool_t bDebug = kTRUE;

void Init (const TString kInputFile);
double GetIonMass(const int A, const int Z);
