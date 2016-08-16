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

const Int_t CHAR_LEN = 1000;

/*Init(){{{*/
void Init (const TString kInputFile){

    cout<<"============================================================"<<endl;
    cout<<"&& Initializing Parameters from "<<kInputFile<<" ..."<<endl;
    int i,j,k;
    vector<TString> inputdata;
    /*Read INPUTfile{{{*/
    FILE* INPUTfile;
    INPUTfile=fopen(kInputFile.Data(),"r");
    char buf[CHAR_LEN];
    char data[CHAR_LEN];
    while ( fgets(buf,CHAR_LEN,INPUTfile) )
    {
        i=0;
        while ( buf[i]==' '|| buf[i]=='\t' )
        {
            i++;
        }
        if ( buf[i]!='#' )
        {
            j=0;
            while ( buf[i]!='#' && buf[i]!='\0' && buf[i]!='\t' && buf[i]!='\n' )
            {
                if( buf[i]!=' ' && buf[i]!='\t' && buf[i]!='\n' )
                    data[j]=buf[i];
                i++; j++;
            }
            data[j]='\0';
            while ( data[--j]==' ' || data[j]=='\t'  || data[j]=='\n' )
            {
                //remove space or tab at the end of data
                data[j]='\0';
            }
            inputdata.push_back(data);
        }
        //else it's comment, skipped
    }
    
    fclose(INPUTfile);
    /*}}}*/

    /*Set Global Value{{{*/
    k=0;
    A=atoi(inputdata[k++]);
    Z=atoi(inputdata[k++]);
    particle_flag=atoi(inputdata[k++]);
    momentum_ele=atof(inputdata[k++]);
    momentum_ion=atof(inputdata[k++]);
    number_of_events=atoi(inputdata[k++]);
    FileNo=atoi(inputdata[k++]);
    config= inputdata[k++];
    model= inputdata[k++];
    cdxs_max=atof(inputdata[k++]);
    bLUND=atoi(inputdata[k++]);
    bXSMode=atoi(inputdata[k++]);
    Output_FileName= inputdata[k++];
    bDebug=atoi(inputdata[k++]);
  
    //A=atoi(inputdata[k++].c_str());
    //Z=atoi(inputdata[k++].c_str());
    //particle_flag=atoi(inputdata[k++].c_str());
    //momentum_ele=atof(inputdata[k++].c_str());
    //momentum_ion=atof(inputdata[k++].c_str());
    //FileNo=atoi(inputdata[k++].c_str());
    //number_of_events=atoi(inputdata[k++].c_str());
    //config= inputdata[k++].c_str();
    //model= inputdata[k++].c_str();
    //cdxs_max=atof(inputdata[k++].c_str());
    //bLUND=atoi(inputdata[k++].c_str());
    //bXSMode=atoi(inputdata[k++].c_str());
    //Output_FileName= inputdata[k++].c_str();
    //bDebug=atoi(inputdata[k++].c_str());
    /*}}}*/

    cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;
    cout<<"---- A = "<<A<<endl;
    cout<<"---- Z = "<<Z<<endl;
    cout<<"---- Partile(pi+/-: +/-1,K+/-: +/-2 ) = ";
    if(particle_flag==1) cout<<"Pi+"<<endl;
    else if(particle_flag==-1) cout<<"Pi-"<<endl;
    else if(particle_flag==2) cout<<"K+"<<endl;
    else if(particle_flag==-2) cout<<"K-"<<endl;
    else cout<<"*** ERROR, a wrong particle flag in your input-file!!!"<<endl;
    cout<<"---- P_e = "<<momentum_ele<<" GeV"<<endl;
    cout<<"---- P_A = "<<momentum_ion<<" GeV"<<endl;
    cout<<"---- #Events = "<<number_of_events<<endl;
    cout<<"---- File# = "<< FileNo;
    if(FileNo==0) cout<<" (from command line, e.g. ./GetSIDIS input.data FileNo)"<<endl;
    else cout<<endl;
    cout<<"---- Configs = "<<config.Data()<<endl;
    cout<<"---- Model = "<<model.Data()<<endl;
    cout<<"---- Save to LUND? = "<<bLUND<<endl;
    cout<<"---- Save in XSMode? = "<<bXSMode<<endl;
    cout<<"---- Rename files to (*.LUND, *_0.root)= "<<Output_FileName;
    if(Output_FileName=="NONE") cout <<" (use the default name)"<<endl;
    else cout<<endl;
    cout<<"---- Debug? = "<<bDebug<<endl;
    cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<endl;

    inputdata.clear();
    cerr<<"&& Initialization done!"<<endl;
    cout<<"============================================================"<<endl;
}
/*}}}*/
