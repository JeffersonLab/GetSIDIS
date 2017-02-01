/////////////////////////////////////////////////////////////////////////////
//      Semi-Inclusive Deep Inelstic Scattering Cross Section Model        //
// Note:                                                                   //
//   This model is developed by Xin in his "collider" code. I extracted    //
//   the cross section parts and coded it into a C++ class which can be    //
//   easily embeded by other programs.                                     //
//                                                                         //
// ------------                                                            //
// In this version, I used only LO PDF and remove s, sbar and g  --11/30/16//
/////////////////////////////////////////////////////////////////////////////
/*C/C++ Includes{{{*/
#include <stdio.h>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <map>
#include <cmath>
//#include <experimental/string_view>
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
//#include "LHAPDF/LHAPDF.h"
//#include <TMatrix.h>
/*}}}*/
#include "cteqpdf.h"
#include "eps09.h"

using namespace std;
//using namespace LHAPDF;

//char *LHAPDF_Dir = std::getenv("LHAPDF");
const double DEG=180./3.1415926;
const double PI=3.1415926;
const double EXP = exp(1.);
const double GeV2_to_nbarn = 0.3894 * 1e6; //GeV^2 to nbarn
//ofstream outlog("pdf_check.dat");
class SIDIS
{
    public:
        SIDIS(TString kModel){/*{{{*/
            fModel = kModel;
            if(fModel=="EPS09"){
                //1->LO need CTEQ6L1,  2->NLO need CTEQ6.1M, 
                //fOrder = 2;
                fOrder = 1;
                SetEPS09();
            }
            else if(fModel=="CTEQPDF"){
                //3->free L0 CTEQ6L1 PDF, 4->free NL0 CTEQ6.1M PDF
                //fOrder = 4;
                fOrder = 3;
                SetCTEQ();
            }
            //else if(fModel=="LHAPDF") {
                //SetLHAPDF();
                //fOrder = 0;
            //}
            else{
                cerr<<"*** ERROR, I don't understand the XS model (not EPS09) :"<<fModel.Data()<<endl;
                exit(-2);
            }
        }/*}}}*/

        virtual ~SIDIS(){/*{{{*/
            if(fModel=="CTEQPDF"||fModel=="EPS09") 
                cteq_pdf_free(fPDF);
        }/*}}}*/

        /*void Print(){{{*/
        void Print(){
            cerr<<"===================================================================================================="<<endl;	
            cerr<<" --- Useage:"<<endl;
            cerr<<"      1, SetKin(momentum_ele, momentum_ion, P_ele, Th_ele, Ph_ele, P_had, Th_had, Ph_had, A, Z, ptc_flg)"<<endl;
            cerr<<"           where momentum_* and P_* are in GeV/c, and Theta_* and Phi_* are in Deg,"<<endl;
            cerr<<"                  pct_flag = 1->Pion+, 2->Pion-, 3->Kaon+, 4->Kaon-;"<<endl;
            cerr<<"      2, CalcXS() --> Calculated dxs_hp and dxs_hm; "<<endl; 
            cerr<<"      3, GetXS_HP() --> Get SIDIS XS of Pion+ or Kaon+ ,"<<endl;
            cerr<<"         GetXS_HM() --> Get SIDIS XS of Pion- or Kaon- ,"<<endl;
            cerr<<"         GetDilute_HP() --> Get Dilution of Pion+ or Kaon+ ,"<<endl;
            cerr<<"         GetDilute_HM() --> Get Dilution of Pion- or Kaon- ,"<<endl;
            cerr<<"      4, Other Output Quantities can be calculated after running CalcXS: "<<endl;
            cerr<<"           Q2, W, Wp, x, y, z, pt, nu, s, gamma, epsilon, jacoF "<<endl;
            cerr<<"===================================================================================================="<<endl;	

        }
        /*}}}*/

        /*SetLHAPDF{{{*/
 /*       void SetLHAPDF(){*/
            //const int SUBSET = 1;
            ////const string NAME = "CT10nlo";
            //const string NAME = "cteq6m";

            //TString LHAPDF_path=Form("%s/share/lhapdf",LHAPDF_Dir);
            //setPDFPath(LHAPDF_path.Data());
            //LHAPDF::initPDFSet(NAME, LHAPDF::LHPDF, SUBSET);
        /*}*/
        /*}}}*/

        void SetCTEQ(int mode){/*{{{*/
            //set ./cteq-pdf-1.0.4/Cteq6Pdf-2008.txt for details
            fPDF = cteq_pdf_alloc_id(mode);
        }/*}}}*/
        
        void SetCTEQ(){/*{{{*/
            int mode = 0;
            if(fOrder==1||fOrder==3)
                mode = 4; //CTEQ6L1, see JHEP04 (2009) 065
            else if(fOrder==2||fOrder==4)
                mode = 200; //CTEQ6.1M, see JHEP04 (2009) 065

             //set ./cteq-pdf-1.0.4/Cteq6Pdf-2008.txt for details
            fPDF = cteq_pdf_alloc_id(mode);
        }/*}}}*/

        void SetEPS09(){/*{{{*/
            fErrSet = 1;
            SetCTEQ();
        }/*}}}*/

        void SetCTEQOrder(int kOrder){/*{{{*/
            fOrder = kOrder;

            //Reinitialize:
            if(fModel=="EPS09"){
                if(fOrder!=1 && fOrder!=2)
                    fOrder = 2;//default is 2 if neither 1 or 2
                SetEPS09();
            }
            else if(fModel=="CTEQPDF"){
                if(fOrder!=3 && fOrder!=4)
                    fOrder = 4;//default is 4 if neither 3 or 4
                SetCTEQ();
            }
            else{
                cerr<<"*** ERROR, I don't understand the XS model (not EPS09) :"<<fModel.Data()<<endl;
                exit(-2);
            }
        }/*}}}*/

        /*void SetKin( double mom_beam_ele, double mom_beam_ion,...){{{*/
        void SetKin(double mom_beam_ele,double mom_beam_ion,                        // GeV          GeV   
                double p_ele, double th_ele, double ph_ele,         // GeV/c         DEG            DEG 
                double p_had, double th_had, double ph_had,         // GeV/c         DEG            DEG 
                double mass_target,int iA, int iZ,  int ptcl_flag){

            /*Define{{{*/
            fA = iA;
            fZ = iZ;
            double mom_ion;
            double energy_ion;
            double mom_ini_ele; 
            double energy_ini_ele;

            particle_flag = ptcl_flag;
            double mass_e = 0.511e-3; // electron mass in GeV
            //double mass_p = 0.93827;
            double mass_pi = 0.13957;
            //double mass_pi0 = 0.1349766;
            double mass_kaon = 0.493667;
            //double mass_n = 0.939566;

            ////initialize
            ///////////////////////////////////////////////////////////////////////////////////
            //--- I have question here: -- Z. Ye 07/27
            //  1) to get x, use the ion_mass, proton-mass or the average mass target_mass/A?
            //  2) then what is the momentum? Is it Mom_ion or Mom_ion/A?
            //  I will use the average values temprately 
            ///////////////////////////////////////////////////////////////////////////////////
            mass_target /= fA;
            mom_beam_ion /= fA;

            if (particle_flag == 1){
                mass_hadron = mass_pi;
            }else if (particle_flag == -1){
                mass_hadron = mass_pi;
            }else if (particle_flag == 2){
                mass_hadron = mass_kaon;
            }else if (particle_flag == -2){
                mass_hadron = mass_kaon;
            }else{
                cout << "particle_flag is wrong +-1 and +-2" << endl;
            }	

            E0_beam_ele = mom_beam_ele;
            E0_beam_ion = mom_beam_ion;

            // //define the 4-momentum
            //define electron direction as +z assuming a proton/neutron for the ion mass now 
            // approximation for SIDIS
            TLorentzVector *P4_ini_ele = new TLorentzVector(0.,0.,mom_beam_ele,sqrt(mom_beam_ele*mom_beam_ele + mass_e*mass_e));
            TLorentzVector *P4_ini_ion = new TLorentzVector(0.,0.,-mom_beam_ion,sqrt(mom_beam_ion*mom_beam_ion + mass_target*mass_target));
            TLorentzVector *P4_fin_had = new TLorentzVector(0.,0.,0.,1.);
            TLorentzVector *P4_fin_ele = new TLorentzVector(0.,0.,0.,1.);
            TLorentzVector *P4_q = new TLorentzVector(0.,0.,0.,1.);

            TLorentzVector *lrz_P4_q = new TLorentzVector(0.,0.,0.,1.);
            TLorentzVector *lrz_P4_h = new TLorentzVector(0.,0.,0.,1.);
            TLorentzVector *lrz_P4_ef = new TLorentzVector(0.,0.,0.,1.);
            //TLorentzVector *lrz_P4_ion = new TLorentzVector(0.,0.,0.,1.);
            //TLorentzVector *lrz_P4_ei = new TLorentzVector(0.,0.,0.,1.);

            double vn=mom_beam_ion/sqrt(mom_beam_ion*mom_beam_ion + mass_target*mass_target);
            TVector3 vnboost(0.,0.,vn);

            TVector3 p3_q,p3_fin_ele,p3_fin_had;
            TVector3 p3_target_spin(0.,1.,0);
            mom_ini_ele = fabs(mom_beam_ele);
            energy_ini_ele = sqrt(mom_beam_ele*mom_beam_ele + mass_e*mass_e);

            mom_ion = fabs(mom_beam_ion);
            energy_ion = sqrt(mom_beam_ion*mom_beam_ion + mass_target*mass_target);
            /*}}}*/

            /*Kinematics Quantities{{{*/
            //For electron
            px_ele = p_ele*sin(th_ele/180.*3.1415926)*cos(ph_ele/180.*3.1415926);
            py_ele = p_ele*sin(th_ele/180.*3.1415926)*sin(ph_ele/180.*3.1415926);
            pz_ele = p_ele*cos(th_ele/180.*3.1415926);
            E_ele = sqrt(p_ele*p_ele+mass_e*mass_e);
            P4_fin_ele->SetPxPyPzE(p_ele*sin(th_ele/180.*3.1415926)*cos(ph_ele/180.*3.1415926),
                    p_ele*sin(th_ele/180.*3.1415926)*sin(ph_ele/180.*3.1415926),
                    p_ele*cos(th_ele/180.*3.1415926)
                    ,sqrt(p_ele*p_ele+mass_e*mass_e));

            //For hadron
            px_had = p_had*sin(th_had/180.*3.1415926)*cos(ph_had/180.*3.1415926);
            py_had = p_had*sin(th_had/180.*3.1415926)*sin(ph_had/180.*3.1415926);
            pz_had = p_had*cos(th_had/180.*3.1415926);
            E_had = sqrt(p_had*p_had+mass_hadron*mass_hadron);
            P4_fin_had->SetPxPyPzE(p_had*sin(th_had/180.*3.1415926)*cos(ph_had/180.*3.1415926),
                    p_had*sin(th_had/180.*3.1415926)*sin(ph_had/180.*3.1415926),
                    p_had*cos(th_had/180.*3.1415926)
                    ,sqrt(p_had*p_had+mass_hadron*mass_hadron));


            *P4_q = *P4_ini_ele - *P4_fin_ele;
            Q2 = - (*P4_q)*(*P4_q);
            W = (*P4_ini_ele + *P4_ini_ion - *P4_fin_ele)*(*P4_ini_ele + *P4_ini_ion - *P4_fin_ele);
            Wp = (*P4_ini_ele + *P4_ini_ion - *P4_fin_ele - *P4_fin_had)*(*P4_ini_ele + *P4_ini_ion - *P4_fin_ele - *P4_fin_had);

            s = (*P4_ini_ele + *P4_ini_ion )*(*P4_ini_ele + *P4_ini_ion);
            nu = (*P4_ini_ion) * (*P4_q)/mass_target;
            x = Q2/(2 * (*P4_ini_ion)*(*P4_q));
            z = ((*P4_ini_ion)*(*P4_fin_had))/((*P4_ini_ion)*(*P4_q));
            y = ((*P4_ini_ion)*(*P4_q))/((*P4_ini_ion)*(*P4_ini_ele));
            W = sqrt(W);
            Wp = sqrt(Wp);

            *lrz_P4_q = *P4_q;
            *lrz_P4_h = *P4_fin_had;
            *lrz_P4_ef = *P4_fin_ele;

            lrz_P4_ef->Boost(vnboost);
            lrz_P4_h->Boost(vnboost); 
            lrz_P4_q->Boost(vnboost);

            gamma = 2*mass_target *x/sqrt(Q2);
            epsilon = (1-y-0.25*gamma*gamma*y*y)/(1-y+0.5*y*y+0.25*gamma*gamma*y*y);

            for(Int_t j=0;j<3;j++){
                p3_q(j)=(*lrz_P4_q)(j);
                p3_fin_ele(j)=(*lrz_P4_ef)(j);
                p3_fin_had(j)=(*lrz_P4_h)(j);
            }

            pt = p3_fin_had.Perp(p3_q);

            theta_q =p3_q.Theta();
            phi_q = p3_q.Phi();
            p3_target_spin.SetXYZ(0.,1.,0.);

            p3_fin_ele.RotateZ(-phi_q);
            p3_fin_ele.RotateY(-theta_q);
            p3_fin_had.RotateZ(-phi_q);
            p3_fin_had.RotateY(-theta_q);
            p3_target_spin.RotateZ(-phi_q);
            p3_target_spin.RotateY(-theta_q);

            phi_s = (Azimuthalphi(p3_target_spin(0),p3_target_spin(1))-Azimuthalphi(p3_fin_ele(0),p3_fin_ele(1)));
            phi_h = (Azimuthalphi(p3_fin_had(0),p3_fin_had(1))-Azimuthalphi(p3_fin_ele(0),p3_fin_ele(1)));
            theta_s = p3_target_spin.Theta();

            mom_ele = P4_fin_ele->P();
            theta_ele = P4_fin_ele->Theta();
            mom_had = P4_fin_had->P();
            theta_had = P4_fin_had->Theta();
            phi_ele = P4_fin_ele->Phi();
            phi_had = P4_fin_had->Phi();
            
            rapidity = 0.5 * log((mom_had + pz_had)/(mom_had - pz_had)  );

            rapidity = 0.5 * log((mom_had + pz_had)/(mom_had - pz_had)  );

            //phi_ele *= 180/3.1415926; if(py_ele<0.) phi_ele+=360;
            //phi_had *= 180/3.1415926; if(py_had<0.) phi_had+=360;
            //cout<<Form("--Electron: phi_gen = %f,  phi=%f", ph_ele, phi_ele)<<endl;
            //cout<<Form("--  Hadron: phi_gen = %f,  phi=%f", ph_had, phi_had)<<endl;

            jacoF = Jacobian(mom_ele,theta_ele,phi_ele,mom_had,theta_had,phi_had,mom_ion,energy_ion,mom_beam_ele);
            /*}}}*/

            delete P4_ini_ele; delete P4_ini_ion;  
            delete P4_fin_had; delete P4_fin_ele; delete P4_q; 
            delete lrz_P4_q; delete lrz_P4_h; delete lrz_P4_ef;
        } 
        /*}}}*/

        /*int CalcXS(){{{*/
        int CalcXS(){
            double bpt_m=4.694;
            double bpt_p=4.661;
            double pt_tmp = 0.0;	
            double dxs_all[3][4];
            double dxs_hp = -1000.0;
            double dxs_hm = -1000.0;

            /*Calculate XS{{{*/
            double mass_p = 0.938272;
            double mass_n = 0.939566;
            double kSinSQ = pow( sin(theta_ele*0.5),2);
            double kCosSQ = pow( cos(theta_ele*0.5),2);

            ////This definition is in the CM frame, i.e. dSigma/dx/dy
            //double xs_p = 2.0*PI*pow(1/137.036,2)*s/Q2/Q2*(1+pow(1-y, 2)) * fF2p * GeV2_to_nbarn;
            //double xs_n = 2.0*PI*pow(1/137.036,2)*s/Q2/Q2*(1+pow(1-y, 2)) * fF2n * GeV2_to_nbarn;
            //double Trans = 2.*mass_p*E0_beam_ele/mom_ele * PI *y; 

            ////This definition is in the lab frame, e.g. dSigma/dE'dOmega, where dE'dOmega = 2.*M*E0/E'*PI*y*dxdy
            double xs_p= 4.0*pow(1/137.036,2) * mom_ele*mom_ele/Q2/Q2
                * (1./nu *kCosSQ + 1/mass_p/x * kSinSQ) * fF2p * GeV2_to_nbarn;
            double xs_n= 4.0*pow(1/137.036,2) * mom_ele*mom_ele/Q2/Q2
                * (1./nu *kCosSQ + 1/mass_n/x * kSinSQ) * fF2n * GeV2_to_nbarn;

            XS_Inclusive = fZ * xs_p + (fA-fZ)*xs_n;//nbarn

            if (pt<0.8){
                //first method 	 
                //bpt_p = 1./(0.2+z*z*0.25);// <pt^2> = 0.2 GeV^2 (quark internal momentum)
                //Disable the TMD feature temparately, 02/01/2017
                bpt_p = 1./(0.2);// <pt^2> = 0.2 GeV^2 (quark internal momentum)
                bpt_m = bpt_p; // <kt^2> = 0.25 GeV^2 (struck quark internal momentum)

                //second method  use the original value
                // bpt_m=4.694; // HERMES parameterization
                // bpt_p=4.661;

                pt_tmp = pt;
                dxs(pt_tmp, bpt_m, bpt_p, &dxs_hp, &dxs_hm,dxs_all[0]);

                dilute_hp = dxs_all[0][0]/dxs_all[0][1];
                dilute_hm = dxs_all[0][2]/dxs_all[0][3];

            }
            else{
                // this part is to generate factor K
                //make sure the DXS is the same at PT= 0.8 GeV
                //calculating the TMD part
                double K[2],dxs_temp[10];
                //bpt_p = 1./(0.2+z*z*0.25);// <pt^2> = 0.2 GeV^2 (quark internal momentum)
                //Disable the TMD feature temparately, 02/01/2017
                bpt_p = 1./(0.2);// <pt^2> = 0.2 GeV^2 (quark internal momentum)
                bpt_m = bpt_p; // <kt^2> = 0.25 GeV^2 (struck quark internal momentum)

                //second method  use the original value
                // bpt_m=4.694; // HERMES parameterization
                // bpt_p=4.661;

                pt_tmp = 0.8;
                dxs(pt_tmp, bpt_m, bpt_p, &dxs_temp[0], &dxs_temp[1],dxs_all[0]);	   

                //calculating the TMD parts using the new vpt_p values
                //bpt_p = 1./(0.25+z*z*0.28);// <pt^2> = 0.25 GeV^2 (quark internal momentum)
                //Disable the TMD feature temparately, 02/01/2017
                bpt_p = 1./(0.25);// <pt^2> = 0.25 GeV^2 (quark internal momentum)
                bpt_m = bpt_p; // <kt^2> = 0.28 GeV^2 (struck quark internal momentum)

                //taking into account the NLO etc
                pt_tmp = 1.0;
                dxs(pt_tmp, bpt_m, bpt_p, &dxs_temp[2], &dxs_temp[3],dxs_all[1]);	   
                dxs_hpt(pt_tmp, &dxs_temp[4],&dxs_temp[5],dxs_all[2]);

                K[0] = (dxs_temp[0]-dxs_temp[2])/dxs_temp[4];
                K[1] = (dxs_temp[1]-dxs_temp[3])/dxs_temp[5];

                if (K[0]<0.) K[0] = 0.;
                if (K[1]<0.) K[1] = 0.;

                if (pt>1.2){
                    pt_tmp = pt;
                    dxs(pt_tmp,bpt_m, bpt_p, &dxs_temp[2], &dxs_temp[3],dxs_all[1]);	   
                    dxs_hpt(pt_tmp,&dxs_temp[4],&dxs_temp[5],dxs_all[2]);

                    dxs_hp = dxs_temp[2] + K[0]*dxs_temp[4];
                    dxs_hm = dxs_temp[3] + K[1]*dxs_temp[5];

                    dilute_hp = (dxs_all[1][0] + K[0]*dxs_all[2][0])/(dxs_all[1][1] + K[0]*dxs_all[2][1]);
                    dilute_hm = (dxs_all[1][2] + K[1]*dxs_all[2][2])/(dxs_all[1][3] + K[1]*dxs_all[2][3]);

                }else{
                    pt_tmp = 1.2;
                    dxs(pt_tmp, bpt_m, bpt_p, &dxs_temp[2], &dxs_temp[3],dxs_all[1]);	   
                    dxs_hpt(pt_tmp,&dxs_temp[4],&dxs_temp[5],dxs_all[2]);

                    dxs_temp[2] = dxs_temp[2] + K[0]*dxs_temp[4];
                    dxs_temp[3] = dxs_temp[3] + K[1]*dxs_temp[5];

                    dxs_temp[4] = (dxs_temp[0]-dxs_temp[2])/(1./0.8/0.8-1./1.2/1.2);
                    dxs_temp[5] = dxs_temp[0]-dxs_temp[4]/0.8/0.8;
                    dxs_temp[6] = (dxs_temp[1]-dxs_temp[3])/(1./0.8/0.8-1./1.2/1.2);
                    dxs_temp[7] = dxs_temp[1]-dxs_temp[6]/0.8/0.8;

                    dxs_hp = dxs_temp[4]/pt/pt + dxs_temp[5];
                    dxs_hm = dxs_temp[6]/pt/pt + dxs_temp[7];

                    dilute_hp = ((dxs_all[1][0] + K[0]*dxs_all[2][0])
                            /(dxs_all[1][1] + K[0]*dxs_all[2][1])+dxs_all[0][0]/dxs_all[0][1])/2.;
                    dilute_hm = ((dxs_all[1][2] + K[1]*dxs_all[2][2])
                            /(dxs_all[1][3] + K[1]*dxs_all[2][3])+dxs_all[0][2]/dxs_all[0][3])/2.;
                }
            }
            /*}}}*/
 
            //Caculate Boer-Mulder TMD Asymmetry (A= dXS*cos2phi/dXS
            double Asym_BM_hp, Asym_BM_hm;
            double Asym_BM_all[4];
            pt_tmp = pt;
            int kErrBM = GetAsym_Cos2Phi(x, y, z, Q2, pt_tmp, phi_h, &Asym_BM_hp, &Asym_BM_hm, &Asym_BM_all[0]);

            /*try to take care of the decay{{{*/
            double decay_part = 0.0;
            if (abs(particle_flag)==1){
                if (theta_had>155./180.*3.1415926){
                    decay_part = exp(7./cos(theta_had)*mass_hadron/(2.6*mom_had*3.0));
                }else if (theta_had<=155./180.*3.1415926&&theta_had>=140./180.*3.1415926){
                    decay_part = exp(4.5/cos(theta_had)*mass_hadron/(2.6*mom_had*3.0));
                }else if (theta_had <140/180.*3.1415926){
                    decay_part = exp(-2.5/sin(theta_had)*mass_hadron/(2.6*mom_had*3.0));
                }
            }else{
                if (theta_had>155./180.*3.1415926){
                    decay_part = exp(7./cos(theta_had)*mass_hadron/(1.24*mom_had*3.0));
                }else if (theta_had<=155./180.*3.1415926&&theta_had>=140./180.*3.1415926){
                    decay_part = exp(4.5/cos(theta_had)*mass_hadron/(1.24*mom_had*3.0));
                }else if (theta_had <140/180.*3.1415926){
                    decay_part = exp(-2.5/sin(theta_had)*mass_hadron/(1.24*mom_had*3.0));
                }
            }
            /*}}}*/

            //cerr<<Form("++++ pt=%f, jacoF=%f, dxs_hp=%f,dxs_hm=%f, decay=%f, dilu_hp=%f, dilu_hm=%f", pt, jacoF,dxs_hp,dxs_hm,decay_part,dilute_hp,dilute_hm)<<endl;

            /*Save&Return{{{*/
            dxs_hp = decay_part * dxs_hp; 
            dxs_hm = decay_part * dxs_hm;

            int err_hp = -1, err_hm=-1;
            if (dxs_hp>0.&& dxs_hp<10000.){
                err_hp = 0;
            }
            else{
                dxs_hp = 0;
                err_hp = 1;
            }

            XS_HP = dxs_hp;
            XS_HM = dxs_hm;

            if (dxs_hm>0.&& dxs_hm<10000.){
                err_hm = 0;
            }
            else{
                dxs_hm = 0;
                err_hm = 1;
            }

            if(err_hp==0&&err_hm==0) 
                return 0;
            else if(err_hp==0&&err_hm==1)
                return 1;
            else if(err_hp==1&&err_hm==0)
                return -1;
            else
                return 2;
            /*}}}*/
        } 
        /*}}}*/

        /*Return Values{{{*/

        double GetXS_Inclusive(){
            return XS_Inclusive;
        }

        double GetXS_HP(){
            return XS_HP;
        }
        double GetXS_HM(){
            return XS_HM;
        }
        double GetDilute_HP(){
            return dilute_hp;
        }
        double GetDilute_HM(){
            return dilute_hm;
        }
        double get_uA(){ return fuA; }//return the average u
        double get_dA(){ return fdA; }//return the average d
        double get_s(){ return fs; }//return the average s
        double get_g(){ return fg; }//return the average g
        double get_ubar(){ return fubar; }//return the average u
        double get_dbar(){ return fdbar; }//return the average d
        double get_sbar(){ return fsbar; }//return the average s
        /*}}}*/
        
        /*double Run_CTEQPDF(double x,double Q2){{{*/          
        void Run_CTEQPDF(double ix,double iQ2)          
        {
            fubar = Get_CTEQPDF(-1, ix, iQ2);
            fdbar = Get_CTEQPDF(-2, ix, iQ2);
            fsbar = Get_CTEQPDF(-3, ix, iQ2);
            fuA = Get_CTEQPDF(1, ix, iQ2);
            fdA = Get_CTEQPDF(2, ix, iQ2);
            fs = Get_CTEQPDF(3, ix, iQ2);
            fg = Get_CTEQPDF(0, ix, iQ2);

            //Temperately set to zero, 02/01/2017
            fs =0.0;
            fsbar=0.0;
            fg = 0.0;
        }
        /*}}}*/

        void RunEPS09(int iOrder, int iErrSet,int iA,int iZ, double ix, double iQ2){/*{{{*/

            double iR_uv = 0.0, iR_dv = 0.0, iR_u = 0.0, iR_d = 0.0, iR_s = 0.0, iR_c = 0.0, iR_b = 0.0, iR_g = 0.0; 
            //Order; //1->LO use CTEQ6L1, 2->NLO use CETQ6.1M
            //ErrSet;//1->central fit, 2,3-> err set#1, 4,5->err set#2, ...,30,31->err set#15
            eps09(iOrder, iErrSet, iA, ix, sqrt(iQ2), iR_uv, iR_dv, iR_u, iR_d, iR_s, iR_c,iR_b,iR_g);

            double ubar = Get_CTEQPDF(-1, ix, iQ2);
            double dbar = Get_CTEQPDF(-2, ix, iQ2);
            double sbar = Get_CTEQPDF(-3, ix, iQ2);
            double u = Get_CTEQPDF(1, ix, iQ2);
            double d = Get_CTEQPDF(2, ix, iQ2);
            double s = Get_CTEQPDF(3, ix, iQ2);
            double uv = u-ubar;//uv = u - usea
            double dv = d-ubar;//dv = d - dsea
            double g = Get_CTEQPDF(0, ix, iQ2);

            fuA = (double)iZ/(double)iA * ( iR_uv * uv + iR_u * ubar)
                + (double)(iA-iZ)/(double)iA * (iR_dv * dv + iR_d * dbar);

            fdA = (double)iZ/(double)iA * ( iR_dv * dv + iR_d * dbar)
                + (double)(iA-iZ)/(double)iA * (iR_uv * uv + iR_u * ubar);

            fubar =  iR_u * ubar;
            fdbar =  iR_d * dbar;

            fs =  iR_s * s;
            fsbar =  iR_s *sbar;

            fg =  iR_g * g;


            //Temperately set to zero, 02/01/2017
            fs =0.0;
            fsbar=0.0;
            fg = 0.0;
        }/*}}}*/
        
        /*GetAsym_Cos2Phi(){{{*/
        //The unpolarized SIDIS XS is already a six-fold differential function of (x,y,z, pT, phi_S, phi_h)
        //while the pT dependence is modeled by Gaussian Ansetz, while the smearing is given by:
        //
        //   //first method 	 
        //        bpt_p = 1./(0.2+z*z*0.25);// 
        //        bpt_m = bpt_p; 
        //   where the intial values in Xin's old code
        //        <pt^2> = 0.2 GeV^2 (quark internal momentum)
        //        <kt^2> = 0.25 GeV^2 (struck quark internal momentum)
        //   Using Anselmino 2014 Fit, JHEP 04 (2014)005
        //        <pt^2> = 0.57+/-0.08 GeV^2 (quark internal momentum)
        //        <kt^2> = 0.12+/-0.01 GeV^2 (struck quark internal momentum)
        //
        //   //or, second method  use the original value
        //        // bpt_m=4.694; // HERMES parameterization
        //        // bpt_p=4.661;
        //pT Broadening has been considered for nuclear targets
        //And the SIDIS cross section is:
        //  *dxs_hp = (*dxs_hp)*bpt_p/PI*exp(-bpt_p*pt_tmp*pt_tmp)*x;
        //  *dxs_hm = (*dxs_hm)*bpt_m/PI*exp(-bpt_m*pt_tmp*pt_tmp)*x;
        //
        //Now, add the Boer-Mulder Amplitude which is the convolution of the BM-TMD and Collins FF.
        //Besides, the Twist-4 effect, e.g. Cahn effect,can also be added 
        //And, we also need to take into account if the target is a nucleus instead of free nucleon.
        //
        int GetAsym_Cos2Phi(double kXb, double kY, double kZ, double kQ2, double kPt, double kPhi, double* kA_hp, double* kA_hm, double *kA_all){
        //Barone et. al., PRD 81, 114026 2010.
            /*f1 and D1 {{{*/
            const double qu=2./3.;
            const double qd=-1./3.;
            const double qs=-1./3.;
            const double kMKtSQ = 0.25;//GeV^2
            const double kMPtSQ = 0.20;//GeV^2
            const double kMass_Proton= 0.938272;

            double kD_fav=0.0, kD_unfav=0.0, kD_s = 0.0, kD_g=0.0;
            double kMass_Hadron=0.0;
            double kPtSQ = pow(kPt,2); 

            double kA1_p_hp=0.0, kA1_n_hp =0.0;
            double kA1_p_hm=0.0, kA1_n_hm =0.0;
            if (fabs(particle_flag)==1) {/*{{{*/
                kMass_Hadron = 0.13957;
                //pion fragmentation functions
                Unpol_FF(1,z,Q2,&kD_fav,&kD_unfav,&kD_s,&kD_g);

                //proton
                kA1_p_hp = qu*qu*fuA*kD_fav + qu*qu*fubar*kD_unfav + qd*qd*fdA*kD_unfav + qd*qd*fdbar*kD_fav + qs*qs*fs*kD_s + qs*qs*fsbar*kD_s;
                kA1_p_hm = qu*qu*fuA*kD_unfav + qu*qu*fubar*kD_fav + qd*qd*fdA*kD_fav + qd*qd*fdbar*kD_unfav + qs*qs*fs*kD_s + qs*qs*fsbar*kD_s;

                //neutron
                //u->d, d->u, ubar->dbar, dbar->ubar
                kA1_n_hp = qu*qu*fdA*kD_fav + qu*qu*fdbar*kD_unfav + qd*qd*fuA*kD_unfav + qd*qd*fubar*kD_fav + qs*qs*fs*kD_s + qs*qs*fsbar*kD_s;
                kA1_n_hm = qu*qu*fdA*kD_unfav + qu*qu*fdbar*kD_fav + qd*qd*fuA*kD_fav + qd*qd*fubar*kD_unfav + qs*qs*fs*kD_s + qs*qs*fsbar*kD_s;
            }else{
                kMass_Hadron = 0.493667;
                //kaon, fragmentation functions
                Unpol_FF(2,z,Q2,&kD_fav,&kD_unfav,&kD_s,&kD_g);

                //proton
                kA1_p_hp = qu*qu*fuA*kD_fav + qu*qu*fubar*kD_unfav + qd*qd*fdA*kD_s + qd*qd*fdbar*kD_s + qs*qs *fs * kD_unfav + qs*qs*fsbar*kD_fav;
                kA1_p_hm = qu*qu*fuA*kD_unfav + qu*qu*fubar*kD_fav + qd*qd*fdA*kD_s + qd*qd*fdbar*kD_s + qs*qs *fs * kD_fav + qs*qs*fsbar*kD_unfav;

                //neutron
                //u->d, d->u, ubar->dbar, dbar->ubar
                kA1_n_hp = qu*qu*fdA*kD_fav   + qu*qu*fdbar*kD_unfav + qd*qd*fuA*kD_s + qd*qd*fubar*kD_s + qs*qs *fs * kD_unfav + qs*qs*fsbar*kD_fav;
                kA1_n_hm = qu*qu*fdA*kD_unfav + qu*qu*fdbar*kD_fav   + qd*qd*fuA*kD_s + qd*qd*fubar*kD_s + qs*qs *fs * kD_fav   + qs*qs*fsbar*kD_unfav;
            }/*}}}*/

            double kA2 = (PI*PI*kMKtSQ*kMPtSQ);
            double kA3 = (pow(kMKtSQ,3)*kMPtSQ*PI*kPtSQ*pow(kZ,2))/(pow(EXP,kPtSQ/(kMPtSQ + kMKtSQ*pow(kZ,2)))*kQ2*pow(kMPtSQ + kMKtSQ*pow(kZ,2),3));
            double kA_p_hp = kA1_p_hp/kA2 * kA3;
            double kA_p_hm = kA1_p_hm/kA2 * kA3;
            double kA_n_hp = kA1_n_hp/kA2 * kA3;
            double kA_n_hm = kA1_n_hm/kA2 * kA3;

            double kC3 = (kMKtSQ*kMPtSQ*PI)/(pow(EXP,kPtSQ/(kMPtSQ + kMKtSQ*pow(kZ,2)))*(kMPtSQ + kMKtSQ*pow(kZ,2)));
            double kC_p_hp = kA1_p_hp/kA2 * kC3;
            double kC_p_hm = kA1_p_hm/kA2 * kC3;
            double kC_n_hp = kA1_n_hp/kA2 * kC3;
            double kC_n_hm = kA1_n_hm/kA2 * kC3;

            //cout<<Form("A: 1 = %f, 2 = %f, 3 = %f, C=%f", kA1_p_hp, kA2, kA3, kC3)<<endl;
            
            /*}}}*/

            /*Boer Mulder TMD and Collins FF{{{*/
            /*BM-TMD: Fitting from Barone et. al., PRD 81, 114026 2010{{{*/
            const double kAu = -0.35;
            const double kAd = 0.90;
            const double kAs = 0.24;
            //const double kAubar = -0.04;
            //const double kAdbar = 0.40;
            //const double kAsbar = -1.0;
            const double kAlphaU = 0.73; 
            const double kAlphaD = 1.08; 
            const double kAlphaSea = 0.79; 
            const double kBeta = 3.46; 
            const double kM1SQ = 0.34; // (GeV/c)^2
            const double kLambdaU = 2.0;
            const double kLambdaD = -1.111;
            const double kLambdaSea = 0.0; //Set to zero so far

            double kEta = sqrt(2.0*EXP) * kMass_Proton/sqrt(kM1SQ);//move this term into the integral (kB3) exp(-kKt*kKt/kMcSQ);
            double kRhoU = kAu*pow(kXb,kAlphaU)*pow((1-kXb),kBeta)
                *pow((kAlphaU+kBeta),(kAlphaU+kBeta))/pow(kAlphaU, kAlphaU)/pow(kBeta, kBeta);
            double kRhoD = kAd*pow(kXb,kAlphaD)*pow((1-kXb),kBeta)
                *pow((kAlphaD+kBeta),(kAlphaD+kBeta))/pow(kAlphaD, kAlphaD)/pow(kBeta, kBeta);
            double kRhoSea = kAs*pow(kXb,kAlphaSea)*pow((1-kXb),kBeta)
                *pow((kAlphaSea+kBeta),(kAlphaSea+kBeta))/pow(kAlphaSea, kAlphaSea)/pow(kBeta, kBeta);
            
            //Now constract h1_perp (BM) for different quark:
            double kU_BM = kLambdaU * kRhoU * kEta * fuA;
            double kD_BM = kLambdaD * kRhoD * kEta * fdA;
            double kUbar_BM = kLambdaSea * kRhoSea * kEta * fubar;
            double kDbar_BM = kLambdaSea * kRhoSea * kEta * fdbar;
            double kS_BM = kLambdaSea * kRhoSea * kEta * fs;
            double kSbar_BM = kLambdaSea * kRhoSea * kEta * fsbar;
            /*}}}*/

            /*Collins FF(){{{*/
            //BM-TMD: Fitting from Barone et. al., PRD 81, 114026 2010.
            const double kAc_fav = 0.44;
            const double kAc_unfav = -1.00;
            const double kGamma = 0.96;
            const double kDelta = 0.01;
            const double kMcSQ = 0.91; //GeV^2/c

            const double kRhoC = pow(kZ,kGamma)*pow((1-kZ),kDelta)
                *pow((kGamma+kDelta),(kGamma+kDelta))/pow(kGamma, kGamma)/pow(kDelta, kDelta);
            const double kEtaC = sqrt(2.0*EXP) * kZ * kMass_Hadron/sqrt(kMcSQ);//moving into the integral in B3: * exp(-kP_perp*kP_perp/kMcSQ);

            double kH_fav = kAc_fav * kRhoC * kEtaC * kD_fav;
            double kH_unfav = kAc_unfav * kRhoC * kEtaC * kD_unfav;
            double kH_s = 0.0; //set to zero temperately
            //double kH_g = 0.0; //set to zero temperately
            /*}}}*/

            double kB1_p_hp=0.0, kB1_n_hp =0.0;
            double kB1_p_hm=0.0, kB1_n_hm =0.0;
            if (fabs(particle_flag)==1) {/*{{{*/
                //proton
                kB1_p_hp = qu*qu*kU_BM*kH_fav + qu*qu*kUbar_BM*kH_unfav + qd*qd*kD_BM*kH_unfav + qd*qd*kDbar_BM*kH_fav + qs*qs*kS_BM*kH_s + qs*qs*kSbar_BM*kH_s;
                kB1_p_hm = qu*qu*kU_BM*kH_unfav + qu*qu*kUbar_BM*kH_fav + qd*qd*kD_BM*kH_fav + qd*qd*kDbar_BM*kH_unfav + qs*qs*kS_BM*kH_s + qs*qs*kSbar_BM*kH_s;

                //neutron
                //u->d, d->u, ubar->dbar, dbar->ubar
                kB1_n_hp = qu*qu*kD_BM*kH_fav + qu*qu*kDbar_BM*kH_unfav + qd*qd*kU_BM*kH_unfav + qd*qd*kUbar_BM*kH_fav + qs*qs*kS_BM*kH_s + qs*qs*kSbar_BM*kH_s;
                kB1_n_hm = qu*qu*kD_BM*kH_unfav + qu*qu*kDbar_BM*kH_fav + qd*qd*kU_BM*kH_fav + qd*qd*kUbar_BM*kH_unfav + qs*qs*kS_BM*kH_s + qs*qs*kSbar_BM*kH_s;
            }else{
                //proton
                kB1_p_hp = qu*qu*kU_BM*kH_fav + qu*qu*kUbar_BM*kH_unfav + qd*qd*kD_BM*kH_s + qd*qd*kDbar_BM*kH_s + qs*qs *kS_BM * kH_unfav + qs*qs*kSbar_BM*kH_fav;
                kB1_p_hm = qu*qu*kU_BM*kH_unfav + qu*qu*kUbar_BM*kH_fav + qd*qd*kD_BM*kH_s + qd*qd*kDbar_BM*kH_s + qs*qs *kS_BM * kH_fav + qs*qs*kSbar_BM*kH_unfav;

                //neutron
                //u->d, d->u, ubar->dbar, dbar->ubar
                kB1_n_hp = qu*qu*kD_BM*kH_fav   + qu*qu*kDbar_BM*kH_unfav + qd*qd*kU_BM*kH_s + qd*qd*kUbar_BM*kH_s + qs*qs *kS_BM * kH_unfav + qs*qs*kSbar_BM*kH_fav;
                kB1_n_hm = qu*qu*kD_BM*kH_unfav + qu*qu*kDbar_BM*kH_fav   + qd*qd*kU_BM*kH_s + qd*qd*kUbar_BM*kH_s + qs*qs *kS_BM * kH_fav   + qs*qs*kSbar_BM*kH_unfav;
            }/*}}}*/

            double kB2 = PI*PI*kMKtSQ*kMPtSQ;
            double kB3 = (pow(kM1SQ,2)*kMcSQ*pow(kMKtSQ,2)*kMPtSQ*PI*sqrt(kPtSQ)
                    *(pow(kMcSQ,2)*kMKtSQ*pow(kMPtSQ,2)*sqrt((pow(kMcSQ + kMPtSQ,2)*kPtSQ*pow(kZ,2))/(pow(kMcSQ,2)*pow(kMPtSQ,2))) 
                        - kM1SQ*(kMKtSQ*pow(kMPtSQ,2)*sqrt(kPtSQ)*pow(kZ,3) 
                            - kMcSQ*kMKtSQ*kMPtSQ*pow(kZ,2)*(-2*sqrt(kPtSQ)*kZ + kMPtSQ*sqrt((pow(kMcSQ + kMPtSQ,2)*kPtSQ*pow(kZ,2))/(pow(kMcSQ,2)*pow(kMPtSQ,2)))) 
                            + pow(kMcSQ,2)*(-(pow(kMPtSQ,2)*sqrt((pow(kMcSQ + kMPtSQ,2)*kPtSQ*pow(kZ,2))/(pow(kMcSQ,2)*pow(kMPtSQ,2)))) 
                                + kMKtSQ*pow(kZ,2)*(sqrt(kPtSQ)*kZ - kMPtSQ*sqrt((pow(kMcSQ + kMPtSQ,2)*kPtSQ*pow(kZ,2))/(pow(kMcSQ,2)*pow(kMPtSQ,2))))))))
                /(pow(EXP,((kM1SQ + kMKtSQ)*(kMcSQ + kMPtSQ)*kPtSQ)/(kMcSQ*kMKtSQ*kMPtSQ + kM1SQ*(kMKtSQ*kMPtSQ*pow(kZ,2) + kMcSQ*(kMPtSQ + kMKtSQ*pow(kZ,2)))))
                        *kMass_Proton*kMass_Hadron*kZ*pow(kMcSQ*kMKtSQ*kMPtSQ + kM1SQ*(kMKtSQ*kMPtSQ*pow(kZ,2) + kMcSQ*(kMPtSQ + kMKtSQ*pow(kZ,2))),3));

            double kB_p_hp = kB1_p_hp / kB2 * kB3;
            double kB_p_hm = kB1_p_hm / kB2 * kB3;
            double kB_n_hp = kB1_n_hp / kB2 * kB3;
            double kB_n_hm = kB1_n_hm / kB2 * kB3;

            //cout<<Form("B: 1 = %f, 2 = %f, 3 = %f", kB1_p_hp, kB2, kB3)<<endl;
            /*}}}*/

            double kNorm = 4.0*(1.0-kY)/(1.0+pow((1.0+kY), 2));

            double kACos2Phi_p_hp = kNorm*(kA_p_hp + kB_p_hp/2.0)/kC_p_hp;
            double kACos2Phi_n_hp = kNorm*(kA_n_hp + kB_n_hp/2.0)/kC_n_hp;
            double kACos2Phi_p_hm = kNorm*(kA_p_hm + kB_p_hm/2.0)/kC_p_hm;
            double kACos2Phi_n_hm = kNorm*(kA_n_hm + kB_n_hm/2.0)/kC_n_hp;

            Asym_Cos2Phi_p_hp = kNorm*(kA_p_hp + kB_p_hp/2.0)/kC_p_hp;
            Asym_Cos2Phi_p_hm = kNorm*(kA_p_hm + kB_p_hm/2.0)/kC_p_hm;
            Asym_Cos2Phi_n_hp = kNorm*(kA_n_hp + kB_n_hp/2.0)/kC_n_hp;
            Asym_Cos2Phi_n_hm = kNorm*(kA_n_hm + kB_n_hm/2.0)/kC_n_hm;
            
            Asym_Cahn_p_hp = kNorm*(kA_p_hp)/kC_p_hp;
            Asym_Cahn_p_hm = kNorm*(kA_p_hm)/kC_p_hm;
            Asym_Cahn_n_hp = kNorm*(kA_n_hp)/kC_n_hp;
            Asym_Cahn_n_hm = kNorm*(kA_n_hm)/kC_n_hm;

            Asym_BM_p_hp = kNorm*(kB_p_hp/2.0)/kC_p_hp;
            Asym_BM_p_hm = kNorm*(kB_p_hm/2.0)/kC_p_hm;
            Asym_BM_n_hp = kNorm*(kB_n_hp/2.0)/kC_n_hp;
            Asym_BM_n_hm = kNorm*(kB_n_hm/2.0)/kC_n_hm;

            *(kA_all+0) = kACos2Phi_p_hp;
            *(kA_all+1) = kACos2Phi_n_hp;
            *(kA_all+2) = kACos2Phi_p_hm;
            *(kA_all+3) = kACos2Phi_n_hm;

            //*kA_hp = (kACos2Phi_p_hp*fZ+kACos2Phi_n_hp*(fA-fZ));
            //*kA_hp = (kACos2Phi_p_hm*fZ+kACos2Phi_n_hm*(fA-fZ));

            *kA_hp = kNorm*( (kA_p_hp*fA+kA_n_hp*(fA-fZ)) + (kB_p_hp*fA+kB_n_hp*(fA-fZ))/2.0 )/( kC_p_hp*fA+kC_n_hp*(fA-fZ) );
            *kA_hm = kNorm*( (kA_p_hm*fA+kA_n_hm*(fA-fZ)) + (kB_p_hm*fA+kB_n_hm*(fA-fZ))/2.0 )/( kC_p_hm*fA+kC_n_hm*(fA-fZ) );

            Asym_Cos2Phi_hp = kNorm*( (kA_p_hp*fA+kA_n_hp*(fA-fZ)) + (kB_p_hp*fA+kB_n_hp*(fA-fZ))/2.0 )/( kC_p_hp*fA+kC_n_hp*(fA-fZ) );
            Asym_Cos2Phi_hm = kNorm*( (kA_p_hm*fA+kA_n_hm*(fA-fZ)) + (kB_p_hm*fA+kB_n_hm*(fA-fZ))/2.0 )/( kC_p_hm*fA+kC_n_hm*(fA-fZ) );

            Asym_Cahn_hp = kNorm*( (kA_p_hp*fA+kA_n_hp*(fA-fZ)) )/( kC_p_hp*fA+kC_n_hp*(fA-fZ) );
            Asym_Cahn_hm = kNorm*( (kA_p_hm*fA+kA_n_hm*(fA-fZ)) )/( kC_p_hm*fA+kC_n_hm*(fA-fZ) );

            Asym_BM_hp = kNorm*( (kB_p_hp*fA+kB_n_hp*(fA-fZ))/2.0 )/( kC_p_hp*fA+kC_n_hp*(fA-fZ) );
            Asym_BM_hm = kNorm*( (kB_p_hm*fA+kB_n_hm*(fA-fZ))/2.0 )/( kC_p_hm*fA+kC_n_hm*(fA-fZ) );

            if(isnan(*kA_hp) || isnan(*kA_hp))
                return -1;
            else
                return 0;
        }
        /*}}}*/

    private:
        /*double Azimuthalphi(double vx, double vy){{{*/
        double Azimuthalphi(double vx, double vy){
            double pmod, phi, cosf;
            pmod=vx*vx+vy*vy;
            if(pmod>0.) {
                pmod=sqrt(pmod);
                cosf=vx/pmod;
            }
            else{
                phi=-100.;
                cosf=10.;
            }
            if(fabs(cosf)<=1.0) phi=acos(cosf);
            if(vy<0.) phi=2*PI-phi;
            return phi;
        }
        /*}}}*/

        /*double Jacobian(double P_ef,double theta_ef,double phi_ef,...){{{*/
        double Jacobian(double P_ef,double theta_ef,double phi_ef,double P_hf,double theta_hf,double phi_hf,double mom_ion,double mom_beam_ion,double P_ei){


            // need to shift the electron phi angle to 0 in order to do calculation for jacobian

            phi_hf -= phi_ef;
            phi_ef = 0.;

            double qbeta = -1.*mom_ion/mom_beam_ion;
            //double qbeta = 0.;
            double qgamma = 1./sqrt(1-qbeta*qbeta);


            double mass_pi = mass_hadron;
            double P_hx,P_hy,P_hz;  
            P_hx = P_hf*sin(theta_hf)*cos(phi_hf);
            P_hy = P_hf*sin(theta_hf)*sin(phi_hf);
            P_hz = qgamma * (P_hf*cos(theta_hf) - qbeta * sqrt(P_hf*P_hf+mass_pi*mass_pi));



            //P_hz = P_hf*cos(theta_hf);

            double qz = qgamma*(P_ei-P_ef*cos(theta_ef)-qbeta*(P_ei-P_ef));

            double alpha,beta,gamma1;
            alpha = atan(P_ef*sin(theta_ef)*cos(phi_ef)/qz);
            beta = atan(-P_ef*sin(theta_ef)*sin(phi_ef)*
                    qz/cos(alpha)/
                    (pow(P_ef*sin(theta_ef)*sin(phi_ef),2) + 
                     pow(qz,2)));
            gamma1 = atan(sin(beta)/tan(alpha));

            double a2,b2,c2,d2,e2,f2;
            a2 = (-sin(gamma1)*cos(alpha)+cos(gamma1)*sin(beta)*sin(alpha));
            b2 = cos(gamma1)*cos(beta);
            c2 = -(sin(gamma1)*sin(alpha)+cos(gamma1)*sin(beta)*cos(alpha));
            d2 = cos(gamma1)*cos(alpha) + sin(gamma1)*sin(beta)*sin(alpha);
            e2 = sin(gamma1)*cos(beta);
            f2 = cos(gamma1)*sin(alpha) - sin(gamma1)*sin(beta)*cos(alpha);

            double a1,b1,c1,d1,e1,f1;
            a1 = cos(alpha)*cos(alpha) + sin(beta)*sin(beta) * sin(alpha)*sin(alpha);
            b1 = cos(beta)* cos(beta);
            c1 = sin(alpha)*sin(alpha) + sin(beta)*sin(beta) * cos(alpha)*cos(alpha);      
            d1 = 2*cos(beta)*sin(beta)*sin(alpha);
            e1 = -2*cos(beta)*sin(beta)*cos(alpha);
            f1 = 2*cos(beta)*cos(beta)*sin(alpha)*cos(alpha);

            double Px1,Py1,Pz1,Px2,Py2,Pz2;
            Px1 = P_hf*sin(theta_hf)*(-sin(phi_hf));
            Py1 = P_hf*sin(theta_hf)*cos(phi_hf);
            Pz1 = 0.;

            Px2 = -P_hf*cos(phi_hf)/tan(theta_hf);
            Py2 = -P_hf*sin(phi_hf)/tan(theta_hf);
            Pz2 = P_hf*qgamma;

            double Px3 = sin(theta_hf)*cos(phi_hf);
            double Py3 = sin(theta_hf)*sin(phi_hf);
            double Pz3 = qgamma*(cos(theta_hf) -qbeta/sqrt(P_hf*P_hf+mass_pi*mass_pi)*P_hf);

            double P_hyp,P_hxp;
            P_hyp = a2*P_hx + b2*P_hy + c2*P_hz;
            P_hxp = d2*P_hx + e2*P_hy + f2*P_hz;


            //cout << a1*P_hx*P_hx + b1*P_hy*P_hy + c1*P_hz*P_hz + d1*P_hx*P_hy + e1 *P_hy*P_hz + f1 * P_hx*P_hz << "\t" << atan(P_hyp/P_hxp) << endl;

            double jaco_F;
            double jaco2; // py/pa
            double jaco1; // px/pb
            double jaco3; // pz/pd
            double jaco4; // pl/pc
            double jaco5; // pm/pe
            double jaco6; // pm/pf
            double jaco7; // pn/pe
            double jaco8; // pn/pf

            double jaco9;  // px/pa
            double jaco10; // py/pb
            double jaco11; // pz/pe
            double jaco12; // pm/pd
            double jaco13; // pn/pd

            jaco9 = (P_ei*(1.-cos(theta_ef))*(P_ei*(mom_beam_ion+mom_ion)-P_ef*(mom_beam_ion+mom_ion*cos(theta_ef)))+P_ef*P_ei*(1-cos(theta_ef))*(mom_beam_ion+mom_ion*cos(theta_ef)))/(P_ei*(mom_beam_ion+mom_ion) - P_ef*(mom_beam_ion+mom_ion*cos(theta_ef)))/(P_ei*(mom_beam_ion+mom_ion) - P_ef*(mom_beam_ion+mom_ion*cos(theta_ef))); // px / pa

            jaco1 = (-P_ef*P_ei*(P_ei*(mom_ion+mom_beam_ion)-P_ef*(mom_beam_ion+mom_ion*cos(theta_ef)))+P_ei*P_ef*(1-cos(theta_ef))*P_ef*mom_ion)/(P_ei*(mom_beam_ion+mom_ion) - P_ef*(mom_beam_ion+mom_ion*cos(theta_ef)))/(P_ei*(mom_beam_ion+mom_ion) - P_ef*(mom_beam_ion+mom_ion*cos(theta_ef))); // p x / pb



            jaco2 = - (mom_beam_ion+mom_ion*cos(theta_ef))/(P_ei*(mom_beam_ion+mom_ion)); // py/pa

            jaco10 = (-1*P_ef*mom_ion)/(P_ei*(mom_beam_ion+mom_ion)); // py / pb

            jaco3 = (P_hf*mom_beam_ion/sqrt(P_hf*P_hf + mass_pi*mass_pi) + mom_ion*cos(theta_hf))/(P_ei*(mom_beam_ion+mom_ion) - P_ef*(mom_beam_ion+mom_ion*cos(theta_ef))); // p z / p d
            jaco11 = mom_ion*P_hf/(P_ei*(mom_beam_ion+mom_ion) - P_ef*(mom_beam_ion+mom_ion*cos(theta_ef))); // p z/ p e

            jaco4 = -1;

            jaco5 = (2*a1*P_hx+d1*P_hy+f1*P_hz)*Px1
                + (2*b1*P_hy+d1*P_hx+e1*P_hz)*Py1; // p pt^2/ p (phi_pi)   or p m / pf

            jaco6 = (2*a1*P_hx+d1*P_hy+f1*P_hz)*Px2
                + (2*b1*P_hy+d1*P_hx+e1*P_hz)*Py2
                + (2*c1*P_hz+e1*P_hy+f1*P_hx)*Pz2; // ppt^2/ p (cos(theta_pi)) or p m / p e

            jaco7 = (P_hxp*(a2*Px1 + b2*Py1) - 
                    P_hyp*(d2*Px1 + e2*Py1) )/
                (P_hxp*P_hxp + P_hyp*P_hyp); // p n / pf

            jaco8 = (P_hxp*(a2*Px2 + b2*Py2 + c2*Pz2) - 
                    P_hyp*(d2*Px2 + e2*Py2 + f2*Pz2) )/
                (P_hxp*P_hxp + P_hyp*P_hyp); // p n / p e



            jaco12 = (2*a1*P_hx+d1*P_hy+f1*P_hz)*Px3
                + (2*b1*P_hy+d1*P_hx+e1*P_hz)*Py3
                + (2*c1*P_hz+e1*P_hy+f1*P_hx)*Pz3; // or p m / p d
            jaco13 = (P_hxp*(a2*Px3 + b2*Py3 + c2*Pz3) - 
                    P_hyp*(d2*Px3 + e2*Py3 + f2*Pz3) )/
                (P_hxp*P_hxp + P_hyp*P_hyp); // p n / p d

            //cout << jaco13 << endl;

            //(pz/pd *(p m / p e * p n / pf - p m / pf * // p n / p e) - pz/pe*(p m / p d * p n / pf - p n / p d * p m / pf))

            jaco_F = jaco4*(jaco2*jaco1-jaco9*jaco10)*(jaco3*(jaco6*jaco7-jaco5*jaco8)-jaco11*(jaco12*jaco7-jaco13*jaco5));
            jaco_F = fabs(jaco_F);

            //cout << jaco_F << "\t" << jaco1 << "\t" << jaco2 << "\t" << jaco3 << "\t" << jaco5 << "\t" << jaco6 << "\t" << jaco7 << "\t" << jaco8 << endl;

            return jaco_F;
        }
        /*}}}*/

        /*void dxs(double x, double y, double Q2, ...){{{*/
        //////////////////////////////////////////////////////////////////////////////////////////////////
        // This subroutine is called to compute new XS for any nuclear A
        // where the PDFs contain the free PDF from CTEQ and the nulcear-medium-modification from  EPS09
        // -- Zhihong Ye, 06/27/2016
        //////////////////////////////////////////////////////////////////////////////////////////////////
        void dxs(double pt_tmp,double bpt_m, double bpt_p, double* dxs_hp, double* dxs_hm,double* dxs_all){
            double qu=2./3.;
            double qd=-1./3.;
            double qs=-1./3.;

            *dxs_hp = 1./137.035/137.035/x/y/Q2 *y*y/(2*(1-epsilon))*(1+gamma*gamma/2./x)*GeV2_to_nbarn;
            *dxs_hm = *dxs_hp;

            *dxs_hp = (*dxs_hp)*bpt_p/PI*exp(-bpt_p*pt_tmp*pt_tmp)*x;
            *dxs_hm = (*dxs_hm)*bpt_m/PI*exp(-bpt_m*pt_tmp*pt_tmp)*x;

            double df_p_hp=0,df_p_hm=0;
            double df_n_hp=0,df_n_hm=0;

            double uquark=0.0,dquark=0.0,squark=0.0,ubarquark=0.0,dbarquark=0.0,sbarquark=0.0;
            //if(fOrder==0){
                //uquark = Get_LHAPDF(1,x,Q2);
                //dquark = Get_LHAPDF(2,x,Q2);
                //squark = Get_LHAPDF(3,x,Q2);
                //ubarquark = Get_LHAPDF(-1,x,Q2);
                //dbarquark = Get_LHAPDF(-2,x,Q2);
                //sbarquark = Get_LHAPDF(-3,x,Q2);
            //}
            //else 
            if(fOrder==1 || fOrder==2){/*{{{*/
                //Calculate medium modified PDFs:
                RunEPS09(fOrder, fErrSet, fA, fZ, x, Q2);

                //uquark = get_uA();
                //ubarquark = get_ubar();
                //dquark = get_dA();
                //dbarquark = get_dbar();
                //squark = get_s();
                //sbarquark = get_sbar();
            }
            else if(fOrder==3 || fOrder==4){
               //free PDF from CTEQ
                Run_CTEQPDF(x, Q2);         
                
                //uquark = Get_CTEQPDF(1,x,Q2);
                //ubarquark = Get_CTEQPDF(-1,x,Q2);
                //dquark = Get_CTEQPDF(2,x,Q2);
                //dbarquark = Get_CTEQPDF(-2,x,Q2);
                //squark = Get_CTEQPDF(3,x,Q2);
                //sbarquark = Get_CTEQPDF(-3,x,Q2);
            }else{
                cerr<<"***, in dxs(....), I don't know the type of fOrder"<<endl;
            }

            uquark = fuA;
            ubarquark = fubar;
            dquark = fdA;
            dbarquark = fdbar;
            squark = fs;
            sbarquark = fsbar;

            ////Test: see the different of nuclear-PDF and free-PDF in the same (Q2,x)
            //if(x>0.&&x<1.0){
            //double u1 = Get_LHAPDF(1,x,Q2);
            //double u2 = Get_CTEQPDF(1,x,Q2);
            ////outlog<<Form("%f   %f   %f   %f    %f",Q2, x, u1, u2, uquark)<<endl;
            ////cout<<Form("Q2=%f, x=%f, LHAPDF=%f, CTEQ=%f, EPS09=%f",Q2, x, u1, u2, uquark)<<endl;
            //}/*}}}*/

            //calculate F2p and F2n for inclusive XS calculations
            fF2p = pow(qu,2) * (uquark+ubarquark) + pow(qd,2)*(dquark+dbarquark) + pow(qs,2)*(squark+sbarquark); 
            fF2n = pow(qd,2) * (uquark+ubarquark) + pow(qu,2)*(dquark+dbarquark) + pow(qs,2)*(squark+sbarquark); //u-->d, d-->u, 

            double D_fav,D_unfav,D_s,D_g;
            if (fabs(particle_flag)==1) {
                //pion fragmentation functions
                Unpol_FF(1,z,Q2,&D_fav,&D_unfav,&D_s,&D_g);

                //proton
                df_p_hp = qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_unfav + qd*qd*dbarquark*D_fav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s;
                df_p_hm = qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_fav + qd*qd*dbarquark*D_unfav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s;

                //neutron
                //u->d, d->u, ubar->dbar, dbar->ubar
                df_n_hp = qu*qu*dquark*D_fav + qu*qu*dbarquark*D_unfav + qd*qd*uquark*D_unfav + qd*qd*ubarquark*D_fav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s;
                df_n_hm = qu*qu*dquark*D_unfav + qu*qu*dbarquark*D_fav + qd*qd*uquark*D_fav + qd*qd*ubarquark*D_unfav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s;
            }else{
                //kaon, fragmentation functions
                Unpol_FF(2,z,Q2,&D_fav,&D_unfav,&D_s,&D_g);

                //proton
                df_p_hp = qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_unfav + qs*qs*sbarquark*D_fav;
                df_p_hm = qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_fav + qs*qs*sbarquark*D_unfav;

                //neutron
                //u->d, d->u, ubar->dbar, dbar->ubar
                df_n_hp = qu*qu*dquark*D_fav   + qu*qu*dbarquark*D_unfav + qd*qd*uquark*D_s + qd*qd*ubarquark*D_s + qs*qs *squark * D_unfav + qs*qs*sbarquark*D_fav;
                df_n_hm = qu*qu*dquark*D_unfav + qu*qu*dbarquark*D_fav   + qd*qd*uquark*D_s + qd*qd*ubarquark*D_s + qs*qs *squark * D_fav   + qs*qs*sbarquark*D_unfav;
            }

            *dxs_hp = *dxs_hp*jacoF;
            *dxs_hm = *dxs_hm*jacoF;

            *(dxs_all+0) = *dxs_hp * df_p_hp/x;
            *(dxs_all+1) = *dxs_hp * df_n_hp/x;
            *(dxs_all+2) = *dxs_hm * df_p_hm/x;
            *(dxs_all+3) = *dxs_hm * df_n_hm/x;

            *dxs_hp *= (df_p_hp*fZ+df_n_hp*(fA-fZ))/x;
            *dxs_hm *= (df_p_hm*fZ+df_n_hm*(fA-fZ))/x;

        }
        /*}}}*/

        /*void dxs_hpt(double pt_tmp, double* dxs_hp,...){{{*/
        void dxs_hpt(double pt_tmp, double* dxs_hp, double* dxs_hm,double* dxs_all){

            double alpha_s = 0.0;
            //if(fOrder==0)
                //alpha_s = alphasPDF(sqrt(Q2));
            //else
                alpha_s = cteq_pdf_evolveas(fPDF, sqrt(Q2) );


            //cout << y << "\t" << Q2 << endl;
            //cout << alpha_s << endl;

            *dxs_hp = 1./137.035/137.035/16./PI/PI/Q2/Q2*y/4./PI*(197.3*197.3/100.*1.0e9)*alpha_s/1000000.;
            *dxs_hm = 1./137.035/137.035/16./PI/PI/Q2/Q2*y/4./PI*(197.3*197.3/100.*1.0e9)*alpha_s/1000000.;

            double qu=2./3.;
            double qd=-1./3.;
            double qs=-1./3.;

            double Pqq, Pqg, Pgq;
            double zp,xp;
            double df_p_hp = 0,df_p_hm = 0,df_n_hp = 0,df_n_hm = 0;

            double uquark=0.0,dquark=0.0,squark=0.0,ubarquark=0.0,dbarquark=0.0,sbarquark=0.0, gluon=0.0;
            double D_fav,D_unfav,D_s,D_g;
            // doing the integral
            for (Int_t i=0;i!=100;i++){
                xp = x + (1-x)/100.*(i+0.5);
                zp =  z / (z + xp*pt_tmp*pt_tmp/z/(1-xp)/Q2);
                if (z/zp>1) continue;

                Pqq = 64*PI/3.*Q2/y/y*((1+(1-y)*(1-y))*((1-xp)*(1-zp) + (1+xp*xp*zp*zp)/(1-xp)/(1-zp)) 
                        + 8*xp*zp*(1-y) //- 4 *sqrt(xp*zp*(1-y)/(1-xp)/(1-zp))*(2-y)*(xp*zp+(1-xp)*(1-zp))*cos(phi_h)
                        //+ 4*xp*zp*(1-y)*cos(2. * phi_h)
                        ) *(1-x)/100. /(xp*pt_tmp*pt_tmp+z*z*(1-xp)*Q2); //13

                Pqg = 64*PI/3.*Q2/y/y*( (1+(1-y)*(1-y))*((1-xp)*zp + (1+xp*xp*(1-zp)*(1-zp))/(1-xp)/zp) 
                        + 8*xp*(1-y)*(1-zp) 
                        //+ 4.*sqrt(xp*(1-y)*(1-zp)/(1-xp)/zp)*(2-y)*(xp*(1-zp)+(1-xp)*zp)*cos(phi_h)
                        //+ 4.*xp*(1-y)*(1-zp)*cos(2.*phi_h)
                        )*(1-x)/100. /(xp*pt_tmp*pt_tmp+z*z*(1-xp)*Q2);  //14

                Pgq = 8.*PI*Q2/y/y*((1+(1-y)*(1-y))*(xp*xp+(1-xp)*(1-xp))*(zp*zp+(1-zp)*(1-zp))/zp/(1-zp) 
                        + 16*xp*(1-xp)*(1-y)
                        //-4.*sqrt(xp*(1-xp)*(1-y)/zp/(1-zp))*(2-y)*(1-2*xp)*(1-2.*zp)*cos(phi_h)
                        //+ 8*xp*(1-xp)*(1-y)*cos(2.*phi_h)
                        )*(1-x)/100. /(xp*pt_tmp*pt_tmp+z*z*(1-xp)*Q2);  //15

                //cout << xp << "\t" << zp << "\t" << z/zp << "\t" << x/xp << "\t" << Pqq << "\t" << Pqg << "\t" << Pgq <<"\t" << (1-x)/100./(xp*pt_tmp*pt_tmp+z*z*(1-xp)*Q2)<<  endl;
                // Pqq = 1.;
                //     Pqg = 0.;
                //     Pgq = 0.;


                //if(fOrder==0){
                    ////free PDF from LHAPDF
                    //uquark = Get_LHAPDF(1,x/xp,Q2);
                    //dquark = Get_LHAPDF(2,x/xp,Q2);
                    //squark = Get_LHAPDF(3,x/xp,Q2);
                    //ubarquark = Get_LHAPDF(-1,x/xp,Q2);
                    //dbarquark = Get_LHAPDF(-2,x/xp,Q2);
                    //sbarquark = Get_LHAPDF(-3,x/xp,Q2);
                    //gluon = Get_LHAPDF(0,x/xp,Q2);
                //}
                //else
                if(fOrder==1 || fOrder==2){
                    //Calculate medium modified PDFs:
                    RunEPS09(fOrder, fErrSet, fA, fZ, x/xp, Q2);

                    //uquark = get_uA();
                    //dquark = get_dA();
                    //squark = get_s();
                    //ubarquark = get_ubar();
                    //dbarquark = get_dbar();
                    //sbarquark = get_sbar();
                    //gluon = get_g();
                }
                else if(fOrder==3 || fOrder==4){
                    //free PDF from CTEQ
                    Run_CTEQPDF(x, Q2);         

                    //uquark = Get_CTEQPDF(1,x/xp,Q2);
                    //dquark = Get_CTEQPDF(2,x/xp,Q2);
                    //squark = Get_CTEQPDF(3,x/xp,Q2);
                    //ubarquark = Get_CTEQPDF(-1,x/xp,Q2);
                    //dbarquark = Get_CTEQPDF(-2,x/xp,Q2);
                    //sbarquark = Get_CTEQPDF(-3,x/xp,Q2);
                    //gluon = Get_CTEQPDF(0,x/xp,Q2);
                }
                else{
                    cerr<<"***, in dxs(....), I don't know the type of fOrder"<<endl;
                }

                uquark = fuA;
                ubarquark = fubar;
                dquark = fdA;
                dbarquark = fdbar;
                squark = fs;
                sbarquark = fsbar;
                gluon = fg;

                if (fabs(particle_flag)==1){
                    //pion fragmentation function
                    Unpol_FF(1,z/zp,Q2,&D_fav,&D_unfav,&D_s,&D_g);
                    //proton
                    df_p_hp += Pqq*(qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_unfav + qd*qd*dbarquark*D_fav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_fav + qu*qu*D_unfav + qd*qd*D_unfav + qd*qd*D_fav + qs*qs*D_s + qs*qs*D_s);
                    df_p_hm += Pqq*(qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_fav + qd*qd*dbarquark*D_unfav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_fav + qd*qd*D_unfav + qs*qs*D_s + qs*qs*D_s);

                    //neutron
                    //u->d, d->u, ubar->dbar, dbar->ubar
                    df_n_hp += Pqq*(qu*qu*dquark*D_fav + qu*qu*dbarquark*D_unfav + qd*qd*uquark*D_unfav + qd*qd*ubarquark*D_fav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s) + Pqg*(qu*qu*dquark + qu*qu*dbarquark + qd*qd*uquark + qd*qd*ubarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_fav + qu*qu*D_unfav + qd*qd*D_unfav + qd*qd*D_fav + qs*qs*D_s + qs*qs*D_s);
                    df_n_hm += Pqq*(qu*qu*dquark*D_unfav + qu*qu*dbarquark*D_fav + qd*qd*uquark*D_fav + qd*qd*ubarquark*D_unfav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s) + Pqg*(qu*qu*dquark + qu*qu*dbarquark + qd*qd*uquark + qd*qd*ubarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_fav + qd*qd*D_unfav + qs*qs*D_s + qs*qs*D_s);
                }else{
                    //kaon fragmentation functions
                    Unpol_FF(2,z/zp,Q2,&D_fav,&D_unfav,&D_s,&D_g);

                    //proton
                    df_p_hp += Pqq*(qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_unfav + qs*qs*sbarquark*D_fav) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_s + qd*qd*D_s + qs*qs*D_fav + qs*qs*D_unfav);
                    df_p_hm += Pqq*(qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_fav + qs*qs*sbarquark*D_unfav) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_s + qd*qd*D_s + qs*qs*D_fav + qs*qs*D_unfav);

                    //neutron
                    df_n_hp += Pqq*(qu*qu*dquark*D_fav + qu*qu*dbarquark*D_unfav + qd*qd*uquark*D_s + qd*qd*ubarquark*D_s + qs*qs *squark * D_unfav + qs*qs*sbarquark*D_fav) + Pqg*(qu*qu*dquark + qu*qu*dbarquark + qd*qd*uquark + qd*qd*ubarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_s + qd*qd*D_s + qs*qs*D_fav + qs*qs*D_unfav);
                    df_n_hm += Pqq*(qu*qu*dquark*D_unfav + qu*qu*dbarquark*D_fav + qd*qd*uquark*D_s + qd*qd*ubarquark*D_s + qs*qs *squark * D_fav + qs*qs*sbarquark*D_unfav) + Pqg*(qu*qu*dquark + qu*qu*dbarquark + qd*qd*uquark + qd*qd*ubarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_s + qd*qd*D_s + qs*qs*D_fav + qs*qs*D_unfav);
                }
            }
            //cout << *dxs_hp << " \t" << df_hp/x << "\t" << jacoF << endl;

            *dxs_hp = *dxs_hp*jacoF;
            *dxs_hm = *dxs_hm*jacoF;

            *(dxs_all+0) = *dxs_hp * df_p_hp/x;
            *(dxs_all+1) = *dxs_hp * df_n_hp/x;
            *(dxs_all+2) = *dxs_hm * df_p_hm/x;
            *(dxs_all+3) = *dxs_hm * df_n_hm/x;

            *dxs_hp = *dxs_hp * (df_p_hp*fZ+df_n_hp*(fA-fZ))/x;
            *dxs_hm = *dxs_hm * (df_p_hm*fZ+df_n_hm*(fA-fZ))/x;
        }
        /*}}}*/

        /*double Unpol_FF(int flag, double z, double Q2, double* D_fav, double* D_unfav, double* D_s, double* D_g){{{*/
        double Unpol_FF(int flag, double z, double Q2, double* D_fav, double* D_unfav, double* D_s, double* D_g){
            //flag ==1 for pion, flag==2 for kaon
            double Q2zero=2.0;
            double lambda=0.227;

            if (z>1) z=1;

            double sv = log(log(Q2/lambda/lambda)/log(Q2zero/lambda/lambda));

            if (flag == 1){//pion
                double N = 1.150 - 1.522*sv + 1.378*pow(sv,2) - 0.527*pow(sv,3);
                double a1 = -0.740 - 1.680*sv + 1.546*pow(sv,2) - 0.596*pow(sv,3);
                double a2 = 1.430 + 0.543*sv - 0.023*pow(sv,2);
                double Ns = 4.250 - 3.147*sv + 0.755*pow(sv,2);
                double a1s = -0.770 -0.573*sv + 0.117*pow(sv,2);
                double a2s = 4.48 + 0.890*sv - 0.138*pow(sv,2);

                double Ng = 5.530 - 9.228*sv + 5.192*sv*sv - 0.966 * sv*sv*sv;
                double a1g = -0.320 + 0.318*sv - 0.561*sv*sv;
                double a2g = 2.7+2.553*sv-0.907*sv*sv;
                double a3g = 0.751*sv+0.496*sv*sv;

                double D_sum = N*pow(z,a1)*pow((1.0-z),a2);
                double D_sum_s = Ns*pow(z,a1s)*pow((1.0-z),a2s);
                double R_D = pow((1.0-z),0.083583)/pow((1.0+z),1.9838);
                *D_fav = D_sum/(1.0+R_D);
                *D_unfav = D_sum/(1.0+1.0/R_D);
                *D_s = D_sum_s/2.0;
                *D_g = Ng * pow(z,a1g) *pow(1-z,a2g)*(1+a3g/z)/2.;

            }else{//kaon
                double N = 0.310 - 0.038*sv - 0.042*pow(sv,2);
                double a1 = -0.980 - 0.260*sv + 0.008*pow(sv,2);
                double a2 = 0.970 + 0.978*sv - 0.229*pow(sv,2);
                double Ns =  1.080 - 0.469*sv + 0.003*pow(sv,2);
                double a1s = -0.820 -0.240*sv - 0.035*pow(sv,2);
                double a2s = 2.550 + 1.026*sv - 0.246*pow(sv,2);

                double Ng = 0.310 - 0.325*sv -0.092*sv*sv;
                double a1g = -0.17 - 0.214*sv-0.184*sv*sv;
                double a2g = 0.89 + 2.185*sv-0.471*sv*sv;
                double a3g = 1.154*sv-0.026*sv*sv;

                double D_sum = N*pow(z,a1)*pow((1.0-z),a2);
                double D_sum_s = Ns*pow(z,a1s)*pow((1.0-z),a2s);
                double R_D = pow((1.0-z),0.083583)/pow((1.0+z),1.9838);
                *D_fav = D_sum/(1.0+R_D);
                *D_unfav = D_sum/(1.0+1.0/R_D);
                *D_s = D_sum_s/2.0;
                *D_g = Ng * pow(z,a1g) *pow(1-z,a2g)*(1+a3g/z)/2.;
            }

            return 0;
        }
        /*}}}*/

        /*double Get_LHAPDF(int iparton,double x,double Q2){{{*/          
/*        double Get_LHAPDF(int iparton,double x,double Q2)          */
        //{
            //double Q = sqrt(Q2);
            //double result = 0;

            /////////////////in LHAPDF->xfx(): 
            //// 1->d, 2->u, 3->s, 4->c, 5->b, 6->t, 21->g 
            ////-1->dbar, -2->ubar,-3->sbar,-4->cbar,-5->bbar,-6->tbar
            ////For valance quarks:
            ////  dv-> d-dbar, u->u-ubar


            //if (iparton==1){
                //// u quark
                //result = xfx(x, Q, 2);
            //}else if (iparton==2){
                //// d quark
                //result = xfx(x, Q, 1);
            //}else if (iparton==-1){
                //// \bar{u} quark
                //result = xfx(x, Q, -2);
            //}else if (iparton==-2){
                //// \bar{d} quark
                //result = xfx(x, Q, -1);
            //}else if (iparton==3){
                //// strange quark
                //result = xfx(x, Q, 3);
            //}else if (iparton==-3){
                //// \bar{s} quark
                //result = xfx(x, Q, -3);
            //}else if (iparton==0){
                //result = xfx(x,Q,0);
            //}

            //return(result);
        /*}*/
        /*}}}*/

        /*double Get_CTEQPDF(int iparton,double x,double Q2){{{*/          
        double Get_CTEQPDF(int iprtn,double ix,double iQ2)          
        {
            double iQ = sqrt(iQ2);

            /*if(kParton=="g") iprtn = 0;*/
            //else if(kParton=="u") iprtn = 1;
            //else if(kParton=="d") iprtn = 2;
            //else if(kParton=="s") iprtn = 3;
            //else if(kParton=="c") iprtn = 4;
            //else if(kParton=="b") iprtn = 5;
            //else if(kParton=="ubar") iprtn =-1;
            //else if(kParton=="dbar") iprtn =-2;
            //else if(kParton=="sbar") iprtn =-3;
            //else if(kParton=="cbar") iprtn =-4;
            //else if(kParton=="bbar") iprtn =-5;
            //else{
            //cerr<<"***EORROR, in un_Get_CTEQPDF(...), unknown parton name = "<<kParton.Data()<<endl;
            /*}*/

            double result = cteq_pdf_evolvepdf(fPDF, iprtn, ix, iQ);
            //cout<<Form("for %d: x=%f, Q2=%f, PDF=%f", iprtn, ix, iQ, result)<<endl;
            return(result*ix);//return xf(x) instead of f(x)
        }
        /*}}}*/

        /*Kinematic Quantties{{{*/
    private:
        double mass_hadron;
        int	particle_flag;

        int fA;
        int fZ;
        cteq_pdf_t *fPDF;
        double fuA;
        double fdA;
        double fs;
        double fg;
        double fubar;
        double fdbar;
        double fsbar;
        double fF2p;
        double fF2n;

        TString fModel;
        int fOrder;
        int fErrSet;

    public:
        double E0_beam_ele;
        double E0_beam_ion;
        double mom_ele;
        double theta_ele; 
        double phi_ele; 
        double mom_had;
        double theta_had;
        double phi_had;

        double theta_q;
        double theta_s;
        double phi_q;
        double phi_s;
        double phi_h;

        double px_ele;
        double py_ele;
        double pz_ele;
        double E_ele;
        double px_had;
        double py_had;
        double pz_had;
        double E_had;

        double Q2;
        double W;
        double Wp;
        double x;
        double y;
        double z;
        double pt;
        double nu;
        double s;
        double rapidity;
        double gamma;
        double epsilon;
        double jacoF;
        double XS_Inclusive;
        double XS_HP;
        double XS_HM;
        double dilute_hp;
        double dilute_hm;

        double Asym_Cos2Phi_p_hp;
        double Asym_BM_p_hp;
        double Asym_Cahn_p_hp;
        double Asym_Cos2Phi_p_hm;
        double Asym_BM_p_hm;
        double Asym_Cahn_p_hm;

        double Asym_Cos2Phi_n_hp;
        double Asym_BM_n_hp;
        double Asym_Cahn_n_hp;
        double Asym_Cos2Phi_n_hm;
        double Asym_BM_n_hm;
        double Asym_Cahn_n_hm;

        double Asym_Cos2Phi_hp;
        double Asym_BM_hp;
        double Asym_Cahn_hp;
        double Asym_Cos2Phi_hm;
        double Asym_BM_hm;
        double Asym_Cahn_hm;
        /*}}}*/
};