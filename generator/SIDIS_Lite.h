/////////////////////////////////////////////////////////////////////////////
//      Semi-Inclusive Deep Inelstic Scattering Cross Section Model        //
// Note:                                                                   //
//   This model is developed by Xin in his "collider" code. I extracted    //
//   the cross section parts and coded it into a C++ class which can be    //
//   easily embeded by other programs.                                     //
//                                                                         //
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
const double GeV2_to_nbarn = 0.3894 * 1e6; //GeV^2 to nbarn
//ofstream outlog("pdf_check.dat");
class SIDIS
{
    public:
        SIDIS(TString kModel){/*{{{*/
            fModel = kModel;
            if(fModel=="EPS09"){
                //1->LO need CTEQ6L1,  2->NLO need CTEQ6.1M, 
                fOrder = 2;
                SetEPS09();
            }
            else if(fModel=="CTEQPDF"){
                //3->free L0 CTEQ6L1 PDF, 4->free NL0 CTEQ6.1M PDF
                fOrder = 4;
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

        /*void SetKin( double mom_beam_ele, double mom_target,...){{{*/
        void SetKin(double kMom_beam_ele,double kMom_beam_ion,                        // GeV          GeV   
                double kP_ele, double kTh_ele, double kPh_ele,         // GeV/c         DEG            DEG 
                double kP_had, double kTh_had, double kPh_had,         // GeV/c         DEG            DEG 
                double kMass_ion,int kA, int kZ,  int kPtcl_flag){

            /*Define{{{*/
            fA = kA;
            fZ = kZ;

            particle_flag = kPtcl_flag;
            double mass_e = 0.511e-3; // electron mass in GeV
            double mass_p = 0.93827;
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
            //double mass_target = kMass_ion/kA;
            double mass_target = mass_p;
            double mom_target = fabs(kMom_beam_ion)/kA;
            double mom_beam_ele = fabs(kMom_beam_ele);
            double energy_beam_ele = sqrt(kMom_beam_ele*kMom_beam_ele + mass_e*mass_e);
            double energy_target= sqrt(mom_target*mom_target + mass_target*mass_target);

            if (abs(particle_flag) == 1){
                fMass_Hadron = mass_pi;
            }else if (abs(particle_flag) == 2){
                fMass_Hadron = mass_kaon;
            }else{
                cout << "particle_flag is wrong +-1 and +-2" << endl;
            }	

            // //define the 4-momentum
            //define electron direction as +z assuming a proton/neutron for the ion mass now 
            // approximation for SIDIS
            TLorentzVector *P4_ini_ele = new TLorentzVector(0.,0.,mom_beam_ele, energy_beam_ele);
            TLorentzVector *P4_ini_ion = new TLorentzVector(0.,0.,-mom_target, energy_target);
            TLorentzVector *P4_fin_had = new TLorentzVector(0.,0.,0.,1.);
            TLorentzVector *P4_fin_ele = new TLorentzVector(0.,0.,0.,1.);
            TLorentzVector *P4_q = new TLorentzVector(0.,0.,0.,1.);

            TLorentzVector *lrz_P4_q = new TLorentzVector(0.,0.,0.,1.);
            TLorentzVector *lrz_P4_h = new TLorentzVector(0.,0.,0.,1.);
            TLorentzVector *lrz_P4_ef = new TLorentzVector(0.,0.,0.,1.);

            double vn=mom_target/energy_target;
            TVector3 vnboost(0.,0.,vn);

            TVector3 p3_q,p3_fin_ele,p3_fin_had;
            TVector3 p3_target_spin(0.,1.,0);
            /*}}}*/

            /*Kinematics Quantities{{{*/
            //For electron
            fPx_ele = kP_ele*sin(kTh_ele/DEG)*cos(kPh_ele/DEG);
            fPy_ele = kP_ele*sin(kTh_ele/DEG)*sin(kPh_ele/DEG);
            fPz_ele = kP_ele*cos(kTh_ele/DEG);
            fE_ele = sqrt(kP_ele*kP_ele+mass_e*mass_e);
            P4_fin_ele->SetPxPyPzE(kP_ele*sin(kTh_ele/DEG)*cos(kPh_ele/DEG),
                    kP_ele*sin(kTh_ele/DEG)*sin(kPh_ele/DEG),
                    kP_ele*cos(kTh_ele/DEG)
                    ,sqrt(kP_ele*kP_ele+mass_e*mass_e));

            //For hadron
            fPx_had = kP_had*sin(kTh_had/DEG)*cos(kPh_had/DEG);
            fPy_had = kP_had*sin(kTh_had/DEG)*sin(kPh_had/DEG);
            fPz_had = kP_had*cos(kTh_had/DEG);
            fE_had = sqrt(kP_had*kP_had+fMass_Hadron*fMass_Hadron);
            P4_fin_had->SetPxPyPzE(kP_had*sin(kTh_had/DEG)*cos(kPh_had/DEG),
                    kP_had*sin(kTh_had/DEG)*sin(kPh_had/DEG),
                    kP_had*cos(kTh_had/DEG)
                    ,sqrt(kP_had*kP_had+fMass_Hadron*fMass_Hadron));

            *P4_q = *P4_ini_ele - *P4_fin_ele;
            fQ2 = - (*P4_q)*(*P4_q);
            fW = (*P4_ini_ele + *P4_ini_ion - *P4_fin_ele)*(*P4_ini_ele + *P4_ini_ion - *P4_fin_ele);
            fWp = (*P4_ini_ele + *P4_ini_ion - *P4_fin_ele - *P4_fin_had)*(*P4_ini_ele + *P4_ini_ion - *P4_fin_ele - *P4_fin_had);

            fS = (*P4_ini_ele + *P4_ini_ion )*(*P4_ini_ele + *P4_ini_ion);
            fNu = (*P4_ini_ion) * (*P4_q)/mass_target;
            fXb = fQ2/(2 * (*P4_ini_ion)*(*P4_q));
            fZ_h = ((*P4_ini_ion)*(*P4_fin_had))/((*P4_ini_ion)*(*P4_q));
            fY = ((*P4_ini_ion)*(*P4_q))/((*P4_ini_ion)*(*P4_ini_ele));
            fW = sqrt(fW);
            fWp = sqrt(fWp);

            *lrz_P4_q = *P4_q;
            *lrz_P4_h = *P4_fin_had;
            *lrz_P4_ef = *P4_fin_ele;

            lrz_P4_ef->Boost(vnboost);
            lrz_P4_h->Boost(vnboost); 
            lrz_P4_q->Boost(vnboost);

            fGamma = 2*mass_target *fXb/sqrt(fQ2);
            fEpsilon = (1-fY-0.25*fGamma*fGamma*fY*fY)/(1-fY+0.5*fY*fY+0.25*fGamma*fGamma*fY*fY);

            for(Int_t j=0;j<3;j++){
                p3_q(j)=(*lrz_P4_q)(j);
                p3_fin_ele(j)=(*lrz_P4_ef)(j);
                p3_fin_had(j)=(*lrz_P4_h)(j);
            }

            fPt = p3_fin_had.Perp(p3_q);

            fTheta_q =p3_q.Theta();
            fPhi_q = p3_q.Phi();
            p3_target_spin.SetXYZ(0.,1.,0.);//transverse case only

            p3_fin_ele.RotateZ(-fPhi_q);
            p3_fin_ele.RotateY(-fTheta_q);
            p3_fin_had.RotateZ(-fPhi_q);
            p3_fin_had.RotateY(-fTheta_q);
            p3_target_spin.RotateZ(-fPhi_q);
            p3_target_spin.RotateY(-fTheta_q);

            fPhi_s = (Azimuthalphi(p3_target_spin(0),p3_target_spin(1))-Azimuthalphi(p3_fin_ele(0),p3_fin_ele(1)));
            fPhi_h = (Azimuthalphi(p3_fin_had(0),p3_fin_had(1))-Azimuthalphi(p3_fin_ele(0),p3_fin_ele(1)));
            fTheta_s = p3_target_spin.Theta();

            fMom_ele = P4_fin_ele->P();
            fTheta_ele = P4_fin_ele->Theta();
            fMom_had = P4_fin_had->P();
            fTheta_had = P4_fin_had->Theta();
            fPhi_ele = P4_fin_ele->Phi();
            fPhi_had = P4_fin_had->Phi();
            
            fRapidity = 0.5 * log((fMom_had + fPz_had)/(fMom_had - fPz_had)  );

            fJacobF = Jacobian(fMom_ele,fTheta_ele,fPhi_ele,fMom_had,fTheta_had,fPhi_had,mom_target,energy_target,mom_beam_ele);
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
            double kSinSQ = pow( sin(fTheta_ele*0.5),2);
            double kCosSQ = pow( cos(fTheta_ele*0.5),2);

            ////This definition is in the CM frame, i.e. dSigma/dx/dy
            //double xs_p = 2.0*PI*pow(1/137.036,2)*s/fQ2/fQ2*(1+pow(1-fY, 2)) * fF2p * GeV2_to_nbarn;
            //double xs_n = 2.0*PI*pow(1/137.036,2)*s/fQ2/fQ2*(1+pow(1-fY, 2)) * fF2n * GeV2_to_nbarn;

            ////This definition is in the lab frame, e.g. dSigma/dE'dOmega, where dE'dOmega = 2.*M*E0/E'*PI*fY*dxdy
            double xs_p= 4.0*pow(1/137.036,2) * fMom_ele*fMom_ele/fQ2/fQ2
                * (1./fNu *kCosSQ + 1/mass_p/fXb * kSinSQ) * fF2p * GeV2_to_nbarn;
            double xs_n= 4.0*pow(1/137.036,2) * fMom_ele*fMom_ele/fQ2/fQ2
                * (1./fNu *kCosSQ + 1/mass_n/fXb * kSinSQ) * fF2n * GeV2_to_nbarn;

            fXS_Inclusive = fZ * xs_p + (fA-fZ)*xs_n;//nbarn

            pt_tmp = fPt;
            if(fPt<0.8){
                //first method 	 
                bpt_p = 1./(0.2+fZ_h*fZ_h*0.25);// <pt^2> = 0.2 GeV^2 (quark internal momentum)
                ////Disable the TMD feature temparately, 02/01/2017
                //bpt_p = 1./(0.2);// <pt^2> = 0.2 GeV^2 (quark internal momentum)
                //bpt_m = bpt_p; // <kt^2> = 0.25 GeV^2 (struck quark internal momentum)

                //second method  use the original value
                // bpt_m=4.694; // HERMES parameterization
                // bpt_p=4.661;

                GetXS(pt_tmp, bpt_m, bpt_p, &dxs_hp, &dxs_hm,dxs_all[0]);

                fDilute_hp = dxs_all[0][0]/dxs_all[0][1];
                fDilute_hm = dxs_all[0][2]/dxs_all[0][3];

            }
            else{
                // this part is to generate factor K
                //make sure the DXS is the same at PT= 0.8 GeV
                //calculating the TMD part
                double K[2],dxs_temp[10];
                bpt_p = 1./(0.2+fZ_h*fZ_h*0.25);// <pt^2> = 0.2 GeV^2 (quark internal momentum)
                ////Disable the TMD feature temparately, 02/01/2017
                //bpt_p = 1./(0.2);// <pt^2> = 0.2 GeV^2 (quark internal momentum)
                //bpt_m = bpt_p; // <kt^2> = 0.25 GeV^2 (struck quark internal momentum)

                //second method  use the original value
                // bpt_m=4.694; // HERMES parameterization
                // bpt_p=4.661;

                pt_tmp = 0.8;
                GetXS(pt_tmp, bpt_m, bpt_p, &dxs_temp[0], &dxs_temp[1],dxs_all[0]);	   

                //calculating the TMD parts using the new vpt_p values
                bpt_p = 1./(0.25+fZ_h*fZ_h*0.28);// <pt^2> = 0.25 GeV^2 (quark internal momentum)
                ////Disable the TMD feature temparately, 02/01/2017
                //bpt_p = 1./(0.25);// <pt^2> = 0.25 GeV^2 (quark internal momentum)
                //bpt_m = bpt_p; // <kt^2> = 0.28 GeV^2 (struck quark internal momentum)

                //taking into account the NLO etc
                pt_tmp = 1.0;
                GetXS(pt_tmp, bpt_m, bpt_p, &dxs_temp[2], &dxs_temp[3],dxs_all[1]);	   
                GetXS_hPt(pt_tmp, &dxs_temp[4],&dxs_temp[5],dxs_all[2]);

                K[0] = (dxs_temp[0]-dxs_temp[2])/dxs_temp[4];
                K[1] = (dxs_temp[1]-dxs_temp[3])/dxs_temp[5];

                if (K[0]<0.) K[0] = 0.;
                if (K[1]<0.) K[1] = 0.;

                if (fPt>1.2){
                    pt_tmp = fPt;
                    GetXS(pt_tmp,bpt_m, bpt_p, &dxs_temp[2], &dxs_temp[3],dxs_all[1]);	   
                    GetXS_hPt(pt_tmp,&dxs_temp[4],&dxs_temp[5],dxs_all[2]);

                    dxs_hp = dxs_temp[2] + K[0]*dxs_temp[4];
                    dxs_hm = dxs_temp[3] + K[1]*dxs_temp[5];

                    fDilute_hp = (dxs_all[1][0] + K[0]*dxs_all[2][0])/(dxs_all[1][1] + K[0]*dxs_all[2][1]);
                    fDilute_hm = (dxs_all[1][2] + K[1]*dxs_all[2][2])/(dxs_all[1][3] + K[1]*dxs_all[2][3]);

                }else{
                    pt_tmp = 1.2;
                    GetXS(pt_tmp, bpt_m, bpt_p, &dxs_temp[2], &dxs_temp[3],dxs_all[1]);	   
                    GetXS_hPt(pt_tmp,&dxs_temp[4],&dxs_temp[5],dxs_all[2]);

                    dxs_temp[2] = dxs_temp[2] + K[0]*dxs_temp[4];
                    dxs_temp[3] = dxs_temp[3] + K[1]*dxs_temp[5];

                    dxs_temp[4] = (dxs_temp[0]-dxs_temp[2])/(1./0.8/0.8-1./1.2/1.2);
                    dxs_temp[5] = dxs_temp[0]-dxs_temp[4]/0.8/0.8;
                    dxs_temp[6] = (dxs_temp[1]-dxs_temp[3])/(1./0.8/0.8-1./1.2/1.2);
                    dxs_temp[7] = dxs_temp[1]-dxs_temp[6]/0.8/0.8;

                    dxs_hp = dxs_temp[4]/fPt/fPt + dxs_temp[5];
                    dxs_hm = dxs_temp[6]/fPt/fPt + dxs_temp[7];

                    fDilute_hp = ((dxs_all[1][0] + K[0]*dxs_all[2][0])
                            /(dxs_all[1][1] + K[0]*dxs_all[2][1])+dxs_all[0][0]/dxs_all[0][1])/2.;
                    fDilute_hm = ((dxs_all[1][2] + K[1]*dxs_all[2][2])
                            /(dxs_all[1][3] + K[1]*dxs_all[2][3])+dxs_all[0][2]/dxs_all[0][3])/2.;
                }
            }
            /*}}}*/

            /*try to take care of the decay{{{*/
            double decay_part = 0.0;
            if (abs(particle_flag)==1){
                if (fTheta_had>155./DEG){
                    decay_part = exp(7./cos(fTheta_had)*fMass_Hadron/(2.6*fMom_had*3.0));
                }else if (fTheta_had<=155./DEG&&fTheta_had>=140./DEG){
                    decay_part = exp(4.5/cos(fTheta_had)*fMass_Hadron/(2.6*fMom_had*3.0));
                }else if (fTheta_had <140/DEG){
                    decay_part = exp(-2.5/sin(fTheta_had)*fMass_Hadron/(2.6*fMom_had*3.0));
                }
            }else{
                if (fTheta_had>155./DEG){
                    decay_part = exp(7./cos(fTheta_had)*fMass_Hadron/(1.24*fMom_had*3.0));
                }else if (fTheta_had<=155./DEG&&fTheta_had>=140./DEG){
                    decay_part = exp(4.5/cos(fTheta_had)*fMass_Hadron/(1.24*fMom_had*3.0));
                }else if (fTheta_had <140/DEG){
                    decay_part = exp(-2.5/sin(fTheta_had)*fMass_Hadron/(1.24*fMom_had*3.0));
                }
            }
            /*}}}*/

            /*Save&Return{{{*/
            dxs_hp *= decay_part; 
            dxs_hm *= decay_part;

            //to avoid some wired behavior in log scale/*{{{*/
            int err_hp = 0, err_hm= 0;
            if((fXS_Inclusive)<1e-34) fXS_Inclusive=1e-34;
            if((dxs_hp)<1e-34) dxs_hp=1e-34;
            if((dxs_hm)<1e-34) dxs_hm=1e-34;
            if((fDilute_hp)<1e-34) fDilute_hp=1e-34;
            if((fDilute_hm)<1e-34) fDilute_hm=1e-34;

            if(isnan(dxs_hp)||isinf(dxs_hp)) {
                dxs_hp=1e-34;
                err_hp = 1;
            }
            if(isnan(dxs_hm)||isinf(dxs_hm)){
                dxs_hm=1e-34;
                err_hm = 1;
            }
            if(isnan(fXS_Inclusive)||isinf(fXS_Inclusive)) fXS_Inclusive=1e-34;
            if(isnan(fDilute_hp)||isinf(fDilute_hp)) fDilute_hp=1e-34;
            if(isnan(fDilute_hm)||isinf(fDilute_hm)) fDilute_hm=1e-34;
            /*}}}*/
            fXS_HP = dxs_hp;
            fXS_HM = dxs_hm;

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
            return fXS_Inclusive;
        }

        double GetXS_HP(){
            return fXS_HP;
        }
        double GetXS_HM(){
            return fXS_HM;
        }
        double GetDilute_HP(){
            return fDilute_hp;
        }
        double GetDilute_HM(){
            return fDilute_hm;
        }
        double get_uA(){ return fuA; }//return the average u
        double get_dA(){ return fdA; }//return the average d
        double get_s(){ return fs; }//return the average s
        double get_g(){ return fg; }//return the average g
        double get_ubar(){ return fubar; }//return the average u
        double get_dbar(){ return fdbar; }//return the average d
        double get_sbar(){ return fsbar; }//return the average s
        
        double get_Dfav(){ return fD_fav; }//return the favor Fragmentation Function
        double get_Dunfav(){ return fD_unfav; }//return the unfavor Fragmentation Function
        double get_Ds(){ return fD_s; }//return the sea Fragmentation Function
        double get_Dg(){ return fD_g; }//return the gluon Fragmentation Function
        /*}}}*/
        
        /*double RunCTEQPDF(double x,double fQ2){{{*/          
        void RunCTEQPDF(double ix,double iQ2)          
        {
            fubar = Get_CTEQPDF(-1, ix, iQ2);
            fdbar = Get_CTEQPDF(-2, ix, iQ2);
            fsbar = Get_CTEQPDF(-3, ix, iQ2);
            fuA = Get_CTEQPDF(1, ix, iQ2);
            fdA = Get_CTEQPDF(2, ix, iQ2);
            fs = Get_CTEQPDF(3, ix, iQ2);
            fg = Get_CTEQPDF(0, ix, iQ2);
        }
        /*}}}*/

        void RunEPS09(int iOrder, int iErrSet,int iA,int iZ, double ix, double iQ2){/*{{{*/

            double iR_uv = 0.0, iR_dv = 0.0, iR_u = 0.0, iR_d = 0.0, iR_s = 0.0, iR_c = 0.0, iR_b = 0.0, iR_g = 0.0; 
            //Order; //1->LO use CTEQ6L1, 2->NLO use CETQ6.1M
            //ErrSet;//1->central fit, 2,3-> err set#1, 4,5->err set#2, ...,30,31->err set#15
            if(iA>2)
                eps09(iOrder, iErrSet, iA, ix, sqrt(iQ2), iR_uv, iR_dv, iR_u, iR_d, iR_s, iR_c,iR_b,iR_g);
            else{
                iR_uv = 1.0; iR_dv = 1.0; iR_u = 1.0; iR_d = 1.0; iR_s = 1.0; iR_c = 1.0; iR_b = 1.0; iR_g = 1.0; 
            }

            double ubar = Get_CTEQPDF(-1, ix, iQ2);
            double dbar = Get_CTEQPDF(-2, ix, iQ2);
            double sbar = Get_CTEQPDF(-3, ix, iQ2);
            double u = Get_CTEQPDF(1, ix, iQ2);
            double d = Get_CTEQPDF(2, ix, iQ2);
            double s = Get_CTEQPDF(3, ix, iQ2);
            double g = Get_CTEQPDF(0, ix, iQ2);

            double uv = u-ubar;//uv = u - usea
            double dv = d-dbar;//dv = d - dsea
            
/*            fuA = (1.0*iZ)/(1.0*iA) * ( iR_uv * uv + iR_u * ubar)*/
                //+ (1.0*(iA-iZ))/(1.0*iA) * (iR_dv * dv + iR_d * dbar);

            //fdA = (1.0*iZ)/(1.0*iA) * ( iR_dv * dv + iR_d * dbar)
                //+ (1.0*(iA-iZ))/(1.0*iA) * (iR_uv * uv + iR_u * ubar);

            fuA = iR_uv * uv + iR_u * ubar;
            fdA = iR_dv * dv + iR_d * dbar;

            fubar =  iR_u * ubar;
            fdbar =  iR_d * dbar;

            fs =  iR_s * s;
            fsbar =  iR_s *sbar;

            fg =  iR_g * g;
        }/*}}}*/

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
        double Jacobian(double P_ef,double theta_ef,double phi_ef,
                        double P_hf,double theta_hf,double phi_hf,
                        double P_tgt,double E_tgt,double P_beam_ele){


            // need to shift the electron phi angle to 0 in order to do calculation for jacobian

            phi_hf -= phi_ef;
            phi_ef = 0.;

            double qbeta = -1.*P_tgt/E_tgt;
            //double qbeta = 0.;
            double qgamma = 1./sqrt(1-qbeta*qbeta);


            double mass_pi = fMass_Hadron;
            double P_hx,P_hy,P_hz;  
            P_hx = P_hf*sin(theta_hf)*cos(phi_hf);
            P_hy = P_hf*sin(theta_hf)*sin(phi_hf);
            P_hz = qgamma * (P_hf*cos(theta_hf) - qbeta * sqrt(P_hf*P_hf+mass_pi*mass_pi));



            //P_hz = P_hf*cos(theta_hf);

            double qz = qgamma*(P_beam_ele-P_ef*cos(theta_ef)-qbeta*(P_beam_ele-P_ef));

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

            double Px1 = P_hf*sin(theta_hf)*(-sin(phi_hf));
            double Py1 = P_hf*sin(theta_hf)*cos(phi_hf);
            //double Pz1 = 0.;

            double Px2 = -P_hf*cos(phi_hf)/tan(theta_hf);
            double Py2 = -P_hf*sin(phi_hf)/tan(theta_hf);
            double Pz2 = P_hf*qgamma;

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

            jaco9 = (P_beam_ele*(1.-cos(theta_ef))*(P_beam_ele*(E_tgt+P_tgt)-P_ef*(E_tgt+P_tgt*cos(theta_ef)))+P_ef*P_beam_ele*(1-cos(theta_ef))*(E_tgt+P_tgt*cos(theta_ef)))/(P_beam_ele*(E_tgt+P_tgt) - P_ef*(E_tgt+P_tgt*cos(theta_ef)))/(P_beam_ele*(E_tgt+P_tgt) - P_ef*(E_tgt+P_tgt*cos(theta_ef))); // px / pa

            jaco1 = (-P_ef*P_beam_ele*(P_beam_ele*(P_tgt+E_tgt)-P_ef*(E_tgt+P_tgt*cos(theta_ef)))+P_beam_ele*P_ef*(1-cos(theta_ef))*P_ef*P_tgt)/(P_beam_ele*(E_tgt+P_tgt) - P_ef*(E_tgt+P_tgt*cos(theta_ef)))/(P_beam_ele*(E_tgt+P_tgt) - P_ef*(E_tgt+P_tgt*cos(theta_ef))); // p x / pb



            jaco2 = - (E_tgt+P_tgt*cos(theta_ef))/(P_beam_ele*(E_tgt+P_tgt)); // py/pa

            jaco10 = (-1*P_ef*P_tgt)/(P_beam_ele*(E_tgt+P_tgt)); // py / pb

            jaco3 = (P_hf*E_tgt/sqrt(P_hf*P_hf + mass_pi*mass_pi) + P_tgt*cos(theta_hf))/(P_beam_ele*(E_tgt+P_tgt) - P_ef*(E_tgt+P_tgt*cos(theta_ef))); // p z / p d
            jaco11 = P_tgt*P_hf/(P_beam_ele*(E_tgt+P_tgt) - P_ef*(E_tgt+P_tgt*cos(theta_ef))); // p z/ p e

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

        /*void GetXS(double x, double y, double Q2, ...){{{*/
        //////////////////////////////////////////////////////////////////////////////////////////////////
        // This subroutine is called to compute new XS for any nuclear A
        // where the PDFs contain the free PDF from CTEQ and the nulcear-medium-modification from  EPS09
        // -- Zhihong Ye, 06/27/2016
        //////////////////////////////////////////////////////////////////////////////////////////////////
        void GetXS(double pt_tmp,double bpt_m, double bpt_p, double* kXS_HP, double* kXS_HM,double* kXS_All){
            double qu=2./3.;
            double qd=-1./3.;
            double qs=-1./3.;

            double dxs_hp_temp = 1./137.035/137.035/fXb/fY/fQ2 *fY*fY/(2*(1-fEpsilon))*(1+fGamma*fGamma/2./fXb)*GeV2_to_nbarn*fJacobF;
            double dxs_hm_temp = dxs_hp_temp;

            ////Temparetly set to zero, 02/10/2017
            //dxs_hp_temp *= fXb;
            //dxs_hm_temp *= fXb; 

            dxs_hp_temp *= bpt_p/PI*exp(-bpt_p*pt_tmp*pt_tmp)*fXb;
            dxs_hm_temp *= bpt_m/PI*exp(-bpt_m*pt_tmp*pt_tmp)*fXb;

            double df_p_hp=0,df_p_hm=0;
            double df_n_hp=0,df_n_hm=0;

            double uquark=0.0,dquark=0.0,squark=0.0,ubarquark=0.0,dbarquark=0.0,sbarquark=0.0;
            //if(fOrder==0){
                //uquark = Get_LHAPDF(1,fXb,fQ2);
                //dquark = Get_LHAPDF(2,fXb,fQ2);
                //squark = Get_LHAPDF(3,fXb,fQ2);
                //ubarquark = Get_LHAPDF(-1,fXb,fQ2);
                //dbarquark = Get_LHAPDF(-2,fXb,fQ2);
                //sbarquark = Get_LHAPDF(-3,fXb,fQ2);
            //}
            //else 
            if(fOrder==1 || fOrder==2){/*{{{*/
                //Calculate medium modified PDFs:
                RunEPS09(fOrder, fErrSet, fA, fZ, fXb, fQ2);
            }
            else if(fOrder==3 || fOrder==4){
               //free PDF from CTEQ
                RunCTEQPDF(fXb, fQ2);         
            }else{
                cerr<<"***, in GetXS(....), I don't know the type of fOrder"<<endl;
            }

            uquark = fuA;
            ubarquark = fubar;
            dquark = fdA;
            dbarquark = fdbar;
            squark = fs;
            sbarquark = fsbar;

            ////Test: see the different of nuclear-PDF and free-PDF in the same (fQ2,fXb)
            //if(fXb>0.&&fXb<1.0){
            //double u1 = Get_LHAPDF(1,fXb,fQ2);
            //double u2 = Get_CTEQPDF(1,fXb,fQ2);
            ////outlog<<Form("%f   %f   %f   %f    %f",fQ2, fXb, u1, u2, uquark)<<endl;
            ////cout<<Form("fQ2=%f, fXb=%f, LHAPDF=%f, CTEQ=%f, EPS09=%f",fQ2, fXb, u1, u2, uquark)<<endl;
            //}/*}}}*/

            //calculate F2p and F2n for inclusive XS calculations
            fF2p = pow(qu,2) * (uquark+ubarquark) + pow(qd,2)*(dquark+dbarquark) + pow(qs,2)*(squark+sbarquark); 
            fF2n = pow(qd,2) * (uquark+ubarquark) + pow(qu,2)*(dquark+dbarquark) + pow(qs,2)*(squark+sbarquark); //u-->d, d-->u, 

            double D_fav,D_unfav,D_s,D_g;
            if (fabs(particle_flag)==1) {
                //pion fragmentation functions
                GetFF_Unpol(1,fZ_h,fQ2,&D_fav,&D_unfav,&D_s,&D_g);

                //proton
                df_p_hp = qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_unfav + qd*qd*dbarquark*D_fav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s;
                df_p_hm = qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_fav + qd*qd*dbarquark*D_unfav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s;

                //neutron
                //u->d, d->u, ubar->dbar, dbar->ubar
                df_n_hp = qu*qu*dquark*D_fav + qu*qu*dbarquark*D_unfav + qd*qd*uquark*D_unfav + qd*qd*ubarquark*D_fav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s;
                df_n_hm = qu*qu*dquark*D_unfav + qu*qu*dbarquark*D_fav + qd*qd*uquark*D_fav + qd*qd*ubarquark*D_unfav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s;
            }else{
                //kaon, fragmentation functions
                GetFF_Unpol(2,fZ_h,fQ2,&D_fav,&D_unfav,&D_s,&D_g);

                //proton
                df_p_hp = qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_unfav + qs*qs*sbarquark*D_fav;
                df_p_hm = qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_fav + qs*qs*sbarquark*D_unfav;

                //neutron
                //u->d, d->u, ubar->dbar, dbar->ubar
                df_n_hp = qu*qu*dquark*D_fav   + qu*qu*dbarquark*D_unfav + qd*qd*uquark*D_s + qd*qd*ubarquark*D_s + qs*qs *squark * D_unfav + qs*qs*sbarquark*D_fav;
                df_n_hm = qu*qu*dquark*D_unfav + qu*qu*dbarquark*D_fav   + qd*qd*uquark*D_s + qd*qd*ubarquark*D_s + qs*qs *squark * D_fav   + qs*qs*sbarquark*D_unfav;
            }
            fD_unfav = D_unfav;
            fD_fav = D_fav;
            fD_s = D_s;
            fD_g = D_g;

            *(kXS_All+0) = dxs_hp_temp * df_p_hp/fXb;
            *(kXS_All+1) = dxs_hp_temp * df_n_hp/fXb;
            *(kXS_All+2) = dxs_hm_temp * df_p_hm/fXb;
            *(kXS_All+3) = dxs_hm_temp * df_n_hm/fXb;

            *kXS_HP = dxs_hp_temp * (df_p_hp*fZ+df_n_hp*(fA-fZ))/fXb;
            *kXS_HM = dxs_hm_temp * (df_p_hm*fZ+df_n_hm*(fA-fZ))/fXb;
        }
        /*}}}*/

        /*void GetXS_hPt(double pt_tmp, double* dxs_hp,...){{{*/
        void GetXS_hPt(double pt_tmp, double* kXS_HP, double* kXS_HM,double* kXS_All){

            double alpha_s = 0.0;
            //if(fOrder==0)
                //alpha_s = alphasPDF(sqrt(fQ2));
            //else
                alpha_s = cteq_pdf_evolveas(fPDF, sqrt(fQ2) );


            //cout << fY << "\t" << fQ2 << endl;
            //cout << alpha_s << endl;

            double dxs_hp_temp = 1./137.035/137.035/16./PI/PI/fQ2/fQ2*fY/4./PI*(197.3*197.3/100.*1.0e9)*alpha_s/1000000.*fJacobF;
            double dxs_hm_temp = dxs_hp_temp;

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
                xp = fXb + (1-fXb)/100.*(i+0.5);
                zp =  fZ_h / (fZ_h + xp*pt_tmp*pt_tmp/fZ_h/(1-xp)/fQ2);
                if (fZ_h/zp>1) continue;

                Pqq = 64*PI/3.*fQ2/fY/fY*((1+(1-fY)*(1-fY))*((1-xp)*(1-zp) + (1+xp*xp*zp*zp)/(1-xp)/(1-zp)) 
                        + 8*xp*zp*(1-fY) //- 4 *sqrt(xp*zp*(1-fY)/(1-xp)/(1-zp))*(2-fY)*(xp*zp+(1-xp)*(1-zp))*cos(fPhi_h)
                        //+ 4*xp*zp*(1-fY)*cos(2. * fPhi_h)
                        ) *(1-fXb)/100. /(xp*pt_tmp*pt_tmp+fZ_h*fZ_h*(1-xp)*fQ2); //13

                Pqg = 64*PI/3.*fQ2/fY/fY*( (1+(1-fY)*(1-fY))*((1-xp)*zp + (1+xp*xp*(1-zp)*(1-zp))/(1-xp)/zp) 
                        + 8*xp*(1-fY)*(1-zp) 
                        //+ 4.*sqrt(xp*(1-fY)*(1-zp)/(1-xp)/zp)*(2-fY)*(xp*(1-zp)+(1-xp)*zp)*cos(fPhi_h)
                        //+ 4.*xp*(1-fY)*(1-zp)*cos(2.*fPhi_h)
                        )*(1-fXb)/100. /(xp*pt_tmp*pt_tmp+fZ_h*fZ_h*(1-xp)*fQ2);  //14

                Pgq = 8.*PI*fQ2/fY/fY*((1+(1-fY)*(1-fY))*(xp*xp+(1-xp)*(1-xp))*(zp*zp+(1-zp)*(1-zp))/zp/(1-zp) 
                        + 16*xp*(1-xp)*(1-fY)
                        //-4.*sqrt(xp*(1-xp)*(1-fY)/zp/(1-zp))*(2-fY)*(1-2*xp)*(1-2.*zp)*cos(fPhi_h)
                        //+ 8*xp*(1-xp)*(1-fY)*cos(2.*fPhi_h)
                        )*(1-fXb)/100. /(xp*pt_tmp*pt_tmp+fZ_h*fZ_h*(1-xp)*fQ2);  //15

                //cout << xp << "\t" << zp << "\t" << fZ_h/zp << "\t" << fXb/xp << "\t" << Pqq << "\t" << Pqg << "\t" << Pgq <<"\t" << (1-fXb)/100./(xp*pt_tmp*pt_tmp+fZ_h*fZ_h*(1-xp)*fQ2)<<  endl;
                // Pqq = 1.;
                //     Pqg = 0.;
                //     Pgq = 0.;


                //if(fOrder==0){
                    ////free PDF from LHAPDF
                    //uquark = Get_LHAPDF(1,fXb/xp,fQ2);
                    //dquark = Get_LHAPDF(2,fXb/xp,fQ2);
                    //squark = Get_LHAPDF(3,fXb/xp,fQ2);
                    //ubarquark = Get_LHAPDF(-1,fXb/xp,fQ2);
                    //dbarquark = Get_LHAPDF(-2,fXb/xp,fQ2);
                    //sbarquark = Get_LHAPDF(-3,fXb/xp,fQ2);
                    //gluon = Get_LHAPDF(0,fXb/xp,fQ2);
                //}
                //else
                if(fOrder==1 || fOrder==2){
                    //Don't call RunEPS09DF(fXb, fQ2) again since now we are only dealing with hPt corrections         
                    uquark = Get_EPS09(fOrder, fErrSet, fA, fZ, 1, fXb, fQ2);
                    dquark = Get_EPS09(fOrder, fErrSet, fA, fZ, 2, fXb, fQ2);
                    squark = Get_EPS09(fOrder, fErrSet, fA, fZ, 3, fXb, fQ2);
                    ubarquark = Get_EPS09(fOrder, fErrSet, fA, fZ,-1, fXb, fQ2);
                    dbarquark = Get_EPS09(fOrder, fErrSet, fA, fZ,-2, fXb, fQ2);
                    sbarquark = Get_EPS09(fOrder, fErrSet, fA, fZ,-3, fXb, fQ2);
                    gluon = Get_EPS09(fOrder, fErrSet, fA, fZ, 0, fXb, fQ2);
                }
                else if(fOrder==3 || fOrder==4){
                    //free PDF from CTEQ
                    //Don't call RunCTEQPDF(fXb, fQ2) again since now we are only dealing with hPt corrections         
                    uquark = Get_CTEQPDF(1,fXb/xp,fQ2);
                    dquark = Get_CTEQPDF(2,fXb/xp,fQ2);
                    squark = Get_CTEQPDF(3,fXb/xp,fQ2);
                    ubarquark = Get_CTEQPDF(-1,fXb/xp,fQ2);
                    dbarquark = Get_CTEQPDF(-2,fXb/xp,fQ2);
                    sbarquark = Get_CTEQPDF(-3,fXb/xp,fQ2);
                    gluon = Get_CTEQPDF(0,fXb/xp,fQ2);
                }
                else{
                    cerr<<"***, in GetXS(....), I don't know the type of fOrder"<<endl;
                }

                if (fabs(particle_flag)==1){
                    //pion fragmentation function
                    GetFF_Unpol(1,fZ_h/zp,fQ2,&D_fav,&D_unfav,&D_s,&D_g);
                    //proton
                    df_p_hp += Pqq*(qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_unfav + qd*qd*dbarquark*D_fav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_fav + qu*qu*D_unfav + qd*qd*D_unfav + qd*qd*D_fav + qs*qs*D_s + qs*qs*D_s);
                    df_p_hm += Pqq*(qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_fav + qd*qd*dbarquark*D_unfav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_fav + qd*qd*D_unfav + qs*qs*D_s + qs*qs*D_s);

                    //neutron
                    //u->d, d->u, ubar->dbar, dbar->ubar
                    df_n_hp += Pqq*(qu*qu*dquark*D_fav + qu*qu*dbarquark*D_unfav + qd*qd*uquark*D_unfav + qd*qd*ubarquark*D_fav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s) + Pqg*(qu*qu*dquark + qu*qu*dbarquark + qd*qd*uquark + qd*qd*ubarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_fav + qu*qu*D_unfav + qd*qd*D_unfav + qd*qd*D_fav + qs*qs*D_s + qs*qs*D_s);
                    df_n_hm += Pqq*(qu*qu*dquark*D_unfav + qu*qu*dbarquark*D_fav + qd*qd*uquark*D_fav + qd*qd*ubarquark*D_unfav + qs*qs*squark*D_s + qs*qs*sbarquark*D_s) + Pqg*(qu*qu*dquark + qu*qu*dbarquark + qd*qd*uquark + qd*qd*ubarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_fav + qd*qd*D_unfav + qs*qs*D_s + qs*qs*D_s);
                }else{
                    //kaon fragmentation functions
                    GetFF_Unpol(2,fZ_h/zp,fQ2,&D_fav,&D_unfav,&D_s,&D_g);

                    //proton
                    df_p_hp += Pqq*(qu*qu*uquark*D_fav + qu*qu*ubarquark*D_unfav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_unfav + qs*qs*sbarquark*D_fav) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_s + qd*qd*D_s + qs*qs*D_fav + qs*qs*D_unfav);
                    df_p_hm += Pqq*(qu*qu*uquark*D_unfav + qu*qu*ubarquark*D_fav + qd*qd*dquark*D_s + qd*qd*dbarquark*D_s + qs*qs *squark * D_fav + qs*qs*sbarquark*D_unfav) + Pqg*(qu*qu*uquark + qu*qu*ubarquark + qd*qd*dquark + qd*qd*dbarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_s + qd*qd*D_s + qs*qs*D_fav + qs*qs*D_unfav);

                    //neutron
                    df_n_hp += Pqq*(qu*qu*dquark*D_fav + qu*qu*dbarquark*D_unfav + qd*qd*uquark*D_s + qd*qd*ubarquark*D_s + qs*qs *squark * D_unfav + qs*qs*sbarquark*D_fav) + Pqg*(qu*qu*dquark + qu*qu*dbarquark + qd*qd*uquark + qd*qd*ubarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_s + qd*qd*D_s + qs*qs*D_fav + qs*qs*D_unfav);
                    df_n_hm += Pqq*(qu*qu*dquark*D_unfav + qu*qu*dbarquark*D_fav + qd*qd*uquark*D_s + qd*qd*ubarquark*D_s + qs*qs *squark * D_fav + qs*qs*sbarquark*D_unfav) + Pqg*(qu*qu*dquark + qu*qu*dbarquark + qd*qd*uquark + qd*qd*ubarquark + qs*qs*squark + qs*qs*sbarquark)*D_g + Pgq*gluon*(qu*qu*D_unfav + qu*qu*D_fav + qd*qd*D_s + qd*qd*D_s + qs*qs*D_fav + qs*qs*D_unfav);
                }
            }
            //cout << *dxs_hp << " \t" << df_hp/fXb << "\t" << fJacobF << endl;

            *(kXS_All+0) = dxs_hp_temp * df_p_hp/fXb;
            *(kXS_All+1) = dxs_hp_temp * df_n_hp/fXb;
            *(kXS_All+2) = dxs_hm_temp * df_p_hm/fXb;
            *(kXS_All+3) = dxs_hm_temp * df_n_hm/fXb;

            *kXS_HP = dxs_hp_temp * (df_p_hp*fZ+df_n_hp*(fA-fZ))/fXb;
            *kXS_HM = dxs_hm_temp * (df_p_hm*fZ+df_n_hm*(fA-fZ))/fXb;
        }
        /*}}}*/

        /*double GetFF_Unpol(int flag, double z, double Q2, double* D_fav, double* D_unfav, double* D_s, double* D_g){{{*/
        double GetFF_Unpol(int flag, double z, double Q2, double* D_fav, double* D_unfav, double* D_s, double* D_g){
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

        /*double Get_LHAPDF(int iparton,double ix,double iQ2){{{*/          
/*        double Get_LHAPDF(int iparton,double ix,double iQ2)          */
        //{
            //double iQ = sqrt(iQ2);
            //double result = 0;

            /////////////////in LHAPDF->xfx(): 
            //// 1->d, 2->u, 3->s, 4->c, 5->b, 6->t, 21->g 
            ////-1->dbar, -2->ubar,-3->sbar,-4->cbar,-5->bbar,-6->tbar
            ////For valance quarks:
            ////  dv-> d-dbar, u->u-ubar


            //if (iparton==1){
                //// u quark
                //result = xfx(ix, iQ, 2);
            //}else if (iparton==2){
                //// d quark
                //result = xfx(ix, iQ, 1);
            //}else if (iparton==-1){
                //// \bar{u} quark
                //result = xfx(ix, iQ, -2);
            //}else if (iparton==-2){
                //// \bar{d} quark
                //result = xfx(ix, iQ, -1);
            //}else if (iparton==3){
                //// strange quark
                //result = xfx(ix, iQ, 3);
            //}else if (iparton==-3){
                //// \bar{s} quark
                //result = xfx(ix, iQ, -3);
            //}else if (iparton==0){
                //result = xfx(ix,iQ,0);
            //}

            //return(result);
        /*}*/
        /*}}}*/

        /*double Get_CTEQPDF(int iparton,double ix,double iQ2){{{*/          
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

        double Get_EPS09(int iOrder, int iErrSet,int iA,int iZ, int iParton, double ix, double iQ2){/*{{{*/

            double iR_uv = 0.0, iR_dv = 0.0, iR_u = 0.0, iR_d = 0.0, iR_s = 0.0, iR_c = 0.0, iR_b = 0.0, iR_g = 0.0; 
            //Order; //1->LO use CTEQ6L1, 2->NLO use CETQ6.1M
            //ErrSet;//1->central fit, 2,3-> err set#1, 4,5->err set#2, ...,30,31->err set#15
            if(iA>2)
                eps09(iOrder, iErrSet, iA, ix, sqrt(iQ2), iR_uv, iR_dv, iR_u, iR_d, iR_s, iR_c,iR_b,iR_g);
            else{
                iR_uv = 1.0; iR_dv = 1.0; iR_u = 1.0; iR_d = 1.0; iR_s = 1.0; iR_c = 1.0; iR_b = 1.0; iR_g = 1.0; 
            }

            if(iParton==1){
                double u = Get_CTEQPDF(1, ix, iQ2);
                double ubar = Get_CTEQPDF(-1, ix, iQ2);
                double uv = u-ubar;//uv = u - usea
                double uA = iR_uv * uv + iR_u * ubar;
                return uA;
            }
            if(iParton==2){
                double d = Get_CTEQPDF(2, ix, iQ2);
                double dbar = Get_CTEQPDF(-2, ix, iQ2);
                double dv = d-dbar;//dv = d - dsea
                double dA = iR_dv * dv + iR_d * dbar;
                return dA;
            }
            if(iParton==3){
                double s = Get_CTEQPDF(3, ix, iQ2);
                return iR_s * s;
            }
            if(iParton==-1){
                double ubar = Get_CTEQPDF(-1, ix, iQ2);
                return iR_u * ubar;
            }
            if(iParton==-2){
                double dbar = Get_CTEQPDF(-2, ix, iQ2);
                return iR_d * dbar;
            }
            if(iParton==-3){
                double sbar = Get_CTEQPDF(-3, ix, iQ2);
                return iR_s * sbar;
            }
            if(iParton==0){
                double g = Get_CTEQPDF(0, ix, iQ2);
                return iR_g * g;
            }
        }/*}}}*/

        /*Kinematic Quantties{{{*/
    private:
        double fMass_Hadron;
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
        double fD_unfav;
        double fD_fav;
        double fD_s;
        double fD_g;

        double fXS_Inclusive;
        double fXS_HP;
        double fXS_HM;
        double fDilute_hp;
        double fDilute_hm;

        TString fModel;
        int fOrder;
        int fErrSet;

    public:
        double fMom_ele;
        double fTheta_ele; 
        double fPhi_ele; 
        double fMom_had;
        double fTheta_had;
        double fPhi_had;

        double fTheta_q;
        double fTheta_s;
        double fPhi_q;
        double fPhi_s;
        double fPhi_h;

        double fPx_ele;
        double fPy_ele;
        double fPz_ele;
        double fE_ele;
        double fPx_had;
        double fPy_had;
        double fPz_had;
        double fE_had;

        double fQ2;
        double fW;
        double fWp;
        double fXb;
        double fY;
        double fZ_h;
        double fPt;
        double fNu;
        double fS;
        double fRapidity;
        double fGamma;
        double fEpsilon;
        double fJacobF;
        /*}}}*/
};
