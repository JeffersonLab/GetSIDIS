/////////////////////////////////////////////////////////////////////////////
//      Semi-Inclusive Deep Inelstic Scattering Cross Section Model        //
// Note:                                                                   //
//   This model is developed by Xin in his "collider" code. I extracted    //
//   the cross section parts and coded it into a C++ class which can be    //
//   easily embeded by other programs.                                     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
#include "Lsidis.h" //Tianbo Liu's SIDIS model

using namespace std;

const double DEG=180./3.1415926;
const double PI=3.1415926;
const double GeV2_to_nbarn = 0.3894 * 1e6; //GeV^2 to nbarn
class SIDIS
{
    public:
        SIDIS(TString kPDFSet, TString kFFSet){/*{{{*/
            fPDFSet= kPDFSet;
            fFFSet= kFFSet;
            fDebug=false;
        }/*}}}*/

        virtual ~SIDIS(){/*{{{*/
        
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
            cerr<<"      4, Other Output Quantities can be calculated after running CalcXS: "<<endl;
            cerr<<"           Q2, W, Wp, x, y, z, pt, nu, s, gamma, epsilon, jacoF "<<endl;
            cerr<<"===================================================================================================="<<endl;	

        }
        /*}}}*/

        void Init(double kMass_ion,int kA, int kZ,  int kPtcl_flag){/*{{{*/
            fA = kA;
            fZ = kZ;
            fMass_Ion=kMass_ion/fA;
            particle_flag = kPtcl_flag;
 
            double mass_pi = 0.13957;
            //double mass_pi0 = 0.1349766;
            double mass_kaon = 0.493667;
            if (abs(particle_flag) == 1){
                fMass_Hadron = mass_pi;
            }else if (abs(particle_flag) == 2){
                fMass_Hadron = mass_kaon;
            }else{
                cout << "particle_flag is wrong +-1 and +-2" << endl;
            }	

            sidis_hp.SetNucleus(fZ, fA-fZ);
            sidis_hm.SetNucleus(fZ, fA-fZ);
            if(particle_flag==1){//pion
                sidis_hp.SetHadron("pi+");
                sidis_hp.CheckHadron();
                sidis_hm.SetHadron("pi-");
                sidis_hm.CheckHadron();
            }
            if(particle_flag==2){//kaon
                sidis_hp.SetHadron("K+");
                sidis_hm.SetHadron("K-");
            }
           
            //Define PDF sets
<<<<<<< HEAD
            sidis_hp.SetPDFset(fPDFSet);
            sidis_hm.SetPDFset(fPDFSet);

            //Define Fragmentation Function sets
            sidis_hp.SetFFset(fFFSet);
            sidis_hm.SetFFset(fFFSet);
=======
            TString pdfSet="";
            const TString pdfEPPS="EPPS16nlo_CT14nlo_";
            const TString pdfNNPDF="nNNPDF20_nlo_as_0118_";
            const TString nnPDF = pdfNNPDF;
            if (fA==2 && fZ==1)
                pdfSet=nnPDF+"D2";
            else if (fA==4 && fZ==2)
                pdfSet=nnPDF+"He4";
            else if (fA==6 && fZ==3)
                pdfSet=nnPDF+"Li6";
            else if (fA==9 && fZ==4)
                pdfSet=nnPDF+"Be9";
            else if (fA==12 && fZ==6)
                pdfSet=nnPDF+"C12";
            else if (fA==14 && fZ==7)
                pdfSet=nnPDF+"N14";
            else if (fA==16 && fZ==8)
                pdfSet=nnPDF+"Li9";
            else if (fA==27 && fZ==13)
                pdfSet=nnPDF+"Al27";
            else if (fA==40 && fZ==20)
                pdfSet=nnPDF+"Ca40";
            else if (fA==56 && fZ==26)
                pdfSet=nnPDF+"Fe56";
            else if (fA==64 && fZ==29)
                pdfSet=nnPDF+"Cu64";
            else
                pdfSet="CJ15lo";

            sidis_hp.SetPDFset(pdfSet);
            sidis_hm.SetPDFset(pdfSet);

            //Define Fragmentation Function sets
            sidis_hp.SetFFset("DSSFFlo");
            sidis_hm.SetFFset("DSSFFlo");
>>>>>>> personal/master
        }/*}}}*/

        /*void SetKin( double mom_beam_ele, double mom_target,...){{{*/
        void SetKin(double kMom_beam_ele,double kMom_beam_ion,                        // GeV          GeV   
                double kP_ele, double kTh_ele, double kPh_ele,         // GeV/c         DEG            DEG 
                double kP_had, double kTh_had, double kPh_had         // GeV/c         DEG            DEG 
                ){

            /*Define{{{*/
            double mass_e = 0.511e-3; // electron mass in GeV
            double mass_p = 0.93827;
            //double mass_n = 0.939566;

            ////initialize
            ///////////////////////////////////////////////////////////////////////////////////
            //--- I have question here: -- Z. Ye 07/27
            //  1) to get x, use the ion_mass, proton-mass or the average mass target_mass/A?
            //  2) then what is the momentum? Is it Mom_ion or Mom_ion/A?
            //  I will use the average values temprately 
            ///////////////////////////////////////////////////////////////////////////////////
            //double mass_target = fMass_Ion/fA;
            double mass_target = mass_p;
            double mom_target = fabs(kMom_beam_ion)/fA;
            double mom_beam_ele = fabs(kMom_beam_ele);
            double energy_beam_ele = sqrt(kMom_beam_ele*kMom_beam_ele + mass_e*mass_e);
            double energy_target= sqrt(mom_target*mom_target + mass_target*mass_target);


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
            fPhysical = p3_q * p3_fin_had;

            fJacobF = Jacobian(fMom_ele,fTheta_ele,fPhi_ele,fMom_had,fTheta_had,fPhi_had,mom_target,energy_target,mom_beam_ele);
            /*}}}*/

            sidis_hp.SetInitialState(*P4_ini_ele, *P4_ini_ion);
            sidis_hp.SetFinalState(*P4_fin_ele, *P4_fin_had);
            sidis_hp.CalculateVariables();

            sidis_hm.SetInitialState(*P4_ini_ele, *P4_ini_ion);
            sidis_hm.SetFinalState(*P4_fin_ele, *P4_fin_had);
            sidis_hm.CalculateVariables();

            if(fDebug){
                TLorentzVector lp = sidis_hp.GetLorentzVector("lp");
                TLorentzVector Ph = sidis_hp.GetLorentzVector("Ph");
                cout<<Form("Electron: Px=%f (%f), Py=%f (%f), Pz=%f (%f), P=%f (%f)", fPx_ele,lp.X(), fPy_ele,lp.Y(), fPz_ele,lp.Z(),fMom_ele, lp.P())<<endl;
                cout<<Form("  Hadron: Px=%f (%f), Py=%f (%f), Pz=%f (%f), P=%f (%f)", fPx_had,Ph.X(), fPy_had,Ph.Y(), fPz_had,Ph.Z(),fMom_had, Ph.P())<<endl;
            }

            delete P4_ini_ele; delete P4_ini_ion;  
            delete P4_fin_had; delete P4_fin_ele; delete P4_q; 
            delete lrz_P4_q; delete lrz_P4_h; delete lrz_P4_ef;
        } 
        /*}}}*/

        /*int CalcXS(){{{*/
        int CalcXS(){
            /*SIDIS Model from Lsidis*/
            if(particle_flag==2){//kaon
                sidis_hp.ChangeTMDpars(0.604, 0.131);
                sidis_hm.ChangeTMDpars(0.604, 0.131);
            }  

            sidis_hp.CalculateFinalState();
            fXS_HP = sidis_hp.dsigma(0)*GeV2_to_nbarn*fJacobF;
            fIsPhy_HP = sidis_hp.CheckIsPhy();

            sidis_hm.CalculateFinalState();
            fXS_HM = sidis_hm.dsigma(0)*GeV2_to_nbarn*fJacobF;
            fIsPhy_HM = sidis_hm.CheckIsPhy();

            if(fIsPhy_HP && fIsPhy_HM)
                return 0;
            else
                return -1;
        } 
        /*}}}*/

        /*Return Values{{{*/
        double GetXS_HP(){
            return fXS_HP;
        }
        double GetXS_HM(){
            return fXS_HM;
        }
        Int_t IsPhy_HP(){
            return fIsPhy_HP;
        }
        Int_t IsPhy_HM(){
            return fIsPhy_HM;
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

        /*Kinematic Quantties{{{*/
    private:
        double fMass_Hadron;
        double fMass_Ion;
        int	particle_flag;
        int fA;
        int fZ;
        
        double fXS_HP;
        double fXS_HM;

        TString fPDFSet;
        TString fFFSet;

        Lsidis sidis_hp;
        Lsidis sidis_hm;
        Int_t fIsPhy_HP;
        Int_t fIsPhy_HM;

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
        double fPhysical;
        double fGamma;
        double fEpsilon;
        double fJacobF;
        bool fDebug;
        /*}}}*/
};
