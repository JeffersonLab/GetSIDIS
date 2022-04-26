#ifndef __CLAS12_FASTMC__
#define __CLAS12_FASTMC__

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
#include <TH1D.h>
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
#include <TRandom3.h>
/*}}}*/

using namespace std;

#include "constants.h"
#include "conf.h"

#define R2D 57.29578
#define D2R 
//all functions here//
int accep_func(int Id, double P, double Theta, double Phi0, double *out);
void smear_fun(double Beta, double P, double Theta, double Phi0, double *out);
void Charge_Type(int Id, double *out);
double findBdL(double P, double E, double Theta);
void EC_XYZ_UVW(double x, double y, double z, double *out);
void type(double Charge, double Theta, int *out);
int clas12(int Id, double P, double Theta, double Phi0, double *out);

int    InBend, OutBend, Neutral, Fwd, Central;

//###########################################################################################################################################################################################################
int accep_func(int Id, double P, double Theta, double Phi0, double *out)/*{{{*/
{
    fCnfg_num = 0 ;
    fTTmax = 3375.;        fTor_Cur =  2250.;      fThf_lim = 40.;
    fThacc1= 5.;           fThacc2= 17.;           fThacc3=35.;            fThacc4= 0.25 ;
    fThacc5= 0.2;          fThacc6=1.1;            fThacc01=5.;
    fThmom1= 0.6;          fThmom2= 4.;            fThmom3=45.;            fThmom4= 6.;
    fPhimom1=25.;          fPhimom2= 5.;           fPhimom3=0.15;          fPhimom4=0.089;         fPhimom5= 0.0416667;
    fFoMT1S1=0.20678E-3;   fFoMT1S2= -0.57466E-5 ; fFoMT1S3=0.38913E-6;    fFoMT2S1=0.26652E-2;    fFoMT2S2=-0.80631E-4;   fFoMT2S3= 0.87599E-5;
    fFoInMoS1= 0.005;      fFoInMoS2= 0.0 ;        fFoOuMoS1=0.007 ;       fFoOuMoS2=0.0 ;         fFoNeMoS1=0.050;        fFoNeMoS2= 0.0;
    fFoInThS1=3.684E-05;   fFoInThS2= 0.00162613;  fFoOuThS1= 0.0286;      fFoOuThS2= 0.00;        fFoNeThS1=0.0286;       fFoNeThS2=0.00;
    fFoA0FiS1=0.3125E-2;   fFoA1FiS1= 0.262E-4;    fFoA2FiS1= -0.350E-6 ;  fFoA0FiS2=0.542E-3;     fFoA1FiS2=-0.531E-5;    fFoA2FiS2=0.886E-7;
    fFoNeFiS1=0.0573;      fFoNeFiS2= 0.0;         fCeInMoS1= 0.065;       fCeInMoS2= 0.0138;      fCeOuMoS1=0.065;        fCeOuMoS2= 0.0138 ;
    fCeNeMoS1= 0.050 ;     fCeNeMoS2=0.0;          fCeA0ThS2= 0.4236E-3;   fCeA1ThS2= -0.4451E-5;  fCeA2ThS2= 0.247E-7 ;   fCeInThS1= 0.8776E-2;
    fCeInThS2= 0.00000;    fCeOuThS1= 0.0286;      fCeOuThS2=0.00000;      fCeNeThS1=0.0286;       fCeNeThS2=0.00000;      fCeInFiS1=0.0069;



    double cLight(0.299792458), Gap(2.0);
    double Theta_min_cut(5.0), E_el_min(0.5), P_pi_min(0.2);
    double Theta_IC_min(2.5), Theta_IC_max(4.5), coil_CALO_width(0.0);	// Chgt de Theta_IC_max de 4.5 à 6.5deg
    double EC_r(510.32), EC_Th(0.436332);
    double ylow(-182.974), yhi(189.956), tgrho(1.95325), dl(100.), rc0(120.);

    double weight(0.), dPhi(0.);
    double Phis(0.), Pnorm(0.), Thcut(0.), xcut(0.);
    double Charge(0.), Mass(0.), E(0.), BdL(0.), defl_Phi(0.);
    double sinth(0.), costh(0.), sinph(0.), cosph(0.);
    double xh(0.), yh(0.), zh(0.), uh(0.), vh(0.), wh(0.), xi(0.), yi(0.);
    double rcoil(0.), hcoil(0.), dcoil(0.);
    double r(0.), Th(0.), expon(0.); 


    //Initialization
    out[0] = 0.; out[1] = 0.;
    double out1[2] = {0,0};
    double out2[5]; for(int i=0; i<5; i++) out2[i] = 0.;
    bool   Pcutl;

    Charge_Type(Id, out1);
    Charge = out1[0];
    Mass   = out1[1];
    E      = sqrt(P*P + Mass*Mass);
    if(Charge == -999 || Mass == -999) 
    {
        return -1;
    }

    //Take into account phi deflection due to solenoid magnet
    if(TMath::Abs(P*sin(Theta*TORAD)) > 0)
    {
        BdL      = findBdL(P, E, Theta*TORAD);
        defl_Phi = TODEG*cLight*BdL/(P*sin(Theta*TORAD)); // BAPTISTE a enlever pour virer la deflection qui impactt la distribution phi Vs Theta des electrons
        //defl_Phi = 0.0;
    }
    else
    {
        defl_Phi = 0.0;
    }

    double Phi00  = Phi0 + TMath::Sign(1.,Charge)*defl_Phi;
    double dPhi00 = Phi00 - (int)Phi00;
    if(Charge != 0)
    {
        Phis = (((int)Phi00 + 30 + 60000)% 60) - 30 + dPhi00;
        //cout<<"Phi "<<Phi0<<"  "<<Phi00<<"  "<<Phis<<endl;
    }
    else
    {
        Phis = Phi0;
    }


    //since our acceptance functions care about sagitta and not momentum, per 
    //sector, we will form a "normalized momentum", pnorm
    //cout<<"fTor_Cur="<<fTor_Cur<<endl;
    Pnorm = P*fTTmax/TMath::Abs(fTor_Cur);

    //initialize values
    weight = 0.;
    dPhi   = 0.;

    if(Theta < Theta_min_cut && Id != PHOT_ID) return -1;
    if(Id == ELEC_ID && P < E_el_min)          return -1;
    //if(Id == ELEC_ID && P > E_el_min)   cout<<"e- in accept; "<<P<<endl;
    if((TMath::Abs(Id) == PION_ID || TMath::Abs(Id) == KAON_ID) && P < P_pi_min) return -1;

    //negatives 
    if(Id == ELEC_ID || Id == -PION_ID || Id == -KAON_ID || Id == -PROT_ID)
        //if(Id == ELEC_ID)
    {    
        //Define acceptance cuts for electrons
        if(Theta > Theta_min_cut && Fwd)
        {
            if(OutBend)
            {	
                Thcut = fThacc01 - 25.*powf((1.-P*1500./TMath::Abs(fTor_Cur)/5.),24);
                Pnorm = P*fTTmax/TMath::Abs(fTor_Cur);
                Pcutl = (Pnorm > (fThmom1+fThmom2*powf((Theta/fThmom3),fThmom4))*TMath::Abs(fTor_Cur)/fTTmax);
                if(Theta > Thcut && Pcutl)
                {
                    dPhi = fPhimom1*powf(sin((Theta-fPhimom2)*TORAD),fPhimom3)*powf((P*1500./TMath::Abs(fTor_Cur)/5.),fPhimom5);

                    if(TMath::Abs(Phis) < dPhi)
                    {
                        weight=1.;
                    }
                    else
                    {
                        weight=0.;
                    }
                }
                else
                {
                    weight = 0.;
                    dPhi   = 0.;
                }
            }
            else if(InBend)
            {
                Thcut = fThacc1 + fThacc2/((P+0.05)*fTTmax/TMath::Abs(fTor_Cur));
                if(Theta > Thcut)
                {
                    expon = fThacc4*pow(Pnorm,fThacc5);
                    dPhi  = fThacc3*pow(sin((Theta-Thcut)*TORAD),expon); 
                    if(Id == -PION_ID) dPhi *= fThacc6;
                    if(TMath::Abs(Phis) < dPhi)
                    {
                        weight=1.;
                    }
                    else
                    {
                        weight=0.;
                    }
                }
                else
                {
                    dPhi   = 0.;
                    weight = 0.;
                }
            }
        }
        else if(Theta > (fThf_lim+Gap) && ((Id == -PION_ID) || (Id == -KAON_ID) || (Id == -PROT_ID)))
        {
            weight =  1.;
            dPhi   = 30.;
        }
        dPhi /= 30.;
    }//End of simulation for negatives

    //Positively charged particles
    else if(Id == PROT_ID || Id == PION_ID || Id == KAON_ID || Id == 45)
    {
        //minimum momentum cut
        if(P > 0.4 && Fwd && Theta > Theta_min_cut)
        {
            if(InBend)
            {
                Thcut = 5.-17.*powf((1.-P*1500./TMath::Abs(fTor_Cur)/5.),24);
                Pnorm = P*fTTmax/TMath::Abs(fTor_Cur);
                Pcutl = (Pnorm > (0.6 + 4.*powf((Theta/45.),6))*TMath::Abs(fTor_Cur)/fTTmax);
                if(Theta > Thcut && Pcutl)
                {
                    dPhi = 25.*powf(sin((Theta-5.)*TORAD),0.15)*pow((P*1500./TMath::Abs(fTor_Cur)/5.),(1./24.));
                }
                else
                {
                    dPhi = 0.;
                }
            }
            else if(OutBend)
            {
                dPhi = 28.*cos((Theta-37.)*TORAD)*pow((P*fTTmax/TMath::Abs(fTor_Cur)/5.),(1./64.));
            }
            if(TMath::Abs(Phis) < dPhi) weight=1.;
        }
        else if(Theta > (fThf_lim+Gap) && Theta < 100. && (P > (0.25/sin(Theta*TORAD))) && P < 2.)
        {
            dPhi   = 30.;
            weight = 1.;
        }
        else
        {
            dPhi   = 0.;
            weight = 0.;
        }
        dPhi /= 30.;
    }//End of simulation for positives

    //Photons
    else if(Id == PHOT_ID)
    {
        if(P < 0.1) return -1;


        // Central calorimeter: 
        if(Theta > (fThf_lim+Gap) && Theta < 100.)
        {
            //Comment out these two  line if you want central EC            
            dPhi   = 0.0;
            weight = 0.0;
            //Uncomment these two lines if you want central EC    
            //d_phi=28.;
            //weight=1.;
        }
        // Forward EC
        else if(Theta < fThf_lim && Theta > Theta_IC_max)
            //else if(Theta < fThf_lim && Theta > 7)  
        {
            //      cout<<"EC photon!!!!!!!!!!!!!";
            if(P < 0.2) return -1; //check this

            sinth = sin(Theta*TORAD);
            costh = cos(Theta*TORAD);
            sinph = sin(Phis*TORAD);
            cosph = cos(Phis*TORAD);
            r     = (EC_r*(sin(EC_Th)*sin(EC_Th) + cos(EC_Th)) + dl)/(costh + sinth*cosph*sin(EC_Th));
            xh    = r*cosph*sinth;
            yh    = r*sinph*sinth;
            zh    = r*costh - dl;

            EC_XYZ_UVW(xh,yh,zh,out2);
            uh = out2[0];
            vh = out2[1];
            wh = out2[2];
            xi = out2[3];
            yi = out2[4];

            rcoil = dl*costh + sqrt(rc0*rc0 - (dl*sinth)*(dl*sinth));
            hcoil = rcoil*sinth*cosph;
            dcoil = hcoil*sin(0.5236)/cos(0.5236) - TMath::Abs(rcoil*sinth*sinph);

            if(yi > ylow+10. && yi < (yhi-5.))
            {
                xcut = (yi-ylow)/tgrho;
                if(TMath::Abs(xi) < (xcut-5.))
                {
                    weight = 1;
                }
                else if(TMath::Abs(xi) > xcut && dcoil < coil_CALO_width)
                {
                    weight = 1;
                }
                dPhi = 1.-5./(TMath::Abs(r*sinph*sinth) + 8.);
            }

        }

        // Only IC part
        else if(Theta > Theta_IC_min && Theta <= Theta_IC_max)
        {//cout<<"IC photon3!!!!!!!!!!!!!"<<endl;
            //weight =  1.0;
            //dPhi   = 30.0;
            weight =  0.0;
            dPhi   = 0.0;
        }

        else if(Theta > 0. && Theta <= Theta_IC_min)  // simple else avant
        { //cout<<"IC photon!!!!!!!!!!!!!"<<endl;
            weight = 0.;
            dPhi   = 0.;
        }

        else 
        {//cout<<"IC photon2!!!!!!!!!!!!!"<<endl;
            weight = 0.;
            dPhi   = 0.;
        }


    }


    else
    {
        // No acceptance and smearing for this particle
        weight = 0.;
        dPhi   = 0.;
    }

    out[0] = weight;
    out[1] = dPhi;

    //  cout << "out[0]" << out[0] << endl;

    return 0;
}/*}}}*/
//###########################################################################################################################################################################################################

//###########################################################################################################################################################################################################
void smear_fun(double Beta, double P, double Theta, double Phi0, double *out)/*{{{*/
{


    /*
       c---- ----------------------------------------------------------
       c      Smearing momentum and theta 
       c      !common for fastmc and clasev packages! 
       c      April 15 2005
       c---- ----------------------------------------------------------
       c      May 23 2005. K.Mikhailov: add new parameters to conf-files
       c      for the momentum resolution depedence on theta
       c---- ----------------------------------------------------------
       c      May 28 2005. K.Mikhailov: add new parameters to conf-files
       c      for the forward phi-resolution depedence on theta
       c      Change the coge according Gail's phi-resolution study:
       c sigma_1 = 0.0142/sin(theta) in units of mrad (where theta is the polar angle)
       c sigma_1=sigma_1/sin(12^0)=sigma_1/0.2079=0.0142/0.2079/sin(theta)=
       c         0.0683/sin(theta)[mrad]=0.683E-4/sin(theta)[rad]
       c sigma_2 = a0 + a1*theta + a2*theta^2   where theta is in degrees
       c       a0 = 530.6 mrad-MeV; a1= -5.38 mrad-MeV/deg; a2=0.0898 mrad-MeV/deg2
       c       a0 = 0.5306 rad-MeV; a1= -0.00538 rad-MeV/deg; a2=0.0000898 rad-MeV/deg2
       c       a0 = 0.0005306 rad-GeV; a1= -0.00000538 rad-GeV/deg; a2=0.0000000898 rad-GeV/deg2
       c       a0 = 0.5306E-3 rad-GeV; a1= -0.538E-5 rad-GeV/deg; a2=0.898E-7 rad-GeV/deg2
       c---- ----------------------------------------------------------
       c      May 31 2005. K.Mikhailov: add new parameters to conf-files
       c      The new results for the forward phi resolution (vertex+svt+dc)
       c      sigma_1 = a0/theta + a1 + a2*theta where theta is in degrees
       c	    a0 = 3.125; a1 = 0.0262; a2 = -0.00035
       c	    a0 = 0.3125E-2 rad; a1 = 0.262E-4 rad; a2 = -0.35E-6 rad
       c
       c      sigma_2 = a0 + a1*theta + a2*theta^2 where theta is in degrees
       c           a0 = 541.86; a1 = -5.307; a2 = 0.0886
       c           a0 = 0.542E-3 rad-GeV; a1 = -0.531E-5 rad-GeV; a2 = 0.886E-7rad-Gev
       c---------------------------------------------------------------
       */
    double Pout, Thout, Phiout;
    double Thr, Phir;
    double Pnorm;
    // int    InBend, OutBend, Neutral, Fwd, Central;  
    double sig1p, sig2p, sig1th, sig2th, sig1fi, sig2fi;
    double Pperp, dPperp, Pperpout;
    double Erg, IC_Erg_Res, IC_Ang_Res;
    double Theta_IC_max(4.5), Theta_IC_min(2.5);
    double P_IC_Erg[4] = {0.0299700, -0.0076880, 0.00106800, -4.63600e-05};
    double P_IC_Ang[4] = {0.0572505, -0.0144679, 0.00204576, -9.26138e-05};  

    int   ithree(3);
    float_t ran[3];
    //  rnorml_(ran,ithree);
    //ran[0]=.5; // MICHEL
    //ran[1]=.5;
    //ran[2]=.5;
    ran[0] = gRandom->Gaus(0,1);
    ran[1] = gRandom->Gaus(0,1);
    ran[2] = gRandom->Gaus(0,1);


    //TRandom rndm(0); //MUST BE INITIALIZED SOMEWHERE IN THE FIRST PLACE
    //double rran[3], sig[3];
    //rndm.Gaus(rran[0], sig[0]);
    //rndm.Gaus(rran[1], sig[1]);
    //rndm.Gaus(rran[2], sig[2]);

    // since our acceptance functions care about sagitta and not momentum, per se,
    // we will form a "normalized momentum", pnorm
    //cmac this normalized momentum is only used for the fwd. tracker

    Pnorm = P*fTTmax/TMath::Abs(fTor_Cur);
    Thr   = Theta*TORAD;
    Phir  = Phi0*TORAD;

    // select which of the 6 functions we want
    if(InBend && Fwd)
    {

        //ckm 14-apr-2005 (fit on all thetas), theta, phi resolutions:
        sig1p = (fTTmax/fTor_Cur)*(fFoMT1S1 + fFoMT1S2*Theta + fFoMT1S3*(Theta)*(Theta));
        sig2p = (fTTmax/fTor_Cur)*(fFoMT2S1 + fFoMT2S2*Theta + fFoMT2S3*(Theta)*(Theta));

        //ckm   sigma_theta is at 5 degrees(as upper limit)
        sig1th = fFoInThS1;
        sig2th = fFoInThS2;

        sig1fi = (fFoA0FiS1/(Thr*TODEG)) + fFoA1FiS1 + fFoA2FiS1*(Thr*TODEG); 		//Theta in degrees
        sig2fi = fFoA0FiS2 + fFoA1FiS2*Thr*TODEG + fFoA2FiS2*(Thr*TODEG)*(Thr*TODEG);  	//Theta in degrees

        Pout   = P    + ((double)ran[0])*P*sqrt((sig1p*P)*(sig1p*P) + (sig2p/Beta)*(sig2p/Beta));
        Thout  = Thr  + ((double)ran[1])*sqrt((sig1th)*(sig1th)    + (sig2th/P/Beta)*(sig2th/P/Beta));
        Phiout = /*Phir +*/ ((double)ran[2])*sqrt((sig1fi)*(sig1fi)    + (sig2fi/P/Beta)*(sig2fi/P/Beta)); 
    }
    else if(Central && !Neutral)
    {
        sig1p = fCeInMoS1;
        sig2p = fCeInMoS2;
        //ckm May 31, 2005: sigma_2 = a0 + a1*theta + a2*theta^2
        sig1th   = fCeInThS1*(sin(Thr))*(sin(Thr));
        sig2th   = fCeInThS2*(sin(Thr))*(sin(Thr));
        //cmac         sig1th=CeInThS1*(sin(thr))**2
        //cmac         sig2th=CeA0ThS2 + CeA1ThS2*(thr*r2d) + CeA2ThS2*(thr*r2d)**2

        sig1fi   = fCeInFiS1;
        sig2fi   = fCeInFiS2;
        Pperp    = P*sin(Thr);
        //cmac         dpperp=p*sqrt((CeInMoS1*p)**2+(CeInMoS2/beta)**2)
        dPperp   = Pperp*sqrt((fCeInMoS1*P)*(fCeInMoS1*P) + (fCeInMoS2/Beta)*(fCeInMoS2/Beta));
        Pperpout = Pperp + ((double)ran[0])*dPperp;
        Thout    = Thr   + ((double)ran[1])*sqrt((sig1th)*(sig1th) + (sig2th/P/Beta)*(sig2th/P/Beta));
        Pout     = Pperpout/sin(Thout);
        Phiout   = /*Phir +*/ ((double)ran[2])*sqrt((sig1fi)*(sig1fi) + (sig2fi/P/Beta)*(sig2fi/P/Beta));

    }
    else if(OutBend && Fwd)
    {
        //ckm 14-apr-2005 (fit on all thetas), theta, phi resolutions:
        sig1p = fTTmax/fTor_Cur*(fFoMT1S1 + fFoMT1S2*Theta + fFoMT1S3*(Theta)*(Theta));
        sig2p = fTTmax/fTor_Cur*(fFoMT2S1 + fFoMT2S2*Theta + fFoMT2S3*(Theta)*(Theta));
        //ckm   sigma_theta is at 5 degrees(as upper limit)
        sig1th = fFoInThS1;
        sig2th = fFoInThS2;

        //km May 31, 2005
        sig1fi = fFoA0FiS1/(Thr*TODEG) + fFoA1FiS1 + fFoA2FiS1*(Thr*TODEG);  //Theta in degrees
        sig2fi = fFoA0FiS2 + fFoA1FiS2*Thr*TODEG + fFoA2FiS2*(Thr*TODEG)*(Thr*TODEG); //Theta in degrees

        Pout   = P    + ((double)ran[0])*P*sqrt((sig1p*P)*(sig1p*P) + (sig2p/Beta)*(sig2p/Beta));
        Thout  = Thr  + ((double)ran[1])*sqrt((sig1th)*(sig1th)    + (sig2th/P/Beta)*(sig2th/P/Beta));
        Phiout = /*Phir +*/ ((double)ran[2])*sqrt((sig1fi)*(sig1fi)    + (sig2fi/P/Beta)*(sig2fi/P/Beta));

    }
    else if(Neutral && Fwd)
    {//cout<<"COUCOU4!!!!!!!!"<<endl;
        if(Theta > Theta_IC_max)
        {
            // Neutral in FWD EC
            Pout   = P    + ((double)ran[0])*sqrt(P)*fFoNeMoS1;
            //      Pout   = P    + ((double)ran[0])*sqrt(P)*.012;
            Thout  = Thr  + ((double)ran[1])*fFoNeThS1;
            //      Thout  = Thr  + ((double)ran[1])/sqrt(P)*.02*TORAD;
            Phiout = /*Phir +*/ ((double)ran[2])*fFoNeFiS1;											// CHgt BAPTISTE
        }
        else if(Theta > Theta_IC_min && Theta < Theta_IC_max)
        {
            // Neutral in IC 
            if(Beta > 0 && Beta < 1)
            {
                Erg = P/Beta;
            }
            else
            {
                Erg = P;
            }

            // Find IC energy resolution from Rustam parameterization
            //IC_Erg_Res = P_IC_Erg[0]*Erg + P_IC_Erg[1]*(Erg)*(Erg) + P_IC_Erg[2]*pow(Erg,3) + P_IC_Erg[3]*pow(Erg,4);
            //cout<<"IC_Erg_Res 1 ="<<P_IC_Erg[0]*Erg + P_IC_Erg[1]*(Erg)*(Erg) + P_IC_Erg[2]*pow(Erg,3) + P_IC_Erg[3]*pow(Erg,4)<<endl;

            //IC_Erg_Res = Erg*sqrt(pow((0.053/Erg),2)+pow((0.035/sqrt(Erg)),2)+pow(0.017,2));
            //cout<<"IC_Erg_Res 2 ="<<Erg*sqrt(pow((0.053/Erg),2)+pow((0.035/sqrt(Erg)),2)+pow(0.017,2))<<endl;
            //IC_Erg_Res = 0.02*sqrt(Erg);
            IC_Erg_Res=Erg*sqrt((0.02/Erg)*(0.02/Erg) + 0.01*0.01);

            // Find IC angular resolution from Rustam parameterization 
            // Assume that Theta and Phi resolutions are the same
            //      IC_Ang_Res = P_IC_Ang[0] + P_IC_Ang[1]*Erg + P_IC_Ang[2]*(Erg)*(Erg) + P_IC_Ang[3]*pow(Erg,3);
            IC_Ang_Res = 57.2557*0.2/(186*sqrt(Erg));

            // Smear energy and angles with a gaussian
            Pout   = P    + ((double)ran[0])*IC_Erg_Res; 
            Thout  = Thr  + ((double)ran[1])*IC_Ang_Res*TORAD;
            Phiout = /*Phir +*/ ((double)ran[2])*IC_Ang_Res*TORAD;	
        }
    }
    else if(Neutral && Central)
    {//cout<<"COUCOU5!!!!!!!!"<<endl;
        Pout   = P    + ((double)ran[0])*sqrt(P)*fCeNeMoS1;
        Thout  = Thr  + ((double)ran[1])*fCeNeThS1;
        Phiout = /*Phir +*/ ((double)ran[2])*fCeNeFiS1;

    }

    //Convert from radians to degrees
    Thout  *= TODEG;
    Phiout *= TODEG;
    //cout<<"Phiout_correction"<<Phiout<<endl;

    out[0] = Pout;
    out[1] = Thout;
    out[2] = Phiout;



}/*}}}*/
//###############################################################################

//###############################################################################
void Charge_Type(int Id, double *out)/*{{{*/
{
    double Charge, Mass;
    Charge = (double)(TMath::Sign(1,Id));
    if(TMath::Abs(Id) == ELEC_ID)
    {
        Charge = -(double)(TMath::Sign(1,Id));
        Mass   = M_ELEC;
    }
    else if(TMath::Abs(Id) == PION_ID)
    {
        Mass   = M_PION;
    }
    else if(TMath::Abs(Id) ==  KAON_ID)
    {
        Mass   = M_KAON;
    }
    else if(TMath::Abs(Id) == PROT_ID)
    {
        Mass   = M_PROT;
    }
    else if(TMath::Abs(Id) == 45)//?
    {
        Mass   = 0.5;
    }
    else if(Id == PHOT_ID)
    {
        Charge = 0.;
        Mass   = M_PHOT;
    }
    else if(Id == NEUT_ID)
    {
        Charge = 0.;
        Mass   = M_NEUT;
    }
    else
    {
        //cout<<"Unsupported ID ->  "<<Id<<endl;
        Charge = -999.;
        Mass   = -999.;
    }
    out[0] = Charge;
    out[1] = Mass;
}/*}}}*/
//###########################################################################################################################################################################################################

//###########################################################################################################################################################################################################
// A function to calculate BdL of the solenoid 
// assuming uniform magnetic field of solenoid 
// specified in the pCDR of CLAS12
double findBdL(double P, double E, double Theta)/*{{{*/
{
    //Input: momentum, energy and the polar angle of 
    //       charged particle in GeV and rad
    //Output: calculated BdL in T x m 

    double BdL = -999;
    double Time, TimeL, TimeR;
    //B-field, diameter and length of the solenoid from Vlassov plots
    double B(5.0), dBarrel(0.78), lBarrel(0.90), rBarrel;
    //speed of light in m/ns 
    double cLight(0.299792458);
    double v, v_t, v_l;
    double P_t, r_curv;

    //Find the velocity and the components  
    if(E > 0 && P <= E)
    {
        v = cLight*P/E; //beta.c
    }
    else
    {
        cout<<"Error: Bad kinematics"<<endl;
        v = 0.0;
    }
    v_l = v*cos(Theta);
    v_t = v*sin(Theta);
    P_t = P*sin(Theta);

    //Find the radious of the curvature in the solenoid 
    if(TMath::Abs(B)> 0)
    {
        r_curv = P_t/(cLight*B);
    }
    else
    {
        r_curv = 1.0e+15;
    }

    //Find the time it would take to leave the solenoid 
    //through the base of the cylinder
    if(v_l > 0 )
    {
        TimeL = lBarrel/2.0/v_l;
    }
    else
    {
        TimeL = 1.0e+10;
    }

    //Find te time it would take to leave the solenoid
    //through the walls of the cylinder
    rBarrel = dBarrel/2.0;
    if(rBarrel < (2.0*r_curv) && v_t > 0)
    {
        TimeR = 2.0*r_curv*asin(rBarrel/(2.0*r_curv))/v_t;
    }
    else
    {
        TimeR = TimeL;
    }

    //Take the minimum of the two
    Time = min(TimeL, TimeR);

    //Final value
    BdL = B * v_t * Time;

    return BdL;
}/*}}}*/
//###########################################################################################################################################################################################################

//###########################################################################################################################################################################################################
void EC_XYZ_UVW(double x, double y, double z, double *out)/*{{{*/
{
    double EC_Th(0.4363323), ylow(-182.974), yhi(189.956);
    double tgrho(1.95325), sinrho(0.8901256), cosrho(0.455715);
    double EC_Phi, Phi;
    double xi,yi,zi,u,v,w;
    double rot[3][3];

    Phi = atan2(y,x)*57.29578;
    if(Phi < 0.) Phi += 360.;
    Phi += 30.;
    if(Phi >= 360.) Phi -= 360.;
    EC_Phi = (int)(Phi/60.)*1.0471975;

    rot[0][0] = cos(EC_Th)*cos(EC_Phi);
    rot[0][1] = -sin(EC_Phi);
    rot[0][2] = sin(EC_Th)*cos(EC_Phi);
    rot[1][0] = cos(EC_Th)*sin(EC_Phi);
    rot[1][1] = cos(EC_Phi);
    rot[1][2] = sin(EC_Th)*sin(EC_Phi);
    rot[2][0] = -sin(EC_Th);
    rot[2][1] = 0.;
    rot[2][2] = cos(EC_Th);

    yi = x*rot[0][0] + y*rot[1][0] + z*rot[2][0];
    xi = x*rot[0][1] + y*rot[1][1] + z*rot[2][1];
    zi = x*rot[0][2] + y*rot[1][2] + z*rot[2][2] - 510.32;

    u  = (yi-ylow)/sinrho;
    v  = (yhi-ylow)/tgrho-xi+(yhi-yi)/tgrho;
    w  = ((yhi-ylow)/tgrho+xi+(yhi-yi)/tgrho)/2./cosrho;

    out[0] = u;
    out[1] = v;
    out[2] = w;
    out[3] = xi;
    out[4] = yi;
}/*}}}*/
//###########################################################################################################################################################################################################

//###########################################################################################################################################################################################################
void type(double Charge, double Theta, int *out)/*{{{*/
{


    // subroutine which figures out the track type
    // common block keeping the logical values of the track type;      

    //bool InBend, OutBend, Neutral, Fwd, Central;
    //common/track_type/inbend,outbend,neutral,fwd,central
    //#include "read_par_clas12g.inc"
    // initialize logicals to false
    int inbend, outbend, neutral, fwd, central;
    inbend  = 0;
    outbend = 0;
    neutral = 0;
    fwd     = 0;
    central = 0;
    //  6 functions: (inbend,outbend,neutral) track X (fwd,central) det.

    if(Charge == 0) neutral = 1;
    if(TMath::Abs(Charge) == 1)
    {
        if((Charge*fTor_Cur) > 0.) outbend = 1;
        if((Charge*fTor_Cur) < 0.) inbend  = 1;
    }
    if(Theta <  fThf_lim) fwd     = 1;  //Thf_lim max forward angle
    if(Theta >= fThf_lim) central = 1;
    out[0] = inbend;
    out[1] = outbend;
    out[2] = neutral;
    out[3] = fwd;
    out[4] = central;

}/*}}}*/
//###########################################################################################################################################################################################################

//###########################################################################################################################################################################################################
int clas12(int Id, double P, double Theta, double Phi0, double *out)/*{{{*/
{

    double Charge(0.), Mass(0.), E(0.);
    double Beta(1.0);
    double out1[2]={0.,0.};
    int    out2[5]={0,0,0,0,0};
    double out3[2]={0.,0.};
    double out4[5]={0.,0.,0.,0.,0.};
    int    iacc(-1);
    double weight(0.), dphi(0.);
    for(int i=0; i<5; i++) out[i] = 0.;

    // 1st decide if particle charge is +1, -1, 0 or unknown
    Charge_Type(Id, out1);
    Charge = out1[0]; 
    Mass   = out1[1];  //cout<<"Chargecla12="<<Charge<<"Masscla12="<<Mass<<"P="<<P<<endl; //BAPTISTE
    E      = sqrt(P*P + Mass*Mass);
    if(Charge == -999 || Mass == -999)
    {
        return -1;
    }
    Beta = P/sqrt(P*P + Mass*Mass);

    // 2nd, determine which track type we are concerned with
    type(Charge, Theta, out2);
    InBend  = out2[0]; 
    OutBend = out2[1]; 
    Neutral = out2[2];
    Fwd     = out2[3]; 
    Central = out2[4]; 

    // 3rd, call acceptance function with charge, torus current,
    // and p,th,phi as inputs; return is weight (1, accepted; 0, not)
    // and d_phi which is the phi part of the acceptance.  Note, this
    // is the tracking acceptance only, it does NOT contain the additional
    // fiducial volume cuts of Cerenkov counters imposed on electrons.

    // get a 3-dim. random number
    int   ithree(3);
    float_t ran[3];
    //  rnorml_(ran,ithree);
    //ran[0]=.5; // MICHEL
    //ran[1]=.5;
    //ran[2]=.5;
    ran[0] = gRandom->Gaus(0,1);
    ran[1] = gRandom->Gaus(0,1);
    ran[2] = gRandom->Gaus(0,1);


    if(TMath::Abs(Charge) == 1 || Charge == 0)
    {
        iacc = accep_func(Id, P, Theta, Phi0, out3);
        if(iacc != -1)
        {
            weight = out3[0];
            dphi   = out3[1];
        }
    }

    // 4th, if particle is accepted (weight=1.) call smearing routine
    // with input momenta and output momenta
    if(iacc != -1 && weight == 1)
    {
        smear_fun(Beta, P, Theta, Phi0, out4);

        // Now set p,th,phi to smeared values
        for(int i=0; i<3; i++) out[i] = out4[i];
        out[3] = weight;
        out[4] = dphi;
    }

    return iacc;
}/*}}}*/
//###########################################################################################################################################################################################################
#endif
