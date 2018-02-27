/****************************************************************************

		 	EPPS16.cpp

 A C++ interface for the scale dependent nuclear modifications

   		R_f^A(x,Q) = f_A(x,Q)/f_p(x,Q) 

 where f_A is the distribution of the parton flavour f for a PROTON in a
 nucleus A, and f_p is the corresponding parton distribution in the free proton.
  
 When using this interface, please refer to:
  
 K.J. Eskola, P. Paakkinen, H. Paukkunen and C.A. Salgado,
 "EPPS16: Nuclear parton distributions with LHC data."
 Published in EPJC
 Eprint: arXiv:1612.05741

 Questions & comments to:
   hannu.paukkunen@jyu.fi
   petja.paakkinen@jyu.fi
   kari.eskola@jyu.fi
   carlos.salgado@usc.es
 
 ***************************************************************************
 Instructions:

 For given input values of

     order (integer): dummy variable, retained for
                      compatibility with EPS09 routine
     
     pset (integer):
            1     = central fit
            2,3   = error sets S{+1}, S{-1}
            4,5   = error sets S{+2}, S{-2}
            ...   ...
            40,41 = error sets {S+20}, {S-20}

     A (integer): atomic number
     x (double) : Bjorken-x
     Q (double) : scale in GeV

 the function (returning void)

   EPPS16(order, pset, A, x, Q, ruv, rdv, ru, rd, rs, rc, rb, rg)

 returns the bound proton nuclear corrections R_f^A(x,Q) (double) for
	
	ruv = up valence
	rdv = down valence
	ru  = up sea
	rd  = down sea
	rs  = strange
	rc  = charm
	rb  = bottom
	rg  = gluons

 The nuclear corrections for bound neutrons can be obtained
 by the isospin symmetry, e.g. the total up quark distribution
 per nucleon in a nucleus A with Z protons is

  u_A(x,Q) =    Z/A * [ruv*uV_p(x,Q) + ru*uSea_p(x,Q)] +
            (A-Z)/A * [rdv*dV_p(x,Q) + rd*dSea_p(x,Q)]

 Note that the parametrization should only be applied at the
 kinematical domain

             1e-7 <= x <= 1
              1.3 <= Q <= 10000 GeV.

 Outside these boundaries the code stops and
 throws an error message.

 The data used by the program for required order
 and atomic number A, are stored in separate files

   NLO: EPPS16NLOR_A

 which are assumed to be located in the working directory.

 The error bands for absolute cross-sections and for
 their nuclear ratios should be computed as explained
 in Secs. 2.5 and 4 of arXiv:0902.4154 [hep-ph]. For
 the absolute cross sections, both the errors in the
 free-proton PDFs f_p(x,Q) and the errors in
 the modifications R_f^A(x,Q) should be accounted for.
 For the nuclear ratios, it is usually sufficient to 
 account only for the errors in the modifications R_f^A(x,Q).

*********************************************************
*********************************************************/

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstdlib>

using namespace std;

void luovi(const double* f, const double* arg, int mmm, double z, double &sum);

// *************************************************************

 void EPPS16(const int    &order, 
             const int    &pset, 
             const int    &A, 
             const double &xin,
             const double &Qin,
             double       &ruv,
             double       &rdv,
             double       &ru,
             double       &rd,
             double       &rs,
             double       &rc,
             double       &rb,
             double       &rg)
{

     int    const xsteps=80;
     int    const Qsteps=30;
     int    const polorder=4;

     int    charmflag, bottomflag;

     double const x_i=0.0000001;
     double const Q2min=1.690;
     double const Q2max=100000000.0;

     double x, Q2, realQ, LSTEP, n_x, arg[5],fu[5],fg[5];
     char filenimi[50];   
     double dummydouble, res, result[9];
     int Qpoint, xpoint;

     static int Alast=-1;
     static double allvalues[42][9][51][81];

     charmflag  = 0;
     bottomflag = 0;

// *********************************************
// Stop if the set specifications are wrong ones
// *********************************************

  if (pset < 1 || pset > 41){
      cout << "In EPPS16: Wrong set!" << " " << pset << endl;
      cout << "Central set: pset = 1" << endl;
      cout << "Error sets : pset = 2...41" << endl;
      std::exit(EXIT_FAILURE);
  }

      Q2 = Qin*Qin;
       x = xin;

// *******************************
// Freeze x if it's < 10E-6 or > 1
// *******************************

      if (x<x_i){
      cout << "In EPPS16: x<1E-7" << endl;
      std::exit(EXIT_FAILURE);
      }
      if (x > 1){
      cout << "In EPPS16: x>1" << endl;
      std::exit(EXIT_FAILURE);
      }

// ************************************
// Freeze Q^2 if it's < 1.69 or > 10E+6
// ************************************

      if (Q2 < Q2min){
      cout << "In EPPS16: Q<1.3GeV" << endl;
      std::exit(EXIT_FAILURE);
      }
      if (Q2 > Q2max){
      cout << "In EPPS16: Q>10000GeV" << endl;
      std::exit(EXIT_FAILURE);
      }

// If the set specifications have been changed, read the tables again

 if (A != Alast){
  sprintf(filenimi,"./grid/EPPS16NLOR_%d", A);
  
  fstream myfile;

  myfile.open (filenimi, ios::in);
  if (!myfile.is_open()){
  cout << "In EPPS16: Missing file " << filenimi << endl;
  exit(EXIT_FAILURE);
  }

  for (int setnumber=1;setnumber<42;setnumber++){

   for (int k=0;k<(Qsteps+1);k++){
   myfile >> dummydouble;

    for (int t=0;t<xsteps;t++){

     for (int flavour=1;flavour<9;flavour++){
     myfile >> allvalues[setnumber][flavour][k][t];
     }

    } // end xloop

   }  // end Qloop

  } 

 myfile.close();

 Alast     = A;
 }

// Find out the position in the loglog Q^2-grid

      realQ  = Qsteps * (log(log(Q2)/log(Q2min)))/
                        (log(log(Q2max)/log(Q2min)));
      Qpoint = realQ;

      if (Qpoint == 0){
         Qpoint = 1;
      }
      if (Qpoint > (Qsteps-2)){
      Qpoint = Qsteps-2;
      }

// Find the position in the x-grid

      LSTEP = (0.0 - (log(1.0/x_i) + 2.0*(1.0-x_i)))/(1.0*xsteps);

      n_x  = ((log(1.0/x)+2.0*(1-x))-(log(1.0/x_i)+2.0*(1.0-x_i)))/LSTEP;
      xpoint = n_x;

// *********************
// Interpolate the grids 
// *********************

  for (int flavour=1;flavour<9;flavour++){

      if (flavour > 2 && flavour < 8){
      if (xpoint == 0){
      xpoint = 1;
      }
      else if (xpoint > (xsteps-6)){
      xpoint = xsteps - 6;
      }
      }

      else{
      if (xpoint == 0){
      xpoint = 1;
      }
      else if (xpoint > (xsteps-4)){
      xpoint = xsteps - 4;
      }
      }


   for (int k=1;k<5;k++){
   arg[k] = xpoint-2+k;
   }

// For charm, don't use the Q=1.3 for interpolation

      if (flavour == 6 && Qpoint == 1){
      Qpoint = 2;
      charmflag = 1;
      }

// For bottom, only use points Q>4.75 for interpolation

      if (flavour == 7 && Qpoint < 17 && Qpoint > 1){
      bottomflag = Qpoint;
      Qpoint = 17;
      }


   for (int j=1;j<5;j++){
   fu[1] = allvalues[pset][flavour][Qpoint-2+j][xpoint-1];
   fu[2] = allvalues[pset][flavour][Qpoint-2+j][xpoint+0];
   fu[3] = allvalues[pset][flavour][Qpoint-2+j][xpoint+1];
   fu[4] = allvalues[pset][flavour][Qpoint-2+j][xpoint+2];

   luovi(fu,arg,polorder,n_x,res);

   fg[j] = res;
   }

   for (int k=1;k<5;k++){
   arg[k] = Qpoint-2+k;
   }

   luovi(fg,arg,polorder,realQ,res);

   result[flavour] = res;

// For charm, put back the original Qpoint

      if (charmflag == 1){
      Qpoint = 1;
      charmflag = 0;
      }

// For bottom, put back the original Qpoint

      if (bottomflag > 1){
      Qpoint = bottomflag;
      bottomflag = 0;
      }

  }

      ruv = result[1];
      rdv = result[2];
      ru  = result[3];
      rd  = result[4];
      rs  = result[5];
      rc  = result[6];
      rb  = result[7];
      rg  = result[8];

// Put bottom to zero below the mass threshold

      if (Qin < 4.75){
      rb  = 0.0;
      }

}

//
// Modified version of Cern Library
// interpolation routine E100
//

void luovi(const double* f, const double* arg, int mmm, double z, double &sum)
{
    double *cof = new double[mmm+1];
    int jndex, index;
    
    int mm = mmm < 20 ? mmm : 20;
    int m = mm - 1;
    for (int i=1; i<=mm; i++) cof[i] = f[i];
    for (int i=1; i<=m; i++) {
        for (int j=i; j<=m; j++){
            jndex = mm - j;
            index = jndex + i;
            cof[index] = (cof[index]-cof[index-1])/(arg[index]-arg[jndex]);
        }
    }
    sum = cof[mm];
    for (int i=1; i<=m; i++){
        index = mm - i;
        sum = (z-arg[index])*sum + cof[index];
    }
    delete [] cof;
}
