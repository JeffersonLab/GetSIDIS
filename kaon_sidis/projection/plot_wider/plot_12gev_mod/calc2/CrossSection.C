#include<math.h>
#include<stdlib.h>
#include "vegas_C.h"
#include "iostream.h"
#include "target.h"

double s0 ;
//double Q2min=1.0,W2min=5.29,s0=21.5,xmin=0.05,xmax=0.55,ymin=0.34,ymax=0.9,zmin=0.3,zmax=0.7,phtmin=0.0,phtmax=1.6;
//double Q2min=1.0,W2min=25.0,s0=300.0,xmin=0.003,xmax=0.21,ymin=0.1,ymax=0.9,zmin=0.2,zmax=0.7,phtmin=0.0,phtmax=1.6;
//double Q2min=1.0,W2min=10.0,s0=51.7,xmin=0.023,xmax=0.4,ymin=0.1,ymax=0.85,zmin=0.2,zmax=0.7,phtmin=0.0,phtmax=1.6;

long ncall=5000;
int init=0,itmx=5,nprn=5;
double x0,Q20,z0,pht0;
int hadron0;

#define PI 3.1415927
#define ALPHA (1.0/137)
double Mh = 0.1396; 


// #ifdef PROTON
// double	Mn=0.9383;    //the mass of a proton
// int u=1,d=2,s=3,ubar=-1,dbar=-2,sbar=-3;
// #else
double Mn=0.9396;     //the mass of a neutron
int u=2,d=1,s=3,ubar=-2,dbar=-1,sbar=-3;
// #endif

long idum = 13213241;


double Cro_Sec_Unp(double x,double Q2,double z,double pht,double s,int hadron)
{
  //cout << u << endl;

  double Cro_Sec_Unp_kt(double[],double);
	int ndim=2;
	double tgral,sd,chi2a,regn[5]={0,0,0,5.0,2*PI};
	x0=x;
	Q20=Q2;
	z0=z;
	pht0=pht;
	hadron0=hadron;
	s0 = s;
	vegas(regn,ndim,Cro_Sec_Unp_kt,init,ncall,itmx,nprn,&tgral,&sd,&chi2a);
	return(tgral);
}

double Cro_Sec_Tran(double x,double Q2,double z,double pht,double s,int hadron)
{
	double Cro_Sec_Tran_kt(double[],double);
	int ndim=2;
	double tgral,sd,chi2a,regn[5]={0,0,0,5.0,2*PI};
	x0=x;
	Q20=Q2;
	z0=z;
	pht0=pht;
	hadron0=hadron;
	s0 = s;
	vegas(regn,ndim,Cro_Sec_Tran_kt,init,ncall,itmx,nprn,&tgral,&sd,&chi2a);
	return(tgral);
}

double Cro_Sec_Siv(double x,double Q2,double z,double pht,double s,int hadron)
{
	double Cro_Sec_Siv_kt(double[],double);
	int ndim=2;
	double tgral,sd,chi2a,regn[5]={0,0,0,5.0,2*PI};
	x0=x;
	Q20=Q2;
	z0=z;
	pht0=pht;
	hadron0=hadron;
	s0 = s;
	vegas(regn,ndim,Cro_Sec_Siv_kt,init,ncall,itmx,nprn,&tgral,&sd,&chi2a);
	return(tgral);
}

double Cro_Sec_Pret(double x,double Q2,double z,double pht,double s,int hadron)
{
	double Cro_Sec_Pret_kt(double[],double);
	int ndim=2;
	double tgral,sd,chi2a,regn[5]={0,0,0,5.0,2*PI};
	x0=x;
	Q20=Q2;
	z0=z;
	pht0=pht;
	hadron0=hadron;
	s0 = s;
	vegas(regn,ndim,Cro_Sec_Pret_kt,init,ncall,itmx,nprn,&tgral,&sd,&chi2a);
	return(tgral);
}

double Cro_Sec_Unp_kt(double x[],double wgt)
{
	double un_pdf(int,double,double);
	double un_ff(double,double,int,int);
	int iparton;
	double charge,angle_dep,pdf,frag,pt,y0,result=0.0;
	pt=sqrt(pht0*pht0/z0/z0+x[1]*x[1]+2*pht0/z0*x[1]*cos(x[2]));
	angle_dep=1.0;
	y0=Q20/s0/x0;
	for(iparton=-3;iparton<=3;iparton++)
	{
		switch(iparton)
		{
		case -3:
			{
				charge=1.0/9;
				pdf=un_pdf(sbar,x0,pt);
				break;
			}
		case -2:
			{
				charge=1.0/9;
				pdf=un_pdf(dbar,x0,pt);
				break;
			}
		case -1:
			{
				charge=4.0/9;
				pdf=un_pdf(ubar,x0,pt);
				break;
			}
		case 0:
			break;
		case 1:
			{
				charge=4.0/9;
				pdf=un_pdf(u,x0,pt);
				break;
			}
		case 2:
			{
				charge=1.0/9;
				pdf=un_pdf(d,x0,pt);
				break;
			}
		case 3:
			{
				charge=1.0/9;
				pdf=un_pdf(s,x0,pt);
				break;
			}
		}
		frag=un_ff(z0,z0*x[1],iparton,hadron0);
		result=result+2*ALPHA*ALPHA/s0/x0/y0/y0*(1-y0+0.5*y0*y0)
			*charge*2*PI*angle_dep*pht0*x[1]*pdf*frag;
	}
	return(result);
}

double Cro_Sec_Tran_kt(double x[],double wgt)
{
	double trans_pdf(int,double,double);
	double col_ff(double,double,int,int);
	int iparton;
	double charge,angle_dep,pdf,frag,pt,y0,result=0.0;
	pt=sqrt(pht0*pht0/z0/z0+x[1]*x[1]+2*pht0/z0*x[1]*cos(x[2]));
	angle_dep=-1*x[1]*cos(x[2])/Mh;
	y0=Q20/s0/x0;
	for(iparton=1;iparton<=2;iparton++)
	{
		if(iparton==1)
		{
			charge=4.0/9;
			pdf=trans_pdf(u,x0,pt);
		}
		else if(iparton==2)
		{
			charge=1.0/9;
			pdf=trans_pdf(d,x0,pt);
		}
		else
		{
			charge=0;
			pdf=0.0;
		}
		frag=col_ff(z0,z0*x[1],iparton,hadron0);
		result=result+2*ALPHA*ALPHA/s0/x0/y0/y0*(1-y0)
			*charge*2*PI*angle_dep*pht0*x[1]*pdf*frag;
	}
	return(result);
}

double Cro_Sec_Siv_kt(double x[],double wgt)
{
	double siv_pdf(int,double,double);
	double un_ff(double,double,int,int);
	int iparton;
	double charge,angle_dep,pdf,frag,pt,y0,result=0.0;
	pt=sqrt(pht0*pht0/z0/z0+x[1]*x[1]+2*pht0/z0*x[1]*cos(x[2]));
	angle_dep=-1.0*(x[1]*cos(x[2])+pht0/z0)/Mn;
	y0=Q20/s0/x0;
	for(iparton=-3;iparton<=3;iparton++)
	{
		switch(iparton)
		{
		case -3:
			{
				charge=1.0/9;
				pdf=siv_pdf(sbar,x0,pt);
				break;
			}
		case -2:
			{
				charge=1.0/9;
				pdf=siv_pdf(dbar,x0,pt);
				break;
			}
		case -1:
			{
				charge=4.0/9;
				pdf=siv_pdf(ubar,x0,pt);
				break;
			}
		case 0:
			break;
		case 1:
			{
				charge=4.0/9;
				pdf=siv_pdf(u,x0,pt);
				break;
			}
		case 2:
			{
				charge=1.0/9;
				pdf=siv_pdf(d,x0,pt);
				break;
			}
		case 3:
			{
				charge=1.0/9;
				pdf=siv_pdf(s,x0,pt);
				break;
			}
		}
		frag=un_ff(z0,z0*x[1],iparton,hadron0);
		result=result+2*ALPHA*ALPHA/s0/x0/y0/y0*(1-y0+0.5*y0*y0)
			*charge*2*PI*angle_dep*pht0*x[1]*pdf*frag;
	}
	return(result);
}

double Cro_Sec_Pret_kt(double x[],double wgt)
{
	double pretz_pdf(int,double,double);
	double col_ff(double,double,int,int);
	int iparton;
	double charge,angle_dep,pdf,frag,pt,y0,result=0.0;
	pt=sqrt(pht0*pht0/z0/z0+x[1]*x[1]+2*pht0/z0*x[1]*cos(x[2]));
	angle_dep=-1*x[1]*(pht0*pht0/z0/z0*cos(x[2])+2*pht0/z0*x[1]*cos(2*x[2])+x[1]*x[1]*cos(3*x[2]))
		/2/Mn/Mn/Mh;
	y0=Q20/s0/x0;
	for(iparton=1;iparton<=2;iparton++)
	{
		if(iparton==1)
		{
			charge=4.0/9;
			pdf=pretz_pdf(u,x0,pt);
		}
		else if(iparton==2)
		{
			charge=1.0/9;
			pdf=pretz_pdf(d,x0,pt);
		}
		else
		{
			charge=0;
			pdf=0.0;
		}
		frag=col_ff(z0,z0*x[1],iparton,hadron0);
		result=result+2*ALPHA*ALPHA/s0/x0/y0/y0*(1-y0)
			*charge*2*PI*angle_dep*pht0*x[1]*pdf*frag;
	}
	return(result);
}
