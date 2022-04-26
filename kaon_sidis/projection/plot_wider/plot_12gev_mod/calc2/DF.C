#include<math.h>
#define PI 3.1415927

extern double Mn;
static double mq,md;

#include "LHAPDF/LHAPDF.h"
using namespace LHAPDF;

double un_pdf(int iparton,double x,double pt)          //带横动量的非极化分布函数
{
	double Q=1.5,result=0;
	
	if (iparton==1){
	  result = xfx(x, Q, 2)/x;
	}else if (iparton==2){
	  result = xfx(x, Q, 1)/x;
	}else if (iparton==-1){
	  result = xfx(x, Q, -2)/x;
	}else if (iparton==-2){
	  result = xfx(x, Q, -1)/x;
	}


	result*=exp(-pt*pt/0.25)/PI/0.25;
	return(result);
}

double trans_pdf(int iparton,double x,double pt)		//带横动量的transversity
{
	double un_pdf(int,double,double);
	double ws,wv,Q=2.0,result,M_cal;
	int u=1,d=2,ubar=-1,dbar=-2;
	mq=0.33;
	md=0.6;
	M_cal=sqrt((mq*mq+pt*pt)/x+(md*md+pt*pt)/(1-x));
	ws=(x*M_cal+mq)*(x*M_cal+mq)/((x*M_cal+mq)*(x*M_cal+mq)+pt*pt);
	md=0.8;
	M_cal=sqrt((mq*mq+pt*pt)/x+(md*md+pt*pt)/(1-x));
	wv=(x*M_cal+mq)*(x*M_cal+mq)/((x*M_cal+mq)*(x*M_cal+mq)+pt*pt);
	if(iparton==1)
		result=((un_pdf(u,x,pt)-un_pdf(ubar,x,pt))-0.5*(un_pdf(d,x,pt)-un_pdf(dbar,x,pt)))*ws
		-1.0/6*(un_pdf(d,x,pt)-un_pdf(dbar,x,pt))*wv;
	else if(iparton==2)
		result=-1.0/3*(un_pdf(d,x,pt)-un_pdf(dbar,x,pt))*wv;
	else
		result=0.0;
	return(result);
}

double siv_pdf(int iparton,double x,double pt)
{
	double un_pdf(int,double,double);
	double result;
	double u=1,d=2,s=3,ubar=-1,dbar=-2,sbar=-3;
	switch(iparton)
	{
	case -3:result=1*pow(x,0.79)*pow(1-x,3.46)*7.697;break;
	case -2:result=-0.4*pow(x,0.79)*pow(1-x,3.46)*7.697;break;
	case -1:result=0.04*pow(x,0.79)*pow(1-x,3.46)*7.697;break;
	case 0:break;
	case 1:result=0.35*pow(x,0.73)*pow(1-x,3.46)*6.945;break;
	case 2:result=-0.9*pow(x,1.08)*pow(1-x,3.46)*12.07;break;
	case 3:result=-0.24*pow(x,0.79)*pow(1-x,3.46)*7.697;break;
	}
	result=result*-Mn*un_pdf(iparton,x,pt)*2.332/0.5831*exp(-pt*pt/0.34);
	return(result);
}

double pretz_pdf(int iparton,double x,double pt)		//带横动量的pretzelosity
{
	double un_pdf(int,double,double);
	double ws,wv,Q=2.0,result,M_cal;
	int u=1,d=2,ubar=-1,dbar=-2;
	mq=0.33;
	md=0.6;
	M_cal=sqrt((mq*mq+pt*pt)/x+(md*md+pt*pt)/(1-x));
	ws=-2*Mn*Mn/((x*M_cal+mq)*(x*M_cal+mq)+pt*pt);
	md=0.8;
	M_cal=sqrt((mq*mq+pt*pt)/x+(md*md+pt*pt)/(1-x));
	wv=-2*Mn*Mn/((x*M_cal+mq)*(x*M_cal+mq)+pt*pt);
	if(iparton==1)
		result=((un_pdf(u,x,pt)-un_pdf(ubar,x,pt))-0.5*(un_pdf(d,x,pt)-un_pdf(dbar,x,pt)))*ws
		-1.0/6*(un_pdf(d,x,pt)-un_pdf(dbar,x,pt))*wv;
	else if(iparton==2)
		result=-1.0/3*(un_pdf(d,x,pt)-un_pdf(dbar,x,pt))*wv;
	else
		result=0.0;
	return(result);
}
