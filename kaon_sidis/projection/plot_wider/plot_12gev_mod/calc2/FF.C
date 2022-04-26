#include<math.h>
#define PI 3.1415927
extern double Mh;

double un_ff(double z,double kt,int parton,int meson)
{			
	double un_ff_fav,un_ff_unf,value;
	double a1=-1.039,b1=1.241,N1=0.689,a2=-1.805,b2=2.037,N2=0.217;
	int flag;
	un_ff_fav=N1*pow(z,a1)*pow(1-z,b1);
	un_ff_unf=N2*pow(z,a2)*pow(1-z,b2);
	if(meson!=0)
	{
		flag=parton*meson;
		switch(flag)
		{
		case -2:
		case 1:
			{
				value=un_ff_fav;
				break;
			}
		default:{value=un_ff_unf;}
		}
	}
	else
		value=0.5*(un_ff_fav+un_ff_unf);
	value=value*exp(-kt*kt/0.2)/PI/0.2;
	return(value);
}

double col_ff(double z,double kt,int parton,int meson)			
{
	double  col_ff_fav,col_ff_unf,value,un_temp;
	double M=0.9593;
	int flag;
	un_temp=un_ff(z,kt,parton,meson);
	col_ff_fav=0.43*1.057*pow(z,0.96)*pow(1-z,0.01)*2.332*z*Mh/M
		*un_temp*exp(-kt*kt/M/M);
	col_ff_unf=-1.00*1.057*pow(z,0.96)*pow(1-z,0.01)*2.332*z*Mh/M
		*un_temp*exp(-kt*kt/M/M);
	if(meson!=0)
	{
		flag=parton*meson;
		switch(flag)
		{
		case -2:
		case 1:
			{
				value=col_ff_fav;
				break;
			}
		default:{value=col_ff_unf;}
		}
	}
	else
		value=0.5*(col_ff_fav+col_ff_unf);
	return(value);
}

