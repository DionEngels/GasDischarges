#include "disc.h"
#include <cmath>

void calcf(double& f1, double& f2, double z)
{
	if (z<10)
	{
		if ((z>1E-6)||(z<-1E-6))
		{
			const double expz=std::exp(z);
			f1=z/(expz-1);
		}
		else f1=1;
	}
	else
	{
		const double expz=std::exp(-z);
		f1=z*(1/(1-expz)-1);
	}
	f2=f1+z;
}

void calcgh(double& g1, double& g2, double& h, double z)
{
	double zz;

	if (z<10)
	{
		if ((z>1E-6)||(z<-1E-6))
		{
			const double expz=std::exp(z);

			zz=1.0/(expz-1.0);
			zz*=zz;
			g1=zz*((1.0-z)*expz-1.0);
			h=z*z*expz*zz;
		}
		else
		{
			g1=-1.0;
			h=1.0;
		}
		g2=g1+1.0;
	}
	else
	{
		const double expz=std::exp(-z);
		zz=1.0/(expz-1.0);
		zz*=zz;
		g2=zz*(1.0-(1.0+z)*expz);
		h=z*z*expz*zz;
		g1=g2-1.0;
	}
}

void calch(double& h, double z)
{
	double zz;

	if (z<10)
	{
		if ((z>1E-6)||(z<-1E-6))
		{
			const double expz=std::exp(z);
			zz=1.0/(expz-1.0);
			zz*=zz;
			h=z*z*expz*zz;
		}
		else
		{
			h=1.0;
		}
	}
	else
	{
		const double expz=std::exp(-z);
		zz=1.0/(expz-1.0);
		zz*=zz;
		h=z*z*expz*zz;
	}
}

