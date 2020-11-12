#ifndef H_EEDF_H
#define H_EEDF_H

#include "physconst.h"
#include "basics.h"
#include "field.h"
#include "phivar.h"
#include "crosssec.h"
#include <cmath>

// The Maxwell distribution function [(s/m)^3] as a function of energy [J]
inline double eedf_maxwellian(double eps, double T)
{
	const double kT=PhysConst::k_b*T;
	const double factor = std::pow(PhysConst::me/(2*PhysConst::pi*kT),1.5);
	return factor*std::exp(-eps/kT);
}

inline double eedf_1_sqrteps_dP(const PhiVariable& eedf, unsigned i)
{
	const double factor = 4*PhysConst::pi*std::sqrt(2.0)*std::pow(PhysConst::me,-1.5);
	return factor*eedf[i]*eedf.grid().dx_np(i);
} 

inline double eedf_dP(const PhiVariable& eedf, unsigned i)
{
	return std::sqrt(eedf.grid().pos_np(i))*eedf_1_sqrteps_dP(eedf,i);
} 

/** returns the integral \f$ \int f(\epsilon)\sqrt(\epsilon)d\epsilon \f$
 *  where f is the field \a field.
 */
inline double eedf_average(const PhiVariable& eedf, const Field& f)
{
	double integral = 0.0;
	for (unsigned i=0;i<eedf.grid().num_np();++i)
	{
		integral += f[i]*eedf_dP(eedf,i);
	}
	return integral;
};

/** Normalise the \a eedf. The required scale factor is returned.
 */
inline double eedf_renormalize(PhiVariable& eedf)
{
	double prob = 0.0;
	for (unsigned i=0;i<eedf.grid().num_np();++i)
	{
		prob += eedf_dP(eedf,i);
	}

	const double scale_factor = 1.0/prob;
	for (unsigned i=0;i<eedf.size();++i)
	{
		eedf[i] *= scale_factor;
	}
	return scale_factor;
}

/** calculates the difference between the fields \a f1 and \a f2.
 *  where this difference is the norm ||f1* -f2*||_inf 
 *  with f1* f1 rescaled according to f1*[i]= f1[i]/min{f1[i],f2[i]} 
 *  values with min{f1,f2}=0 are skipped.
 */
inline double eedf_calculate_residual(const Field& f1, const Field& f2)
{
	assert(f1.size()==f2.size());
	double res=0.0;
	for (unsigned p=0; p< f1.size();++p)
	{
		const double diff  = std::abs( f1[p] - f2[p] );
		const double scale = std::min(std::abs(f1[p]),std::abs(f2[p]));
		if (scale!=0.0)
			res = std::max( res, diff/scale );
	}
	return res;
}

// The convection coefficient beta_cv as a function of energy [J]
inline double eedf_conv(double eps, const CrossSection& sig_elas, const double m_M)
{
	return -4*sqr(eps)*m_M*sig_elas(eps);
}

// The generalised diffusion coefficient lambda as a function of energy [J]
inline double eedf_diff(double eps, const CrossSection& sig_elas, const double m_M, double E_N, double Tg)
{
	const double sigm = sig_elas(eps);
	return (2.0/3.0)*eps*sqr(PhysConst::e*E_N)/sigm
		+ 4*sigm*m_M*sqr(eps)*PhysConst::k_b*Tg;
}

#endif // H_EEDF_H

