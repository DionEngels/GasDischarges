/** \file
 *
 *  Implementation of the code which has been declared in poisson.h
 *  $Id$
 */
	
#include "poisson.h"
#include "basics.h"
#include "grid.h"
#include "phivar.h"
#include "disc.h"

#include <iostream>
#include <cmath>

PoissonVariable::PoissonVariable( const Grid & grid,
		const BoundaryCondition & bc_w,
		const BoundaryCondition & bc_e)
 : Field( grid.num_np() ),
   lambda_cv( grid.num_ew(), 0.0 ),
   rho( grid.num_np(), 0.0 ),
   m_grid( grid ),
   m_system( grid.num_np() ),
   m_bc_w( bc_w ),
   m_bc_e( bc_e ),
   m_urf( 1.0 ),
   m_flux_dens( grid.num_ew(), 0.0 )
{
}

void PoissonVariable::DiscFluxDens(unsigned ndx_cv, double& DA, double& DB) const
{
	const double dLH = (ndx_cv==0 || ndx_cv==m_grid.num_ew()-1) ? m_grid.del()/2 : m_grid.del();
	DA =  lambda_cv[ndx_cv]/dLH;
	DB = DA;
}

double PoissonVariable::CalculateFluxDens(unsigned ndx_cv) const
{
	const unsigned ndx_L = ndx_cv;
	const unsigned ndx_H = ndx_L+1;
	double DA, DB;
	DiscFluxDens(ndx_cv,DA,DB);
	return DB*(*this)[ndx_L] - DA*((*this)[ndx_H]);
}

void PoissonVariable::UpdateFluxDens()
{
        for( unsigned i=0; i<m_grid.num_ew(); ++i)
        {
                m_flux_dens[i] = CalculateFluxDens(i);
        }
}

void PoissonVariable::Discretize()
{
	// we calculate the coefficients in the *internal* nodal points.
	for( unsigned p=1; p<m_grid.num_np()-1; ++p)
	{
		double DA,DB;
		m_system.ap( p )  = 0;

		unsigned ncv_w = p-1;
		const double areaW = m_grid.area_ew(ncv_w);
		DiscFluxDens(ncv_w,DA,DB);
		m_system.aw( p )  = DB*areaW;
		m_system.ap( p ) += DA*areaW;

		unsigned ncv_e = p;
		const double areaE = m_grid.area_ew(ncv_e);
		DiscFluxDens(ncv_e,DA,DB);
		m_system.ae( p )  = DA*areaE;
		m_system.ap( p ) += DB*areaE;
		
		const double vol = m_grid.vol_np(p);
		m_system.b(p)    = rho[p] * vol;
	}
	// the boundary conditions
	m_bc_w.ModifyCoefs(m_system, true,  m_grid.del()/2.0);
	m_bc_e.ModifyCoefs(m_system, false, m_grid.del()/2.0);
}

void PoissonVariable::ApplyRelaxation()
{
	if (urf()==1.0)
	{
		return;
	}
	for( unsigned p=1; p< size()-1; ++p)
	{
		m_system.ap(p) /= urf();
		m_system.b(p)  += (1.0-urf()) * m_system.ap(p) * (*this)[p];
	}
}

#if 0
double calc_eps_star_contrib(double q, const PhiVariable& dens, unsigned i_cv)
{
	const unsigned iL = i_cv;
	const unsigned iH = i_cv+1;
	double g1, g2, h;
	calcgh(g1, g2, h, dens.Pe_cv()[i_cv]);
	return q*dens.beta_cv[i_cv]*(g2*dens[iL]-g1*dens[iH]);
}

double calc_rho_star_contrib(double q, const PhiVariable& dens, unsigned i)
{
	const unsigned il = i-1;
	const unsigned ih = i;
	double h;
	// upper side
	calch(h, dens.Pe_cv()[ih]);
	const double dx_h = (i==dens.grid().num_np()-2) ? dens.grid().del()/2 : dens.grid().del();
	const double flux_diff_h = -dens.grid().area_ew(ih)*(dens.lambda_cv(ih)/dx_h)*h*(dens[i+1]-dens[i]);
	// lower side
	calch(h, dens.Pe_cv()[il]);
	const double dx_l = (i==1) ? dens.grid().del()/2 : dens.grid().del();
	const double flux_diff_l = -dens.grid().area_ew(il)*(dens.lambda_cv(il)/dx_l)*h*(dens[i]-dens[i-1]);
	return -q*(flux_diff_h-flux_diff_l)/dens.grid().vol_np(i);
}
//
//#define DIAGNOSE_VENTZEK
#if defined VENTZEK && defined DIAGNOSE_VENTZEK
		for (unsigned i=1; i<grid.num_np()-1; ++i)
		{
			const unsigned il = i-1;
			const unsigned ih = i;
			const double dtJint =
				+ eps_star[ih]*E[ih]*grid.area_ew(ih)-
				- eps_star[il]*E[il]*grid.area_ew(il);
			const double dtQ = rho_star[i]*grid.vol_np(i);
			std::cout << "ratio: " << (dtJint+dtQ) << '\t' << (rho[i]-rho_old[i])/dt << std::endl;
		}
#endif

#endif

// Update the field by calculation the discretization coefficients and
// solving the resulting matrix-vector problem. Returns the residual, i.e.
// a measure of how much this solution changed with respect to the previous
// one.
double PoissonVariable::Update(double dt, PhiVariable& ne, PhiVariable& ni)
{
	Discretize();

	if (dt!=0.0)
	{
		AddVentzek(dt, ne, -1*PhysConst::e);
		AddVentzek(dt, ni, +1*PhysConst::e);
	}

	ApplyRelaxation();

	const Field old(*this);

	m_system.solve(*this);

	UpdateFluxDens();

	// calculate the residual.
	double max_val=0, max_diff=0;

	for( unsigned p=0; p<size(); ++p)
	{
    
		const double diff = std::abs( old[p] - (*this)[p] );
		max_diff = std::max( diff, max_diff );

		const double val = std::abs( (*this)[p] );
		max_val = std::max( val, max_val );
	}
	return (max_val==0.0) ? max_diff : (max_diff/max_val);
}

void PoissonVariable::AddVentzek(double dt, PhiVariable& nx, double qx)
{
	return;
	const double qdt_dx=qx*dt/nx.grid().del();
	for (unsigned i=1; i!= nx.grid().num_ew()-1; ++i)
	{
		const unsigned iW = i;
		const unsigned iE = i+1;
		double mumid=nx.beta_cv[i];
		double Dmid=nx.lambda_cv(i);
		double g1,g2,h;
		calcgh(g1,g2,h,nx.Pe_cv()[i]);
		const double cf =qdt_dx*nx.grid().area_ew(i)*mumid*
				(nx[iW]*g1 - nx[iE]*g2);
		system().ae(iW)+=cf;
		system().ap(iW)+=cf;
		system().aw(iE)+=cf;
		system().ap(iE)+=cf;
		const double src=qdt_dx*nx.grid().area_ew(i)*Dmid*(nx[iE]-nx[iW])*h;
		system().b(iW)-=src;
		system().b(iE)+=src;
	}
};

