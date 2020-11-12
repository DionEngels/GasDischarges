/** \file
 *
 *  Boltzmann Solver based on the phi-equation solver in phivar.*
 *
 */
  
#include <iostream>
#include <cassert>
#include <fstream>
#include <cmath>

#include "eedf.h"
#include "phivar.h"
#include "basics.h"
#include "argon.h"

// ********** numerical configuration **********

// lower boundary of the energy domain [J]
const double eps_min = 1e-3*PhysConst::e;

// upper boundary of the energy domain [J]
const double eps_max = 30*PhysConst::e;

// the number of points of the (discretised) energy grid [-]
const unsigned eps_size = 10000;

//maximum number of iterations [-]
const unsigned iter_max = 100;

// the tolerance (the required precision) [-]
const double tolerance = 1e-13;

// ********** physical configuration **********

// electric field strength [V/m]
const double Efield = 0;

// temperature of the background gas [K]
const double Tgas = 1*PhysConst::e/PhysConst::k_b; // 1: 1*eV

// density of the background gas [m^-3]
const double Ngas = 1e23;

// mass of the background gas particles [kg]
const double M = Argon::Mass();

const Argon::Elastic sigma_elas;

// initial guess for the `temperature' T==(2/3k)<eps>) of the electrons
// used for the initialisation of the EEDF with a Maxwellian function.
const double Te_init = Tgas;

// ********** derived auxiliary constants **********

const double E_N = Efield/Ngas;

const double m_M = PhysConst::me/M;

int main( int argc, char **argv )
{
	// Create a Cartesian grid. The coordinate represents the energy value.
	Grid grid( eps_size, eps_min, eps_max, Grid::Cartesian);

	// We use the variable flow_vel to represent the generalized 
	// velocity.
 	const Field flow_vel( grid.num_ew(), 1.0 );

	// Left boundary condition.
	NeumannBndCond bc_interp(0,0);
	// the value at the first internal grid point. This will
	// be updated to honor the normalisation condition.
	double value_guess = 1;

	// Right boundary condition: from integral considerations.
	// the integrated collision integral vanishes, this assumes
	// that no flux escapes at the right-hand side grid boundary.
	const double epsr=grid.pos_np(grid.num_np()-1);
	NeumannBndCond eedf_slope(0.0, eedf_conv(epsr,sigma_elas,m_M)/eedf_diff(epsr,sigma_elas,m_M,E_N,Tgas));

	// The Boltzmann equation is modelled as a Phi-equation
	PhiVariable eedf( grid, flow_vel, bc_interp, eedf_slope );
	//set the under-relaxation factor
	eedf.set_urf(1.0);

	// set the (guessed) starting value for the distribution function
	for( unsigned i=0; i<grid.num_np(); ++i)
	{
		const double eps = grid.pos_np(i);
		eedf[i] = eedf_maxwellian(eps,Te_init);
	}
	value_guess=eedf_renormalize(eedf);

	/*  the convective part. In our PhiVariable class the total flux is
	 *  written as beta*V. We keep flow_vel unity
	 *  and set beta*V by giving beta_cv an appropriate value.
	 *
	 *  Note that the convective coefficient beta_cv is defined on the
	 *  control volume boundary ('ew') grid.
	 */
	eedf.beta_cv = 1.0;
	for (unsigned i=0; i<grid.num_ew(); ++i )
	{	
		//const double eps = grid.pos_ew(i);
		//eedf.beta_cv[i] = ...;
	}

	/*  calculate the diffusion coefficient.
	 *
	 *  Note that the diffusion coefficient is defined on the
	 *  nodal ('np') grid.
	 */
	for( unsigned i=0; i<grid.num_np(); ++i)
	{
		//const double eps = grid.pos_np(i);
		//eedf.lambda[i] = ...;
	}

	// start our iterative procedure.
	unsigned iter = 0;
	double residue = 1.0;
	do 
	{
		++iter;

		/* In this implementation the coefficients and sources are
		 * constants, no updates are needed. In the more general case
		 * this would be the time to update any coefficient or
		 * source term that may have changed because we have just
		 * determined a new guess for the EDF.
		 * For example, this is the case when inelastic collisions
		 * are taken into account:
		 *
	 	for( unsigned i=0; i<grid.num_np(); ++i)
		{	
			eedf.sc[i] = ...
			eedf.sp[i] = ...
		}
		 */
		// override for i==1 (Patankar, page 145)
		const double huge=1e60;
		eedf.sc[1] = huge;
		eedf.sp[1] = -huge/value_guess;

		PhiVariable eedf_old(eedf);
		eedf.Update();
		//renormalize:
		double scale_factor = eedf_renormalize(eedf);
		// we adjust the (Dirichlet) boundary condition at
		// the lower energy boundary such that the
		// normalisation constraint is respected.
		value_guess *= scale_factor;
		const double slope=(std::log(eedf[2])-std::log(eedf[1]))*eedf[1]/grid.del();
		bc_interp.Set(-slope,0);
		const double epsr=grid.pos_np(grid.num_np()-1);
		eedf_slope.Set(0, eedf_conv(epsr,sigma_elas,m_M)/eedf_diff(epsr,sigma_elas,m_M,E_N,Tgas));

		residue = eedf_calculate_residual(eedf,eedf_old);
		std::cerr << "Iteration " << iter
			<< ", residue: " << residue
			<< ", value_guess: " << value_guess << '\n';

	} 
	while ( residue>tolerance && iter<iter_max );

	/// print the resulting EDF to the file eedf.dat
	std::ofstream eedf_stream("eedf.dat");
	grid.write(eedf_stream,eedf);

	// calculate the average electron energy. In case the field
	// is zero and no inelastic processes are present, this should
	// be equal to the average energy of the background gas,
	// (3/2)*kb*Tgas.
	const double en_av = eedf_average(eedf, grid.pos_np());

	// note that in equilibrium (3/2)kT = epsilon (with epsilon in J).
	// we use this here as a defining equation for the temperature, also
	// when we have a non-equilibrium EDF.
	const double Telec = en_av*(2.0/3.0)/PhysConst::k_b;
	std::cerr << "Field strength: " << Efield << std::endl;
	std::cerr << "Gas density:    " << Ngas << std::endl;
	std::cerr << "E/N:            " << Efield/Ngas << 
		 " Vm^2 (" << Efield/Ngas*1e21 << " Td)" << std::endl;
	std::cerr << "Average electron energy: " << en_av << 
		" (" << Telec << " K)" << std::endl;
	std::ofstream maxw("maxwell.dat");
	for (unsigned i=0;i<grid.num_np();++i)
	{
		const double eps=grid.pos_np(i);
		maxw << eps << '\t' << eedf_maxwellian(eps,Telec) << std::endl;
	}
	

	// print diagnostics and exit. We return the value 0 on success,
	// 1 if the simulation did not converge.

	if (residue < tolerance)
	{
		std::cerr << "Simulation converged.\n";
		return 0;
	}
	else
	{
		std::cerr << "Simulation did NOT converge.\n";
		return 1;
	}
}
  
