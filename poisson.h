/** \file
 *  Declarations of the classes that describe the stationary Phi equation
 *  in one dimension and the accompanying Dirichlet and Neumann boundary
 *  conditions.
 *
 *  $Id$
 *
 *  \author Bart Hartgers and Jan van Dijk
 */
	
#ifndef H_POISSON_H
#define H_POISSON_H

#include "basics.h"
#include "3diagsys.h"
#include "bndcond.h"
#include "grid.h"

class PhiVariable;

/** The stationary Phi-variable in one dimension.
 *
 */
class PoissonVariable : public Field
{
  public:
	/** Constructs a Phi-variable.
	 *  Firstly, this constructs the base class by passing the
	 *  number of nodal points as obtained from the \a grid argument.
	 *  Next, this member initializes the Field members that represent
	 *  the various transport and source coefficients. Some of these
	 *  are defined on the nodal grid, others on the ew-grid.
	 *  Also the references to the grid, to the generalised convective
	 *  flux and to the boundary conditions are stored in reference members.
	 *
	 *  The coefficients beta_cv, lambda, sc and sp are all zero-initialised,
	 *  the under-relaxation coefficient is initialised to 1.0 (no relaxation).
	 */
	PoissonVariable( const Grid& grid,
		const BoundaryCondition & l,
		const BoundaryCondition & r);

	// lambda, i.e. the conduction/diffusion coefficient. Defined on the np grid.
	Field lambda_cv;
	// constant source-term.
	Field rho;

	double Update(double dt, PhiVariable& ne, PhiVariable& ni);

	/// return the presently configured  relaxation factor.
	double urf() const { return m_urf; }
	/// set a new value for the relaxation factor.
	void set_urf(double urf) { m_urf=urf; }

	//set all values on the np to val
	const PoissonVariable& operator=(double val)
	{
		for (unsigned i=0; i<m_grid.num_np(); i++)
		{
			(*this)[i] = val;
		}
		return *this;
	}
	TridiagonalSystem& system() { return m_system; }
	const TridiagonalSystem& system() const { return m_system; }

	Field& flux_density() { return m_flux_dens; }
	const Field& flux_density() const { return m_flux_dens; }

	/** Returns a constant reference to the grid on which this Phi-variable
	 *  has been defined.
	 */
	const Grid& grid() const { return m_grid; }
	void DiscFluxDens(unsigned ndx_cv, double& DA, double& DB) const;
	double CalculateFluxDens(unsigned cv_ndx) const;
	void UpdateFluxDens();
  private:
	void AddVentzek(double dt, PhiVariable& nx, double qx);
	void Discretize();
	void ApplyRelaxation();

	/** Constant reference to the grid on which this Phi-variable
	 *  has been defined.
	 */
	const Grid& m_grid;
	// discretization coefficients.
	TridiagonalSystem m_system;
	// left and right boundary conditions.
	const BoundaryCondition& m_bc_w;
	const BoundaryCondition& m_bc_e;

	// relaxation factor.
	double m_urf;
	// flux density
	Field m_flux_dens;
};

#endif // H_POISSON_H

