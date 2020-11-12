/** \file
 *  Declarations of the classes that describe the stationary Phi equation
 *  in one dimension and the accompanying Dirichlet and Neumann boundary
 *  conditions.
 *
 *  $Id$
 *
 *  \author Bart Hartgers and Jan van Dijk
 */
	
#ifndef H_PHIVAR_H
#define H_PHIVAR_H

#include "basics.h"
#include "3diagsys.h"
#include "bndcond.h"
#include "grid.h"

/** The stationary Phi-variable in one dimension.
 *
 */
class PhiVariable : public Field
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
	PhiVariable( const Grid& grid,
		const Field& w_cv,
		const BoundaryCondition & l,
		const BoundaryCondition & r);

	// beta, i.e. the convection coefficient. Defined on the ew-grid.
	Field beta_cv;
	// lambda, i.e. the conduction/diffusion coefficient. Defined on the np grid.
	Field lambda;
	// constant source-term.
	Field sc;
	// linear source-term.
	Field sp;

	// interpolate lambda from the nodal to the ew grid.
	// recicprocal interpolation is used.
	double lambda_cv(unsigned ndx_cv) const;

	double Update();

	/// return the presently configured  relaxation factor.
	double urf() const { return m_urf; }
	/// set a new value for the relaxation factor.
	void set_urf(double urf) { m_urf=urf; }

	//set all values on the np to val
	const PhiVariable& operator=(double val)
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

	// the Peclet numbers on the ew-grid.
	const Field& Pe_cv() const { return m_Pe_cv; }


	/** Returns a constant reference to the grid on which this Phi-variable
	 *  has been defined.
	 */
	const Grid& grid() const { return m_grid; }
	void DiscFluxDens(unsigned ndx_cv, double& DA, double& DB) const;
	double CalculateFluxDens(unsigned cv_ndx) const;
	void UpdateFluxDens();
  private:
	void Discretize();
	void ApplyRelaxation();

	/** Constant reference to the grid on which this Phi-variable
	 *  has been defined.
	 */
	const Grid& m_grid;
	// a reference to the generalised convection velocity
	const Field& m_w_cv;
	// discretization coefficients.
	TridiagonalSystem m_system;
	// left and right boundary conditions.
	const BoundaryCondition& m_bc_w;
	const BoundaryCondition& m_bc_e;

	// relaxation factor.
	double m_urf;
	// flux density
	Field m_flux_dens;
	mutable Field m_Pe_cv;
};

#endif // H_PHIVAR_H

