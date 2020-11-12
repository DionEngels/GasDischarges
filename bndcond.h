#ifndef H_BNDCOND_H
#define H_BNDCOND_H

#include "field.h"
#include "3diagsys.h"

/** Abstract base-class for boundary conditions.
 *  This defines the interface that is shared by all bounday conditions.
 *  The boundary condition can be applied to a Phi-problem by calling
 *  the member ModifyCoefs. This calls the protected member DoModifyCoefs. The latter
 *  member is abstract and must be provided by derived classes.
 *
 *  \author Jan van Dijk
 */
class BoundaryCondition
{
  public:
	/** C++ memory management requires that a class has a virtual destructor
	 *  if it has other virtual members, even if it does nothing otherwise.
	 */
	virtual ~BoundaryCondition(){}
	/** Modifies the coefficients of a tridiagonal system at the indicated
	 *  boundary to reflect this boundary condition.
	 *
	 *  \param sys     A reference to the system matrix that needs to be modified.
	 *  \param low_bnd Indicates the side: true means the lower boundary,
	 *                 false the upper boundary.
	 *  \param dx      The distance [m] between the boundary point and the
	 *                 adjacent internal point.
	 */
	void ModifyCoefs(TridiagonalSystem& sys, bool low_bnd, double dx) const
	{
		const unsigned pbnd = low_bnd ? 0 : sys.size()-1;
		double& ax = low_bnd ? sys.ae(pbnd) : sys.aw(pbnd);
		DoModifyCoefs(sys.b(pbnd), sys.ap(pbnd), ax, dx);
	}
  protected:
	/** This member must be overridden in a derived class to do the
	 *  boundary condition modification in a particular way.
	 *
	 *  \param b      The coefficient b (the source term)
	 *  \param a_p    The discretisation coefficient coefficient a_p
	 *  \param a_side The discretisation coefficient a_w or a_e, depending
	 *                on whether this member is called for the lower or upper side.
	 *  \param dx     The distance [m] between the boundary point and the
	 *                adjacent internal point.
	 */
	virtual void DoModifyCoefs(
			double& b,
			double& a_p, 
			double& a_side, 
			double dx) const=0;
};

/** A class for Dirichlet boundary conditions, meaning that the boundary value
 *  is specified.
 *
 *  The boundary value is passed to the constructor and stored in a private
 *  member variable. This value may be changed later by calling member Set.
 *
 *  The virtual member DoModifyCoefs has been overridden to set the boundary value
 *  and modify the coefficients according to the Dirichlet scheme (see syllabus).
 *
 *  \author Jan van Dijk
 */
class DirichletBndCond : public BoundaryCondition
{
  public:
	/** Constructs a Dirichlet boundary condition object.
	 *  Argument \a bval is the boundary value.
	 */
	DirichletBndCond(double bval) 
	 : BoundaryCondition(), m_bval(bval)
	{}
	/// Set the boundary value to the new value \a bval.
	void Set(double bval)
	{
		m_bval=bval; 
	}
  protected:
	/** Overridden to do the Dirichlet coefficient modification.
	 *  For a description of the parameters, see the base class
	 *  documentation.
	 */
	virtual void DoModifyCoefs(
			double& b,
			double& a_p, 
			double& a_side, 
			double dx) const;
  private:
	double m_bval;
};

/** A class for Neumann boundary conditions, meaning that the normal component
 *  of the gradient at the boundary is specified.
 *
 *  The component of the boundary gradient is supposed to be linearised as
 *  \f$\partial f/\partial{\bf n} = F_0 + f_{bnd}F_1 \f$. The F-coeficients are
 *  passed to the constructor and stored in two private member variable. The
 *  value may be changed later by calling member Set.
 *
 *  The virtual member DoModifyCoefs has been overridden to update the boundary value
 *  and modify the coefficients according to the Neumann scheme (see syllabus).
 *
 *  \author Jan van Dijk
 */
class NeumannBndCond : public BoundaryCondition
{
  public:
	/** Constructs a neumann boundary condition object.
	 *  Arguments \a f0 and f1 parametrise the boundary gradient
	 *  that is assumed to be of the form
	 *
	 *    \f$\left[\frac{df}{dn}\right]_{bnd}=f0 + f_{bnd}f_1.\f$
	 */
	NeumannBndCond(double f0, double f1)
	 : BoundaryCondition(), m_F0(f0), m_F1(f1)
	{}
	/// Set the boundary coefficients  the new valus \a f0 and \a f1.
	void Set(double f0, double f1)
	{
		m_F0=f0;
		m_F1=f1; 
	}
  protected:
	/** Overridden to do the Neumann coefficient modification.
	 *  For a description of the parameters, see the base class
	 *  documentation.
	 */
	virtual void DoModifyCoefs(
			double& b,
			double& a_p, 
			double& a_side, 
			double dx) const;
  private:
	double m_F0;
	double m_F1;
};

class ExternalBndCond : public BoundaryCondition
{
  public:
	ExternalBndCond()
	 : BoundaryCondition()
	{}
  protected:
	virtual void DoModifyCoefs(
			double& b,
			double& a_p, 
			double& a_side, 
			double dx) const
	{
	}
};

#endif // H_BNDCOND_H

