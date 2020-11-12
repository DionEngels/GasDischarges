#ifndef H_3DIAGSYS_H
#define H_3DIAGSYS_H

#include "field.h"

/** Class for managing systems of equations that are governed by a
 *  tridiagonal matrix.
 *
 *  The equation is assumed to be in the form
 *
 *     \f$a_p\Phi_P=\sum_{nb}a_{nb}\Phi_{nb}+b\f$
 *  
 *  The size of the system, which must be passed to the constructor
 *  of this class, is equal to the number of internal points, augmented
 *  with 2. 
 *
 *  \author Jan van Dijk
 */
class TridiagonalSystem
{
  public:
	/** Constructs a tridiagonal system. The argument \a size is the
	 *  number of points of this system.
	 */
	TridiagonalSystem(unsigned size)
	 : m_ap(size), m_aw(size), m_ae(size), m_b(size)
	{
	}
	/// returns a reference to the nodal coefficient at point \a i
	double& ap(unsigned i) {return m_ap[i];}
	/// returns a reference to the west coefficient at point \a i
	double& aw(unsigned i) {return m_aw[i];}
	/// returns a reference to the east coefficient at point \a i
	double& ae(unsigned i) {return m_ae[i];}
	/// returns a reference to the source at point \a i
	double& b(unsigned i) {return m_b[i];}

	/// returns the nodal coefficient at point \a i
	double ap(unsigned i) const {return m_ap[i];}
	/// returns the west coefficient at point \a i
	double aw(unsigned i) const {return m_aw[i];}
	/// returns the east coefficient at point \a i
	double ae(unsigned i) const {return m_ae[i];}
	/// returns the source at point \a i
	double b(unsigned i) const {return m_b[i];}
	/** Solve the tridiagonal matrix-vector problem using the TDMA algorithm.
	 *  See Patankar's textbook pages 52-55 for more information.
	 *  The argument \a f is the solution field.
	 *  Precondition: field.size()==size().
	 */
	void solve(Field& f) const;
	/// Returns the number of points of this system.
	unsigned size() const { return m_ap.size(); }
  private:
	// discretization coefficients.
	Field m_ap;
	Field m_aw;
	Field m_ae;
	Field m_b;
};

#endif // H_3DIAGSYS_H
