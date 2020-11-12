#include "bndcond.h"

void DirichletBndCond::DoModifyCoefs (
			double& b,
			double& a_p, 
			double& a_side, 
			double dx ) const
{
	b = m_bval;
	a_p=1;
}

void NeumannBndCond::DoModifyCoefs (
			double& b,
			double& a_p, 
			double& a_side, 
			double dx ) const
{
	const double d0 = m_F0*dx;
	const double d1 = m_F1*dx;
	a_p = 1-d1;
	a_side=1;
	b = d0;
}


