#include "3diagsys.h"
#include "basics.h"

void TridiagonalSystem::solve(Field& f) const
{
	if (f.size()!=m_ap.size())
	{
		npfFatalError("Solver: field does not have the correct dimensions.");
	}
	Field P( f.size() );
	Field Q( f.size() );

	const int first=0;
	const int last=f.size()-1;

	if (m_ap[first]==0.0)
	{
		npfFatalError("Solver error");
	}
	P[first] = m_ae[first]/m_ap[first];
	Q[first] =  m_b[first]/m_ap[first];
	for( int p=first+1; p <= last ; ++p)
	{
		const double denom = m_ap[p]-m_aw[p]*P[p-1];
		if (denom==0.0)
		{
			 npfFatalError("Solver error");
		}
		P[p] = m_ae[p]/denom;
		Q[p] = (m_b[p]+m_aw[p]*Q[p-1])/denom;
	}

	f[last]=Q[last];
	for( int p=last-1; p>=first; --p)
	{
		f[p] = P[p]*f[p+1] + Q[p];
	}
}

