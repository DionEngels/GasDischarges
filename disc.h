#ifndef H_DISC_H
#define H_DISC_H

#include <cassert>
#include <cmath>
#include <algorithm>

/** f1: GJH eqn 3.18
 *  f2: GJH eqn 3.19
 */
void calcf(double& f1, double& f2, double z);

/** g1: GJH eqn 3.41
 *  g2: GJH eqn 3.42
 *   h: GJH eqn 3.34
 */
void calcgh(double& g1, double& g2, double& h, double z);

/** h: GJH eqn 3.34
 *  \todo JvD - h calculation has been implemented as part of g1,g2 as well. ???
 */
void calch(double& h, double z);

namespace disc {
inline double A(double P)
{
	assert (P>=0);
	return P==0 ? 1 : P/(std::exp(P)-1);
}
inline void AB(const double P, double& A, double& B)
{
	calcf(A,B,P);
return;
	const double APabs = disc::A(std::abs(P));
	A = APabs + std::max(-P,0.0);
	B = APabs + std::max( P,0.0);
}

} // namespace disc


#endif // define H_DISC_H

