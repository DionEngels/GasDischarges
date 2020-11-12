/** \file
 *
 *  Implementation of the code which has been declared in monte_carlo.h
 *  $Id$
 */
	
#include "monte_carlo.h"
#include "basics.h"

EnergyTable::EnergyTable(unsigned nrpoints, double maxenergy)
{
	for (unsigned i=0;i<nrpoints;++i)
	{
		m_x.push_back(i*maxenergy/nrpoints);
		m_y.push_back(0.0);
	}
}

void Swarm::AddParticle(double pos, const npfGeomVector& vel, double Tcol)
{
	Particle newprtcl(m_mass, m_charge, pos, vel, Tcol);
	push_back( newprtcl );
}

void Swarm::RemoveInactiveParticles()
{
	// this could be done a lot neater.
	const unsigned original_number = size();
	for (unsigned i=0;i<original_number;++i)
	{
		const unsigned index = original_number-i-1;
		if (!(*this)[index].Active())
		{
			std::vector<Particle>::iterator itRemove;
			itRemove = begin() + index;
			// now, actually remove this element from the array
			erase(itRemove);
		}
	}
}

void Swarm::CalcRho(const Grid& grid, Field& rho)
{
	const double dx = (grid.pos_np(grid.num_np()-1))/(grid.num_np()-2);
	rho = 0.0;
	for (const_iterator p=begin(); p!=end(); ++p)
	{
		// i) first find out between which nodal points the
		//    particle is located.
		const double pos = p->Position();
		const unsigned j = (int)(pos/dx + 0.5);
		// ii) then we distribute it between the two 
		//     nearest grid cells.
		rho[j  ] += (grid.pos_np(j+1)-pos)* m_weight*Charge()/(dx*dx);
		rho[j+1] +=-(grid.pos_np(j  )-pos)* m_weight*Charge()/(dx*dx);
	}
	// because of the specifics of the grid (nodal point 0 and nrpoints+1
	// are really fake points) we need to do a little correction.
	rho[0]=0.0;
	rho[1] *= 4.0/3.0;
	rho[grid.num_np()-1]=0.0;
	rho[grid.num_np()-2] *= 4.0/3.0;
}

const EnergyTable Swarm::EDF() const
{
	// first determine the maximum energy:
	double maxenergy = 0.0;
	for (const_iterator p=begin(); p!=end(); ++p)
	{
		if (maxenergy < p->Energy())
		{
			maxenergy = p->Energy();
		}
	}
	// then construct an energytable:
	EnergyTable edf(100,maxenergy);
	// finaly bin the particles:
	for (const_iterator p=begin(); p!=end(); ++p)
	{
		const unsigned j = (int)(100.0*p->Energy()/maxenergy);
		if (j<edf.Value().size())
		{
			 edf.Value()[j]++;
		}
	}
	// and return the table:
	return edf;
}
	
void Swarm::WriteEDF(std::ostream& os) const
{
	// print the eedf.
	EnergyTable eedf = EDF();
	for (unsigned i=0; i<eedf.size();++i)
	{
		os << eedf.Energy()[i] << '\t'
			<< eedf.Value()[i] << std::endl;
	}
	os << std::endl << std::endl;
}


const double Swarm::MeanEnergy() const
{
	double sum = 0.0;
	for (const_iterator p=begin(); p!=end(); ++p)
	{
		sum += p->Energy();
	}
	// calculating the size of a list is expensive
	return sum/size();
}


npfGeomVector Swarm::MeanVelocity() const
{
	npfGeomVector v(0.0,3);
	for (const_iterator p=begin(); p!=end(); ++p)
	{
		v += p->Velocity();
	}
	v /= size();
	return v;
}




Field Swarm::MeanEnergyField(const Grid& grid) const
{
	Field energy(grid.num_np());
	const double dx = (grid.pos_np(grid.num_np()-1))/(grid.num_np()-2);
	std::vector<double> nrincell;
	for (unsigned i=0;i<grid.num_np();++i)
	{
		nrincell.push_back(0.0);
	}
	for (unsigned p=0; p<size();++p)
	{
		// i) first find out between which nodal points the
		//    particle is located.
		const unsigned j = (int)((*this)[p].Position()/dx + 0.5);
		// ii) then we distribute its energy between the two 
		//     nearest grid cells.
		energy[j] += (grid.pos_np(j+1)-(*this)[p].Position())*								(*this)[p].Energy()/dx;
		nrincell[j] += (grid.pos_np(j+1)-(*this)[p].Position())/dx;
		energy[j+1] += ((*this)[p].Position()-grid.pos_np(j))*
						(*this)[p].Energy()/dx;
		nrincell[j+1] += ((*this)[p].Position()-grid.pos_np(j))/dx;
	}
	// because of the specifics of the grid (nodal point 0 and nrpoints+1
	// are really fake points) we need to do a little correction.
	energy[0]=0.0;	energy[grid.num_np()-1]=0.0;
	for (unsigned i=1;i<grid.num_np()-1;++i)
	{
		// if there is no particle in the cell, the mean energy
		// is set to zero. This is of course not correct, but how
		// can one do this differently?
		if (nrincell[i]!=0) energy[i] /= nrincell[i];
		else energy[i] = 0.0;
	}
	return energy;
}

void Swarm::WriteMeanEnergyField(const Grid& grid, 
					std::ostream& os) const
{
	Field energies = MeanEnergyField(grid);
	for (unsigned i=0; i<grid.num_np(); i++)
	{
		os << i << '\t' << grid.pos_np(i) << '\t' 
			<<  energies[i]/PhysConst::eV << std::endl;
	}
	os << std::endl << std::endl;
}

Process::Process(const CrossSection* csp, double delta_eps) 
	: m_delta_eps(delta_eps), m_csp(csp), m_nr_called(0)
{
	assert(csp);
}

void Process::Collide(Particle& p1, Particle& p2)
{
	m_nr_called++;
	DoCollide(p1,p2);
}

void Process::RotateVel(npfGeomVector& v, double chi, double psi) const
{
	// length of projection on x-y plane.
	const double rho = sqrt(v[0]*v[0]+v[1]*v[1]);
	// length of complete vector v.
	const double vabs = abs(v);
	// theta is in the range [0,pi).
	const double theta = (rho == 0.0 ) ? ( (v[2] > 0) ? 0 : -PhysConst::pi )
						: PhysConst::pi/2 - atan(v[2]/rho);
	// phi is in the range [-pi,pi).
	const double phi = (v[0] == 0.0 ) ? 0.0 : 2*atan(v[1]/v[0]);
	// find the cosines and sines of the theta and phi.
	const double costheta = cos(theta);
	const double cosphi = cos(phi);
	const double sintheta = sin(theta);
	const double sinphi = sin(phi);
	// find the sines and cosines of chi and psi.
	const double coschi = cos(chi);
	const double cospsi = cos(psi);
	const double sinchi = sin(chi);
	const double sinpsi = sin(psi);
	// do the rotation.
	v[0] = vabs * ( sintheta*cosphi*coschi + 
			sinchi*(costheta*cosphi*cospsi-sinphi*sinpsi) );
	v[1] = vabs * ( sintheta*sinphi*coschi + 
			sinchi*(costheta*sinphi*cospsi+cosphi*sinpsi) );
	v[2] = vabs * ( costheta*coschi - sintheta*sinchi*cospsi );
}

void HardSphereProcess::DoCollide(Particle& p1, Particle& p2) const
{
	// first calculate the center of mass properties.
	const double mtotal = p1.Mass()+p2.Mass();
	const npfGeomVector u_M = ( p1.Mass()*p1.Velocity() +
			p2.Mass()*p2.Velocity() ) / mtotal;
	npfGeomVector u_r = p1.Velocity() - p2.Velocity();
	// then draw two random scattering angles.
	const double chi = 2.0*PhysConst::pi*m_rnd();
	const double psi = std::acos(1.0-2.0*m_rnd());
	// turn the relative velocity over these angles.
	RotateVel(u_r,chi,psi);
	// if the collision is inelastic, it changes the size of
	// the relative velocity. We first check if this is 
	// possible and then do it.
	const double u_r_abssqr = u_r[0]*u_r[0]+u_r[1]*u_r[1]+u_r[2]*u_r[2];
	if (u_r_abssqr*p1.Mass()*p2.Mass()/mtotal < 2.0*m_delta_eps)
	{
		// Don't stop here, just go on with p1 an p2 unchanged. If
		// this does not happen too often, then the mistake we make
		// by ignoring this is small.
		return;
	}
	u_r *= sqrt( 1.0 - 
		m_delta_eps*(2.0*mtotal/(p1.Mass()*p2.Mass()))/u_r_abssqr );
	// now add the center of mass velocity again to obtain
	// the new velocities of the particles.
	p1.Velocity() = u_M + (p2.Mass()/mtotal)*u_r;
	p2.Velocity() = u_M - (p1.Mass()/mtotal)*u_r;
}

void ProcessList::AddProcess(Process* prcss)
{
	push_back(prcss);
}

void ChargeExchangeProcess::DoCollide(Particle& p1, 
					Particle& p2) const
{
	// Do a little check:
	if (!( ((p1.Charge()==0.0)&&(p2.Charge()!=0.0)) ||
	     ((p1.Charge()!=0.0)&&(p2.Charge()==0.0)) ))
	{
		npfFatalError("One particle needs to have a charge.");
	}
	// The charge jumps over from one particle to the next. p1
	// is assumed to be charged, p2 is not. Effectively the 
	// particles exchange their velocities:
	std::swap(p1.Velocity(),p2.Velocity());
}

Random::Random()
{
	initialise();
	m_seed = NewSeed();
	sgenrand(m_seed);
}

Random::~Random()
{
	delete [] m_mt;
}

void Random::initialise()
{
	m_N = 624;
	m_M = 397;
	/* constant vector a */
	m_A = 0x9908b0df;
	/* masks ... maybe these should be inserted directly
	 * where they are used. That could make the code a
	 * bit faster I suppose. */   
	m_upper_mask = 0x80000000;
	m_lower_mask = 0x7fffffff;
	m_tempering_mask_B = 0x9d2c5680;
	m_tempering_mask_C = 0xefc60000;
	/* allocate memory */
	m_mt = new unsigned long[m_N];
	m_nr = 0;
}

void Random::sgenrand(unsigned long seed)
{
	/* setting initial seeds to m_mt[N] using the generator Line 25 of 
	 * Table 1 in [KNUTH 1981, The Art of Computer Programming Vol. 2 
	 * (2nd Ed.), pp102] 
	 * The random number generator is thus bootstrapped by another one. */
 
	m_mt[0]= seed & 0xffffffff;
	for (m_mti=1; m_mti < m_N; ++m_mti) 
	{
		/* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier. */
		/* In the previous versions, MSBs of the seed affect   */
		/* only MSBs of the array mt[].                        */
		m_mt[m_mti] = (1812433253 * 
				(m_mt[m_mti-1] ^ (m_mt[m_mti-1] >> 30)) 
				+ m_mti);
		m_mt[m_mti] &= 0xffffffff;
	}
}

const double Random::operator()()
{
	/* Distribution on [0,1]: */
//	return calcnextI()*(1.0/4294967295.0); // divided by 2^32-1

	/* Distribution on [0,1): */
	return calcnextI()*(1.0/4294967296.0); // divided by 2^32
 
}

unsigned long Random::calcnextI()
{
	/* This routine generates N random numbers at once and then uses
	 * those until it runs out of them ... if so it will generate a
	 * new batch */
	unsigned long y;
	unsigned long mag01[2] = {0x0, m_A};

	if (m_mti >= m_N) 
	{ 
		unsigned long kk;
		for (kk=0; kk < m_N - m_M; ++kk) 
		{
			y = (m_mt[kk] & m_upper_mask) | (m_mt[kk+1] 
							& m_lower_mask);
			m_mt[kk] = m_mt[kk+m_M] ^ (y >> 1) ^ mag01[y & 0x1];
		}
		for ( ; kk < m_N-1; ++kk) 
		{
			y = (m_mt[kk] & m_upper_mask) | (m_mt[kk+1] 
							& m_lower_mask);
			m_mt[kk] = m_mt[kk+(m_M-m_N)] ^ (y >> 1) 	
							^ mag01[y & 0x1];
		}
		y = (m_mt[m_N-1] & m_upper_mask) | (m_mt[0] & m_lower_mask);
		m_mt[m_N-1] = m_mt[m_M-1] ^ (y >> 1) ^ mag01[y & 0x1];
		m_mti = 0;
	}
	y = m_mt[m_mti++];
	y ^= (y >> 11);
	y ^= (y << 7) & m_tempering_mask_B;
	y ^= (y << 15) & m_tempering_mask_C;
	y ^= (y >> 18);
	m_nr++;
	return y;
}

