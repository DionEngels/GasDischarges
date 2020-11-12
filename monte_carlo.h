// File:	monte_carlo.h
// Author:	W.J.M. Brok <wjmb@etpmod.phys.tue.nl>
// Date:	januari 2004
// Version:     Id: monte_carlo.h,v 1.16 2004/08/31 14:22:41 wjmb Exp
// Description:	This file defines and implements some basic functionality
//		for constructing Monte Carlo based models. Several assumptions
//		are made, most important the ones concerning dimension: one
//		dimension in configuration space, one in velocity space.
 
#ifndef H_MONTE_CARLO_H
#define H_MONTE_CARLO_H

#include <valarray>
#include <vector>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <cassert>
#include "basics.h"
#include "grid.h"
#include "crosssec.h"

/** The type of the Cartesian vector we use to represent positions,
 *  velocity et cetera
 */
typedef std::valarray<double> npfGeomVector;

/// Returns the norm of the npfGeomVector \a v
inline double abs(const npfGeomVector& v)
{
	return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]);
}

/// Container in which we store an energy table.
class EnergyTable
{
  public:
	// Constructor: needs the number of points and the maximum energy.
	// Note that the values still need to be defined after contruction !!
	EnergyTable(unsigned nrpoints, double maxenergy);
	// accessors.
	const std::vector<double>& Energy() const { return m_x; }
	const std::vector<double>& Value() const { return m_y; }
	std::vector<double>& Value() { return m_y; }
	// the number of points in our table:
	unsigned size() const { return m_x.size(); }
  private:
	std::vector<double> m_x;
	std::vector<double> m_y;
};

// A particle. In this definition it is assumed that we are using it for
// a 1d place and 3d velocity MC simulation. Note that this particle 
// has its own mass and charge variables: if you have 10000 particles this
// is a gross waste of memory space. Typically one would like to make a
// class `Species' which keeps this information for a particle swarm of 
// alike particles and have each particle have a reference to this Species
// object.
class Particle {
  public:
	// Constructor: you need to create this particle with a mass, charge,
	// initial position and initial velocity.
	Particle(double mass, double charge, double position, const npfGeomVector& velocity,
		 double Tcol)
		: m_mass(mass), m_charge(charge), m_active(true),
		  m_pos(position), m_vel(velocity), m_Tcol(Tcol) {};
	// Constructor which sets all the extrinsic properties of the 
	// particle to zero.
	Particle(double mass, double charge)
		: m_mass(mass), m_charge(charge), m_active(true),
		  m_pos(0.0), m_vel(0.0,3), m_Tcol(0.0) {};
	// Move the particle: calculate new position and velocity for this 
	// particle under the influence of an electric field Efield during
	// a time-interval dt.
	void Move(double Efield, double dt)
	{
		const double force = Efield*m_charge;
		m_pos += dt*(m_vel[0]+force*(dt/(2.0*m_mass)));
		m_vel[0] += force*dt/m_mass;
	}
	// const accessor for the particle's mass.
	const double Mass() const { return m_mass; }
	// const accessor for the particle's charge.
	const double Charge() const { return m_charge; }
	// const accessor for the particle's position.
	const double Position() const { return m_pos; }
	// const accessor for the particle's position.
	double& Position() { return m_pos; }
	// const accessor for the particle's velocity (3d cartesian vector).
	const npfGeomVector& Velocity() const { return m_vel; }
	// non-const accessor for a reference to the velocity.
	npfGeomVector& Velocity() { return m_vel; }
	// time of the next collision.
	const double Tcol() const { return m_Tcol; }
	// reference to the time of the next collision.
	double& Tcol() { return m_Tcol; }
	// flag which tells if the particle is still active.
	const bool Active() const { return m_active; }
	// reference to flag which tells if the particle is still active.
	bool& Active() { return m_active; }
	// calculates the kinetic energy (in Joule) of the particle.
	const double Energy() const
	{ 
		return 0.5*m_mass*( m_vel[0]*m_vel[0] +
					m_vel[1]*m_vel[1] +
				  	m_vel[2]*m_vel[2] );
	}
  private:
	// intrinsic properties:
	double m_mass;
	double m_charge;
	bool m_active;
	// extrinsic properties:
	double m_pos;
	npfGeomVector m_vel;
	double m_Tcol;
};

// Swarm: collection of particles of the same species.
class Swarm : public std::vector<Particle>
{
  public:
	// constructor: mass and charge of the particles contained in this
	// swarm have to be given as arguments, as well as the weight of the
	// particles in the swarm.
	Swarm(double mass, double charge, double weight) 
		: m_mass(mass), m_charge(charge), m_weight(weight) {};
	// adds a particle to the swarm. This needs to know the place and 
	// velocity of the newly-born, as well as the moment in time it will
	// experience its first collision.
	void AddParticle(double pos, const npfGeomVector& vel, double Tcol);
	// const accessor for the particles' mass.
	const double Mass() const { return m_mass; }
	// const accessor for the particles' charge.
	const double Charge() const { return m_charge; }
	// const accessor for the particles' weight (the number of real
	// particles it represents).
	const double Weight() const { return m_weight; }
	// calculates the energy distribution function of the particles in
	// this swarm and returns it in the form of a table:
	const EnergyTable EDF() const;
	// Writes the energy distribution function to ostream.
	void WriteEDF(std::ostream& os) const;
	// calculates the mean energy of the particles in the swarm.
	const double MeanEnergy() const;
	npfGeomVector MeanVelocity() const;
	// Calculates the mean energy of the particles on grid.
	Field MeanEnergyField(const Grid& grid) const;
	// Writes the mean energy of particles on the grid to stream.
	void WriteMeanEnergyField(const Grid& grid, 
				std::ostream& os) const;
	// Removes all the particles in this swarm of which Active() is false.
	void RemoveInactiveParticles();
	// Calculates the space charge density rho of the particles on the 
	// grid.
	void CalcRho(const Grid& grid, Field& rho);
  private:
	// remove particle with index from the swarm.
	void RemoveParticle(unsigned index);
	double m_mass;
	double m_charge;
	double m_weight;
};

/** Declaration of Mersenne Twister MT19937ra random number generator 
 *  class, adapted from http://www.math.keio.ac.jp/matumoto/emt.html */
class Random
{
  public:
	/** Default constructor: when this is called for the first time
	 *  within a project, the seed will be set at 4357. Every
	 *  subsequent creation of a Random object will have its seed
	 *  incremented compared to the previous one. If you want to be
	 *  sure to have on specific seed, use the next constructor */
	Random(); 
	/** Destructor. */
	~Random();
	/** Generate new random number, probably on [0,1),
	 *  have a look in the implemenation to be sure. */
	const double operator()();
	/** Return the first seed used. */
	const unsigned long Seed() const { return m_seed; }
	/** Returns the number of random numbers generated. */
	const unsigned long Nr() const { return m_nr; }
  private:
	/** This function is called for by the default constructor
	 *  in order to provide each instance of the class with a 
	 *  different seed. */
	unsigned long NewSeed() {
		static unsigned long seed=4357; 
		// with this seed a random number is produced which
		// functions as the real seed for Random. Incrementing
		// this seed should provide sufficient independence
		// between subsequent instances of this class.
		return seed++;
	}
	/** Copy constructor: we do not allow coppies of Random */
	Random(const Random& r); 
	/** Assignment operator: we do not allow assignments either. */
	Random& operator=(const Random& r); 
	/** Calculate the next integer in the sequence. */
	unsigned long calcnextI();
	/** Set constants to their standard values */
	void initialise();
	/** initializing the array with a NONZERO seed */
	void sgenrand(unsigned long seed);
	/** Period parameters */  
	unsigned long m_N;
	unsigned long m_M;

	unsigned long m_A;          /* constant vector a */
	unsigned long m_upper_mask; /* most significant w-r bits */
	unsigned long m_lower_mask; /* least significant r bits */

	unsigned long m_tempering_mask_B;
	unsigned long m_tempering_mask_C;
	
	unsigned long m_nr;

	unsigned long *m_mt; 
	unsigned long m_mti; 

	unsigned long m_seed;
};

class Process
{
  public:
	// Constructor: takes as an argument the amount of energy that is
	// lost if the collision is inelastic. For elastic collisions this 
	// should be 0.0.
	Process(const CrossSection* csp, double delta_eps);
	// The way the compiler likes it: a virtual destructor.
	virtual ~Process()
	{
		delete m_csp;
	}
	// Returns the cross section at energy eps (in Joule).
	const double CrossSec(double eps) const { return (*m_csp)(eps); }
	// Collide particle p1 and p2 with eachother.
	void Collide(Particle& p1, Particle& p2);
	// Returns the number of times Collide() has been called.
 	const unsigned NrCalled() const { return m_nr_called; }
  protected:
	// Pure virtual function which collides particle p1 with p2. This 
	// needs to be implemented in a derived call.
	virtual void DoCollide(Particle& p1, Particle& p2) const=0;
	// function which can turn a relative velocity vector v over 
	// an angle phi and theta.
	void RotateVel(npfGeomVector& v, double chi, double psi) const;
	// energy lost in inelastic collision.
	double m_delta_eps;
  private:
	const CrossSection* m_csp;
	unsigned m_nr_called;
};

class HardSphereProcess: public Process
{
  public:
	// Constructor.
	HardSphereProcess(const CrossSection* csp, double delta_eps=0.0) 
		: Process(csp, delta_eps), m_rnd() {};
  protected:
	// Do the actual collision between two particles.
	void DoCollide(Particle& p1, Particle& p2) const;
  private:
	mutable Random m_rnd;
};

class ProcessList: public std::vector<Process*>
{
  public:
	ProcessList() {};
	void AddProcess(Process* prcss);
};

class ChargeExchangeProcess: public Process
{
  public:
	// Constructor.
	ChargeExchangeProcess(const CrossSection* csp) : Process(csp, 0.0) {};
	// Do the actual collision between two particles.
  protected:
	void DoCollide(Particle& p1, Particle& p2) const;
};

// Function to calculate the Electric field at pos. We do not interpolate,
// so there is room for improvement.
inline double CalcEfield(const Grid& grid, 
			const Field& V, 
			double pos)
{
	const double dx = (grid.pos_np(grid.num_np()-1))/(grid.num_np()-2);
	const unsigned j = (int)(pos/dx+0.5);
	return -(V[j+1]-V[j])/(grid.pos_np(j+1)-grid.pos_np(j));
}

#endif // H_MONTE_CARLO_H
