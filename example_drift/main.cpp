// File:	example_drift.cpp
// Author:	W.J.M. Brok <wjmb@etpmod.phys.tue.nl>, Mark Beks
// Date:	januari 2005
// Version:     $Id$
// Description:	together with monte_carlo.h this example represents
//		a Monte Carlo (MC) model of a drift tube
//              This example sends gnuplot command to std::cout
//              to make use of this feature pipe the results to gnuplot
//               $./example_drift| gnuplot
//              or on MS Windows open a command window and use
//              c:\npf-practicum> example_drift |
//              c:\gnuplot\binaries\pgnuplot.exe (assuming this code is built in
//              c:\npf-practicum and gnuplot is installed in c:\gnuplot).
//              Gnuplot may be obtained
//		from ftp://ftp.gnuplot.info/pub/gnuplot/ more information
//              can be found on http://www.gnuplot.info/

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

#include "argon.h"
#include "monte_carlo.h"
#include "phivar.h"

// Some settings and constants:
const double initial_e_nr = 2000;   // init. nr. of simulation electrons.
const double particleweight = 1e+9; // real per simulation particle.
const unsigned nrpoints = 50;       // nodal points inside plasma.
const double Efield = 100;          // V/m Electric field in the drift tube
const double N = 4e21;              // m^-3 (background density).
const double dt = 1e-11;            // seconds, particle path integr. time.
const double simtime = 1e-6;        // seconds, complete simulation time.
const double dt_Iurf = 1e-6;        // seconds, current relaxation time.
const double dt_out = 5e-9;         // seconds, between writing output.
const double gridsize = 1.0;        // m the size of the drift chamber

// The main() function: this is where you can see the outline of the
// algorithm that is implemented in this code.
int main() {
  // Grid to keep track of the physical position of the particles
  Grid grid(nrpoints, 0.0, gridsize, Grid::Cartesian);
  // We have no flow when solving the poisson equation.
  Field flow_vel(grid.num_ew(), 0.0);
  // We define charge densities on a grid, but we don't need to solve
  // a phi-equation for them: no need for boundary conditions here,
  // however, something needs to be defined.
  DirichletBndCond rho_bc(0.0);

  // A swarm is a collection of particles of the same kind. We have
  // a swarm for the electrons only. For now we just
  // construct it, but only fill it with particles later, when we
  // know more about processes and crosssections.
  Swarm electrons(PhysConst::me, -PhysConst::e, particleweight);

  // We construct processes. The first argument to create a process is a
  // function which returns the cross section as a function of the
  // energy. These can be found in the file argon.h. The second argument
  // is optional and specifies the energy the electron looses in the
  // collision. If this is omitted it will be set to 0.0;
  HardSphereProcess e_elastic(new Argon::Elastic);
  HardSphereProcess e_inelastic(new Argon::Inelastic, 12.0 * PhysConst::eV);
  HardSphereProcess e_ionisation(new Argon::Ionisation, 15.8 * PhysConst::eV);
  // For convenience later: pointers to the processes in which electrons
  // are involved are collected in an object of type ProcessList. This
  // will allow us to do something like e_proclist[0] to get `e_elastic'.
  ProcessList e_proclist;
  e_proclist.AddProcess(&e_elastic);
  e_proclist.AddProcess(&e_inelastic);
  e_proclist.AddProcess(&e_ionisation);
  // Null-collision: in order to set a collision time, it is convenient
  // if the total (of all processes) cross section is independent of the
  // electron/ion energy. We accomplish this by finding the maximum cross
  // section and using this to determine the collision time. Later we
  // correct for this overestimation by allowing for a null-collision.
  double e_max_crosssec = 0.0;
  for (unsigned i = 0; i < 1000; ++i) {
    const double energy = 0.1 * i * PhysConst::eV;
    // electrons:
    const double e_sum = e_elastic.CrossSec(energy) +
                         e_inelastic.CrossSec(energy) +
                         e_ionisation.CrossSec(energy);
    if (e_sum > e_max_crosssec)
      e_max_crosssec = e_sum;
  }
  // for safety, in case we just missed the maximum.
  e_max_crosssec *= 1.04;
  // The particles in the swarms need to be initialised with a
  // position (1 dimensional) and a velocity (3 dimensional). We
  // do this with random numbers. A generator is made like this:
  Random rnd;
  // We add ...
  for (unsigned p = 0; p < initial_e_nr; ++p) {
    // An electron at 1cm from the anode
    const double place = 1.0 - 1e-2 * gridsize;
    const double eps = 4.0 * PhysConst::eV;
    npfGeomVector speed(-0.5, 3); // with random speeds
    speed[0] += rnd();
    speed[1] += rnd();
    speed[2] += rnd();
    speed *= std::sqrt(2.0 * eps / electrons.Mass());
    // set collision times
    const double dtcol =
        -(1.0 / (N * e_max_crosssec * abs(speed))) * std::log(rnd());
    electrons.AddParticle(place, speed, dtcol);
  }
  // At this point one could check if dt is small enough !!

  // ... every nrout iterations we want to update output.
  const unsigned nrout = (int)(dt_out / dt);

  // give some commands common for all plots (logscale and yrange)
  std::cout << "set logscale y" << std::endl;
  std::cout << "set yrange [" << particleweight * 25 << ":"
            << initial_e_nr * particleweight * 25 << "]" << std::endl;
  // we keep track of the positions of the electrons with rho_e
  PhiVariable rho_e(grid, flow_vel, rho_bc, rho_bc);
  // const double inv_charge= - 1.0/ PhysConst::e;
  // running average of the drift velocity
  double vdrift = 0.0;
  double relaxation = 1e-4;
  // Okay all is set, let's roll: finally the real PIC loop:
  for (unsigned long ti = 0; ti < simtime / dt; ++ti) {
    // the present time:
    const double t = ti * dt;

    // Move all the particles and check for collisions,
    // either with the wall or with the background gas.
    unsigned nrelectrons = electrons.size();
    for (unsigned p = 0; p < nrelectrons; ++p) {
      // 4e) move the particle
      electrons[p].Move(Efield, dt);
      // 5e) did it collide with the wall?
      if (electrons[p].Position() < 0.0) {
        // the electron is lost.
        electrons[p].Active() = false;
      }
      if (electrons[p].Position() > gridsize) {
        // the electron is lost.
        electrons[p].Active() = false;
      }
      // 6e) determine if it needs to collide with the
      //     background gas.
      if (electrons[p].Active() && (electrons[p].Tcol() <= t)) {
        // i) the particle is going to collide, but
        //    which process is going to happen?
        double rndsc = rnd() * e_max_crosssec;
        double sum = 0.0;
        unsigned process = 0;
        for (unsigned i = 0; i < e_proclist.size(); ++i) {
          sum += e_proclist[i]->CrossSec(electrons[p].Energy());
          if (rndsc > sum)
            process++;
        }
        if (process < e_proclist.size()) {
          // Make a dummy particle (which is
          // just standing still).
          Particle d(Argon::Mass(), 0.0);
          e_proclist[process]->Collide(electrons[p], d);
          if (process == 2) {
            npfGeomVector speed(0.0, 3);
            electrons.AddParticle(electrons[p].Position(), speed, t);
          }
        }
        // ii) set the particles new time of
        //     collision:
        //     wjmb, BUG: vabs==0? rnd()==0?
        double vabs = abs(electrons[p].Velocity());
        double dtcol = -(1.0 / (N * e_max_crosssec * vabs)) * std::log(rnd());
        electrons[p].Tcol() += dtcol;
      }
    }
    electrons.RemoveInactiveParticles();
    // get the average velocity in the x-direction
    double vz = electrons.MeanVelocity()[0];
    // apply relaxation to the running average
    vdrift = (1 - relaxation) * vdrift + relaxation * vz;
    // write some data to std::cerr so this is shown on screen:
    if (!(ti % nrout)) {
      std::cerr << std::setw(20) << t << std::setw(20) << vz << std::setw(20)
                << vdrift << std::setw(20) << electrons.size() << std::endl;
      // find out where the electrons are  and plot them
      electrons.CalcRho(grid, rho_e);
      // define an extra field to plot the results
      Field inv_charge(rho_e.size());
      for (unsigned i = 0; i < inv_charge.size(); i++) {
        inv_charge[i] = -rho_e[i] / PhysConst::e;
      }
      // define a stream for double to string conversion of the time
      std::ostringstream time_str;
      time_str << t;
      // plot the charge density
      grid.plot(inv_charge, "Position (cm)", "Particle density (m^-3)",
                "Time = " + time_str.str() + " s");
    }
    // A little check: we might as well stop if there are no
    // electrons or no ions left.
    if (electrons.size() == 0) {
      npfFatalError("No particles left anymore");
    }
  }
  return 0;
}
