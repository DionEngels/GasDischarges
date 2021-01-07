// File:	pic.cpp
// Author:	W.J.M. Brok <wjmb@etpmod.phys.tue.nl>
// Date:	januari 2004
// Version:     Id: pic.cpp,v 1.3 2004/10/30 20:04:21 wjmb Exp
// Description:	together with phivar.h and monte_carlo.h this example
//		implements a very basic one-dimensional Particle In Cell
//		(PIC) model including a Monte Carlo (MC) treatment of the
//		collisions. This example sends gnuplot command to std::cout
//              to make use of this feature pipe the results to gnuplot
//               $./example_pic | gnuplot
//              or on MS Windows open a command window and use
//              c:\npf-practicum> example_pic | c:\gnuplot\binaries\pgnuplot.exe
//              (assuming this code is built in c:\npf-practicum and
//              gnuplot is installed in c:\gnuplot). Gnuplot may be obtained
//		from ftp://ftp.gnuplot.info/pub/gnuplot/ more information
//              can be found on http://www.gnuplot.info/
#include <fstream>
#include <iostream>
#include <sstream>

#include "argon.h"
#include "monte_carlo.h"
#include "phivar.h"
#include "poisson.h"

// Some settings and constants:
const double initial_e_nr = 40;      // init. nr. of simulation electrons.
const double initial_i_nr = 40;      // init. nr. of simulation ions.
const double particleweight = 1e+10; // real per simulation particle.
const double gridsize = 0.01;        // meter between electrodes.
const unsigned nrpoints = 200;       // nodal points inside plasma.
const double Vs = 2000.0;            // voltage of the power supply.
const double R = 100;                // Ohm, series resistor.
const double N = 1e23;               // m^-3 (background density).
const double sec_emiss_coeff = 0.07; // electrons per incident ion.
const double dt = 5e-13;             // seconds, particle path integr. time.
const double simtime = 1e-6;         // seconds, complete simulation time.
const double dt_Iurf = 1e-7;         // seconds, current relaxation time.
const double dt_out = 1e-9;          // seconds, between writing output.
const double mi = Argon::Mass();     // kg.

// A function to write some data: this is implemented at the bottom of
// this file.
void WritePoisson(const Grid &grid, const Field &poisson, const Field &rho_e,
                  const Field &rho_i, std::ostream &os);

// Ambipolar diffusion coefficient (implented towards the end of this file).
double Ambi_Diff_Coeff(double Te, double Th);

// The main() function: this is where you can see the outline of the
// algorithm that is implemented in this code.
int main() {
  // A grid with nrpoints, extending from 0.0 to gridsize.
  Grid grid(nrpoints, 0.0, gridsize, Grid::Cartesian);

  // We have no flow when solving the poisson equation.
  Field flow_vel(grid.num_ew(), 0.0);

  // Densities of eletrons and ions, they are required by the member function
  // Update() of class PoissonVariable.
  DirichletBndCond density_bc(0.0);
  PhiVariable n_e(grid, flow_vel, density_bc, density_bc);
  PhiVariable n_i(grid, flow_vel, density_bc, density_bc);

  // Boundary conditions for the Poisson equation: note, we use a
  // voltage source, and might need a current source or couple this
  // to a circuit equation in order to get reasonable behaviour.
  DirichletBndCond poisson_left(-Vs);
  DirichletBndCond poisson_right(0.0);

  // We define charge densities on a grid, but we don't need to solve
  // a phi-equation for them: no need for boundary conditions here,
  // however, something needs to be defined.
  DirichletBndCond rho_bc(0.0);

  // The different things defined on the grid. Again, mind that we only
  // solve the equation for the potential.
  PoissonVariable poisson(grid, poisson_left, poisson_right);
  Field rho_e(grid.num_np());
  Field rho_i(grid.num_np());

  // A swarm is a collection of particles of the same kind. We have
  // a swarm for the electrons and for the argon ions. For now we just
  // construct them, but only fill them with particles later, when we
  // know more about processes and crosssections.
  Swarm electrons(PhysConst::me, -PhysConst::e, particleweight);
  Swarm ions(mi, PhysConst::e, particleweight);

  // We construct processes. The first argument to create a process is a
  // function which returns the cross section as a function of the
  // energy. These can be found in the file argon.h. The second argument
  // is optional and specifies the energy the electron looses in the
  // collision. If this is omitted it will be set to 0.0;
  HardSphereProcess e_elastic(new Argon::Elastic);
  HardSphereProcess e_inelastic(new Argon::Inelastic, 12.0 * PhysConst::eV);
  HardSphereProcess e_ionisation(new Argon::Ionisation, 15.8 * PhysConst::eV);
  ChargeExchangeProcess i_chargeexchange(new Argon::ChargeExchange);
  ChargeExchangeProcess i_elastic(new Argon::IonElastic);

  // For convenience later: pointers to the processes in which electrons
  // are involved are collected in an object of type ProcessList. This
  // will allow us to do something like e_proclist[0] to get `e_elastic'.
  ProcessList e_proclist;
  e_proclist.AddProcess(&e_elastic);
  e_proclist.AddProcess(&e_inelastic);
  e_proclist.AddProcess(&e_ionisation);

  ProcessList i_proclist;
  i_proclist.AddProcess(&i_chargeexchange);
  i_proclist.AddProcess(&i_elastic);

  // Null-collision: in order to set a collision time, it is convenient
  // if the total (of all processes) cross section is independent of the
  // electron/ion energy. We accomplish this by finding the maximum cross
  // section and using this to determine the collision time. Later we
  // correct for this overestimation by allowing for a null-collision.
  double e_max_crosssec = 0.0;
  double i_max_crosssec = 0.0;
  for (unsigned i = 0; i < 1000; ++i) {
    const double energy = 0.1 * i * PhysConst::eV;
    // electrons:
    const double e_sum = e_elastic.CrossSec(energy) +
                         e_inelastic.CrossSec(energy) +
                         e_ionisation.CrossSec(energy);
    if (e_sum > e_max_crosssec)
      e_max_crosssec = e_sum;
    // ions:
    const double i_sum =
        i_chargeexchange.CrossSec(energy) + i_elastic.CrossSec(energy);
    if (i_sum > i_max_crosssec)
      i_max_crosssec = i_sum;
  }
  // for safety, in case we just missed the maximum.
  e_max_crosssec *= 1.04;
  i_max_crosssec *= 1.04;

  // The particles in the swarms need to be initialised with a
  // position (1 dimensional) and a velocity (3 dimensional). We
  // do this with random numbers. A generator is made like this:
  Random rnd;

  // we calculate the current.
  double I = 0.0;

  // We add ...
  for (unsigned p = 0; p < initial_e_nr; ++p) {
    // An electron.
    const double place = gridsize * rnd();
    const double eps = 4.0 * PhysConst::eV;
    npfGeomVector speed(0.0, 3);
    speed[0] = rnd();
    speed[1] = rnd();
    speed[2] = rnd();
    speed *= std::sqrt(2.0 * eps / electrons.Mass());
    const double dtcol =
        -(1.0 / (N * e_max_crosssec * abs(speed))) * std::log(rnd());
    electrons.AddParticle(place, speed, dtcol);
  }
  for (unsigned p = 0; p < initial_i_nr; ++p) {
    // An ion.
    const double place = gridsize * rnd();
    const double eps = 0.5 * PhysConst::eV;
    npfGeomVector speed(0.0, 3);
    speed[0] = rnd();
    speed[1] = rnd();
    speed[2] = rnd();
    speed *= std::sqrt(2.0 * eps / ions.Mass());
    const double dtcol =
        -(1.0 / (N * i_max_crosssec * abs(speed))) * std::log(rnd());
    ions.AddParticle(place, speed, dtcol);
  }
  // At this point one could check if dt is small enough !!

  // ... every nrout iterations we want to update output.
  const unsigned nrout = (int)(dt_out / dt);

  // Okay, all is set, let's roll: finally the real PIC loop:
  for (unsigned long ti = 0; ti < simtime / dt; ++ti) {
    // the present time:
    const double t = ti * dt;

    // we will count the charges coming in on the two electrodes
    // during this timestep.
    double Qleft = 0.0;
    double Qright = 0.0;

    // 1) Let's interpolate the charges to a grid:
    electrons.CalcRho(grid, rho_e);
    ions.CalcRho(grid, rho_i);

    // 2) With this information we can solve the Poisson
    //    equation:
    for (unsigned i = 0; i < grid.num_np(); i++) {
      poisson.rho[i] = (rho_e[i] + rho_i[i]);
      n_e[i] = rho_e[i] / PhysConst::e;
      n_i[i] = rho_i[i] / PhysConst::e;
    }

    for (unsigned i = 0; i < grid.num_ew(); i++) {
      poisson.lambda_cv[i] = PhysConst::epsilon0;
    }

    poisson.Update(dt, n_e, n_i);

    // Move all the particles and check for collisions,
    // either with the wall or with the background gas.
    unsigned nrelectrons = electrons.size();
    for (unsigned p = 0; p < nrelectrons; ++p) {

      // 3e) determine the local field at the particle's
      //    position:
      const double E = CalcEfield(grid, poisson, electrons[p].Position());
      // 4e) move the particle
      electrons[p].Move(E, dt);
      // 5e) did it collide with the wall?
      if (electrons[p].Position() < 0.0) {
        Qleft += electrons.Charge() * electrons.Weight();
        // the electron is lost.
        electrons[p].Active() = false;
      }
      if (electrons[p].Position() > gridsize) {
        Qright += electrons.Charge() * electrons.Weight();
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
            ions.AddParticle(electrons[p].Position(), speed, t);
          }
        }
        // ii) set the particles new time of
        //     collision. BUG: vabs==0?
        double vabs = abs(electrons[p].Velocity());
        double dtcol =
            -(1.0 / (N * e_max_crosssec * vabs)) * std::log(1.0 - rnd());
        electrons[p].Tcol() += dtcol;
      }
    }
    electrons.RemoveInactiveParticles();

    // do the same thing for the ions.
    unsigned nrions = ions.size();
    for (unsigned p = 0; p < nrions; ++p) {
      // 3i) determine the local field at the particle's
      //    position:
      const double E = CalcEfield(grid, poisson, ions[p].Position());
      // 4i) move the particle
      ions[p].Move(E, dt);
      // 5i) did it collide with the wall?
      if (ions[p].Position() < 0.0) {
        Qleft += ions.Charge() * ions.Weight();
        // the ion is lost.
        ions[p].Active() = false;
        // however, we might have created an electron:
        if (rnd() < sec_emiss_coeff) {
          npfGeomVector speed(0.0, 3);
          double place = 0.0;
          electrons.AddParticle(place, speed, t);
          Qleft -= electrons.Charge() * electrons.Weight();
        }
      }
      if (ions[p].Position() > gridsize) {
        Qright += ions.Charge() * ions.Weight();
        // the ion is lost.
        ions[p].Active() = false;
        // however, we might have created an electron:
        if (rnd() < sec_emiss_coeff) {
          npfGeomVector speed(0.0, 3);
          double place = gridsize;
          electrons.AddParticle(place, speed, t);
          Qright -= electrons.Charge() * electrons.Weight();
        }
      }
      // 6i) determine if it needs to collide with the
      //     background gas.
      if (ions[p].Active() && (ions[p].Tcol() <= t)) {
        // i) the particle is going to collide, but
        //    which process is going to happen?
        double rndsc = rnd() * i_max_crosssec;
        double sum = 0.0;
        unsigned process = 0;
        for (unsigned i = 0; i < i_proclist.size(); ++i) {
          sum += i_proclist[i]->CrossSec(ions[p].Energy());
          if (rndsc > sum)
            process++;
        }
        if (process < i_proclist.size()) {
          // Make a dummy particle (which is
          // just standing still).
          Particle d(Argon::Mass(), 0.0);
          i_proclist[process]->Collide(ions[p], d);
        }
        // ii) set the particles new time of
        //     collision. BUG: vabs==0?
        double vabs = abs(ions[p].Velocity());
        double dtcol =
            -(1.0 / (N * i_max_crosssec * vabs)) * std::log(1.0 - rnd());
        ions[p].Tcol() += dtcol;
      }
    }
    ions.RemoveInactiveParticles();

    // the current (which we try to relax a bit).
    const double urf = dt / dt_Iurf;
    I = (1.0 - urf) * I + urf * (-Qright) / dt;
    // and the resulting voltage across the discharge.
    double Vdischarge = Vs - R * I;
    // let's build in a limit.
    if (Vdischarge < 0.0)
      Vdischarge = 0.0;
    // set the voltage.
    poisson_left.Set(-Vdischarge);

    // some diagnostics:
    const double ne_average = electrons.size() * electrons.Weight() / gridsize;
    const double wp = std::sqrt((ne_average * PhysConst::e * PhysConst::e) /
                                (PhysConst::epsilon0 * PhysConst::me));
    const double kTe = 2.0 * electrons.MeanEnergy() / 3.0;
    const double Ld = std::sqrt((PhysConst::epsilon0 * kTe) /
                                (ne_average * PhysConst::e * PhysConst::e));
    const double v_mean =
        std::sqrt(2.0 * electrons.MeanEnergy() / PhysConst::me);
    const double nsigmav = N * e_max_crosssec * v_mean;

    // write some data:
    if (!(ti % nrout))
      std::cerr << ti << '\t' << t << '\t' << electrons.size() << '\t'
                << ions.size() << '\t' << I << '\t' << Vdischarge << '\t' << wp
                << '\t' << Ld << '\t' << nsigmav << '\t' << v_mean * dt << '\t'
                << std::endl;

    // A little check: we might as well stop if there are no
    // electrons or no ions left.
    if (electrons.size() == 0 || ions.size() == 0) {
      npfFatalError("No particles left anymore");
    }
  }

  // convert time (double) to string (via a stream)
  std::ostringstream time_str;
  time_str << simtime;
  // plot
  grid.plot(poisson, rho_e + rho_i, "Position (cm)", "Potential (V)",
            "Charge density (C m^-3)", "time = " + time_str.str() + " s");

  // final diagnostics, assuming equilibrium conditions at the end of
  // simulation
  double av_E = electrons.MeanEnergy();        // obtain average electron energy
  double av_v = abs(electrons.MeanVelocity()); // obtain average velocity vector
  double sigma_tot =
      e_elastic.CrossSec(av_E) + e_inelastic.CrossSec(av_E) +
      e_ionisation.CrossSec(
          av_E); // total cross section at average electron energy
  double l_mfp =
      1.0 / (N * sigma_tot);       // mean free path for e-neutral collisions
  double tau_col = l_mfp / (av_v); // collision time

  std::cout << "Mean free path: " << l_mfp << "\tCollision time: " << tau_col
            << std::endl;

  return 0;
}

void WritePoisson(const Grid &grid, const Field &poisson, const Field &rho_e,
                  const Field &rho_i, std::ostream &os) {
  // print the solution of the poisson equation
  for (unsigned i = 0; i < grid.num_np(); i++) {
    os << i << '\t' << grid.pos_np(i) << '\t' << rho_e[i] << '\t' << rho_i[i]
       << '\t' << poisson[i] << std::endl;
  }
  os << std::endl << std::endl;
}

double Ambi_Diff_Coeff(double Te, double Th) {
  /// \todo Static members are not nice
  static Argon::IonElastic ar_ion_elas;
  // thermal Velocity
  const double vth =
      std::sqrt(PhysConst::pi * Argon::Mass() * Th / (8.0 * PhysConst::k_b));
  const double AmbiFact = (1.0 + Te / Th);
  // ion-atom (induced dipole) cross secion for argon.
  const double n0_ia_crs = ar_ion_elas(Th * PhysConst::k_b);
  return vth * AmbiFact * PhysConst::k_b / (N * Argon::Mass() * n0_ia_crs);
}
