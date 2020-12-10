#include <fstream>
#include <iostream>
#include <sstream>


#include "argon.h"
#include "disc.h"
#include "lutable.h"
#include "phivar.h"
#include "poisson.h"


// Cartesian gradient on ew-grid from np-data
void gradient_np_cv(const Grid &grid, const Field &f, Field &gradf) {
  const unsigned N = grid.num_ew();

  // i=0: forward difference
  gradf[0] = (f[1] - f[0]) / (grid.del() / 2);
  // 1 < i < N-1: central difference
  for (unsigned i = 1; i < N - 1; ++i) {
    gradf[i] = (f[i + 1] - f[i]) / grid.del();
  }
  // i=N-1: backward difference
  gradf[N - 1] = (f[N] - f[N - 1]) / (grid.del() / 2);
}

// physical configuration:

const double gridsize = 30e-2;   // [m]  the electrode separation
const double pressure = 1 * 133; // [Pa] the gas pressure
const double Tgas = 300;         // [K]  the gas temperature

const double dens_ei_init = 1e12; // [m^-3] initial electron and ion density

// numerical configuration:

const double dt_min = 1e-15;     // [s] lowest allowed timestep
const unsigned gridpoints = 100; // the number of grid points
double tolerance =
    5e-2; // the allowed relative change of a variables after a time step
const unsigned iter_log =
    250; // number of iterations after which info is written to screen

#define DC_SIMULATION
#ifdef DC_SIMULATION
const double dt_max = 1e-8; // [s] highest allowed timestep
const double simtime = 10;  // [s] complete simulation time.
const double dt_out = 1e-5; // [s] time between writing output.

double Vwest(double t) {
  const double Vampl = -400; // V
  return Vampl;
}

#else // set up an AC simulation

const double freq = 1 * 1e6; // [Hz]
const unsigned Nminstepspercycle = 50;
const unsigned Ncycles = 10000;
const unsigned NoutPerCycle = 1000;

double Vwest(double t) {
  const double Vampl = -400;                     // [V]
  const double omega = 2 * PhysConst::pi * freq; // [Hz*rad]
  return Vampl * std::cos(omega * t);
}

const double period = 1.0 / freq;                                 // [s]
const double dt_max = std::min(1e-8, period / Nminstepspercycle); // [s]
const double simtime = period * Ncycles; // [s], complete simulation time.
const double dt_out =
    period / (NoutPerCycle * freq); // [s], between writing output.

#endif // DC_SIMULATION

// derived physical quantities
const double eps_i = (3.0 / 2.0) * PhysConst::k_b * Tgas;  // [J]
const double gasdens = pressure / (PhysConst::k_b * Tgas); // [m^-3]

// gas properties:
//
// The code below implemnts Helium properties

const double mgas = PhysConst::AMU * 4;

// the mean electron energy as a function of the reduced field
double eps_e(double E_N) {
  static mdLookupTable lut("lut/epsilon.dat", 1e-21, PhysConst::e);
  return lut.lookup(E_N);
}

// the electron mobility as a function of the reduced field
double mu_e(double E_N) {
  static mdLookupTable lut("lut/muN.dat", PhysConst::e, 1);
  return -lut.lookup(eps_e(E_N)) / gasdens;
}

// the electron diffusion coefficient as a function of the reduced field
// Note: this uses the Einstein relation D_e = mu_e*kT_e/q_e = -(2/3)eps_e/e.
double D_e(double E_N) {
  const double kT_e = (2.0 / 3.0) * eps_e(E_N);
  const double q_e = -PhysConst::e;
  return mu_e(E_N) * kT_e / q_e;
}

// the helium ion mobility in helium.
// See: H.W.  Ellis et al., Atomic Data & Nucl Data Tab, 17,177 (1976)
double mu_i(double E_N) {
  const double xm = 1e-21;
  const double ym = 2.69e19 / 0.01;
  // the table contains muN(E/N) in the units xm, ym given above.
  static mdLookupTable lut("lut/mup-ion.dat", xm, ym);
  const double muN = (lut.lookup(E_N)) / gasdens;
  return muN;
}

// the He ion diffusion coefficient as a function of the reduced field
// Note: this uses the Einstein relation D_p = mu_p*kT_p/q_p.
//       We assume T_p=Tgas. See Ellis (1976) for the definition of
//       an effective temperature that is better used instead.
double D_i(double E_N) {
  const double q_p = +PhysConst::e;
  const double kT_p = PhysConst::k_b * Tgas;
  return mu_i(E_N) * kT_p / q_p;
}

// the ionisation rate coefficient [m^3/s] as a function of the reduced field
double Kion(double E_N) {
  // the table lists K(eps), with eps in eV.
  static mdLookupTable lut("lut/Kion.dat", PhysConst::e, 1);
  return lut.lookup(eps_e(E_N));
}

// helper function.
//
// The arguments are:
//  * a phi variable that we wish to update,
//  * the old values of that field
//  * the highest relative change that we allow.
//
// The function returns false if the tolerance is exceeded,
// true otherwise. The caller can decide what to do when the
// tolerance is exceeded. (In the present program, we undo all
// changes to the variables, reduce the time step and try again.)
bool Solve(PhiVariable &phivar, const Field &old, double tolerance) {
  phivar.Update();
  for (unsigned i = 1; i < phivar.grid().num_np() - 1; ++i) {
    const double diff = fabs(phivar[i] - old[i]);
    const double rel_diff = diff / fabs(phivar[i]);
    if (rel_diff > tolerance) {
      return false;
    }
  }
  return true;
}

double Vth_e_E(double E_N) {
  const double kT_e = (2.0 / 3.0) * eps_e(E_N);
  const double Vth = std::sqrt(8 * kT_e / (PhysConst::pi * PhysConst::me));
  // std::cout << Vth << std::endl;
  return Vth;
}

double Vth_i() {
  const double kT_i = (2.0 / 3.0) * eps_i;
  const double Vth = std::sqrt(8 * kT_i / (PhysConst::pi * mgas));
  // std::cout << Vth << std::endl;
  return Vth;
}

// The main() function: this is where you can see the outline of the
// algorithm that is implemented in this code.
int main() {
  // Grid to keep track of the physical position of the particles
  Grid grid(gridpoints, 0.0, gridsize, Grid::Cartesian);
  DirichletBndCond hom_dirichlet_bc(0.0);

  Field J_cv(grid.num_ew(), 0.0);
  Field sion(grid.num_np(), 0.0);

  DirichletBndCond V_west_bc(0.0);
  PoissonVariable V(grid, V_west_bc, hom_dirichlet_bc);
  V.lambda_cv = PhysConst::epsilon0;
  Field E(grid.num_ew(), 0.0);

  ExternalBndCond neumann_bc_i_west;
  ExternalBndCond neumann_bc_i_east;
  PhiVariable dens_i(grid, E, neumann_bc_i_west, neumann_bc_i_east);

  const double gam = 0.2;
  ExternalBndCond neumann_bc_e_west;
  ExternalBndCond neumann_bc_e_east;
  PhiVariable dens_e(grid, E, neumann_bc_e_west, neumann_bc_e_east);

  // to this file we will write mean values of some
  // quantities at the time intervals given by dt_out.
  std::ofstream meanvals("meanvals.dat");

  // set initial values:

  dens_i = dens_ei_init;
  dens_e = dens_ei_init;

  // the initial time step
  double dt = dt_min;
  // the initial time
  double t = 0;
  // the time we want to write files
  double tout = 0;
  // the iteration number (number of successful time steps)
  long iter = 0;
  // true if we begin the simulation, or after undoing a time step
  bool restarting = true;

  // now simulation the plasma until the simulation time has exired.
  //
  while (t < simtime) {
    double mean_ne = 0, mean_ni = 0;
    for (unsigned i = 0; i < grid.num_np(); ++i) {
      mean_ne += dens_e[i] * grid.vol_np(i);
      mean_ni += dens_i[i] * grid.vol_np(i);
    }
    mean_ne /= grid.volume();
    mean_ni /= grid.volume();

    if ((iter % iter_log) == 0) {
      std::cout << "iter: " << iter << ", t=" << t << ", dt=" << dt
                << ", <ni>=" << mean_ni << ", <ne>=" << mean_ne << std::endl;
    }

    // store the old values of V, ne and ne.
    const Field V_old(V);
    const Field dens_e_old(dens_e);
    const Field dens_i_old(dens_i);

    // Calculate the potential on the western cathode and
    // inform the Dirichlet boundary condition class about it.
    const double Vw = Vwest(t);
    V_west_bc.Set(Vw);

    // note: V.lambda_cv = PhysConst::epsilon0;
    // This does not have to be set every iteeration,
    // since it is constant. Make sure it has been set
    // when the V-object was initialised, though!
    for (unsigned i = 1; i < grid.num_np() - 1; ++i) {
      V.rho[i] = PhysConst::e * (dens_i[i] - dens_e[i]);
    }
    // now we calculate the potential on the new time step
    // was pass ne and ni. When we are starting (or recovering
    // from a failed time step), we pass 0 for the time step;
    // this suppresses the Ventzek correction that relies
    // on the ne and ni properties being correctly calculated.
    V.Update(restarting ? 0.0 : dt, dens_e, dens_i);

    // now update E=-grad(V)
    gradient_np_cv(grid, V, E);
    for (unsigned i = 0; i < grid.num_ew(); ++i) {
      E[i] *= -1;
    }

    // now update the mobilities, which are defined on
    // the control volume boundary (cv) points.
    for (unsigned i = 0; i < grid.num_ew(); ++i) {
      const double E_N = std::abs(E[i] / gasdens);
      dens_e.beta_cv[i] = mu_e(E_N);
      dens_i.beta_cv[i] = mu_i(E_N);
    }

    // next we update the diffusion coefficients
    // and ionisation source terms, all defined
    // on the nodal points.
    for (unsigned i = 0; i < grid.num_np(); ++i) {
      // we require the reduced field on the
      // nodal points. Unfortunately, E is
      // defined on the cv points so we do
      // a linear interpolation. At the boundary
      // points, the vb and nodal points coincide
      // so we can use the E value directly.
      double E_nod;
      if (i == 0) {
        E_nod = E[0];
      } else if (i == grid.num_np() - 1) {
        E_nod = E[grid.num_np() - 2];
      } else {
        E_nod = 0.5 * (E[i - 1] + E[i]);
      }
      const double E_N = std::abs(E_nod / gasdens);

      dens_e.lambda[i] = D_e(E_N);
      dens_i.lambda[i] = D_i(E_N);

      // now we calculate the ionisation source term,
      // Sion = ngas*ne*Kion. In LFA simulations, we
      // use ne_lfa rather than the normal, real, ne.
      // This way, we mimick the Townsend model for
      // the ionisation rate (alpha coefficient).
      // NOTE: the code below is only an approximation
      //       of the Townsend model. The code assumes
      //       n_e \approx Flux_e/V_drift_e
      double ne_lfa;
      if (i == 0) {
        const double ge = dens_e.flux_density()[0];
        const double me = dens_e.beta_cv[0];
        ne_lfa = std::abs(ge / (me * E_nod));
      } else if (i == grid.num_np() - 1) {
        const double ge = dens_e.flux_density()[grid.num_np() - 2];
        const double me = dens_e.beta_cv[grid.num_np() - 2];
        ne_lfa = std::abs(ge / (me * E_nod));
      } else {
        const double ge =
            0.5 * (dens_e.flux_density()[i - 1] + dens_e.flux_density()[i]);
        const double me = 0.5 * (dens_e.beta_cv[i - 1] + dens_e.beta_cv[i]);
        ne_lfa = std::abs(ge / (me * E_nod));
      }
      sion[i] = Kion(E_N) * gasdens * ne_lfa;
    }
    // now copy the sources into the
    // source-members of ne, ni. We do this
    // only for the interior points: the boundary
    for (unsigned i = 1; i < grid.num_np() - 1; ++i) {
      dens_e.sc[i] = sion[i];
      dens_i.sc[i] = sion[i];

      // the time-derivative part.
      // note: the volume is missing
      // here, since that multiplication
      // is done by the phi variable:
      // that expects sc, sp to be source
      // *densities*.
      const double ap0 = 1.0 / dt;

      dens_e.sp[i] = -ap0;
      dens_e.sc[i] += ap0 * dens_e[i];

      dens_i.sp[i] = -ap0;
      dens_i.sc[i] += ap0 * dens_i[i];
    }

    // boundary conditions at the western side:
    //
    // first we calculate the reduced field
    // and the thermal ion and electron velocities.
    const double E_N_w = std::abs(E[0] / gasdens);
    double Ve_w = 0.25 * Vth_e_E(E_N_w);
    double Vi_w = 0.25 * Vth_i();

    // if the field points to the wall, we
    // also take into account the ion drift,
    // otherwise we take into account the
    // electron drift.
    if (E[0] < 0) {
      Vi_w += -dens_i.beta_cv[0] * E[0];
    } else {
      Ve_w += -dens_e.beta_cv[0] * E[0];
    }
    // calculate the ion flux at the wall,
    // using the old ion density field.
    const double flux_i_w = dens_i.CalculateFluxDens(0);
    // if the ion flux points to the wall,
    // we take into account an electron wall source
    // due to secondary emission.
    const double sec_el_w = flux_i_w < 0 ? -flux_i_w * gam : 0;

    // now discretize the ion boundary conditions at the
    // west wall. We set the 'fluid flux' equal to the
    // wall flux due to thermal motion (and possibly drift).
    dens_i.DiscFluxDens(0, dens_i.system().ae(0), dens_i.system().ap(0));
    dens_i.system().ap(0) += Vi_w;

    // we do the same for the electrons, but in this case
    // we also have to take into account the secondary
    // emission source.
    dens_e.DiscFluxDens(0, dens_e.system().ae(0), dens_e.system().ap(0));
    dens_e.system().ap(0) += Ve_w;
    dens_e.system().b(0) = +sec_el_w;

    // boundary conditions at the eastern side:
    //
    // This is similar the procedure at the western
    // wall, only the orientation of the wall is different
    // so (some) minus signs change.

    const unsigned ndx_e_cv = grid.num_ew() - 1;
    const double E_N_e = std::abs(E[ndx_e_cv] / gasdens);
    double Ve_e = 0.25 * Vth_e_E(E_N_e);
    double Vi_e = 0.25 * Vth_i();
    if (E[ndx_e_cv] > 0) {
      Vi_e += dens_i.beta_cv[ndx_e_cv] * E[ndx_e_cv];
    } else {
      Ve_e += dens_e.beta_cv[ndx_e_cv] * E[ndx_e_cv];
    }
    const unsigned ndx_e_np = grid.num_np() - 1;

    const double flux_i_e = dens_i.CalculateFluxDens(ndx_e_cv);
    const double sec_el_e = flux_i_e > 0 ? flux_i_e * gam : 0;
    dens_i.DiscFluxDens(ndx_e_cv, dens_i.system().ap(ndx_e_np),
                        dens_i.system().aw(ndx_e_np));
    dens_i.system().ap(ndx_e_np) += Vi_e;

    dens_e.DiscFluxDens(ndx_e_cv, dens_e.system().ap(ndx_e_np),
                        dens_e.system().aw(ndx_e_np));
    dens_e.system().ap(ndx_e_np) += Ve_e;
    dens_e.system().b(ndx_e_np) = +sec_el_e;

    // this concludes the preparation of the
    // density updates.
    //
    // Now we are going to advance ne and ni in time...

    // move ne to the new time step
    if (!Solve(dens_e, dens_e_old, tolerance)) {
      // if false was returned, the change was
      // bigger than allowed. In this case, if
      // we can still reduce the time step, we
      // undo the time step for all variables.
      // This is one of the rare occasions where
      // a 'goto' statement is useful.
      if (dt > dt_min) {
        std::cout << "Timestep " << dt << " too big for ne" << std::endl;
        goto undo_step;
      }
    }

    // move ni to the new time step, The procedure
    // is the same as for the electrons.
    if (!Solve(dens_i, dens_i_old, tolerance)) {
      if (dt > dt_min) {
        std::cout << "Timestep " << dt << " too big for ni" << std::endl;
        goto undo_step;
      }
    }

    // good. if we reach this point, all updates were
    // succesful. Now it is time to do some postprocessing,
    // calculate a new time step, et cetera. After that,
    // we can continue wih the next time step.

    // Calculate the electric current density
    for (unsigned i = 0; i < grid.num_ew(); ++i) {
      J_cv[i] =
          PhysConst::e * (+dens_i.flux_density()[i] - dens_e.flux_density()[i]);
    }

    // write the output files if t crossed the
    // value of tout:
    if (t >= tout) {
      // std::cout << "Writing output files." << std::endl;

      // append <ne> and <ni> to the meanvals file
      meanvals << t << '\t' << mean_ni << '\t' << mean_ne << std::endl;

      std::ofstream rfs("rho.dat");
      grid.write(rfs, V.rho);

      std::ofstream Vfs("V.dat");
      grid.write(Vfs, V);

      std::ofstream Efs("E.dat");
      grid.write(Efs, E);

      std::ofstream efs("ne.dat");
      grid.write(efs, dens_e);

      std::ofstream flux_efs("flux_ne.dat");
      grid.write(flux_efs, dens_e.flux_density());

      std::ofstream flux_ifs("flux_ni.dat");
      grid.write(flux_ifs, dens_i.flux_density());

      std::ofstream J_cv_fs("J_cv.dat");
      grid.write(J_cv_fs, J_cv);

      std::ofstream ifs("ni.dat");
      grid.write(ifs, dens_i);

      std::ofstream sfs("sion.dat");
      grid.write(sfs, sion);

      // we update tout so the files will
      // be written again 'dt_out' later.
      tout += dt_out;
    }
    // try to increase dt...
    dt *= 1.1;
    // but make sure it does not exceed dt_max.
    dt = std::min(dt, dt_max);
    // we finished a regular time step.
    // this means that the next iteration
    // we are not restarting a failed one
    // (or do the first one). Also, this is the
    // time to increase the iteration counter
    // and set t to the (now) new time value.
    restarting = false;
    ++iter;
    t += dt;
    // and we continue the while-loop.
    continue;

    // when a time step fails, we 'goto'
    // this point in the file. We simply
    // restore the old values of the fields
    // V, ne and ni (which were saved at the
    // beginning of the time step), reduce the
    // time step if possible, and set restarting
    // to true. Then we go back to the beginning
    // of the while block and try again...
  undo_step:
    // std::cout << "Undoin timestep, new dt = " << dt << std::endl;
    for (unsigned i = 0; i < grid.num_np(); ++i) {
      V[i] = V_old[i];
      dens_e[i] = dens_e_old[i];
      dens_i[i] = dens_i_old[i];
    }
    dt = std::max(dt_min, dt / 10.0);
    restarting = true;
  }
  std::cout << "Simulation finished; t =" << t << std::endl;
  return 0;
}
