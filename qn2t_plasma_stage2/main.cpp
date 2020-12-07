#include <cmath>
#include <iostream>

#include "argon.h"
#include "phivar.h"
#include "physconst.h"

// The physical radius of the grid in m.
const double grid_size = 0.0125;
// The initial value for Te
const double Te_init = 12160.3;
// The maximum of the initialization function of Th
const double Th_max = 400.0;
// The minimum of the initialization function of Th
const double Th_min = 300.0;
// The background particle density
const double n0_dens = 4e22;
// The ion density at the wall. Supposed to be low to
// represent recombination.
const double ion_wall_dens = 1e12;
// The maximum of the initialization function of the ion density.
const double ion_dens_max = 1e18;
// The minimum of the initialization function of the ion density.
const double ion_dens_min = 1e12;
// The average density of the filling gas.
const double ave_dens = n0_dens;
// The number of grid points
const int grid_points = 16;
// Power density calculation
const double grid_length = 1.0; // The length of the rod
const double W = 100;           // Wattage
const double power_density = W / (M_PI * sqr(grid_size) * grid_length);
// argon / electron properties
const double k_rate = 1e-15;
const double G = 6.0;
const double q = 0.5;
const double sigma_ea = 2e-20;
// settings
const double urf_ion = 0.95;
const double urf_Te = 0.99;

// B.Broks, 9-1-04:
// A function which computes your local buffer density based on an average
// density of the plasma gas. This avearge density is the argument of your
// function call. It will require the input of following fields: An electron
// temperature Field called Te. A heavy particle temperature Field called Th. An
// electron density Field called ion_dens. A heavy particle density Field
// neutr_dens. It will modify the contents of neutr_dens so the neutr_dens
// becomes the needed buffer density.
void calculate_buffer_density(double ave_dens, const Grid &grid,
                              PhiVariable &neutr_dens,
                              const PhiVariable &ion_dens,
                              const PhiVariable &Te, const PhiVariable &Th) {
  // Step zero: the volume of the grid
  // Iterate over the whole gird, adding the local volume tho the running total.
  double grid_volume = 0.0;
  unsigned i = 0;
  for (i = 0; i < grid.num_np(); ++i)
    grid_volume += grid.vol_np(i);

  // Step one: compute the point where we have the highest.
  // nonbuffer density. While we are at it, we also compute many
  // of our total particles we have used.

  double highest_local_nonbuffer_pressure = 0.0;
  // The amount of particles used.
  double total_particles_used = 0.0;
  for (i = 0; i < grid.num_np(); ++i) {
    // p=n_e k_b T_e (elec) +n_e k_b T_h (ion)
    double nonbuffer_press_local = PhysConst::k_b * Te[i] * ion_dens[i] +
                                   PhysConst::k_b * Th[i] * ion_dens[i];
    total_particles_used += ion_dens[i] * grid.vol_np(i);
    if (nonbuffer_press_local > highest_local_nonbuffer_pressure)
      highest_local_nonbuffer_pressure = nonbuffer_press_local;
  }
  // Okay, now we know what the maximum pressure is.
  // Step two: We will now assign particles to each point to remove
  // the pressure difference.

  for (i = 0; i < grid.num_np(); ++i) {
    double nonbuffer_press_local = PhysConst::k_b * Te[i] * ion_dens[i] +
                                   PhysConst::k_b * Th[i] * ion_dens[i];
    double pressure_to_compensate =
        highest_local_nonbuffer_pressure - nonbuffer_press_local;
    double particles_used_to_compensate =
        pressure_to_compensate / (PhysConst::k_b * Th[i]) * grid.vol_np(i);
    // Again, this costs us particles.
    total_particles_used += particles_used_to_compensate;
    neutr_dens[i] = pressure_to_compensate / (PhysConst::k_b * Th[i]);
  }
  // Final step:Distribute the rest.
  // Compute the average value of 1/T_h. The particles are distributed
  // proportional to 1/Th.
  double one_over_t_ave = 0;
  for (i = 0; i < grid.num_np(); ++i)
    one_over_t_ave += grid.vol_np(i) / grid_volume / Th[i];
  // Ditribute the paricles.
  for (i = 0; i < grid.num_np(); ++i)
    neutr_dens[i] += 1.0 / one_over_t_ave / Th[i] *
                     (ave_dens - total_particles_used / grid_volume);
}

double kion(double Te) {
  return k_rate * pow(Te, q) * exp(-Argon::E_ion() / (PhysConst::k_b * Te));
}

double krec(double Te) {
  return pow((PhysConst::h_planck /
              sqrt(2 * M_PI * PhysConst::me * PhysConst::k_b * Te)),
             3) *
         k_rate * pow(Te, q) / (2 * G);
}

double ambipolar_diffusion_coefficient(double n0, double Te, double Th) {
  double sig_ia = 7e-19;
  double tau_ia =
      sqrt(M_PI * Argon::Mass() / (8 * PhysConst::k_b * Th)) / (sig_ia * n0);

  return PhysConst::k_b * Th * tau_ia / Argon::Mass() * (1 + Te / Th);
}

double e_heat_cond(double ne, double na, double Te) {
  return PhysConst::k_b / sigma_ea * ne / na *
         sqrt(8 * PhysConst::k_b * Te / (M_PI * PhysConst::me));
}

double e_heavy_heat_trans(double ne, double na, double Te) {
  double v_ea = ne * na * sigma_ea *
                sqrt(8 * PhysConst::k_b * Te / (M_PI * PhysConst::me));
  return v_ea * 3 * PhysConst::me / Argon::Mass() * PhysConst::k_b;
}

int main(int argc, char **argv) {
  double S_plus_constant;
  double S_plus_prop;
  double S_eh;

  double res_ion;
  double res_Te;
  int iterations = 0;
  // The grid. The first argument is the amount of grid points, the
  // second argument is the position of the left edge, and the third
  // argument is the position of the right edge.
  Grid grid(grid_points, 0.0, grid_size, Grid::Cylindrical);
  // the velocity is defined on the cv boundaries. The value is zero
  Field flow_vel(grid.num_ew(), 0.0);

  // Te. We have two Neumann conditions.
  NeumannBndCond Te_left(0, 0);
  NeumannBndCond Te_right(0, 0);
  PhiVariable Te(grid, flow_vel, Te_left, Te_right);

  // Th. We cool the right wall to 300 K.
  NeumannBndCond Th_left(0, 0);
  DirichletBndCond Th_right(Th_min);
  PhiVariable Th(grid, flow_vel, Th_left, Th_right);

  // Ion Density. We have recombination at the right wall, represented by the
  // low ion density value.
  NeumannBndCond ion_dens_left(0, 0);
  DirichletBndCond ion_dens_right(ion_wall_dens);
  PhiVariable ion_dens(grid, flow_vel, ion_dens_left, ion_dens_right);

  // neutr_dens
  NeumannBndCond neutr_dens_left(0, 0);
  NeumannBndCond neutr_dens_right(0, 0);
  PhiVariable neutr_dens(grid, flow_vel, neutr_dens_left, neutr_dens_right);

  // initialization
  for (unsigned i = 0; i < grid.num_np(); ++i) {
    Te[i] = Te_init;
    Th[i] = Th_max - sqr(grid.pos_np(i) / grid.pos_np(grid.num_np() - 1)) *
                         (Th_max - Th_min);
  }
  ion_dens = ion_dens_max;
  // set underrelaxation
  Te.set_urf(urf_Te);
  ion_dens.set_urf(urf_ion);

  // Solve
  do {
    // update
    calculate_buffer_density(ave_dens, grid, neutr_dens, ion_dens, Te, Th);
    for (unsigned i = 0; i < grid.num_np(); ++i) {
      // calculate source terms
      S_plus_constant = neutr_dens[i] * ion_dens[i] * kion(Te[i]);
      S_plus_prop = pow(ion_dens[i], 2) * krec(Te[i]);
      S_eh = e_heavy_heat_trans(ion_dens[i], neutr_dens[i], Te[i]);

      // update ion density
      ion_dens.lambda[i] =
          ambipolar_diffusion_coefficient(neutr_dens[i], Te[i], Th[i]);
      ion_dens.sc[i] = S_plus_constant;
      ion_dens.sp[i] = -S_plus_prop;

      // update Te
      Te.lambda[i] = e_heat_cond(ion_dens[i], neutr_dens[i], Te[i]);
      Te.sc[i] =
          power_density -
          (S_plus_constant - S_plus_prop * ion_dens[i]) * Argon::E_ion() -
          S_eh * Th[i];
      Te.sp[i] = S_eh;
    }

    // solve
    res_ion = ion_dens.Update();
    res_Te = Te.Update();

    iterations++;
    std::cout << "Iteration " << iterations << "\t Residue ion " << res_ion
              << "\t Residue Te: " << res_Te << std::endl;
  } while (std::max(res_Te, res_ion) > 1e-8);

  std::cout << "Ion density centre: " << ion_dens[0] << std::endl;
  std::cout << "Neutral density centre: " << neutr_dens[0] << std::endl;
  std::cout << "Te centre: " << Te[0] << std::endl;

  // plot
  grid.plot(ion_dens, "position (m)", "density (m^{-3})",
            "Exercise 5.34: Ion Density");
  getchar();
  grid.plot(neutr_dens, "position (m)", "density (m^{-3})",
            "Exercise 5.34: Neutral Density");
  getchar();
  grid.plot(Te, "position (m)", "temperature (K)",
            "Exercise 5.34: Electron temperature");

  // Ready
  return 0;
}