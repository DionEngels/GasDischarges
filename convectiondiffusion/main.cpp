#include <fstream>
#include <iostream>

#include "phivar.h"

int main(int argc, char **argv) {
  const unsigned N = 100;       // Gridsize
  const double xa = 0.0;        // x at western point
  const double xb = 1.0;        // x at eastern point
  const double T_e = 300.0;     // T at eastern point
  const double T_w = 400.0;     // T at western point
  const double V = 10.0;        // Velocity V, does not exists so set to 0
  const double lambda = 1.0;    // heat conductivivity
  const double beta_cv = 250.0; // Heat capacity

  // Define Cartesian grid using grid class
  Grid grid(N, xa, xb, Grid::Cartesian);

  // Define Field of generalized velocity for all control vlume faces
  Field w(grid.num_ew(), V);

  // Define the two BCs. East has Neumann with derivative 0
  // West has constant T: Tr
  DirichletBndCond temp_left(T_e);
  DirichletBndCond temp_right(T_w);

  // Define PhiVariable class called temp based on grid and field w, with BCs
  // temp_left and temp_right
  PhiVariable temp(grid, w, temp_left, temp_right);

  // Set physical variables for the PhiVariable, such as heat capacity, source
  // terms and heat conductivity.
  temp.beta_cv = beta_cv; // Set value for heat capacity
  // For all nodes, set heat conducitivty and source terms
  for (unsigned i = 0; i < grid.num_np(); ++i) {
    temp.lambda[i] = lambda; // heat conductivity
  }

  // Solve
  temp.Update();

  // Write result to output
  grid.write(std::cout, temp);
  // Write result to T.dat
  std::ofstream ofs_T("output/T.dat");
  grid.write(ofs_T, temp);
  // Write flux density to heat_flux_density.dat
  std::ofstream ofs_flux("output/heat_flux_density.dat");
  grid.write(ofs_flux, temp.flux_density());

  // Plot result
  grid.plot(temp, "position (m)", "temperature (K)", "Exercise 4.18");

  return 0;
}