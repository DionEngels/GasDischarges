#include <fstream>
#include <iostream>

#include "phivar.h"

int main(int argc, char **argv) {
  const unsigned N = 16;      // Gridsize
  const double xa = 0.0;      // x at western point
  const double xb = 10.0;     // x at eastern point
  const double Tr = 400.0;    // T at eastern point
  const double V = 0.0;       // Velocity V, does not exists so set to 0
  const double lambda = 1.0;  // heat conductivivity
  const double Sc = 10.0;     // Source term (constant)
  const double Sp = 0.0;      // Source term (proportonial)
  const double beta_cv = 1.0; // Heat capacity

  // Define Cartesian grid using grid class
  Grid grid(N, xa, xb, Grid::Cartesian);

  // Define Field of generalized velocity for all control vlume faces
  Field w(grid.num_ew(), V);

  // Define the two BCs. East has Neumann with derivative 0
  // West has constant T: Tr
  NeumannBndCond temp_left(0, 0);
  DirichletBndCond temp_right(Tr);

  // Define PhiVariable class called temp based on grid and field w, with BCs
  // temp_left and temp_right
  PhiVariable temp(grid, w, temp_left, temp_right);

  // Set physical variables for the PhiVariable, such as heat capacity, source
  // terms and heat conductivity.
  temp.beta_cv = beta_cv; // Set value for heat capacity
  // For all nodes, set heat conducitivty and source terms
  for (unsigned i = 0; i < grid.num_np(); ++i) {
    temp.lambda[i] = lambda; // heat conductivity
    temp.sc[i] = Sc;         // constant source
    temp.sp[i] = Sp;         // proportional source
  }

  // Solve
  temp.Update();

  // analytical
  Field T_analytical(grid.num_np());
  for (unsigned i = 0; i < grid.num_np(); ++i) {
    T_analytical[i] = -2.5 * sqr(grid.pos_np(i)) + 900;
  }

  // Write result to output
  grid.write(std::cout, temp);
  // Write result to T.dat
  std::ofstream ofs_T("output/T.dat");
  grid.write(ofs_T, temp);
  // Write flux density to heat_flux_density.dat
  std::ofstream ofs_flux("output/heat_flux_density.dat");
  grid.write(ofs_flux, temp.flux_density());

  // Write analytical to output
  grid.write(std::cout, T_analytical);

  // Plot result
  grid.plot(temp, T_analytical, "position (m)", "temperature (K)",
            "temperature (K)", "Exercise 4.16");

  return 0;
}
