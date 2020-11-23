#include "phivar.h"
#include <cmath>
#include <fstream>
#include <iostream>

int main(int argc, char **argv) {
  const unsigned N = 16;
  const double x_e = 0.0;
  const double x_w = 0.001;
  const double T_w = 300.0;
  const double V = 0.0;
  const double lambda = 0.318;
  const double Sc = 1.02 * pow(10, 10);
  const double Sp = 0.0;
  const double k_b = 1.380658e-23;
  const double beta_cv = 2.5 * 5 * pow(10, 23) * k_b;

  Grid grid(N, x_e, x_w, Grid::Cylindrical);

  Field w(grid.num_ew(), V);

  NeumannBndCond temp_east(0, 0);
  DirichletBndCond temp_west(T_w);

  PhiVariable temp(grid, w, temp_east, temp_west);

  temp.beta_cv = beta_cv;
  for (unsigned i = 0; i < grid.num_np(); ++i) {
    temp.lambda[i] = lambda;
    temp.sc[i] = Sc;
    temp.sp[i] = Sp;
  }

  temp.Update();

  grid.write(std::cerr, temp);
  grid.plot(temp, "position (m)", "temperature (K)", "Exercise 5.3");

  return 0;
}
