#include "phivar.h"
#include <cmath>
#include <fstream>
#include <iostream>

// declare function
double heavy_heat_conductivity(double T_h);

int main(int argc, char **argv) {
  const unsigned N = 16;
  const double x_e = 0.0;
  const double x_w = 0.001;
  const double T_w = 300.0;
  const double V = 0.0;
  double T_central = 50000.0;
  double lambda;
  const double Sc = 1.02 * pow(10, 10);
  const double Sp = 0.0;
  const double k_b = 1.380658e-23;
  const double beta_cv = 2.5 * 5e23 * k_b;
  const int n_iterations = 10;

  Grid grid(N, x_e, x_w, Grid::Cylindrical);
  Field w(grid.num_ew(), V);

  Grid it_grid(n_iterations - 1, 0.0, n_iterations, Grid::Cartesian);
  Field it_field(n_iterations + 1);
  it_field[0] = T_central;

  NeumannBndCond temp_east(0, 0);
  DirichletBndCond temp_west(T_w);

  PhiVariable temp(grid, w, temp_east, temp_west);

  temp.beta_cv = beta_cv;
  for (unsigned i = 0; i < grid.num_np(); ++i) {
    temp.sc[i] = Sc;
    temp.sp[i] = Sp;
  }

  for (int i = 1; i < n_iterations + 1; i++) {
    lambda = heavy_heat_conductivity(T_central); // update lambda
    for (unsigned i = 0; i < grid.num_np(); ++i) {
      temp.lambda[i] = lambda; // update PhiVariable
    }
    temp.Update();       // solve
    T_central = temp[0]; // set new T_central
    std::cout << "Iteration " << i << "\t Current T_central " << T_central
              << std::endl;
    it_field[i] = T_central;
  }

  it_grid.plot(it_field, "iterations", "central temperature (K)",
               "Starting value: 50000 K");

  return 0;
}

// implement function
double heavy_heat_conductivity(double T_h) {
  const double k_b = 1.380658e-23;
  const double m_h = 6.623e-26;
  const double s_aa = 6.33e-20;

  return sqrt(8.0 * k_b * T_h / (M_PI * m_h)) * sqrt(2) / s_aa * k_b;
}