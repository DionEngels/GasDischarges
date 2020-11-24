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
  const double Sc = 1.02 * pow(10, 10);
  const double Sp = 0.0;
  const double k_b = 1.380658e-23;
  const double beta_cv = 2.5 * 5e23 * k_b;
  double residue;
  double lambda;
  int iterations = 0;

  // intial conditions
  double T_central = 2000.0;

  // define grid
  Grid grid(N, x_e, x_w, Grid::Cylindrical);
  Field w(grid.num_ew(), V);

  // BC
  NeumannBndCond temp_east(0, 0);
  DirichletBndCond temp_west(T_w);

  // set initial variables
  PhiVariable temp(grid, w, temp_east, temp_west);
  temp.beta_cv = beta_cv;
  for (unsigned i = 0; i < grid.num_np(); ++i) {
    temp.sc[i] = Sc;
    temp.sp[i] = Sp;
  }

  // solve
  do {
    lambda = heavy_heat_conductivity(T_central); // update lambda
    for (unsigned i = 0; i < grid.num_np(); ++i) {
      temp.lambda[i] = lambda; // update PhiVariable
    }
    residue = temp.Update(); // solve
    T_central = temp[0];     // set new T_central
    iterations++;
    std::cout << "Iteration " << iterations << "\t Current residue " << residue
              << std::endl;
  } while (residue > 10e-6);

  grid.plot(temp, "position (m)", "temperature (K)", "Exercise 5.8");
  return 0;
}

// implement function
double heavy_heat_conductivity(double T_h) {
  const double k_b = 1.380658e-23;
  const double m_h = 6.623e-26;
  const double s_aa = 6.33e-20;

  return sqrt(8.0 * k_b * T_h / (M_PI * m_h)) * sqrt(2) / s_aa * k_b;
}