#include "phivar.h"
#include <cmath>
#include <fstream>
#include <iostream>

// declare function
double heavy_heat_conductivity(double T_h);

int main(int argc, char **argv) {
  const unsigned N = 7;
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
  std::vector<double> T_vector;

  // intial conditions
  double T_central_init = 2000.0; // save for later
  double T_central = T_central_init;
  double underrelaxation = 1.01;

  // define grid
  Grid grid(N, x_e, x_w, Grid::Cylindrical);
  Field w(grid.num_ew(), V);

  // BC
  NeumannBndCond temp_east(0, 0);
  DirichletBndCond temp_west(T_w);

  // set initial variables
  PhiVariable temp(grid, w, temp_east, temp_west);
  temp.set_urf(underrelaxation);
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
    T_vector.push_back(T_central);
    iterations++;
    std::cout << "Iteration " << iterations << "\t Current residue " << residue
              << std::endl;
  } while (residue > 10e-6);

  std::cout << "Final T_central: " << T_central << " K" << std::endl;

  // create field and grid to plot
  Grid it_grid(iterations - 1, 0.0, iterations, Grid::Cartesian);
  Field it_field(iterations + 1);
  it_field[0] = T_central_init;
  int index = 1;
  for (std::vector<double>::iterator it = T_vector.begin();
       it != T_vector.end(); ++it) {
    it_field[index] = *it;
    index++;
  }

  grid.plot(temp, "position (m)", "temperature (K)", "Exercise 5.11");
  it_grid.plot(it_field, "iterations", "central temperature (K)",
               "Underrelaxation: " + std::to_string(underrelaxation));
  return 0;
}

// implement function
double heavy_heat_conductivity(double T_h) {
  const double k_b = 1.380658e-23;
  const double m_h = 6.623e-26;
  const double s_aa = 6.33e-20;

  return sqrt(8.0 * k_b * T_h / (M_PI * m_h)) * sqrt(2) / s_aa * k_b;
}