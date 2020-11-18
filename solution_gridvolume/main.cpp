#include "field.h"
#include "grid.h"
#include <cmath>
#include <iostream>

int main() {
  unsigned size = 5;
  double r = 1.0e-2;

  Grid sph_grid(size, 0.0e-2, r, Grid::Spherical);
  double vol_disc = sph_grid.volume();
  double vol_anal = 4.0 / 3.0 * M_PI * pow(r, 3);
  double diff = std::abs(vol_disc - vol_anal);
  double diff_percentage = diff / vol_anal * 100;

  std::cout << "discrete volume\t\tanalytical volume\terror\t\terror (%)"
            << std::endl;
  std::cout << vol_disc << "\t\t" << vol_anal << "\t\t" << diff << "\t"
            << diff_percentage << std::endl;

  return 0;
}
