#include "field.h"
#include "grid.h"
#include <cmath>
#include <iostream>

int main() {
  unsigned size = 4;

  Grid sph_grid(size, 0.0e-2, 1.0e-2, Grid::Spherical);
  Field rsquared(size + 1);
  for (unsigned i = 0; i < size + 1; i++) {
    rsquared[i] = sqr(sph_grid.pos_ew(i));
  }
  std::cout << "write of field" << std::endl;
  rsquared.write(std::cout);
  std::cout << "write of grid" << std::endl;
  sph_grid.write(std::cout, rsquared);
  sph_grid.plot(rsquared, "r (m)", "r^2", "Exercise 4.4");
  return 0;
}
