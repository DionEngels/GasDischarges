#include "3diagsys.h"
#include "field.h"
#include "grid.h"
#include <cmath>
#include <iostream>

int main() {
  unsigned size = 8;
  double r = 1;
  double T_e = 100;
  // double T_w = 200;
  double Q = 10;
  double lambda = 0.02;

  Grid grid_rod(size - 2, 0.0, r, Grid::Cartesian);
  Field field_T = grid_rod.pos_np();
  TridiagonalSystem system_T(size);

  // West BC:
  system_T.ap(0) = 1;
  system_T.ae(0) = 1;
  system_T.b(0) = 0;

  // East BC:
  system_T.ap(size - 1) = 1;
  system_T.aw(size - 1) = 0;
  system_T.b(size - 1) = T_e;

  // Interior points
  for (unsigned i = 1; i < size - 1; ++i) {
    system_T.aw(i) = lambda * grid_rod.area_ew(i - 1) / grid_rod.dx_np(i);
    system_T.ae(i) = lambda * grid_rod.area_ew(i) / grid_rod.dx_np(i);
    system_T.ap(i) = system_T.aw(i) + system_T.ae(i);
    system_T.b(i) = Q * grid_rod.vol_np(i);
  }

  system_T.solve(field_T);

  grid_rod.write(std::cout, field_T);
  grid_rod.plot(field_T, "x (m)", "T (K)", "Exercise 4.10");
  return 0;
}
