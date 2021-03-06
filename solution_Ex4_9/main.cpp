#include "3diagsys.h"
#include "bndcond.h"
#include "field.h"
#include "grid.h"
#include <cmath>
#include <iostream>

int main() {
  unsigned size = 10;
  double r = 1;
  // double T_e = 100;
  double T_w = 200;
  double Q = 10;
  double lambda = 0.02;

  Grid grid_rod(size - 2, 0.0, r, Grid::Cartesian);
  Field field_T = grid_rod.pos_np();
  TridiagonalSystem system_T(size);

  // West BC:
  DirichletBndCond bc_west(T_w);
  bc_west.ModifyCoefs(system_T, false, grid_rod.del());

  // East BC:
  NeumannBndCond bc_east(0, 0);
  bc_east.ModifyCoefs(system_T, true, grid_rod.del());
  // Interior points
  for (unsigned i = 1; i < size - 1; ++i) {
    system_T.aw(i) = lambda * grid_rod.area_ew(i - 1) / grid_rod.dx_np(i);
    system_T.ae(i) = lambda * grid_rod.area_ew(i) / grid_rod.dx_np(i);
    system_T.ap(i) = system_T.aw(i) + system_T.ae(i);
    system_T.b(i) = Q * grid_rod.vol_np(i);
  }

  system_T.solve(field_T);

  grid_rod.write(std::cout, field_T);
  grid_rod.plot(field_T, "x (m)", "T (K)", "Exercise 4.11 Neumann");
  return 0;
}