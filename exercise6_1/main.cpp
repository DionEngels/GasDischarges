#include "3diagsys.h"
#include "field.h"
#include "grid.h"
#include "physconst.h"
#include <iostream>

int main() {
  // 1. Configuration

  unsigned N = 20;                // number of points;
  const double L = 0.01;          // grid size
  const double del = L / (N - 1); // distance between adjacent grid points
  const double S = 1e-4;          // source density
  const double e0 = PhysConst::epsilon0; // e0
  const double e_r = 11.68;

  const double V_w = 0;   // value at west boundary.
  const double V_e = 100; // value at east boundary.

  // 2. Set up the variables

  Field dens(N);
  Grid grid(N - 2, 0, L, Grid::Cartesian);
  TridiagonalSystem sys(N);

  // 3. Discretise the equations

  // West BC:
  sys.ap(0) = 1;
  sys.aw(0) = 0;
  sys.b(0) = V_w;

  // East BC:
  sys.ap(N - 1) = 1;
  sys.ae(N - 1) = 0;
  sys.b(N - 1) = V_e;

  // Interior points (left)
  for (unsigned i = 1; i < N / 2; ++i) {
    sys.ap(i) = +2 * e0 / (del * del);
    sys.aw(i) = +1 * e0 / (del * del);
    sys.ae(i) = +1 * e0 / (del * del);
    sys.b(i) = S;
  }

  // Interior points (right)
  for (unsigned i = N / 2; i < N - 1; ++i) {
    sys.ap(i) = +2 * e0 * e_r / (del * del);
    sys.aw(i) = +1 * e0 * e_r / (del * del);
    sys.ae(i) = +1 * e0 * e_r / (del * del);
    sys.b(i) = S;
  }

  // 4. Solve the systeem of equations
  sys.solve(dens);

  grid.plot(dens, "position (m)", "Voltage (V)", "Exercise 6.2");

  return 0;
}