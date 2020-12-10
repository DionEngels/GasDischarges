#include "3diagsys.h"
#include "cmath"
#include "field.h"
#include "grid.h"
#include "physconst.h"
#include <iostream>

int main() {
  // 1. Configuration

  unsigned N = 100;               // number of points;
  const double L = 0.01;          // grid size
  const double del = L / (N - 1); // distance between adjacent grid points
  const double e0 = PhysConst::epsilon0; // e0

  const double V_w = 0;   // value at west boundary.
  const double V_e = 100; // value at east boundary.
  double n_i;
  double n_e;
  double S; // source density
  double S_max = 0;
  double S_min = 0;

  // 2. Set up the variables

  Field dens(N);
  Grid grid(N - 2, 0, L, Grid::Cartesian);
  TridiagonalSystem sys(N);

  // Space charge
  const double e_mob = 300;
  double drift_vel = e_mob * V_e / L;
  double drift_time = 1e-10; // 1 ns
  double drift_dist = drift_time * drift_vel;

  // 3. Discretise the equations

  // West BC:
  sys.ap(0) = 1;
  sys.aw(0) = 0;
  sys.b(0) = V_w;

  // East BC:
  sys.ap(N - 1) = 1;
  sys.ae(N - 1) = 0;
  sys.b(N - 1) = V_e;

  // Interior points
  for (unsigned i = 1; i < N - 1; ++i) {
    sys.ap(i) = +2 * e0 / (del * del);
    sys.aw(i) = +1 * e0 / (del * del);
    sys.ae(i) = +1 * e0 / (del * del);
    n_e =
        1e18 * exp(-0.5 * sqr((grid.pos_np(i) - drift_dist - 0.003) / 0.0005));
    n_i = 1e18 * exp(-0.5 * sqr((grid.pos_np(i) - 0.003) / 0.0005));
    S = PhysConst::e * (n_i - n_e);
    sys.b(i) = S;
  }

  // 4. Solve the systeem of equations
  sys.solve(dens);

  // output

  std::cout << "Drift distance: " << drift_dist << " m" << std::endl;
  std::cout << "Node position \tCharge density" << std::endl;
  for (unsigned i = 1; i < N - 1; ++i) {
    std::cout << grid.pos_np(i) << "\t" << sys.b(i) << std::endl;
    if (sys.b(i) > S_max) {
      S_max = sys.b(i);
    } else if (sys.b(i) < S_min) {
      S_min = sys.b(i);
    }
  }
  std::cout << "Maximum charge:\t" << S_max << " C/m^{-3}" << std::endl;
  std::cout << "Minimum charge:\t" << S_min << " C/m^{-3}" << std::endl;

  grid.plot(dens, "position (m)", "Voltage (V)",
            "Exercise 6.4 (Timestep = 10^{-10})");

  return 0;
}