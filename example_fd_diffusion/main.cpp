#include "../3diagsys.h"
#include "../field.h"
#include <iostream>

int main() {
  // 1. Configuration

  unsigned N = 3;                 // number of points;
  const double L = 1;             // grid size
  const double del = L / (N - 1); // distance between adjacent grid points
  const double S = 1;             // source density
  const double D = 1;             // diffusion coefficient

  const double dens_w = 2; // value at west boundary.
  const double dens_e = 1; // value at east boundary.

  // 2. Set up the variables

  Field dens(N);
  TridiagonalSystem sys(N);

  // 3. Discretise the equations

  // West BC:
  sys.ap(0) = 1;
  sys.b(0) = 0;
  sys.aw(0) = 1;

  // East BC:
  sys.ap(N - 1) = 1;
  sys.ae(N - 1) = 1;
  sys.b(N - 1) = 0;

  // Interior points: -Dd^2n/dx^2 = S.
  for (unsigned i = 1; i < N - 1; ++i) {
    sys.ap(i) = +2 * D / (del * del);
    sys.aw(i) = +1 * D / (del * del);
    sys.ae(i) = +1 * D / (del * del);
    sys.b(i) = S;
  }

  // 4. Solve the systeem of equations
  sys.solve(dens);

  // 5. write position, the density, the analytical density and the error.
  std::cout << "pos\tres\tanalytic\terror" << std::endl;

  for (unsigned i = 0; i < N; ++i) {
    const double x = del * i;
    const double dens_an =
        0.5 * (S / D) * x * (L - x) + (dens_e - dens_w) * x / L + dens_w;
    const double error = dens[i] - dens_an;
    std::cout << x << '\t' << dens[i] << '\t' << dens_an << "\t \t" << error
              << std::endl;
  }

  return 0;
}
