#include <cmath>
#include <iostream>

#include "argon.h"
#include "basics.h"
#include "physconst.h"

// The physical radius of the grid in m.
const double grid_size = 0.0125;
// The maximum of the initialization function of Te
const double Te_max = 12000;
// The maximum of the initialization function of Th
const double Th_max = 400.0;
// The background particle density
const double n0_dens = 4e22;

double ambipolar_diffusion_coefficient(double n0, double Te, double Th) {
  double sig_ia = 7e-19;
  double tau_ia =
      sqrt(M_PI * Argon::Mass() / (8 * PhysConst::k_b * Th)) / (sig_ia * n0);

  return PhysConst::k_b * Th * tau_ia / Argon::Mass() * (1 + Te / Th);
}

double Te_from_kion(double na, double Te, double Th) {
  const double k_rate = 1e-15;
  const double q = 0.5;

  return -Argon::E_ion() /
         (PhysConst::k_b * log(ambipolar_diffusion_coefficient(na, Te, Th) /
                               (sqr(grid_size) * na * k_rate * pow(Te, q))));
}

int main(int argc, char **argv) {
  double residue;
  int iterations = 0;
  double Te = Te_max;
  double Te_old = Te_max;

  // Solve
  do {
    Te = Te_from_kion(n0_dens, Te, Th_max);
    residue = std::abs(Te - Te_old);
    Te_old = Te;
    iterations++;
    std::cout << "Iteration " << iterations << "\t Current residue " << residue
              << std::endl;
  } while (residue > 1e-8);

  std::cout << "Te: " << Te << std::endl;

  // Ready
  return 0;
}