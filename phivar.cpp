/** \file
 *
 *  Implementation of the code which has been declared in phivar.h
 *  $Id$
 */

#include "phivar.h"
#include "basics.h"
#include "grid.h"

#include <cmath>
#include <iostream>

PhiVariable::PhiVariable(const Grid &grid, const Field &w_cv,
                         const BoundaryCondition &bc_w,
                         const BoundaryCondition &bc_e)
    : Field(grid.num_np()), beta_cv(grid.num_ew(), 0.0),
      lambda(grid.num_np(), 0.0), sc(grid.num_np(), 0.0),
      sp(grid.num_np(), 0.0), m_grid(grid), m_w_cv(w_cv),
      m_system(grid.num_np()), m_bc_w(bc_w), m_bc_e(bc_e), m_urf(1.0),
      m_flux_dens(grid.num_ew(), 0.0), m_Pe_cv(grid.num_ew(), 0.0) {}

namespace disc {
inline double A(double P) {
  if (P >= 0) {

  } else {
    std::cerr << "P >=0 assertion failed" << std::endl;
    std::cerr << "P = " << P << std::endl;
    abort();
  }
#define EXPONENTIAL_SCHEME
#ifdef EXPONENTIAL_SCHEME
  return P == 0 ? 1 : P / (std::exp(P) - 1);
#else
  return 1 - P / 2;
#endif
}
inline void AB(const double P, double &A, double &B) {
  const double APabs = disc::A(std::abs(P));
  A = APabs + std::max(-P, 0.0);
  B = APabs + std::max(P, 0.0);
}

} // namespace disc

void PhiVariable::DiscFluxDens(unsigned ndx_cv, double &DA, double &DB) const {
  const double dLH = (ndx_cv == 0 || ndx_cv == m_grid.num_ew() - 1)
                         ? m_grid.del() / 2
                         : m_grid.del();
  const double D = lambda_cv(ndx_cv) / dLH;
  const double F = beta_cv[ndx_cv] * m_w_cv[ndx_cv];
  const double Pe = F / D;
  m_Pe_cv[ndx_cv] = Pe;
  disc::AB(Pe, DA, DB);
  DA *= D;
  DB *= D;
  // std::cout << "D: " << kL << ", " << kH << std::endl;
}

double PhiVariable::CalculateFluxDens(unsigned ndx_cv) const {
  const unsigned ndx_L = ndx_cv;
  const unsigned ndx_H = ndx_L + 1;
  double DA, DB;
  DiscFluxDens(ndx_cv, DA, DB);
  return DB * (*this)[ndx_L] - DA * ((*this)[ndx_H]);
}

void PhiVariable::UpdateFluxDens() {
  for (unsigned i = 0; i < m_grid.num_ew(); ++i) {
    m_flux_dens[i] = CalculateFluxDens(i);
  }
}

double PhiVariable::lambda_cv(unsigned ndx_cv) const {
  if (ndx_cv == 0) {
    return lambda[0];
  }
  if (ndx_cv == m_grid.num_ew() - 1) {
    return lambda[m_grid.num_np() - 1];
  }
  const double gL = lambda[ndx_cv];
  const double gH = lambda[ndx_cv + 1];
  return 2 * gL * gH / (gL + gH);
}

void PhiVariable::Discretize() {
  // we calculate the coefficients in the *internal* nodal points.
  for (unsigned p = 1; p < m_grid.num_np() - 1; ++p) {
    double DA, DB;
    m_system.ap(p) = 0;

    unsigned ncv_w = p - 1;
    const double areaW = m_grid.area_ew(ncv_w);
    DiscFluxDens(ncv_w, DA, DB);
    m_system.aw(p) = DB * areaW;
    m_system.ap(p) += DA * areaW;

    unsigned ncv_e = p;
    const double areaE = m_grid.area_ew(ncv_e);
    DiscFluxDens(ncv_e, DA, DB);
    m_system.ae(p) = DA * areaE;
    m_system.ap(p) += DB * areaE;

    const double vol = m_grid.vol_np(p);
    m_system.ap(p) -= sp[p] * vol;
    m_system.b(p) = sc[p] * vol;
  }
  // the boundary conditions
  m_bc_w.ModifyCoefs(m_system, true, m_grid.del() / 2.0);
  m_bc_e.ModifyCoefs(m_system, false, m_grid.del() / 2.0);
}

void PhiVariable::ApplyRelaxation() {
  if (urf() == 1.0) {
    return;
  }
  for (unsigned p = 1; p < size() - 1; ++p) {
    m_system.ap(p) /= urf();
    m_system.b(p) += (1.0 - urf()) * m_system.ap(p) * (*this)[p];
  }
}

// Update the field by calculation the discretization coefficients and
// solving the resulting matrix-vector problem. Returns the residual, i.e.
// a measure of how much this solution changed with respect to the previous
// one.
double PhiVariable::Update() {
  Discretize();

  ApplyRelaxation();

  const Field old(*this);

  m_system.solve(*this);

  UpdateFluxDens();

  // calculate the residual.
  double max_val = 0, max_diff = 0;

  for (unsigned p = 0; p < size(); ++p) {

    const double diff = std::abs(old[p] - (*this)[p]);
    max_diff = std::max(diff, max_diff);

    const double val = std::abs((*this)[p]);
    max_val = std::max(val, max_val);
  }
  return (max_val == 0.0) ? max_diff : (max_diff / max_val);
}
