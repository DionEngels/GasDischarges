// File:	argon.h
// Author:	W.J.M. Brok <wjmb@etpmod.phys.tue.nl>
// Date:	januari 2004
// Version:     Id: argon.h,v 1.8 2004/04/16 06:59:49 wjmb Exp
// Description: Various functions relating to argon data. The electron impact
//		cross section data is taken from:
//
//		  @article{Phe1999:1,
//		   author="A. V. Phelps and Z. Lj. Petrovi{\'c}",
//		   title="Review Article: Cold-Cathode discharges and breakdown
//		          in argon: surface and gas phase production of
//		          secondary electrons",
//		   journal="Plasma Sources Sci. Technol."i,
//		   year="1999",
//		   volume="8",
//		   pages="R21-R44"}
//
// Note:	If eps > 100 eV, the 100 eV crosssection is returned.

#ifndef H_ARGON_H
#define H_ARGON_H

#include "crosssec.h"
#include "physconst.h"
#include <cmath>

namespace Argon {

// Mass in kilogram.
double Mass() { return 39.948 * PhysConst::AMU; }

// ionization energy in Joule
double E_ion() { return 15.759 * PhysConst::e; }

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
// Now some implementations.
// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

// Crosssection (m^2) for elastic "e + Ar -> e + Ar" collision as a
// function of the electron energy in Joule.
class Elastic : public CrossSection {
public:
  Elastic() : CrossSection(0.0) {}

protected:
  double get_above_th(double eps) const { return (1e-19); }
};

// Crosssection (m^2) for inelastic "e + Ar -> e + Ar*" collision as
// a function of the electron energy in Joule.
class Inelastic : public CrossSection {
public:
  Inelastic() : CrossSection(5.0 * PhysConst::eV) {}

protected:
  double get_above_th(double eps) const { return 200e-22; }
};

// Crosssection (m^2) for ionisation "e + Ar -> 2e + Ar+" collision as
// a function of the electron energy in Joule.
class Ionisation : public CrossSection {
public:
  Ionisation() : CrossSection(15.8 * PhysConst::eV) {}

protected:
  double get_above_th(double eps) const {
    using std::pow;
    const double eV = 1.6022e-19;
    eps /= eV; // the lines below want eps in eV.
    if (eps > 100.0)
      eps = 100.0;
    double sigma = 0.0;
    if (eps >= 15.8) {
      sigma = 970.0 * (eps - 15.8) / pow((70.0 + eps), 2.0) +
              0.06 * pow((eps - 15.8), 2.0) * exp(-eps / 9.0);
    }
    return (1e-20 * sigma);
  }
};

// Crosssection (m^2) for charge exchange "Ar+ + Ar -> Ar + Ar+"
// collision.
class ChargeExchange : public CrossSection {
public:
  ChargeExchange() : CrossSection(0.0) {}

protected:
  double get_above_th(double eps) const {
    const double eV = 1.6022e-19;
    eps /= eV; // the lines below want eps in eV.
    if (eps > 100.0)
      eps = 100.0;
    double sigma;
    if (eps < 4.0) {
      sigma = -2.95e-19 * std::sqrt(eps) + 10.65e-19;
    } else {
      sigma = 2.0e-19 + 5.5e-19 / (std::sqrt(eps) + 1.0e-30);
    }
    return sigma;
  }
};

// Crosssection (m^2) for elastic "Ar+ + Ar -> Ar+ + Ar"
// collision.
class IonElastic : public CrossSection {
public:
  IonElastic() : CrossSection(0.0) {}

protected:
  double get_above_th(double eps) const {
    const double eV = 1.6022e-19;
    eps /= eV; // the lines below want eps in eV.
    if (eps > 100.0)
      eps = 100.0;
    double sigma;
    if (eps < 4.0) {
      sigma = -2.0e-19 * std::sqrt(eps) + 7.8e-19;
    } else {
      sigma = 1.8e-19 + 4.0e-19 / (std::sqrt(eps) + 1.0e-30);
    }
    return sigma;
  }
};

} // namespace Argon

#endif // H_ARGON_H
