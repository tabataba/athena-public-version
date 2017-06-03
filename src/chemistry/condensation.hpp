#ifndef CONDENSATION_HPP_
#define CONDENSATION_HPP_

// Athena++ headers
#include "../athena.hpp"

Real GasGasSolidNH4SH(Reaction const& rc, Real const prim[NHYDRO], Real time);
Real GasCloudIdeal(Reaction const& rc, Real const prim[NHYDRO], Real time);
Real LiquidSolidIdeal(Reaction const& rc, Real const prim[NHYDRO], Real time);

#endif
