// C++ header
#include <cmath>

// Athena++ headers
#include "condensation.hpp"
#include "../math_funcs.hpp"
#include "reaction.hpp"
#include "sat_vapor_pres.hpp"

Real GasGasSolidNH4SH(Reaction const& rc, Real const prim[NHYDRO], Real time)
{
  Real rate;
  Real xnh3 = prim[rc.reactor[0]];
  Real xh2s = prim[rc.reactor[1]];
  Real xnh4sh = prim[rc.reactor[2]];
  Real s = SatVaporPresNH4SHLewis(prim[IT])/_sqr(prim[IPR]);
  if (s > 1.) return -xnh4sh;
  Real xt = 1.;
  for (int g = NGAS; g < NCOMP; ++g) xt -= prim[g];
  rate = (xnh3+xh2s-4.*s*xt-sqrt(_sqr(xnh3-xh2s)+4.*s*(xt-2.*xnh3)*(xt-2.*xh2s)))/(2.*(1.-4.*s));
  if (rate > 0.) rate = _min(_min(rate, xnh3), xh2s);
  if (rate < 0.) rate = _max(rate, -xnh4sh);
  return rate;
}

Real GasCloudIdeal(Reaction const& rc, Real const prim[NHYDRO], Real time)
{
  Real rate;
  Real tr = rc.coeff[0],
       pr = rc.coeff[1],
       beta = rc.coeff[2],
       gamma = rc.coeff[3];
  Real x1 = prim[rc.reactor[0]];
  Real xc = prim[rc.reactor[1]];
  Real temp = prim[IT];
  Real s = SatVaporPresIdeal(temp/tr,pr,beta,gamma)/prim[IPR];
  if (s > 1.) return -xc;
  Real xt = 1.;
  for (int g = NGAS; g < NCOMP; ++g) xt -= prim[g];
  rate = (x1-s*xt)/(1.-s);
  //std::cout << rc.reactor[0] << " " << rc.reactor[1] << " " << rate << " " << s << std::endl;
  if (rate > 0.) rate = _min(rate, x1);
  if (rate < 0.) rate = _max(rate, -xc);
  return rate;
}

Real LiquidSolidIdeal(Reaction const& rc, Real const prim[NHYDRO], Real time)
{
  Real rate;
  Real tr = rc.coeff[0];
  if (prim[IT] < tr)
    rate = prim[rc.reactor[0]];
  else if (prim[IT] > tr)
    rate = - prim[rc.reactor[1]];
  return rate;
}
