// C++ header
#include <cmath>

// Athena++ headers
#include "condensation.hpp"
#include "../athena_arrays.hpp"
#include "../math_funcs.hpp"
#include "../mesh/mesh.hpp"
#include "reaction.hpp"
#include "sat_vapor_pres.hpp"

void GasGasSolidNH4SH(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, int i1, int i2, AthenaArray<Real>& rate, int r)
{
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = i1; i <= i2; ++i) {
        Real xnh3 = prim(rc.reactor[0],k,j,i);
        Real xh2s = prim(rc.reactor[1],k,j,i);
        Real xnh4sh = prim(rc.reactor[2],k,j,i);
        Real s = SatVaporPresNH4SHLewis(prim(IT,k,j,i))/_sqr(prim(IPR,k,j,i));
        Real xt = 1.;
        for (int g = NGAS; g < NCOMP; ++g) xt -= prim(g,k,j,i);
        rate(r,k,j,i) = (xnh3+xh2s-4.*s*xt-sqrt(_sqr(xnh3-xh2s)+4.*s*(xt-2.*xnh3)*(xt-2.*xh2s)))/(2.*(1.-4.*s));
        if (rate(r,k,j,i) < 0.) rate(r,k,j,i) = _max(rate(r,k,j,i), -xnh4sh);
      }
}

void GasCloudIdeal(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, int i1, int i2, AthenaArray<Real>& rate, int r)
{
  Real tr = rc.coeff[0],
       pr = rc.coeff[1],
       beta = rc.coeff[2],
       gamma = rc.coeff[3];

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = i1; i <= i2; ++i) {
        Real x1 = prim(rc.reactor[0],k,j,i);
        Real xc = prim(rc.reactor[1],k,j,i);
        Real temp = prim(IT,k,j,i);
        Real s = SatVaporPresIdeal(temp/tr,pr,beta,gamma)/prim(IPR,k,j,i);
        Real xt = 1.;
        for (int g = NGAS; g < NCOMP; ++g) xt -= prim(g,k,j,i);
        rate(r,k,j,i) = (x1-s*xt)/(1.-s);
        if (rate(r,k,j,i) < 0.) rate(r,k,j,i) = _max(rate(r,k,j,i), -xc);
      }
}

void LiquidSolidIdeal(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, int i1, int i2, AthenaArray<Real>& rate, int r)
{
  Real tr = rc.coeff[0];

  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = i1; i <= i2; ++i) {
        if (prim(IT,k,j,i) < tr)
          rate(r,k,j,i) = prim(rc.reactor[0],k,j,i);
        else if (prim(IT,k,j,i) > tr)
          rate(r,k,j,i) = - prim(rc.reactor[1],k,j,i);
      }
}
