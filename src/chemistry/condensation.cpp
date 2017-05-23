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
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate)
{
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real xnh3 = prim(rc.reactor[0],k,j,i);
        Real xh2s = prim(rc.reactor[1],k,j,i);
        Real xnh4sh = prim(rc.reactor[2],k,j,i);
        Real svp = SatVaporPresNH4SHLewis(prim(IT,k,j,i));
        Real xt = 0.;
        for (int g = 0; g < NGAS; ++g) xt += prim(g,k,j,i);
        rate(rc.id,k,j,i) = ((xnh3+xh2s)-sqrt(_sqr(xnh3-xh2s)+4.*svp*(xt-2.*xnh3)*(xt-2.*xh2s)))/(2.*(1.-4.*svp));
        if (rate(rc.id,k,j,i) < 0.) rate(rc.id,k,j,i) = _max(rate(rc.id,k,j,i), -xnh4sh);
      }
}

void GasCloudIdeal(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate)
{
}

void LiquidSolidIdeal(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate)
{
}

void GasLiquidH2O(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate)
{
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real xh2og = prim(rc.reactor[0],k,j,i);
        Real xh2ol = prim(rc.reactor[1],k,j,i);
        Real svp = SatVaporPresH2OIdeal(prim(IT,k,j,i));
        Real xt = 0.;
        for (int g = 0; g < NGAS; ++g) xt += prim(g,k,j,i);
        rate(rc.id,k,j,i) = xh2og - xt*svp;
        if (rate(rc.id,k,j,i) < 0.) rate(rc.id,k,j,i) = _max(rate(rc.id,k,j,i), -xh2ol);
      }
}

void GasSolidH2O(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate)
{
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
      }
}

void LiquidSolidH2O(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate)
{
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
      }
}

void GasLiquidNH3(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate)
{
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
      }
}

void GasSolidNH3(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate)
{
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
      }
}

void LiquidSolidNH3(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate)
{
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
      }
}
