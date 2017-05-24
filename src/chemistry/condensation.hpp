#ifndef CONDENSATION_HPP_
#define CONDENSATION_HPP_

// Athena++ headers
#include "../athena.hpp"

class MeshBlock;
template<typename T> class AthenaArray;

void GasGasSolidNH4SH(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, int i, AthenaArray<Real>& rate, int r);

void GasCloudIdeal(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, int i, AthenaArray<Real>& rate, int r);
void LiquidSolidIdeal(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, int i, AthenaArray<Real>& rate, int r);

#endif
