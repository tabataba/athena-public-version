#ifndef CONDENSATION_HPP_
#define CONDENSATION_HPP_

// Athena++ headers
#include "../athena.hpp"

class MeshBlock;
template<typename T> class AthenaArray;

void GasGasSolidNH4SH(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate);

void GasCloudIdeal(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate);
void LiquidSolidIdeal(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate);

void GasLiquidH2O(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate);
void GasSolidH2O(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate);
void LiquidSolidH2O(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate);

void GasLiquidNH3(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate);
void GasSolidNH3(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate);
void LiquidSolidNH3(MeshBlock *pmb, Real const time, Reaction const& rc,
  AthenaArray<Real> const& prim, AthenaArray<Real>& rate);

#endif
