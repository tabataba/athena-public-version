//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file adiabatic_hydro.cpp
//  \brief implements functions in class HeterogeneousHydro for adiabatic hydrodynamics`

// C/C++ headers
#include <cmath>   // sqrt()
#include <cfloat>  // FLT_MIN
#include <cstdlib>  // atof

// Athena++ headers
#include "eos.hpp"
#include "../hydro/hydro.hpp"
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../math_funcs.hpp"
#include "../misc.hpp"

// HeterogeneousHydro constructor

HeterogeneousHydro::HeterogeneousHydro(MeshBlock *pmb, ParameterInput *pin):
  EquationOfState(pmb)
{
  std::stringstream msg;

  std::vector<std::string> str;
  // read gamma, size should be NCOMP
  SplitString(pin->GetString("hydro", "gamma"), str);
  if (str.size() != NCOMP) {
    msg << "### FATAL ERROR in heterogeneous_hydro.cpp::HeterogeneousHydro: number of gases in 'gamma' does not equal NGAS" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  for (int i = 0; i < NCOMP; ++i)
    kappa_[i] = atof(str[i].c_str()) - 1.;

  // read cv, size should be NCOMP, unit is [J/g], convert to [J/kg]
  SplitString(pin->GetString("hydro", "cv"), str);
  if (str.size() != NCOMP) {
    msg << "### FATAL ERROR in heterogeneous_hydro.cpp::HeterogeneousHydro: number of species in 'cv' does not equal NCOMP" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  for (int i = 0; i < NCOMP; ++i)
    cv_[i] = atof(str[i].c_str())*1000.;

  // read latent, size should be NCOMP - NGAS, unit is [J/g], convert to [J/kg]
  if (NCOMP - NGAS > 0) {
    SplitString(pin->GetString("hydro", "latent"), str);
    if (str.size() != NCOMP - NGAS) {
      msg << "### FATAL ERROR in heterogeneous_hydro.cpp::HeterogeneousHydro: number of clouds in 'latent' does not equal " << NCOMP - NGAS << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
  }
  for (int i = 0; i < NGAS; ++i)
    latent_[i] = 0.;
  for (int i = NGAS; i < NCOMP; ++i)
    latent_[i] = atof(str[i].c_str())*1000.;

  //density_floor_  = pin->GetOrAddReal("hydro","dfloor",(1024*(FLT_MIN)));
  density_floor_  = 0.;
  pressure_floor_ = pin->GetOrAddReal("hydro","pfloor",(1024*(FLT_MIN)));
}

// destructor

HeterogeneousHydro::~HeterogeneousHydro()
{
}

//----------------------------------------------------------------------------------------
// \!fn void HeterogeneousHydro::ConservedToPrimitive(AthenaArray<Real> &cons,
//           const AthenaArray<Real> &prim_old, const FaceField &b,
//           AthenaArray<Real> &prim, AthenaArray<Real> &bcc, Coordinates *pco,
//           int is, int ie, int js, int je, int ks, int ke)
// \brief Converts conserved into primitive variables in adiabatic hydro.

void HeterogeneousHydro::ConservedToPrimitive(AthenaArray<Real> &cons,
  const AthenaArray<Real> &prim_old, const FaceField &b, AthenaArray<Real> &prim,
  AthenaArray<Real> &bcc, Coordinates *pco, int is, int ie, int js, int je, int ks, int ke)
{
  int nthreads = pmy_block_->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) num_threads(nthreads)
{
  for (int k=ks; k<=ke; ++k){
#pragma omp for schedule(dynamic)
  for (int j=js; j<=je; ++j){
#pragma simd
    for (int i=is; i<=ie; ++i){
      Real& m1 = cons(IM1,k,j,i);
      Real& m2 = cons(IM2,k,j,i);
      Real& m3 = cons(IM3,k,j,i);
      Real rho = 0., rck = 0., rc = 0., ohr;

      // apply density floor, without changing momentum or energy
      for (int c = 0; c < NCOMP; ++c) {
        cons(c,k,j,i) = _max(cons(c,k,j,i), density_floor_);
        rck += cons(c,k,j,i)*cv_[c]*kappa_[c];
        rc  += cons(c,k,j,i)*cv_[c];
        rho += cons(c,k,j,i);
      }
      ohr = 1./rho;

      // molar mixing ratio
      for (int c = 0; c < NCOMP; ++c)
        prim(c,k,j,i) = cons(c,k,j,i)*cv_[c]*kappa_[c]/rck;

      // subtract cloud components when calculate pressure
      for (int c = NGAS; c < NCOMP; ++c)
        rck -= cons(c,k,j,i)*cv_[c]*kappa_[c];

      // apply pressure floor, correct total energy
      prim(IPR,k,j,i) = rck/rc*(cons(IEN,k,j,i)-0.5*ohr*(_sqr(m1)+_sqr(m2)+_sqr(m3)));
      cons(IEN,k,j,i) = (prim(IPR,k,j,i) > pressure_floor_) ?  cons(IEN,k,j,i) :
        (pressure_floor_/rck*rc + 0.5*ohr*(_sqr(m1)+_sqr(m2)+_sqr(m3)));
      prim(IPR,k,j,i) = _max(prim(IPR,k,j,i), pressure_floor_);

      // temperature
      prim(IT,k,j,i) = prim(IPR,k,j,i)/rck;

      // velocity
      prim(IVX,k,j,i) = ohr*m1;
      prim(IVY,k,j,i) = ohr*m2;
      prim(IVZ,k,j,i) = ohr*m3;
    }
  }}
}

  return;
}


//----------------------------------------------------------------------------------------
// \!fn void HeterogeneousHydro::PrimitiveToConserved(const AthenaArray<Real> &prim,
//           const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
//           int is, int ie, int js, int je, int ks, int ke);
// \brief Converts primitive variables into conservative variables

void HeterogeneousHydro::PrimitiveToConserved(const AthenaArray<Real> &prim,
     const AthenaArray<Real> &bc, AthenaArray<Real> &cons, Coordinates *pco,
     int is, int ie, int js, int je, int ks, int ke)
{
  int nthreads = pmy_block_->pmy_mesh->GetNumMeshThreads();
#pragma omp parallel default(shared) num_threads(nthreads)
{
  Real mixr[NGAS];

  for (int k=ks; k<=ke; ++k){
#pragma omp for schedule(dynamic)
  for (int j=js; j<=je; ++j){
#pragma simd
    for (int i=is; i<=ie; ++i){
      Real const& v1 = prim(IVX,k,j,i);
      Real const& v2 = prim(IVY,k,j,i);
      Real const& v3 = prim(IVZ,k,j,i);

      // total gas mixing ratios
      Real xt = 1.;
      for (int c = NGAS; c < NCOMP; ++c)
        xt -= prim(c,k,j,i);

      // gas density
      Real x1 = 0., rho = 0., rc = 0.;
      for (int c = 1; c < NGAS; ++c) {
        cons(c,k,j,i) = prim(c,k,j,i)/xt*prim(IPR,k,j,i)/(cv_[c]*kappa_[c]*prim(IT,k,j,i));
        x1 += prim(c,k,j,i);
        rho += cons(c,k,j,i);
        rc += cons(c,k,j,i)*cv_[c];
      }
      cons(0,k,j,i) = (xt-x1)/xt*prim(IPR,k,j,i)/(cv_[0]*kappa_[0]*prim(IT,k,j,i));
      rho += cons(0,k,j,i);
      rc += cons(0,k,j,i)*cv_[0];

      // cloud density
      for (int c = NGAS; c < NCOMP; ++c) {
        cons(c,k,j,i) = prim(c,k,j,i)*cv_[0]*kappa_[0]/((xt-x1)*cv_[c]*kappa_[c])*cons(0,k,j,i);
        rho += cons(c,k,j,i);
        rc += cons(c,k,j,i)*cv_[c];
      }

      // momentum
      cons(IM1,k,j,i) = rho*v1;
      cons(IM2,k,j,i) = rho*v2;
      cons(IM3,k,j,i) = rho*v3;

      // total energy
      cons(IEN,k,j,i) = prim(IT,k,j,i)*rc + 0.5*rho*(_sqr(v1)+_sqr(v2)+_sqr(v3));
    }
  }}
}
  return;
}

//----------------------------------------------------------------------------------------
// \!fn Real HeterogeneousHydro::SoundSpeed(Real const prim[])
// \brief returns adiabatic sound speed given vector of primitive variables

Real HeterogeneousHydro::SoundSpeed(Real const prim[])
{
  Real r1, x1 = 0., rho = 0., rck = 0., rc = 0.;
  for (int c = 1; c < NCOMP; ++c) {
    r1 = prim[c]*prim[IPR]/(cv_[c]*kappa_[c]*prim[IT]);
    x1  += prim[c];
    rho += r1;
    rc  += r1*cv_[c];
    rck += r1*cv_[c]*kappa_[c];
  }
  r1 = (1.-x1)*prim[IPR]/(cv_[0]*kappa_[0]*prim[IT]);
  rho += r1;
  rc  += r1*cv_[0];
  rck += r1*cv_[0]*kappa_[0];
  return sqrt((rck/rc+1)*prim[IPR]/rho);
}
