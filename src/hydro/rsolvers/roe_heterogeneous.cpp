//! \file  roe_shallow_water.cpp
//  \brief Roe's linearized Riemann solver for shallow water model

// C/C++ headers
#include <cmath>      // sqrt()
#include <iostream>
#include <algorithm>  // min()

// Athena++ headers
#include "../hydro.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../eos/eos.hpp"
#include "../../math_funcs.hpp"
#include "../../globals.hpp"

inline Real Enthalpy(AthenaArray<Real> &w, int i, Real x1, Real const cv[], Real const latent[])
{
  Real cp = (cv[0] + Globals::Rgas)*x1, lt = 0., ke;
  for (int n = 1; n < NGAS; ++n)
    cp += (cv[n] + Globals::Rgas)*w(n,i);
  for (int n = NGAS; n < NCOMP; ++n) {
    cp += cv[n]*w(n,i);
    lt += latent[n]*w(n,i);
  }
  ke = 0.5*(_sqr(w(IVX,i)) + _sqr(w(IVY,i)) + _sqr(w(IVZ,i)));
  return cp*w(IT,i) + ke + lt;
}

inline Real TotalDensity(AthenaArray<Real> &w, int i, Real const mu[], Real rhon[], Real *x1)
{
  // total gas mixing ratio
  Real xt = 1., rho = 0.;
  for (int n = NGAS; n < NCOMP; ++n)
    xt -= w(n,i);
  *x1 = xt;
  for (int n = 1; n < NGAS; ++n) {
    rhon[n] = w(n,i)/xt*w(IPR,i)*mu[n]/(Globals::Rgas*w(IT,i));
    *x1 -= w(n,i);
    rho += rhon[n];
  }
  rhon[0] = (*x1)/xt*w(IPR,i)*mu[0]/(Globals::Rgas*w(IT,i));
  rho += rhon[0];
  for (int n = NGAS; n < NCOMP; ++n) {
    rhon[n] = w(n,i)/(*x1) * mu[n]/mu[0];
    rho += rhon[n];
  }
  return rho;
}

void Hydro::RiemannSolver(int const k, int const j, int const il, int const iu,
    int const ivx, AthenaArray<Real> const& bx, AthenaArray<Real> &wl,
    AthenaArray<Real> &wr, AthenaArray<Real> &flx)
{
  HeterogeneousHydro *peos = pmy_block->peos;
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  Real wave[3][NHYDRO], speed[3];
  Real ubar, vbar, wbar, kbar, cbar, hbar, hl, hr, rhobar, x1l, x1r;
  Real qbar[NCOMP], rholn[NCOMP], rhorn[NCOMP];
  Real alpha[NHYDRO];
  Real rhol, rhor, sqrhol, sqrhor, isqrho;
  Real du, dv, dw, dp;

  for (int i = il; i <= iu; ++i) {
    rhol = TotalDensity(wl, i, peos->mu_, rholn, &x1l);
    rhor = TotalDensity(wr, i, peos->mu_, rhorn, &x1r);
    sqrhol = sqrt(rhol);
    sqrhor = sqrt(rhor);
    isqrho = 1./(sqrhol + sqrhor);
    rhobar = sqrhol*sqrhor;

    // u,v,w
    ubar = isqrho*(wl(ivx,i)*sqrhol + wr(ivx,i)*sqrhor);
    vbar = isqrho*(wl(ivy,i)*sqrhol + wr(ivy,i)*sqrhor);
    wbar = isqrho*(wl(ivz,i)*sqrhol + wr(ivz,i)*sqrhor);

    hl = Enthalpy(wl, i, x1l, peos->cv_, peos->latent_);
    hr = Enthalpy(wr, i, x1r, peos->cv_, peos->latent_);
    hbar = isqrho*(hl*sqrhol + hr*sqrhor);

    // qbar, molar mixing ratio
    qbar[0] = isqrho*(x1l*sqrhol + x1r*sqrhor);
    for (int n = 1; n < NCOMP; ++n)
      qbar[n] = isqrho*(wl(n,i)*sqrhol + wr(n,i)*sqrhor);

    // kbar
    Real cv = 0.;
    for (int n = 0; n < NGAS; ++n)
      cv += peos->cv_[n]*qbar[n];
    Real xt = 1.;
    for (int n = NGAS; n < NCOMP; ++n) {
      xt -= qbar[n];
      cv += peos->cv_[n]*qbar[n];
    }
    kbar = xt*Globals::Rgas/cv;

    // cbar
    Real latent = 0.;
    for (int n = NGAS; n < NCOMP; ++n)
      latent += peos->latent_[n]*qbar[n];
    cbar = sqrt(kbar*(hbar - 0.5*(_sqr(ubar) + _sqr(vbar) + _sqr(wbar)) - latent));

    // molar mixing ratio to mass mixing ratio
    Real qsum = 0.;
    for (int n = 0; n < NCOMP; ++n)
      qsum += peos->mu_[n]*qbar[n];
    for (int n = 0; n < NCOMP; ++n)
      qbar[n] /= qsum;

    // primitive variable difference
    du = wr(ivx,i) - wl(ivx,i);
    dv = wr(ivy,i) - wl(ivy,i);
    dw = wr(ivz,i) - wl(ivz,i);
    dp = wr(IPR,i) - wl(IPR,i);

    // coefficients of eigenvectors
    alpha[ivx] = - 0.5*rhobar/cbar*du + 0.5*dp/_sqr(cbar);
    alpha[IPR] = + 0.5*rhobar/cbar*du + 0.5*dp/_sqr(cbar);
    alpha[ivy] = rhobar*dv;
    alpha[ivz] = rhobar*dw;
    for (int n = 0; n < NCOMP; ++n)
      alpha[n] = rhorn[n] - rholn[n] - dp/_sqr(cbar)*qbar[n];

    // wave 1, u-c
    for (int n = 0; n < NCOMP; ++n)
      wave[0][n] = alpha[ivx]*qbar[n];
    wave[0][ivx] = alpha[ivx]*(ubar - cbar);
    wave[0][ivy] = alpha[ivx]*vbar;
    wave[0][ivz] = alpha[ivx]*wbar;
    wave[0][IPR] = alpha[ivx]*(hbar - ubar*cbar);

    // wave 3, u+c
    for (int n = 0; n < NCOMP; ++n)
      wave[2][n] = alpha[IPR]*qbar[n];
    wave[2][ivx] = alpha[IPR]*(ubar + cbar);
    wave[2][ivy] = alpha[IPR]*vbar;
    wave[2][ivz] = alpha[IPR]*wbar;
    wave[2][IPR] = alpha[IPR]*(hbar + ubar*cbar);

    // wave 2, u
    Real alpha0 = 0.;
    for (int n = 0; n < NCOMP; ++n) {
      wave[1][n] = alpha[n];
      alpha0 += alpha[n];
    }
    wave[1][ivx] = ubar*alpha0;
    wave[1][ivy] = vbar*alpha0 + alpha[ivy];
    wave[1][ivz] = wbar*alpha0 + alpha[ivz];
    wave[1][IPR] = rhobar*(hr - hl - ubar*du) + alpha0*hbar - dp;

    // speed
    speed[0] = fabs(ubar - cbar);
    speed[1] = fabs(ubar);
    speed[2] = fabs(ubar + cbar);

    // flux
    for (int n = 0; n < NCOMP; ++n)
      flx(n,i) = 0.5*(rholn[n]*wl(ivx,i) + rhorn[n]*wr(ivx,i));
    flx(ivx,i) = 0.5*(rhol*_sqr(wl(ivx,i)) + wl(IPR,i)
                    + rhor*_sqr(wr(ivx,i)) + wr(IPR,i));
    flx(ivy,i) = 0.5*(rhol*wl(ivx,i)*wl(ivy,i)
                    + rhor*wr(ivx,i)*wr(ivy,i));
    flx(ivz,i) = 0.5*(rhol*wl(ivx,i)*wl(ivz,i)
                    + rhor*wr(ivx,i)*wr(ivz,i));
    flx(IPR,i) = 0.5*(rhol*hl*wl(ivx,i)
                    + rhor*hr*wr(ivx,i));
    for (int r = 0; r < 3; ++r)
      for (int n = 0; n < NHYDRO; ++n)
        flx(n,i) -= 0.5*speed[r]*wave[r][n];
  }
}
