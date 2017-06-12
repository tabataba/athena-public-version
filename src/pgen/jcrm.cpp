// C++ headers
#include <sstream>
#include <stdexcept>
#include <cstdlib>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../parameter_input.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../hydro/srcterms/hydro_srcterms.hpp"
#include "../mesh/mesh.hpp"
#include "../chemistry/reaction.hpp"
#include "../chemistry/molecule.hpp"
#include "../chemistry/condensation.hpp"
#include "../utils/utils.hpp"
#include "../globals.hpp"
#include "../misc.hpp"

Real mutop, mubot, cpbot, cptop, xbot[NCOMP], ttop, grav;

void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  for (int k = ks; k <= ke; ++k)
    for (int j = 1; j <= NGHOST; ++j)
      for (int i = is; i <= ie; ++i) {
        prim(IVX,k,js-j,i) = prim(IVX,k,js-j,i);
        prim(IVY,k,js-j,i) = -prim(IVY,k,js-j,i);
        // adiabatic projection
        prim(IT,k,js-j,i)  = prim(IT,k,js+j-1,i) + mubot*grav/cpbot*(2*j-1)*pco->dx2f(js);
        prim(IPR,k,js-j,i) = prim(IPR,k,js+j-1,i)
          *pow(prim(IT,k,js-j,i)/prim(IT,k,js+j-1,i), cpbot/Globals::Rgas);
        for (int n = 1; n < NCOMP; ++n)
          prim(n,k,js-j,i) = xbot[n];
      }
}

void ProjectPressureOuterX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  for (int k = ks; k <= ke; ++k)
    for (int j = 1; j <= NGHOST; ++j)
      for (int i = is; i <= ie; ++i) {
        prim(IVX,k,je+j,i) = prim(IVX,k,je-j+1,i);
        prim(IVY,k,je+j,i) = -prim(IVY,k,je-j+1,i);
        // isothermal projection
        prim(IT,k,je+j,i) = ttop;
        prim(IPR,k,je+j,i) = prim(IPR,k,je-j+1,i)
          *exp(-mutop*grav/(Globals::Rgas*ttop)*(2*j-1)*pco->dx2f(j));
        for (int n = 1; n < NCOMP; ++n)
          prim(n,k,je+j,i) = 0.;
        // adiabatic projection
        //prim(IT,k,je+j,i)  = prim(IT,k,je-j+1,i) - mutop*grav/cptop*(2*j-1)*pco->dx2f(je);
        //prim(IPR,k,je+j,i) = prim(IPR,k,je-j+1,i)
        //  *pow(prim(IT,k,je+j,i)/prim(IT,k,je-j+1,i), cpbot/Globals::Rgas);
        //for (int n = 1; n < NCOMP; ++n)
        //  prim(n,k,je+j,i) = 0.;
      }
}

/*void JovianAtmMicrophysics(MeshBlock *pmb, const Real time, const Real dt, const int step,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  if (step == 1) return;
  Reaction const& r1 = pmb->prg->GetReaction("r1");
  Reaction const& r2 = pmb->prg->GetReaction("r2");
  EquationOfState *peos = pmb->peos;

  Real rate, prim[NHYDRO];
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real xt = 1.;
        for (int n = NGAS; n < NCOMP; ++n) {
          prim[n] = prim(n,k,j,i);
          xt -= prim[n];
        }
        for (int n = 0; n < NGAS; ++n)
          prim[n] = prim(n,k,j,i);
        for (int n = NCOMP; n < NHYDRO; ++n)
          prim[n] = prim(n,k,j,i);
        Real rho1 = prim[IPR]/(xt*Globals::Rgas*prim[IT]);
        rate = GasCloudIdeal(r1, prim, time);
        for (int n = 0; n < 2; ++n)
          cons(r1.reactor[n],k,j,i) += rho1*rate*r1.measure[n]*peos->GetMu(r1.reactor[n]);
        rate = GasCloudIdeal(r2, prim, time);
        for (int n = 0; n < 2; ++n)
          cons(r2.reactor[n],k,j,i) += rho1*rate*r2.measure[n]*peos->GetMu(r2.reactor[n]);
      }
}*/


void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserBoundaryFunction(INNER_X2, ProjectPressureInnerX2);
  EnrollUserBoundaryFunction(OUTER_X2, ProjectPressureOuterX2);
  //EnrollUserExplicitSourceFunction(JovianAtmMicrophysics);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  std::stringstream msg;
  // Step 1. define molecules
  pmol = new Molecule(pin);
  if (Molecule::ntotal != NCOMP) {
    msg << "### FATAL ERROR in ProblemGenerator: Molecule::ntotal != NMOL" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // Step 2. define reactions
  prg = new ReactionGroup(this, "condensation");
  prg->AddReaction(pin, "chemistry", "r1", pmol);
  prg->AddReaction(pin, "chemistry", "r2", pmol);

  // debug
  //Real aa[NHYDRO] = {280.,5.E-3,0.,0.,0.,0.,0.,0.,7.E5};
  //prg->EquilibrateUV(aa, peos);
  //for (int n = 0; n < NHYDRO; ++n) std::cout << aa[n] << std::endl;
  //exit(1.);

  // Step 3. define boundary conditions
  std::vector<std::string> xbot_str;
  Real ptop = pin->GetReal("problem", "ptop");
  Real pref = pin->GetReal("problem", "pref");
  Real tref = pin->GetReal("problem", "tref");
  SplitString(pin->GetString("problem", "xbot"), xbot_str);
  Real prim[NHYDRO];

  ttop = tref;
  grav = - phydro->psrc->GetG2();
  for (size_t n = 0; n < xbot_str.size(); ++n)
    xbot[n + 1] = atof(xbot_str[n].c_str());
  for (int n = NGAS; n < NCOMP; ++n) xbot[n] = 0.;
  for (int n = 1; n < NCOMP; ++n)
    prim[n] = xbot[n];
  mubot = peos->Mass(prim);
  cpbot = peos->HeatCapacityP(prim);
  for (int n = 1; n < NCOMP; ++n)
    prim[n] = 0.;
  mutop = peos->Mass(prim);
  cptop = peos->HeatCapacityP(prim);

  // Step 4. construct adiabatic T-P profile from top down
  // 4.1) T-P profile outside of this MeshBlock
  Real dz = pcoord->dx2f(je), ztop = pmy_mesh->mesh_size.x2max;
  Real mu = mutop, cp = peos->HeatCapacityP(prim);
  prim[IPR] = ptop;
  prim[IT ] = ttop;
  for (; ztop > pcoord->x2v(je); ztop -= dz) {
    if (prim[IPR] < pref) {  // isothermal atm
      prim[IPR] *= exp(mu*grav*dz/(Globals::Rgas*tref));
    } else {  // adiabatic atm
      prim[IPR] *= pow(1. + mu*grav*dz/(cp*prim[IT]), cp/Globals::Rgas);
      prim[IT ] += mu*grav*dz/cp;
    }
    prim[IVX] = prim[IVY] = prim[IVZ] = 0.;

    // remove condensates and update cp, mu
    for (int n = 1; n < NCOMP; ++n) prim[n] = xbot[n];
    int status = prg->EquilibrateTP(prim);
    if (status == 0) { // not converged
      msg << "### FATAL ERROR in ProblemGenerator: EquilibrateTP does not converge."
          << std::endl;
      msg << "prim after: ";
      for (int n = 0; n < NHYDRO; ++n)
        msg << prim[n] << " ";
      throw std::runtime_error(msg.str().c_str());
    }
    for (int n = NGAS; n < NCOMP; ++n) prim[n] = 0.;
    cp = peos->HeatCapacityP(prim);
    mu = peos->Mass(prim);
  }
  ztop += dz;

  // 4.2) T-P profile of this MeshBlock and add pertubation
  for (int j = je; j >= js; --j) {
    if (j == je) dz = ztop - pcoord->x2v(j);
    else dz = pcoord->dx2f(j);

    if (prim[IPR] < pref) { // isothermal atm
      prim[IPR] *= exp(mu*grav*dz/(Globals::Rgas*tref));
    } else { // adiabatic atm
      prim[IPR] *= pow(1. + mu*grav*dz/(cp*prim[IT]), cp/Globals::Rgas);
      prim[IT ] += mu*grav*dz/cp;
    }
    prim[IVX] = prim[IVY] = prim[IVZ] = 0.;

    // remove condensates and update cp, mu
    for (int n = 1; n < NCOMP; ++n) prim[n] = xbot[n];
    int status = prg->EquilibrateTP(prim);
    if (status == 0) { // not converged
      msg << "### FATAL ERROR in ProblemGenerator: EquilibrateTP does not converge."
          << std::endl;
      msg << "prim after: ";
      for (int n = 0; n < NHYDRO; ++n)
        msg << prim[n] << " ";
      msg << std::endl;
      throw std::runtime_error(msg.str().c_str());
    }
    for (int n = NGAS; n < NCOMP; ++n) prim[n] = 0.;
    cp = peos->HeatCapacityP(prim);
    mu = peos->Mass(prim);

    // copy to primitive variables
    for (int n = 0; n < NHYDRO; ++n)
      phydro->w(n,ks,j,is) = prim[n];
  }

  // 4.3) propagate T-P profile to the whole domain, add noise
  // pertubation wavelength
  Real kx = 20.*(PI)/(pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min);
  Real ky = 20.*(PI)/(pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min);
  long int iseed = -1;
  for (int n = 0; n < NHYDRO; ++n)
    for (int k = ks; k <= ke; ++k)
      for (int j = js; j <= je; ++j)
        for (int i = is; i <= ie; ++i) {
          phydro->w(n,k,j,i) = phydro->w(n,ks,j,is);
          if (n == IVY)
            phydro->w(n,k,j,i) = 100.*(ran2(&iseed)-0.5)
              *(1.0+cos(kx*pcoord->x1v(i)))*(1.0+cos(ky*pcoord->x2v(j)))/4.0;
        }

  // Step 5. convert primitive variables to conserved variables
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}

int ReactionGroup::EquilibrateUV(Real prim[], EquationOfState *peos) const
{
  int iter = 1, max_iter = 10;
  Real cv = peos->HeatCapacityV(prim);
  Real dt, rate;
  bool reacted = false;
  Reaction const& r1 = rts_[0];
  Reaction const& r2 = rts_[1];

  do {
    Real xti = 1.;
    for (int n = NGAS; n < NCOMP; ++n)
      xti -= prim[n];
    Real latent = 0.;

    // H2O - H2O(s)
    rate = GasCloudIdeal(r1, prim, peos->latent_[r1.reactor[1]]/cv);
    //std::cout << "rate = " << rate << std::endl;
    if (fabs(rate) > TINY_NUMBER) {
      for (int n = 0; n < 2; ++n)
        prim[r1.reactor[n]] += rate*r1.measure[n];
      latent += 0.5*rate*peos->latent_[r1.reactor[1]]*r1.measure[1];
      reacted = true;
    }

    // NH3 - NH3(s)
    rate = GasCloudIdeal(r2, prim, peos->latent_[r1.reactor[1]]/cv);
    //std::cout << "rate = " << rate << std::endl;
    if (fabs(rate) > TINY_NUMBER) {
      for (int n = 0; n < 2; ++n)
        prim[r2.reactor[n]] += rate*r2.measure[n];
      latent += 0.5*rate*peos->latent_[r2.reactor[1]]*r2.measure[1];
      reacted = true;
    }

    if (!reacted) return 1; // not reacted
    Real xtf = 1.;
    for (int n = NGAS; n < NCOMP; ++n)
      xtf -= prim[n];
    dt = -latent/cv;
    //std::cout << "dt = " << dt << std::endl;
    prim[IPR] *= xtf/xti*(prim[IT]+dt)/prim[IT];
    prim[IT]  += dt;
    iter++;
  } while (fabs(dt) > 0.1 && iter < max_iter);

  if (iter == max_iter) return 0; // not converged
  else return 2; // converged and reacted
}

int ReactionGroup::EquilibrateTP(Real prim[]) const
{
  int iter = 1, max_iter = 10;
  Real rate, max_rate;
  bool reacted = false;
  Reaction const& r1 = rts_[0];
  Reaction const& r2 = rts_[1];

  do {
    // H2O - H2O(s)
    rate = GasCloudIdeal(r1, prim);
    if (fabs(rate) > TINY_NUMBER) {
      for (int n = 0; n < 2; ++n)
        prim[r1.reactor[n]] += rate*r1.measure[n];
      reacted = true;
    }
    max_rate = fabs(rate);

    // NH3 - NH3(s)
    rate = GasCloudIdeal(r2, prim);
    if (fabs(rate) > TINY_NUMBER) {
      for (int n = 0; n < 2; ++n)
        prim[r2.reactor[n]] += rate*r2.measure[n];
      reacted = true;
    }

    if (!reacted) return 1; // not reacted
    max_rate = std::max(max_rate, fabs(rate));
  } while (max_rate > TINY_NUMBER && iter < max_iter);

  if (iter == max_iter) return 0; // not converged
  else return 2; // converged and reacted
}
