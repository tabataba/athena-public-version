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

Real mutop, mubot, cpbot, cptop, xbot[NCOMP], ttop, grav, crate;

void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  for (int k = ks; k <= ke; ++k)
    for (int j = 1; j <= NGHOST; ++j)
      for (int i = is; i <= ie; ++i) {
        prim(IVX,k,js-j,i) = prim(IVX,k,js+j-1,i);
        prim(IVY,k,js-j,i) = -prim(IVY,k,js+j-1,i);
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
          *exp(-mutop*grav/(Globals::Rgas*ttop)*(2*j-1)*pco->dx2f(je));
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

void ProjectPressureInnerX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  for (int k = 1; k <= NGHOST; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        prim(IVX,ks-k,j,i) = prim(IVX,ks+k-1,j,i);
        prim(IVY,ks-k,j,i) = prim(IVY,ks+k-1,j,i);
        prim(IVZ,ks-k,j,i) = -prim(IVZ,ks+k-1,j,i);
        // adiabatic projection
        prim(IT,ks-k,j,i)  = prim(IT,ks+k-1,j,i) + mubot*grav/cpbot*(2*k-1)*pco->dx3f(ks);
        prim(IPR,ks-k,j,i) = prim(IPR,ks+k-1,j,i)
          *pow(prim(IT,ks-k,j,i)/prim(IT,ks+k-1,j,i), cpbot/Globals::Rgas);
        for (int n = 1; n < NCOMP; ++n)
          prim(n,ks-k,j,i) = xbot[n];
      }
}

void ProjectPressureOuterX3(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  for (int k = 1; k <= NGHOST; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        prim(IVX,ke+k,j,i) = prim(IVX,ke-k+1,j,i);
        prim(IVY,ke+k,j,i) = prim(IVY,ke-k+1,j,i);
        prim(IVZ,ke+k,j,i) = -prim(IVZ,ke-k+1,j,i);
        // isothermal projection
        prim(IT,ke+k,j,i) = ttop;
        prim(IPR,ke+k,j,i) = prim(IPR,ke-k+1,j,i)
          *exp(-mutop*grav/(Globals::Rgas*ttop)*(2*k-1)*pco->dx3f(ke));
        for (int n = 1; n < NCOMP; ++n)
          prim(n,ke+k,j,i) = 0.;
      }
}

void ConstantCoolingRate(MeshBlock *pmb, const Real time, const Real dt, const int step,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Coordinates *pcoord = pmb->pcoord;

  if (pmb->pmy_mesh->mesh_size.nx3 == 1) { // 2D problem
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j) 
        if (pcoord->x2v(j) > -35.E3)
          for (int i = pmb->is; i <= pmb->ie; ++i) {
            cons(IEN,k,j,i) -= 2.5*prim(IPR,k,j,i)*crate/prim(IT,k,j,i)*dt;
          }
  } else {  // 3D problem
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      if (pcoord->x3v(k) > -35.E3)
        for (int j = pmb->js; j <= pmb->je; ++j) 
          for (int i = pmb->is; i <= pmb->ie; ++i) {
            cons(IEN,k,j,i) -= 2.5*prim(IPR,k,j,i)*crate/prim(IT,k,j,i)*dt;
          }
  }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  if (mesh_size.nx3 == 1) { // 2D problem
    EnrollUserBoundaryFunction(INNER_X2, ProjectPressureInnerX2);
    EnrollUserBoundaryFunction(OUTER_X2, ProjectPressureOuterX2);
  } else {  // 3D problem
    EnrollUserBoundaryFunction(INNER_X3, ProjectPressureInnerX3);
    EnrollUserBoundaryFunction(OUTER_X3, ProjectPressureOuterX3);
  }
  EnrollUserExplicitSourceFunction(ConstantCoolingRate);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  std::stringstream msg;
  // Step 1. define molecules
  pmol = new Molecule(pin);
  if (pmol->TotalNumber() != NCOMP) {
    msg << "### FATAL ERROR in ProblemGenerator: pmol->TotalNumber() != NMOL" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // Step 2. define reactions
  prg = new ReactionGroup(this, "condensation");
  prg->AddReaction(pin, "chemistry", "r1", pmol);
  prg->AddReaction(pin, "chemistry", "r2", pmol);

  // Step 3. define boundary conditions
  std::vector<std::string> xbot_str;
  Real ptop = pin->GetReal("problem", "ptop");
  Real pref = pin->GetReal("problem", "pref");
  Real tref = pin->GetReal("problem", "tref");
  SplitString(pin->GetString("problem", "xbot"), xbot_str);
  Real prim[NHYDRO];

  ttop = tref;
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

  // Step 4. construct 1D adiabatic T-P profile from top down
  // 4.1) T-P profile outside of this MeshBlock
  Real dz, ztop, zmbtop;
  int nx3 = pmy_mesh->mesh_size.nx3;
  grav = nx3 == 1 ? -phydro->psrc->GetG2() : -phydro->psrc->GetG3();

  if (nx3 == 1) { // 2D
    dz = pcoord->dx2f(je);
    ztop = pmy_mesh->mesh_size.x2max;
    zmbtop = pcoord->x2v(je);
  } else {  // 3D
    dz = pcoord->dx3f(ke);
    ztop = pmy_mesh->mesh_size.x3max;
    zmbtop = pcoord->x3v(ke);
  }
  Real mu = mutop, cp = peos->HeatCapacityP(prim);
  prim[IPR] = ptop;
  prim[IT ] = ttop;
  for (; ztop > zmbtop; ztop -= dz) {
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
  int start, end;
  if (nx3 == 1) { // 2D
    start = js;
    end = je;
  } else {
    start = ks;
    end = ke;
  }
  for (int m = end; m >= start; --m) {
    if (m == end) dz = ztop - zmbtop;
    else dz = nx3 == 1 ? pcoord->dx2f(m) : pcoord->dx3f(m);

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
    if (nx3 == 1) { // 2D
      for (int n = 0; n < NHYDRO; ++n)
        phydro->w(n,ks,m,is) = prim[n];
    } else {  // 3D
      for (int n = 0; n < NHYDRO; ++n)
        phydro->w(n,m,js,is) = prim[n];
    }
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
          if (nx3 == 1) { // 2D
            phydro->w(n,k,j,i) = phydro->w(n,ks,j,is);
          } else {  // 3D
            phydro->w(n,k,j,i) = phydro->w(n,k,js,is);
          }
          if (n == IVY)
            phydro->w(n,k,j,i) = 1.*(ran2(&iseed)-0.5)
              *(1.0+cos(kx*pcoord->x1v(i)))*(1.0+cos(ky*pcoord->x2v(j)))/4.0;
        }

  // Step 5. Cooling rate
  crate = pin->GetReal("problem", "crate");

  // Step 6. convert primitive variables to conserved variables
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
