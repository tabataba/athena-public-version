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
#include "../utils/utils.hpp"
#include "../globals.hpp"
#include "../misc.hpp"

Real mu_top, mu_bot, cp_bot, t_top, grav, xbot[NCOMP];

void ProjectPressureInnerX2(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
     FaceField &b, Real time, Real dt, int is, int ie, int js, int je, int ks, int ke)
{
  for (int k = ks; k <= ke; ++k)
    for (int j = 1; j <= NGHOST; ++j)
      for (int i = is; i <= ie; ++i) {
        prim(IVX,k,js-j,i) = prim(IVX,k,js-j,i);
        prim(IVY,k,js-j,i) = -prim(IVY,k,js-j,i);
        //Real den = prim(IPR,k,js+j-1,i)*mu_bot/(Globals::Rgas*prim(IT,k,js+j-1,i));
        //prim(IPR,k,js-j,i) = prim(IPR,k,js+j-1,i) - den*grav*(2*j-1)*pco->dx2f(j);
        //prim(IT,k,js-j,i) = prim(IT,k,js+j-1,i)*pow(prim(IPR,k,js-j,i)/prim(IPR,k,js+j-1,i), rcp);
        //prim(IT,k,js-j,i) = prim(IPR,k,js-j,i)/prim(IPR,k,js+j-1,i)*prim(IT,k,js+j-1,i);
        // adiabatic projection
        prim(IT,k,js-j,i)  = prim(IT,k,js+j-1,i) - mu_bot*grav/cp_bot*(2*j-1)*pco->dx2f(j);
        prim(IPR,k,js-j,i) = prim(IPR,k,js+j-1,i)
          *pow(prim(IT,k,js-j,i)/prim(IT,k,js+j-1,i), cp_bot/Globals::Rgas);
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
        //Real den = prim(IPR,k,je-j+1,i)*mu_top/(Globals::Rgas*prim(IT,k,je-j+1,i));
        //prim(IPR,k,je+j,i) = prim(IPR,k,je-j+1,i) + den*grav*(2*j-1)*pco->dx2f(j);
        // isothermal projection
        prim(IT,k,je+j,i) = t_top;
        prim(IPR,k,je+j,i) = prim(IPR,k,je-j+1,i)
          *exp(mu_top*grav/(Globals::Rgas*t_top)*(2*j-1)*pco->dx2f(j));
        for (int n = 1; n < NCOMP; ++n)
          prim(n,k,je+j,i) = 0.;
      }
}

void JovianAtmMicrophysics(MeshBlock *pmb, const Real time, const Real dt, const int step,
  const AthenaArray<Real> &prim, const AthenaArray<Real> &bcc, AthenaArray<Real> &cons)
{
  Reaction const& r1 = pmb->prg->GetReaction("r1");
  Reaction const& r2 = pmb->prg->GetReaction("r2");
  EquationOfState *peos = pmb->peos;

  Real prim1[NCOMP], cons1[NCOMP];
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real xt0 = 1.;
        for (int n = 0; n < NCOMP; ++n) {
          prim1[n] = prim(n,k,j,i);
          if (n >= NGAS) xt0 -= prim1[n];
        }
        Real rate, u0 = peos->Energy(prim1);
        rate = GasCloudIdeal(r1, prim1, time);
        for (int n = 0; n < 2; ++n)
          prim1[r1.reactor[n]] += rate*r1.measure[n];
        rate = GasCloudIdeal(r2, prim1, time);
        for (int n = 0; n < 2; ++n)
          prim1[r2.reactor[n]] += rate*r2.measure[n];
        peos->EquilibrateUV(prim1, u0, xt0, cons1);
        for (int n = 0; n < NCOMP; ++n) cons(n,k,j,i) = cons1[n];
      }
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserBoundaryFunction(INNER_X2, ProjectPressureInnerX2);
  EnrollUserBoundaryFunction(OUTER_X2, ProjectPressureOuterX2);
  //EnrollUserExplicitSourceFunction(JovianAtmMicrophysics);
}

void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  std::stringstream msg;
  // molecules
  pmol = new Molecule(pin);
  if (Molecule::ntotal != NCOMP) {
    msg << "### FATAL ERROR in ProblemGenerator: Molecule::ntotal != NMOL" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  //Molecule *p = pmol;
  //while (p != NULL) {
  //  std::cout << *p << std::endl;
  //  p = p->next;
  //}

  // reactions
  prg = new ReactionGroup(this, "condensation");
  prg->AddReaction(pin, "chemistry", "r1", pmol);
  prg->AddReaction(pin, "chemistry", "r2", pmol);
  Reaction const& r1 = prg->GetReaction("r1");
  Reaction const& r2 = prg->GetReaction("r2");

  // pertubation wavelength
  Real kx = 20.*(PI)/(pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min);
  Real ky = 20.*(PI)/(pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min);
  long int iseed = -1;

  // boundary condition
  std::vector<std::string> xbot_str;
  SplitString(pin->GetString("problem", "xbot"), xbot_str);
  xbot[1] = atof(xbot_str[0].c_str());
  xbot[2] = atof(xbot_str[1].c_str());
  for (int i = NGAS; i < NCOMP; ++i)
    xbot[i] = 0.;

  // T-p profile
  Real ptop = pin->GetReal("problem", "ptop");
  Real pref = pin->GetReal("problem", "pref");
  Real tref = pin->GetReal("problem", "tref");
  grav = phydro->psrc->GetG2();
  t_top = tref;

  Real prim1[NCOMP];
  AthenaArray<Real>& prim = phydro->w;
  for (int k = ke; k >= ks; --k)
    for (int j = je; j >= js; --j)
      for (int i = ie; i >= is; --i) {
        for (int n = 0; n < NCOMP; ++n) prim1[n] = prim(n,k,j,i);
        Real dz = pcoord->x2v(j+1) - pcoord->x2v(j);
        Real cp = peos->HeatCapacityP(prim1);
        Real mu = peos->Mass(prim1);
        prim(1,k,j,i)   = xbot[1];
        prim(2,k,j,i)   = xbot[2];
        if (j == je) {
          prim(IT,k,j,i)  = tref;
          prim(IPR,k,j,i) = ptop;
          mu_top = mu;
        } else if (prim(IPR,k,j+1,i) < pref) { // isothermal atm
          prim(IT,k,j,i)  = tref;
          prim(IPR,k,j,i) = prim(IPR,k,j+1,i)*exp(-mu*grav*dz/(Globals::Rgas*tref));
        } else { // adiabatic atm
          prim(IT,k,j,i)  = prim(IT,k,j+1,i) - mu*grav/cp*dz;
          prim(IPR,k,j,i) = prim(IPR,k,j+1,i)*pow(prim(IT,k,j,i)/prim(IT,k,j+1,i), cp/Globals::Rgas);
        }
        if (j == js) {
          mu_bot = mu;
          cp_bot = cp;
        }
        prim(IVX,k,j,i) = prim(IVZ,k,j,i) = 0.;
        prim(IVY,k,j,i) = 0.05*ran2(&iseed)*(1.0+cos(kx*pcoord->x1v(i)))*(1.0+cos(ky*pcoord->x2v(j)))/4.0;

        // remove condensates
        prim1[1] = prim(1,k,j,i);
        prim1[2] = prim(2,k,j,i);
        prim1[IT] = prim(IT,k,j,i);
        prim1[IPR] = prim(IPR,k,j,i);
        Real rate = GasCloudIdeal(r1, prim1, 0.);
        for (int n = 0; n < 2; ++n)
          prim1[r1.reactor[n]] += rate*r1.measure[n];
        rate = GasCloudIdeal(r2, prim1, 0.);
        for (int n = 0; n < 2; ++n)
          prim1[r2.reactor[n]] += rate*r2.measure[n];
        for (int n = 1; n < NCOMP; ++n) prim(n,k,j,i) = prim1[n];
      }
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
