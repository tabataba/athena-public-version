// ! \file thermodynamics.cpp
//   \brief thermodynamic model

// C++ headers
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../chemistry/molecule.hpp"
#include "../chemistry/reaction.hpp"
#include "../chemistry/condensation.hpp"
#include "../radiation/absorber.hpp"
#include "../coordinates/coordinates.hpp"
#include "../hydro/hydro.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"

Real LogPressureCoordinate(Real x, RegionSize rs)
{
  return pow(rs.x1min,1.-x)*pow(rs.x1max,x);
}

void Mesh::InitUserMeshData(ParameterInput *pin)
{
  EnrollUserMeshGenerator(X1DIR, LogPressureCoordinate);
}

// \brief Problem generator for shallow water model
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  std::stringstream msg;

  // Step 1. define molecules
  pmol = new Molecule(pin);
  if (pmol->TotalNumber() != NCOMP) {
    msg << "### FATAL ERROR in ProblemGenerator: pmol->TotalNumber() != NMOL" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  Molecule *p = pmol;
  Real tmols = 0.;
  Real mixr[NCOMP];
  for (int i = 0; i < NCOMP; ++i) {
    Real solar = pin->GetOrAddReal("chemistry", p->myname + ".solar", 0.);
    Real enrich = pin->GetOrAddReal("chemistry", p->myname + ".enrich", 0.);
    mixr[i] = solar*enrich;
    tmols += mixr[i];
    p = p->next;
  }

  // Step 2. define reactions
  prg = new ReactionGroup(this, "condensation");
  prg->AddReaction(pin, "chemistry", "r1", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r2", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r3", pmol, &LiquidSolidIdeal);
  prg->AddReaction(pin, "chemistry", "r4", pmol, &GasGasSolidNH4SH);
  prg->AddReaction(pin, "chemistry", "r5", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r6", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r7", pmol, &LiquidSolidIdeal);

  // Step 3. define absorbers
  std::string kfile = pin->GetString("radiation", "kcoeff_file");
  pabs = new HitranAbsorber("CH4", pmol);
  pabs->LoadCoefficient(kfile);
  pabs->AddAbsorber(HitranAbsorber("C2H2", pmol))
      ->LoadCoefficient(kfile);
  pabs->AddAbsorber(HitranAbsorber("C2H4", pmol))
      ->LoadCoefficient(kfile);
  pabs->AddAbsorber(HitranAbsorber("C2H6", pmol))
      ->LoadCoefficient(kfile);

  // Step 4. set primitive variables
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        phydro->w(IT,k,j,i) = 132.;
        for (int c = 1; c < NCOMP; ++c)
          phydro->w(c,k,j,i) = mixr[c]/tmols;
        phydro->w(IVX,k,j,i) = 0.;
        phydro->w(IVY,k,j,i) = 0.;
        phydro->w(IVZ,k,j,i) = 0.;
        phydro->w(IPR,k,j,i) = pcoord->x1v(i);
      }

  /*Real norm, time = 0.;
  for (int i = is; i <= ie; ++i) {
    std::cout << "i = " << i << std::endl;
    int ncycle = 0;
    do {
      norm = prg->EvolveOneTimeStep(phydro->w, time, 1., i, i, js, je, ks, ke);
      ncycle++;
    } while (norm > 1.E-20 and ncycle < 2000);
    std::cout << "total cycles = " << ncycle << std::endl;
  }

  // test primitive to conservative conversion
  for (int n = 1; n < NCOMP; ++n)
    phydro->w(n,ks,js,is) = n/100.;
  for (int n = 0; n < NHYDRO; ++n)
    std::cout << phydro->w(n,ks,js,is) << std::endl;
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, is, js, js, ks, ks);
  peos->ConservedToPrimitive(phydro->u, phydro->w, pfield->b, phydro->w, pfield->bcc, pcoord, is, is, js, js, ks, ks);
  for (int n = 0; n < NHYDRO; ++n)
    std::cout << phydro->w(n,ks,js,is) << std::endl;
  */

  // convert to conservative variables
  // bcc is cell-centered magnetic fields, it is only a place holder here
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
