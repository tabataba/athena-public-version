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

  pmol = new Molecule(pin);
  if (Molecule::ntotal != NCOMP) {
    msg << "### FATAL ERROR in ProblemGenerator: Molecule::ntotal != NMOL" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  prg = new ReactionGroup(this, "condensation");
  prg->AddReaction(pin, "chemistry", "r1", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r2", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r3", pmol, &LiquidSolidIdeal);
  prg->AddReaction(pin, "chemistry", "r4", pmol, &GasGasSolidNH4SH);
  prg->AddReaction(pin, "chemistry", "r5", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r6", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r7", pmol, &LiquidSolidIdeal);

  Real mixr[NCOMP];

  // read solar abundance and enrichment factor
  Molecule *p = pmol;
  Real tmols = 0.;
  for (int i = 0; i < NCOMP; ++i) {
    Real solar = pin->GetOrAddReal("chemistry", p->name + ".solar", 0.);
    Real enrich = pin->GetOrAddReal("chemistry", p->name + ".enrich", 0.);
    mixr[i] = solar*enrich;
    tmols += mixr[i];
    p = p->next;
  }

  // set primitive variables
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

  Real norm, time = 0.;
  for (int i = is; i <= ie; ++i) {
    std::cout << "i = " << i << std::endl;
    int ncycle = 0;
    do {
      norm = prg->EvolveOneTimeStep(phydro->w, time, 1., i, i, js, je, ks, ke);
      ncycle++;
    } while (norm > 1.E-20 and ncycle < 2000);
    std::cout << "total cycles = " << ncycle << std::endl;
  }

  // convert to conservative variables
  // bcc is cell-centered magnetic fields, it is only a place holder here
  peos->PrimitiveToConserved(phydro->w, pfield->bcc, phydro->u, pcoord, is, ie, js, je, ks, ke);
}
