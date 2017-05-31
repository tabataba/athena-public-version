// ! \file thermodynamics.cpp
//   \brief thermodynamic model

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../chemistry/molecule.hpp"
#include "../chemistry/reaction.hpp"
#include "../chemistry/condensation.hpp"

// \brief Problem generator for shallow water model
void MeshBlock::ProblemGenerator(ParameterInput *pin)
{
  pmol = new Molecule(pin);
  prg = new ReactionGroup(this, "condensation");

  prg->AddReaction(pin, "chemistry", "r1", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r2", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r3", pmol, &LiquidSolidIdeal);
  prg->AddReaction(pin, "chemistry", "r4", pmol, &GasGasSolidNH4SH);
  prg->AddReaction(pin, "chemistry", "r5", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r6", pmol, &GasCloudIdeal);
  prg->AddReaction(pin, "chemistry", "r7", pmol, &LiquidSolidIdeal);

  prg->SetReactionRateArray();
}
