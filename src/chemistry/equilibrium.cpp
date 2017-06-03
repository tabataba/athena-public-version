// C++ header
#include <stdexcept>
#include <sstream>
#include <cmath>

// Athena++ headers
#include "../math_funcs.hpp"
#include "reaction.hpp"

struct EntropySolver {
  EquationOfState *peos;
  ReactionGroup *prg;
  Real prim[];
  Real entropy;
  Real precision;
  Real operator()(Real temp) {
    prim[IT] = temp;
    Real time;
    //prg->EvolveToEquilibrium(prim, time, precision);
    return peos->Entropy(prim) - entropy;
  }
};

/*void EquilibriumSP(Real prim[], EquationOfState *peos, ReactionGroup *prg, 
  Real entropy, Real pres, Real precision)
{
  prim[IPR] = pres;
  EntropySolver solver;
  solver.peos = peos;
  solver.prg  = prg;
  solver.prim = prim;
  solver.entropy = entropy;
  solver.precision = 1.E-20;

  std::stringstream msg;
  int err = _root(prim[IT]/2.,prim[IT]*2, precision, &temp, solver);
  if (error) {
    msg << "### FATAL ERROR in EquilibriumSP: does not converge" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
}*/
