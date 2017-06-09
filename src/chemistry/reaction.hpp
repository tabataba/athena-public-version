#ifndef REACTION_HPP_
#define REACTION_HPP_

// C++ header
#include <string>
#include <vector>
#include <iosfwd>

// Athena++ classes headers
#include "../athena.hpp"
#include "../eos/eos.hpp"

#define NREACTOR 4

class MeshBlock;
class Molecule;
class ParameterInput;

struct Reaction {
  std::string name, tag, comment;
  int reactor[NREACTOR];
  Real measure[NREACTOR];
  std::vector<Real> coeff;

  Reaction();
  ~Reaction();
  Reaction(Reaction const& other);
  Reaction& operator=(Reaction const& other);

  void SetFromString(std::string str, Molecule *pmol, std::string tag_ = "");
};

std::ostream& operator<<(std::ostream &os, Reaction const& rc);

Real NullReaction(Reaction const& rc, Real const prim[], Real time);

class ReactionGroup {
public:
  ReactionGroup(MeshBlock *pmb, std::string _name);
  ~ReactionGroup();

  // data
  MeshBlock* pmy_block;
  std::string name;
  ReactionGroup *prev, *next;

  // functions
  ReactionGroup* AddReactionGroup(std::string name);
  ReactionGroup* AddReaction(ParameterInput *pin, std::string block,
    std::string tag, Molecule *pmol, ReactionFunc_t pfunc_ = NullReaction);
  ReactionGroup* GetReactionGroup(std::string name);
  Reaction& GetReaction(std::string tag);
  Reaction const& GetReaction(std::string tag) const;
  void SetReactionFunction(std::string tag, ReactionFunc_t pfunc);
  void CalculateReactionRates(std::vector<Real>& rates, Real time,
    AthenaArray<Real> const& prim, int i, int j, int k);
  Real EvolveOneTimeStep(AthenaArray<Real>& prim, Real& time, Real dtmax, 
    int is, int ie, int js, int je, int ks, int ke);

protected:
  std::vector<ReactionFunc_t> fns_;
  std::vector<Reaction> rts_;
  std::vector<Real> nrate_;
};

Real GasGasSolidNH4SH(Reaction const& rc, Real const prim[], Real time);
Real GasCloudIdeal(Reaction const& rc, Real const prim[], Real time);
Real LiquidSolidIdeal(Reaction const& rc, Real const prim[], Real time);

#endif
