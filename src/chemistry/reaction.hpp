#ifndef REACTION_HPP_
#define REACTION_HPP_

// C++ header
#include <string>
#include <vector>
#include <iosfwd>

// Athena++ classes headers
#include "../athena.hpp"

#define NREACTOR 4

class MeshBlock;
class Molecule;
class ParameterInput;

struct Reaction {
  std::string name, tag, comment;
  int reactor[NREACTOR];
  Real measure[NREACTOR];
  ReactionFunc_t pfunc;
  std::vector<Real> coeff;

  Reaction();
  ~Reaction();
  Reaction(Reaction const& other);
  Reaction& operator=(Reaction const& other);

  void SetFromString(std::string str, Molecule *pmol, std::string tag_ = "");
};

std::ostream& operator<<(std::ostream &os, Reaction const& rc);

class ReactionGroup {
public:
  ReactionGroup(MeshBlock *pmb, std::string _name);
  ~ReactionGroup();

  // data
  MeshBlock* pmy_block;
  std::string name;
  ReactionGroup *prev, *next;
  AthenaArray<Real> rate;

  // functions
  ReactionGroup* AddReactionGroup(MeshBlock *pmb, std::string name);
  ReactionGroup* AddReaction(ParameterInput *pin, std::string block,
    std::string tag_, Molecule *pmol, ReactionFunc_t pfunc_);
  std::vector<Reaction>& GetReactions(std::string name);
  std::vector<Reaction> const& GetReactions(std::string name) const;
  void SetReactionRateArray();

protected:
  std::vector<Reaction> q;
};

#endif
