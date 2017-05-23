#ifndef REACTION_HPP_
#define REACTION_HPP_

// C++ header
#include <string>
#include <vector>
#include <iosfwd>

// Athena++ classes headers
#include "../athena.hpp"

class MeshBlock;
class Molecule;

struct Reaction {
  int id;
  std::string name;
  std::string comment;
  int reactor[NCOMPONENT];
  Real measure[NCOMPONENT];
  ReactionFunc_t pfunc;
  std::vector<Real> coeff;

  Reaction();
  ~Reaction();
  Reaction(Reaction const& other);
  Reaction& operator=(Reaction const& other);

  void SetFromString(std::string str, Molecule *pmol);
};

std::ostream& operator<<(std::ostream &os, Reaction const& rt);

class ReactionGroup {
public:
  ReactionGroup(MeshBlock *pmb, std::string _name);
  ~ReactionGroup();

  // data
  MeshBlock* pmy_block;
  std::string name;
  ReactionGroup *prev, *next;
  std::vector<Reaction> r;

  // functions
  ReactionGroup* AddReactionGroup(MeshBlock *pmb, std::string name);
  std::vector<Reaction>& GetReactions(std::string name);
  std::vector<Reaction> const& GetReactions(std::string name) const;
  void SetAllReactionIds();
};

#endif
