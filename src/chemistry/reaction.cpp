// C++ headers
#include <sstream>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <cstdlib>

// Athena++ headers
#include "molecule.hpp"
#include "reaction.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"

Reaction::Reaction()
{
  for (int i = 0; i < NREACTORS; ++i) {
    reactor[i] = -1;
    measure[i] = 0.;
  }
  for (int i = 0; i < NCOEFFS; ++i)
    coeff[i] = 0.;
  pfunc = NULL;
}

Reaction::~Reaction() {}

Reaction::Reaction(Reaction const& other)
{
  if (this == &other) return;
  *this = other;
}

Reaction& Reaction::operator=(Reaction const& other)
{
  name = other.name;
  comment = other.comment;
  for (int i = 0; i < NREACTORS; ++i) {
    reactor[i] = other.reactor[i];
    measure[i] = other.measure[i];
  }
  for (int i = 0; i < NCOEFFS; ++i)
    coeff[i] = other.coeff[i];
  pfunc = other.pfunc;
  return *this;
}

void Reaction::SetFromString(std::string str, Molecule *pmol)
{
  std::istringstream ss(str);
  std::string token, type;
  bool in_reagent = true;
  bool in_product = false;
  bool in_coeff = false;
  bool in_comment = false;
  int nreagent = 0;
  int ntotal = 0;
  int ncoeff = 0;
  std::string molecule[NREACTORS];

  while (ss.good()) {
    ss >> token;
    if (token == "--" || token == "->") {
      type = token;
      in_reagent = false;
      in_product = true;
      in_coeff = false;
      in_comment = false;
      ss >> token;
    } 
    if (token == "&") {
      in_reagent = false;
      in_product = false;
      in_coeff = true;
      in_comment = false;
      ss >> token;
    }
    if (token == "!") {
      in_reagent = false;
      in_product = false;
      in_coeff = false;
      in_comment = true;
      //ss >> token;
    }
    if (token == "+") continue;
    if (in_reagent) {
      size_t idigit = token.find_first_not_of("0123456789.");
      if (idigit == 0) {
        measure[ntotal] = -1.;
        reactor[ntotal] = pmol->GetMoleculeId(token);
        molecule[ntotal] = token;
      } else {
        std::string digit = token.substr(0, idigit);
        std::string symbol = token.substr(idigit);
        measure[ntotal] = -atof(digit.c_str());
        reactor[ntotal] = pmol->GetMoleculeId(symbol);
        molecule[ntotal] = symbol;
      }
      nreagent++;
      ntotal++;
    }
    if (in_product) {
      size_t idigit = token.find_first_not_of("0123456789.");
      if (idigit == 0) {
        measure[ntotal] = 1.;
        reactor[ntotal] = pmol->GetMoleculeId(token);
        molecule[ntotal] = token;
      } else {
        std::string digit = token.substr(0, idigit);
        std::string symbol = token.substr(idigit);
        measure[ntotal] = atof(digit.c_str());
        reactor[ntotal] = pmol->GetMoleculeId(symbol);
        molecule[ntotal] = symbol;
      }
      ntotal++;
    }
    if (in_coeff) {
      coeff[ncoeff] = atof(token.c_str());
      ncoeff++;
    }
    if (in_comment) {
      char line[256];
      ss.getline(line, 256);
      comment = line;
      int first = comment.find_first_not_of(" \t"),
          last = comment.find_last_not_of(" \t");
      comment = comment.substr(first, (last - first + 1));
    }
  }
  name = molecule[0] + " ";
  for (int i = 1; i < nreagent; i++)
    name += "+ " + molecule[i] + " ";
  name += type + " " + molecule[nreagent] + " ";
  for (int i = nreagent + 1; i < ntotal; i++)
    name += "+ " + molecule[i] + " ";
  name.erase(name.size() - 1);
}

std::ostream& operator<<(std::ostream &os, Reaction const& rt)
{
  os << rt.name;
  for (int i = 0; i < NCOEFFS; ++i)
    os << std::setw(12) << rt.coeff[i];
  os << " ! " << rt.comment;
  return os;
}

ReactionGroup::ReactionGroup(MeshBlock *pmb, std::string _name):
  pmy_block(pmb), name(_name)
{
  prev = NULL;
  next = NULL;
}

ReactionGroup::~ReactionGroup()
{
  if (prev != NULL) prev->next = next;
  if (next != NULL) next->prev = prev;
}

// functions
ReactionGroup* ReactionGroup::AddReactionGroup(MeshBlock *pmb, std::string name)
{
  std::stringstream msg;
  ReactionGroup *p = this;
  if (p == NULL) {
    msg << "### FATAL ERROR in AddReactionGroup: ReactionGroup is empty, use new ReactionGroup instead" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  while (p->next != NULL) p = p->next;
  p->next = new ReactionGroup(pmb, name);
  p->next->prev = p;
  p->next->next = NULL;

  return p->next;
}

std::vector<Reaction>& ReactionGroup::GetReactions(std::string name)
{
  std::stringstream msg;
  ReactionGroup *p = this;

  while ((p != NULL) && (p->name != name)) p = p->next;
  if (p == NULL) {
    msg << "### FATAL ERROR in GetReaction : ReactionGroup " << name << " not found" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  return p->r;
}

std::vector<Reaction> const& ReactionGroup::GetReactions(std::string name) const
{
  std::stringstream msg;
  ReactionGroup const *p = this;

  while ((p != NULL) && (p->name != name)) p = p->next;
  if (p == NULL) {
    msg << "### FATAL ERROR in GetReaction : ReactionGroup " << name << " not found" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  return p->r;
}
