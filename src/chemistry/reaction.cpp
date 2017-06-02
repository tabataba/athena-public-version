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
#include "../parameter_input.hpp"

Reaction::Reaction()
{
  for (int i = 0; i < NREACTOR; ++i) {
    reactor[i] = -1;
    measure[i] = 0.;
  }
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
  tag = other.tag;
  comment = other.comment;
  for (int i = 0; i < NREACTOR; ++i) {
    reactor[i] = other.reactor[i];
    measure[i] = other.measure[i];
  }
  coeff = other.coeff;
  return *this;
}

void Reaction::SetFromString(std::string str, Molecule *pmol, std::string tag_)
{
  std::stringstream msg;
  std::istringstream ss(str);
  std::string token, type;
  bool in_reagent = true;
  bool in_product = false;
  bool in_coeff = false;
  bool in_comment = false;
  int nreagent = 0;
  int ntotal = 0;
  std::string molecule[NREACTOR];

  // clear coefficients
  coeff.clear();

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
      if (ntotal >= NREACTOR) {
        msg << "### FATAL ERROR in Reaction:SetFromString, number of reagents exceeds NREACTOR" << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
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
      if (ntotal >= NREACTOR) {
        msg << "### FATAL ERROR in Reaction:SetFromString, number of reagents exceeds NREACTOR" << std::endl;
        throw std::runtime_error(msg.str().c_str());
      }
    }
    if (in_coeff)
      coeff.push_back(atof(token.c_str()));
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
  tag = tag_;

  // reset quantities other than ntotal
  for (int i = ntotal; i < NREACTOR; ++i) {
    reactor[i] = -1;
    measure[i] = 0.;
  }
}


std::ostream& operator<<(std::ostream &os, Reaction const& rc)
{
  os << rc.name;
  for (size_t i = 0; i < rc.coeff.size(); ++i)
    os << std::setw(12) << rc.coeff[i];
  if (rc.comment != "")
    os << " ! " << rc.comment;
  if (rc.tag != "")
    os << ", " << rc.tag;
  return os;
}

ReactionGroup::ReactionGroup(MeshBlock *pmb, std::string _name):
  pmy_block(pmb), name(_name)
{
  prev = NULL;
  next = NULL;
  initialized_ = false;

  // this is dummpy allocation of memory for reaction rate array
  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);
  rate.NewAthenaArray(ncells3,ncells2,ncells1);
}

ReactionGroup::~ReactionGroup()
{
  rate.DeleteAthenaArray();
  if (prev != NULL) prev->next = next;
  if (next != NULL) next->prev = prev;
}

// functions
ReactionGroup* ReactionGroup::AddReactionGroup(std::string name)
{
  std::stringstream msg;
  ReactionGroup *p = this;
  if (p == NULL) {
    msg << "### FATAL ERROR in AddReactionGroup: ReactionGroup is empty, use new ReactionGroup instead" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  while (p->next != NULL) p = p->next;
  p->next = new ReactionGroup(pmy_block, name);
  p->next->prev = p;
  p->next->next = NULL;

  return p->next;
}

ReactionGroup* ReactionGroup::AddReaction(ParameterInput *pin, std::string block,
  std::string tag, Molecule *pmol, ReactionFunc_t pfunc)
{
  std::stringstream msg;
  if (initialized_) {
    msg << "### FATAL ERROR in AddReaction: ReactionGroup has been initialized. No more reaction is allowd." << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  std::string rname = pin->GetString(block, tag);
  Reaction rc;
  rc.SetFromString(rname, pmol, tag);
  fns_.push_back(pfunc);
  rts_.push_back(rc);
}

ReactionGroup* ReactionGroup::GetReactionGroup(std::string name)
{
  std::stringstream msg;
  ReactionGroup *p = this;

  while ((p != NULL) && (p->name != name)) p = p->next;
  if (p == NULL) {
    msg << "### FATAL ERROR in GetReactionGroup : ReactionGroup " << name << " not found" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  return p;
}

int ReactionGroup::FindReactionId(std::string tag) const
{
  for (size_t r = 0; r < rts_.size(); ++r)
    if (rts_[r].tag == tag)
      return r;
  std::stringstream msg;
  msg << "### FATAL ERROR in FindReactionId: " << name << " not found" << std::endl;
  throw std::runtime_error(msg.str().c_str());
}

void ReactionGroup::SetReactionFunction(int i, ReactionFunc_t pfunc)
{
  fns_[i] = pfunc;
}

void ReactionGroup::InitReactionRateArray()
{
  rate.DeleteAthenaArray();
  // Allocate memory for reaction rate array
  int ncells1 = pmy_block->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmy_block->block_size.nx2 > 1) ncells2 = pmy_block->block_size.nx2 + 2*(NGHOST);
  if (pmy_block->block_size.nx3 > 1) ncells3 = pmy_block->block_size.nx3 + 2*(NGHOST);
  rate.NewAthenaArray(rts_.size(),ncells3,ncells2,ncells1);
  initialized_ = true;
}

void ReactionGroup::CalculateReactionRates(AthenaArray<Real> const& prim, int i1, int
i2, Real time)
{
  for (size_t r = 0; r < fns_.size(); ++r)
    if (fns_[r] != NULL)
      fns_[r](pmy_block, time, rts_[r], prim, i1, i2, rate, r);
}
