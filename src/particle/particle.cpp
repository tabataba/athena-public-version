// C++ headers
#include <sstream>
#include <stdexcept>
#include <iostream>

// Athena++ headers
#include "particle.hpp"
//#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../mesh/mesh.hpp"
#include "../coordinates/coordinates.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
MPI_Datatype MPI_PARTICLE;
#endif

Particle::Particle() :
  time(0.), x1(0.), x2(0.), x3(0.),
  v1(0.), v2(0.), v3(0.)
{
#if NREAL_PARTICLE_DATA > 0
  for (int i = 0; i < NREAL_PARTICLE_DATA; ++i)
    rdata[i] = 0.;
#endif

#if NINT_PARTICLE_DATA > 0
  for (int i = 0; i < NINT_PARTICLE_DATA; ++i)
    idata[i] = 0;
#endif
}

Particle::~Particle() {}

Particle::Particle(Particle const& other)
{
  if (this == &other) return;
  *this = other;
}

Particle& Particle::operator=(Particle const& other)
{
  time = other.time;
  x1 = other.x1;
  x2 = other.x2;
  x3 = other.x3;
  v1 = other.v1;
  v2 = other.v2;
  v3 = other.v3;
#if NREAL_PARTICLE_DATA > 0
  for (int i = 0; i < NREAL_PARTICLE_DATA; ++i)
    rdata[i] = other.rdata[i];
#endif

#if NINT_PARTICLE_DATA > 0
  for (int i = 0; i < NINT_PARTICLE_DATA; ++i)
    idata[i] = other.idata[i];
#endif
  return *this;
}

std::ostream& operator<<(std::ostream &os, Particle const& pt)
{
  os << "time: "<< pt.time << std::endl
     << "x1: " << pt.x1 << " v1: " << pt.v1 << std::endl
     << "x2: " << pt.x2 << " v2: " << pt.v2 << std::endl
     << "x3: " << pt.x3 << " v3: " << pt.v3 << std::endl;
#if NREAL_PARTICLE_DATA > 0
  os << "rdata: ";
  for (int i = 0; i < NREAL_PARTICLE_DATA; ++i)
    os << pt.rdata[i] << " ";
  os << std::endl;
#endif

#if NINT_PARTICLE_DATA > 0
  os << "idata: ";
  for (int i = 0; i < NREAL_PARTICLE_DATA; ++i)
    os << pt.idata[i] << " ";
#endif
  return os;
}

int ParticleGroup::ntotal = 0;

// constructor, initializes data structure and parameters
ParticleGroup::ParticleGroup(MeshBlock *pmb, std::string _name):
  pmy_block(pmb), name(_name)
{
  prev = NULL;
  next = NULL;
  particle_fn_ = pmb->pmy_mesh->particle_fn_;

  int ncells1 = pmb->block_size.nx1 + 2*(NGHOST);
  int ncells2 = 1, ncells3 = 1;
  if (pmb->block_size.nx2 > 1) ncells2 = pmb->block_size.nx2 + 2*(NGHOST);
  if (pmb->block_size.nx3 > 1) ncells3 = pmb->block_size.nx3 + 2*(NGHOST);

  coordinates_ = new Real [ncells3 + ncells2 + ncells1];

  for (int k = 0; k < ncells3; ++k)
    coordinates_[k] = pmb->pcoord->x3v(k);
  for (int j = 0; j < ncells2; ++j)
    coordinates_[ncells3 + j] = pmb->pcoord->x2v(j);
  for (int i = 0; i < ncells1; ++i)
    coordinates_[ncells3 + ncells2 + i] = pmb->pcoord->x1v(i);

  lengths_[0] = ncells3;
  lengths_[1] = ncells2;
  lengths_[2] = ncells1;

#ifdef MPI_PARALLEL
  if (ntotal == 0) {
    int counts[2] = {7 + NREAL_PARTICLE_DATA, NINT_PARTICLE_DATA};
    MPI_Datatype types[2] = {MPI_ATHENA_REAL, MPI_INT};
    MPI_Aint disps[2];

    Particle p;

    MPI_Address(&p.time, disps);
    #if NINT_PARTICLE_DATA > 0
      MPI_Address(&p.idata, disps + 1);
    #endif

    disps[1] -= disps[0];
    disps[0] = 0;

    if (NINT_PARTICLE_DATA > 0)
      MPI_Type_struct(2, counts, disps, types, &MPI_PARTICLE);
    else
      MPI_Type_contiguous(counts[0], types[0], &MPI_PARTICLE);
    MPI_Type_commit(&MPI_PARTICLE);
  }
#endif

  ntotal++;
}

// destructor
ParticleGroup::~ParticleGroup()
{
  delete[] coordinates_;
  if (prev != NULL) prev->next = next;
  if (next != NULL) next->prev = prev;
  ntotal--;

#ifdef MPI_PARALLEL
  if (ntotal == 0)
    MPI_Type_free(&MPI_PARTICLE);
#endif
}

// functions
ParticleGroup* ParticleGroup::AddParticleGroup(std::string name)
{
  std::stringstream msg;
  ParticleGroup *p = this;
  if (p == NULL) {
    msg << "### FATAL ERROR in AddParticleGroup: ParticleGroup is empty, use new ParticleGroup instead" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  while (p->next != NULL) p = p->next;
  p->next = new ParticleGroup(pmy_block, name);
  p->next->prev = p;
  p->next->next = NULL;

  return p->next;
}

std::vector<Particle>& ParticleGroup::GetParticle(std::string name)
{
  std::stringstream msg;
  ParticleGroup *p = this;

  while ((p != NULL) && (p->name != name)) p = p->next;
  if (p == NULL) {
    msg << "### FATAL ERROR in GetParticle : ParticleGroup " << name << " not found" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  return p->q;
}

std::vector<Particle> const& ParticleGroup::GetParticle(std::string name) const
{
  std::stringstream msg;
  ParticleGroup const *p = this;

  while ((p != NULL) && (p->name != name)) p = p->next;
  if (p == NULL) {
    msg << "### FATAL ERROR in GetParticle : ParticleGroup " << name << " not found" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  return p->q;
}
