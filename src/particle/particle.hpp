#ifndef PARTICLE_HPP
#define PARTICLE_HPP

// C++ headers
#include <vector>
#include <string>

// Athena++ classes headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"

class MeshBlock;
class ParticleTableOutput;

struct Particle {
  Real time, x1, x2, x3;
  Real v1, v2, v3;
  Real rdata[NREAL_PARTICLE_DATA];
  int  idata[NINT_PARTICLE_DATA];
};

class ParticleGroup {
  friend class ParticleTableOutput;
public:
  ParticleGroup(MeshBlock *pmb, std::string _name);
  ~ParticleGroup();
  
  ParticleGroup *prev, *next;

  // functions
  ParticleGroup* AddParticleGroup(MeshBlock *pmb, std::string name);
  std::vector<Particle>& GetParticle(std::string name);
  std::vector<Particle> const& GetParticle(std::string name) const;
  void InterpolateToCellCenter(AthenaArray<Real> &u, Coordinates *pcoord) const;
  void InterpolateToCellFace(AthenaArray<Real> &u, Coordinates *pcoord) const;

protected:
  MeshBlock* pmy_block;
  std::vector<Particle> data;
  std::string name;
};

#endif
