#ifndef DISORT_WRAPPER_HPP
#define DISORT_WRAPPER_HPP

// C/C++ header
#include <vector>

// Athena++ header
#include "../athena.hpp"

struct disort_state;
struct disort_output;
struct disort_bc;
class ParameterInput;
template<typename T> class AthenaArray;

class DisortWrapper {
public:
  DisortWrapper(ParameterInput *pin, int nlayer);
  ~DisortWrapper();

  // data
  std::vector<disort_bc> bc;

  // functions
  void Radiance(AthenaArray<Real>& rad, AthenaArray<Real> const& tau,
    AthenaArray<Real> const& temp, AthenaArray<Real> const& ssalb,
    AthenaArray<Real> const& pmom, Real const *wlo, Real const *whi, int nwave);

private:
  disort_state  ds;
  disort_output ds_out;
};

#endif
