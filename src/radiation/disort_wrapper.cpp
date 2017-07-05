#ifdef USE_DISORT

// third party software DISORT headers
extern "C" {
  #include "cdisort.h"
}

// C++ headers
#include <stdexcept>

// Athena++ headers
#include "../parameter_input.hpp"
#include "../athena_arrays.hpp"
#include "radtran_wrapper.hpp"

DisortWrapper::DisortWrapper(ParameterInput *pin, int nlayer)
{
  ds.nlyr = nlayer;
  ds.nstr = pin->GetInteger("nstr");
  ds.nmom = pin->GetInteger("nmom");
  ds.nphase = pin->GetInteger("nphase");

  ds.ntau = pin->GetOrAddInteger("ntau", 0);
  ds.numu = pin->GetOrAddInteger("numu", 0);
  ds.nphi = pin->GetOrAddInteger("nphi", 0);
  ds.accur = pin->GetOrAddReal("accur", 0.);

  ds.bc.btemp = pin->GetReal("btemp");
  ds.bc.ttemp = pin->GetReal("ttemp");
  ds.bc.fbeam = pin->GetReal("fbeam");
  ds.bc.umu0  = pin->GetReal("umu0");
  ds.bc.phi0  = pin->GetReal("phi0");
  ds.bc.fisot = pin->GetReal("fisot");
  ds.bc.albedo = pin->GetReal("albedo");
  ds.bc.temis = pin->GetReal("temis");

  ds.flag.ibcnd = pin->GetOrAddInteger("ibcnd", 0);
  ds.flag.usrtau = pin->GetOrAddInteger("usrtau", 0);
  ds.flag.usrang = pin->GetOrAddInteger("usrang", 0);
  ds.flag.lamber = pin->GetOrAddInteger("lamber", 1);
  ds.flag.planck = pin->GetOrAddInteger("planck", 0);
  ds.flag.spher = pin->GetOrAddInteger("spher", 0);
  ds.flag.onlyfl = pin->GetOrAddInteger("onlyfl", 1);
  ds.flag.quiet = pin->GetOrAddInteger("quiet", 1);
  ds.flag.intensity_correction = pin->GetOrAddInteger("intensity_correction", 1);
  ds.flag.old_intensity_correction = pin->GetOrAddInteger("old_intensity_correction", 0);
  ds.flag.general_source = pin->GetOrAddInteger("general_source", 0);
  ds.flag.output_uum = pin->GetOrAddInteger("output_uum", 0);

  for (int i = 0; i < 5; ++i) ds.flag.prnt[i] = 0;
  
  c_disort_state_alloc(&ds);
  c_disort_out_alloc(&ds, &ds_out);
}

~DisortWrapper::DisortWrapper()
{
  c_disort_state_free(&ds);
  c_disort_out_free(&ds, &ds_out);
}

void DisortWrapper::Radiance(AthenaArray<Real>& rad, AthenaArray<Real> const& tau,
  AthenaArray<Real> const& temp, AthenaArray<Real> const& ssalb,
  AthenaArray<Real> const& pmom, Real const *wlo, Real const *whi, int nwave)
{
  std::stringstream msg;
  if (ds.flag.ibcnd != 0) {
    msg << "### FATAL ERROR in DisortWrapper::Radiance: ";
    msg << "ibcnd != 0" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // check dimension consistency
  int nlevel = ds.nlyr + 1;
  if (tau.GetDim1() != nlevel || tau.GetDim2() != nwave) {
    msg << "### FATAL ERROR in DisortWrapper::Radiance: ";
    msg << "tau(1) != nlevel or tau(2) != nwave)" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (temp.GetDim1() != nlevel || temp.GetDim2() != nwave) {
    msg << "### FATAL ERROR in DisortWrapper::Radiance: ";
    msg << "temp(1) != nlevel or temp(2) != nwave)" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (ssalb.GetDim1() != nlevel || ssalb.GetDim2() != nwave) {
    msg << "### FATAL ERROR in DisortWrapper::Radiance: ";
    msg << "ssalb(1) != nlevel or ssalb(2) != nwave)" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (pmom.GetDim1() != ds.nmom + 1 || pmom.GetDim2() != nlevel || pmom.GetDim3() != nwave) {
    msg << "### FATAL ERROR in DisortWrapper::Radiance: ";
    msg << "pmom(1) != ds.nmom + 1 or pmom(2) != nlevel or pmom(3) != nwave)" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }
  if (rad.GetDim1() != nwave || rad.GetDim2() != nlevel) {
    msg << "### FATAL ERROR in DisortWrapper::Radiance: ";
    msg << "rad(1) != nwave or rad(2) != nlevel)" << std::endl;
    throw std::runtime_error(msg.str().c_str());
  }

  // run disort
  for (int i = 0; i < nwave; ++i) {
    // boundary condition
    ds.bc = bc[i];

    // source function
    if (ds.flag.planck) {
      for (int j = 0; j < nlevel; ++j)
        ds.temper[j] = temp(i,j);
    }
    ds.wvnmlo = wlo[i];
    ds.wvnmhi = whi[i];

    // absorption
    for (int j = 0; j < nlevel - 1; ++j)
      ds.dtauc[j] = tau(i,j);

    // single scatering albedo
    for (int j = 0; j < nlevel; +j)
      ds.ssalb[j] = ssalb(i,j);

    // Legendre coefficients
    for (int j = 0; j < nlevel; ++j)
      for (int k = 0; k <= ds.nmom; ++k)
        ds.pmom[j*(ds.nmom + 1) + k] = pmom(i,j,k);

    c_disort(&ds, &ds_out);

    if (ds.flag.onlyfl) {
      for (int j = 0; j < nlevel; ++j)
        rad(j,i) = - ds_out.rad[j].rfldir - ds_out.rad[j].rfldn + ds_out.rad[j].flup;
    } else {
      for (int j = 0; j < nlevel; ++j)
        rad(j,i) = ds_out.uu[j];
    }
  }
}

#endif // USE_DISORT
