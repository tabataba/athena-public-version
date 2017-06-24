// C/C++ headers
#include <iostream>
#include <algorithm>

// Athena++ headers
#include "../chemistry/molecule.hpp"
#include "../math_funcs.hpp"
#include "../athena_arrays.hpp"
#include "../globals.hpp"
#include "absorber.hpp"

// External library headers
#ifdef NETCDFOUTPUT
  #include <netcdf.h>
#endif

std::ostream& operator<<(std::ostream &os, HitranAbsorber const& ab)
{
  os << (Absorber)(ab);
  for (int i = 0; i < 3; ++i) {
    os << "Axis " << i << " = " << ab.len_[i] << std::endl;
    if (i == 0) {
      os << "Minimum value = " << *std::min_element(ab.axis_.begin(), ab.axis_.begin() + ab.len_[i]) << std::endl;
      os << "Maximum value = " << *std::max_element(ab.axis_.begin(), ab.axis_.begin() + ab.len_[i]) << std::endl;
    } else {
      os << "Minimum value = " << *std::min_element(ab.axis_.begin() + ab.len_[i - 1], 
              ab.axis_.begin() + ab.len_[i]) << std::endl;
      os << "Maximum value = " << *std::max_element(ab.axis_.begin() + ab.len_[i - 1], 
              ab.axis_.begin() + ab.len_[i]) << std::endl;
    }
  }
  os << "No. of kcoeff = " << ab.kcoeff_.size() << std::endl;
  os << "Minimum value = " << *std::min_element(ab.kcoeff_.begin(), ab.kcoeff_.end()) << std::endl;
  os << "Maximum value = " << *std::max_element(ab.kcoeff_.begin(), ab.kcoeff_.end()) << std::endl;
  return os;
}

double HitranAbsorber::RefTemp_(double pres) const {
  int nlevel = refatm_.GetDim1();
  int jl = -1, ju = nlevel, jm;
  // make sure pressure is in ascending order
  while (ju - jl > 1) {
    jm = (ju + jl) >> 1;
    if (pres < refatm_(IPR,jm)) ju = jm;
    else jl = jm;
  }

  // prevent interpolation problem at boundary
  if (jl == -1) jl = 0;
  if (ju == nlevel) ju = nlevel - 1;
  if (jl == ju) return refatm_(IT,jl);

  //assert(jl >= 0 && ju < nlevel);
  double result = log(refatm_(IPR,jl)/pres)*log(refatm_(IT,ju)) 
    + log(pres/refatm_(IPR,ju))*log(refatm_(IT,jl));
  result = exp(result/log(refatm_(IPR,jl)/refatm_(IPR,ju)));
  return result;
}

void HitranAbsorber::LoadCoefficient(std::string fname)
{
#ifdef NETCDFOUTPUT
  int fileid, dimid, varid;
  char tname[80];
  nc_open(fname.c_str(), NC_NETCDF4, &fileid);

  nc_inq_dimid(fileid, "Wavenumber", &dimid);
  nc_inq_dimlen(fileid, dimid, len_);
  nc_inq_dimid(fileid, "Pressure", &dimid);
  nc_inq_dimlen(fileid, dimid, len_ + 1);
  strcpy(tname, "T_");
  strcat(tname, myname.c_str());
  nc_inq_dimid(fileid, tname, &dimid);
  nc_inq_dimlen(fileid, dimid, len_ + 2);

  axis_.resize(len_[0] + len_[1] + len_[2]);

  nc_inq_varid(fileid, "Wavenumber", &varid);
  nc_get_var_double(fileid, varid, axis_.data());
  nc_inq_varid(fileid, "Pressure", &varid);
  nc_get_var_double(fileid, varid, axis_.data() + len_[0]);
  nc_inq_varid(fileid, tname, &varid);
  nc_get_var_double(fileid, varid, axis_.data() + len_[0] + len_[1]);

  double *temp = new double[len_[1]];
  nc_inq_varid(fileid, "Temperature", &varid);
  nc_get_var_double(fileid, varid, temp);

  refatm_.NewAthenaArray(NCOMP, len_[1]);
  for (int i = 0; i < len_[1]; i++) {
    refatm_(IPR,i) = axis_[len_[0] + i];
    refatm_(IT,i)  = temp[i];
  }

  kcoeff_.resize(len_[0]*len_[1]*len_[2]);
  nc_inq_varid(fileid, myname.c_str(), &varid);
  nc_get_var_double(fileid, varid, kcoeff_.data());
  nc_close(fileid);
  delete[] temp;
#endif
}

double HitranAbsorber::AbsorptionCoefficient(double wave, Real const prim[]) const
{
  // first axis is wavenumber, second is pressure, third is temperature anomaly
  double val, coor[3] = {wave, prim[IPR], prim[IT] - RefTemp_(prim[IPR])};
  _interpn(&val, coor, kcoeff_.data(), axis_.data(), len_, 3);
  double dens = prim[IPR]/(Globals::Rgas*prim[IT]);
  return 1.E-3*exp(val)*dens*prim[id_[0]]; // ln(m*2/kmol) -> 1/m
}
