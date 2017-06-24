#ifndef ABSORBER_HPP_
#define ABSORBER_HPP_

// C++ header
#include <string>
#include <vector>
#include <iosfwd>

// Athena++ headers
#include "../athena.hpp"

/**@file
 * @brief This file contains declaration of Absorber
*
* **Author** : Cheng Li, California Institute of Technology <br>
* **Contact** : cli@gps.caltech.edu <br>
* **Revision history** :
* - June 21 2016, start documenting this file
* - July 28 2016, merge scatter into absorber
* - June 24 2017, adapt to Athena++ framework
*/

/** @brief AbsorpLookup stores a lookup table of absorption coefficients
*/

class Molecule;
template<typename T> class AthenaArray;

class Absorber {
  friend std::ostream& operator<<(std::ostream &os, Absorber const& ab);
public:
  // data
  std::string myname;
  Absorber *prev, *next;
  
  // functions
  Absorber(std::string name = "");
  Absorber(std::string name, Molecule *pmol, std::string mols = "");
  virtual ~Absorber();
  template<typename A> Absorber* AddAbsorber(A const& a);
  virtual void SaveCoefficient(std::string fname) const {}
  virtual void LoadCoefficient(std::string fname) {}
  virtual double AbsorptionCoefficient(double wave, Real const prim[]) const { return 0.; }
  virtual double SingleScateringAlbedo(double wave, Real const prim[]) const { return 0.; }
  virtual void Momentum(double wave, Real const prim[], double *pp, int np) const {}

protected:
  int id_[NCOMP];  /**< id of dependent molecules */
};

template<typename A>
inline Absorber* Absorber::AddAbsorber(A const& a) {
  A* pa = new A(a);
  Absorber *p = this;
  while (p->next != NULL) p = p->next;
  p->next = pa;
  p->next->prev = p;
  p->next->next = NULL;
  return p->next;
}

class HitranAbsorber: public Absorber {
  friend std::ostream& operator<<(std::ostream &os, HitranAbsorber const& ab);
  //friend HitranAbsorber const& MakeCKAbsorber<>(HitranAbsorber const& albl,
  //  int const *ck_axis, double const *ck_wave, int nbins);
public:
  HitranAbsorber() {}
  HitranAbsorber(std::string name, Molecule *pmol) : Absorber(name, pmol, name) {}
  virtual ~HitranAbsorber() {}
  void LoadCoefficient(std::string fname);
  double AbsorptionCoefficient(double wave, Real const prim[]) const;

protected:
  int len_[3];                  /**< length of interpolation axis */
  std::vector<double> axis_;    /**< interpolation axis */
  std::vector<double> kcoeff_;  /**< absorption coefficient */
  AthenaArray<Real>   refatm_;  /**< reference atmosphere */
  double RefTemp_(double pres) const;
};

class N2N2CIA: public Absorber {
public:
  N2N2CIA(Molecule *pmol) : Absorber("N2_N2", pmol, "N2") {}
  virtual ~N2N2CIA() {}
  void LoadCoefficient(std::string fname);
  double AbsorptionCoefficient(double wave, Real const prim[]) const;

protected:
  int rt_len_[2];
  int fd1_len_[2];
  int fd2_len_[2];
  std::vector<double> rt_axis_;
  std::vector<double> fd1_axis_;
  std::vector<double> fd2_axis_;
  std::vector<double> rt_;
  std::vector<double> fd1_;
  std::vector<double> fd2_;
};

class O2O2CIA: public Absorber {
public:
  O2O2CIA(Molecule *pmol) : Absorber("O2_O2", pmol, "O2") {}
  virtual ~O2O2CIA() {}
  void LoadCoefficient(std::string fname);
  double AbsorptionCoefficient(double wave, Real const prim[]) const;

protected:
  int fd_len_[2];
  int a1dg_x3sg00_len_[2];
  int a1dg_x3sg10_len_[2];
  int ab_len_[2];
  int other_len_[2];
  std::vector<double> fd_axis_;
  std::vector<double> a1dg_x3sg00_axis_;
  std::vector<double> a1dg_x3sg10_axis_;
  std::vector<double> ab_axis_;
  std::vector<double> other_axis_;
  std::vector<double> fd_;
  std::vector<double> a1dg_x3sg00_;
  std::vector<double> a1dg_x3sg10_;
  std::vector<double> ab_;
  std::vector<double> other_;
};

/** JunoCIA calculates and stores the CIA absorption coefficients */
class JunoCIA: public Absorber {
public:
  JunoCIA(Molecule *pmol) : Absorber("CIA", pmol, "H2 He CH4") {}
  double AbsorptionCoefficient(double wave, Real const prim[]) const;
};

/** JunoNH3 calculates and stores the NH3 gas absorption coefficients */
class JunoNH3Hanley: public Absorber {
public:
  JunoNH3Hanley(Molecule *pmol) : Absorber("NH3", pmol, "H2 He NH3 H2O") {}
  double AbsorptionCoefficient(double wave, Real const prim[]) const;
};

/** JunoNH3 calculates and stores the NH3 gas absorption coefficients */
class JunoNH3Bellotti: public Absorber {
public:
  JunoNH3Bellotti(Molecule *pmol) : Absorber("NH3", pmol, "H2 He NH3 H2O") {}
  double AbsorptionCoefficient(double wave, Real const prim[]) const;
};

/** JunoPH3 calculates and stores the PH3 gas absorption coefficients */
class JunoPH3Hoffman: public Absorber {
public:
  JunoPH3Hoffman(Molecule *pmol) : Absorber("PH3", pmol, "H2 He PH3") {}
  double AbsorptionCoefficient(double wave, Real const prim[]) const;
};

/** JunoH2S calculates and stores the H2S gas absorption coefficients */
class JunoH2S: public Absorber {
public:
  JunoH2S(Molecule *pmol) : Absorber("H2S", pmol, "H2 He H2S") {}
  double AbsorptionCoefficient(double wave, Real const prim[]) const;
};

/** JunoH2O calculates and stores the H2O gas absorption coefficients */
class JunoH2O: public Absorber {
public:
  JunoH2O(Molecule *pmol) : Absorber("H2O", pmol, "H2 He H2O") {}
  double AbsorptionCoefficient(double wave, Real const prim[]) const;
};

class JunoFreeFree: public Absorber {
public:
  JunoFreeFree() : Absorber("free-free") {}
  double AbsorptionCoefficient(double wave, Real const prim[]) const;
};

#endif
