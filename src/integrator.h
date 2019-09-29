/// \file
///
/// Base class for integrators
#pragma once

#include "action.h"
#include "expansion.h"
#include "lattice.h"

// THIS DOES NOT COMPILE
//................................................................................
//................................................................................
//................................................................................
// namespace Integrators {
// template<typename T, std::size_t ord>
// class Integrator {
//   private:
//   double stepsize;

//   public:
//   virtual void update(Force const & force, Lattice<Expansion<T, ord>> &
//   gauge) const = 0;
//   virtual Integrator~() {};

//   void setStepSize(double stepsize_) { stepsize = stepsize_; }
//   inline double getStepSize() const { return stepsize; }
// };
// }
//................................................................................
//................................................................................
//................................................................................

namespace IntegratorsNoExp {

//  ------------------------------------------------------------------------
/// Base class for Integrators
//  ------------------------------------------------------------------------
template <typename Latticetype, typename Forcefieldtype, typename Forcetype,
          typename Randomtype, typename Generator>
class Integrator {
 private:
  double stepsize;

 public:
  // template<typename Generator>
  virtual void update(Forcetype const &force, Latticetype &lat,
                      Generator &gen) = 0;
  virtual ~Integrator(){};

  void setStepSize(double sz) { stepsize = sz; }
  double getStepSize() const { return stepsize; }
};
}
