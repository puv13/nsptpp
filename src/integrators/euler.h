/// \file
///
/// \brief simple euler integration

#pragma once

#include "expansion.h"
#include "integrator.h"

// namespace Integrators {
//   template<typename T, std::size_t ord>
//   class Euler : public Integrator<T, ord> {
//   public:
//     void update(Force const & force, Lattice<Expansion<T, ord>> & gauge)
//     const {
//       static_assert( ord >= 2, "must have at least 2 orders in expansion
//       ((noise enters at 1st order -- counting from 0))" );
//       Expansion<T, ord> res;
//       forcefield = force(gauge);
//       // F = eps * \Nabla S + sqrt(eps) * \eta
//       // "forcefield" holds, after this transform: exp(-F) or -F, depending
//       on whether force is group or algebra valued.
//       std::transform(forcefield.begin(), forcefield.end(),
//       forcefield.begin(),
// 		     [&stepsize, &force](T const & elem) {
// 		       auto tmp = - stepsize * elem;
// 		       // noise is of order 1:
// 		       tmp[1] -= std::sqrt(stepsize) * random();
// 		       if( not force.isGroupValued() ) {
// 			 tmp = exp(tmp);
// 		       }
// 		     });
//       // multiply "forcefield" with gauge:
//       std::transform(forcefield.begin(), forcefield.end(), gauge.begin(),
//       gauge.begin(),
// 		     [](T const & f, T const & u) {
// 		       return f*u;
// 		     });
//     }
//   private:
//     // is only a private member for performance reasons (avoids unnecessary
//     allocations)
//     Lattice<Expansion<T, ord> forcefield;
//   };
// }

namespace IntegratorsNoExp {

//  ------------------------------------------------------------------------
/// Euler Integrator
//  ------------------------------------------------------------------------
template <typename Latticetype, typename Forcefieldtype, typename Forcetype,
          typename Randomtype, typename Generator>
class EulerInt : public Integrator<Latticetype, Forcefieldtype, Forcetype,
                                   Randomtype, Generator> {
 private:
  /// Member forcefield to avoid new forcefield allocation in every update step
  Forcefieldtype forcefield;
  double beta = 1.0;

 public:
  void setBeta(double const &bt) { beta = bt; }
  double getBeta() const { return beta; }
  ~EulerInt() = default;
  // template<typename Generator>
  void update(Forcetype const &force, Latticetype &lat, Generator &gen)
  {

    using FFscalar = decltype(forcefield[0]);
    using LTscalar = decltype(lat[0]);

    Randomtype random;


    // Pre-compute drift force terms
    // -------------------------------------------------------------
    forcefield = force(lat);
    std::transform(
      forcefield.begin(), forcefield.end(), forcefield.begin(),
      [&](FFscalar &elem) {

        std::complex<double> II(0.0, 1.0);

        auto tmp = -1. * (this->getStepSize() * this->getBeta() * II * elem +
                          II * sqrt(this->getStepSize()) * random(gen));

        return exp(tmp);

      });
    // Update
    // -------------------------------------------------------------
    std::transform(forcefield.begin(), forcefield.end(), lat.begin(),
                   lat.begin(),
                   [](FFscalar const &f, LTscalar const &l) { return f * l; });
  }
};
}
