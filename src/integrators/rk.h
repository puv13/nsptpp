/// \file
///
/// \brief Runge-Kutta  Integrator

#pragma once

#include "../integrator.h"
#include "expansion.h"
#include <sstream>


namespace IntegratorsNoExp {

//  ------------------------------------------------------------------------
/// Runge-Kutta Integrator of "Bali"-type ( PRD87.094517)
//  ------------------------------------------------------------------------
template <typename Latticetype, typename Forcefieldtype, typename Casimirtype,
          typename Forcetype, typename Randomtype, typename Generator>
class RKInt : public Integrator<Latticetype, Forcefieldtype, Forcetype,
                                Randomtype, Generator> {
 private:
  /// Member forcefield to avoid new forcefield allocation in
  /// every update step
  Forcefieldtype forcefield;
  Latticetype copy;
  double beta = 1.0;
  Casimirtype CA;


 public:
  void setBeta(double const &bt) { beta = bt; }
  double getBeta() const { return beta; }
  void setCA(Casimirtype const &c) { CA = c; }
  Casimirtype getCA() const { return CA; }

  ~RKInt() = default;
  void update(Forcetype const &force, Latticetype &lat, Generator &gen)
  {

    // using namespace Helpers;

    using FFscalar = decltype(forcefield[0]);
    using LTscalar = decltype(lat[0]);

    Randomtype random;

    // Constants (named following "Bali et al., PRD87.094517")
    // -------------------------------------------------------------
    double rt2 = sqrt(2.);
    double k1 = (-3. + 2. * rt2) * 0.5;
    double k2 = (-2. + rt2) * 0.5;
    double k5 = (5. - 3. * rt2) / 12.;

    const auto eps = this->getStepSize();


    // -------------------------------------------------------------
    // Tentative update;
    // -------------------------------------------------------------
    copy = Latticetype(lat);
    // Get current state of random number generator
    auto state = gen.get_state();

    // Pre-compute drift force terms
    // -------------------------------------------------------------
    forcefield = force(lat);
    std::transform(forcefield.begin(), forcefield.end(), forcefield.begin(),
                   [&](FFscalar &elem) {

                     std::complex<double> II(0.0, 1.0);
                     auto tmp = II * (eps * k1 * elem * beta +
                                      k2 * sqrt(eps) * random(gen));

                     return exp(tmp);

                   });

    // Update
    // -------------------------------------------------------------
    std::transform(forcefield.begin(), forcefield.end(), lat.begin(),
                   copy.begin(),
                   [](FFscalar const &f, LTscalar const &l) { return f * l; });


    // -------------------------------------------------------------
    // Final update;
    // -------------------------------------------------------------

    // Restore random number generator
    gen.set_state(state);

    forcefield = force(copy);
    std::transform(forcefield.begin(), forcefield.end(), forcefield.begin(),
                   [&](FFscalar &elem) {

                     std::complex<double> II(0.0, 1.0);
                     auto tmp = -II * (eps * elem * beta +
                                       k5 * beta * eps * eps * (CA * elem) +
                                       sqrt(eps) * random(gen));

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
