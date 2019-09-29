/// \file
///
/// \brief Implementation of the force for QCD-like theories
/// (the gauge group is left as a free template parameter)

#pragma once

#include "action.h"
#include "expansion.h"

namespace qcd_internal {
/// Force implementation. Returns \f$\nabla S\f$.
template <typename group, typename algebra, std::size_t dimensions,
          std::size_t Nc>
class WilsonGaugeForce : public Force<LinkLattice<group, dimensions>,
                                      LinkLattice<algebra, dimensions>> {
 public:
  /// returns  \f$\nabla S = (U - U^+) - 1/N_c trace(U - U^+) \f$
  LinkLattice<algebra, dimensions>
  operator()(LinkLattice<group, dimensions> const &gauge) const
  {
    static_assert(dimensions > 0,
                  "less than 1 space time dimension is nonsense.");
    constexpr double oneOverNc = 1. / static_cast<double>(Nc);

    LinkLattice<algebra, dimensions> res(
      gauge.dimensionsArray()); // TODO: construct with correct size?
    auto resit = res.begin();
    group staple_up, staple_down;
    for (auto it = gauge.begin(); it != gauge.end(); ++it, ++resit) {
      for (std::size_t mu = 0; mu < dimensions; ++mu) {
        for (std::size_t nu = 0; nu < dimensions; ++nu) {
          // 	.---<---.
          // 	|       |
          // 	|  up   ^
          // 	v       |
          // 	X       o
          // 	^       |
          // 	| down  |
          // 	|       v
          // 	.---<---.

          const auto startingpt =
            it; // just for sanity check. may be removed later.
          if (mu == nu)
            continue;
          // TODO: replace by summing over nu > mu.
          staple_up = (*it)[nu]; // X upwards
          staple_up =
            staple_up * (*(it.advanceInDirection(
                          nu, 1)))[mu]; // go one up, get right-pointing link
          staple_down = (*(it.advanceInDirection(
            nu, -2)))[nu]; // go two down, get upwards link
          staple_down =
            staple_down * (*(it))[mu]; // stay here, get right-pointing link
          staple_down =
            staple_down *
            (*(it.advanceInDirection(mu, 1)))[nu]; // go right, get upwards
                                                   // link. staple_down complete
                                                   // up to dagger.
          staple_up =
            staple_up * (*(it.advanceInDirection(
                          nu, 1)))[nu]; // complete staple_up with upwards link

          // sanity check.
          assert(it.advanceInDirection(mu, -1) == startingpt);

          algebra tmp =
            static_cast<algebra>(staple_up) +
            static_cast<algebra>(staple_down); // tmp is dagger of U_P
          tmp = dagger(tmp) - tmp;
          decltype(tmp) retr(real(trace(tmp)) * oneOverNc);
          tmp -= retr;

          (*resit)[mu] = static_cast<algebra>(tmp);
          // subtract re(trace)..
        }
      }
    }
    return res;
  }
  ~WilsonGaugeForce(){};
  constexpr bool isGroupValued() const { return false; }
};
}
template <typename group, std::size_t dimensions, std::size_t orders>
using WilsonGaugeForceExpansion =
  qcd_internal::WilsonGaugeForce<Expansion<group, orders>,
                                 Expansion<typename group::algebra, orders>,
                                 dimensions, group::algebra::Nc>;

template <typename group, std::size_t dimensions>
using WilsonGaugeForce =
  qcd_internal::WilsonGaugeForce<group, typename group::algebra, dimensions,
                                 group::algebra::Nc>;
