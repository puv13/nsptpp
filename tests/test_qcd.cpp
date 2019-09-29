#include "actions/qcd.h"
#include "gaugegroups/su3.h"
#include "lattice.h"
#include <gtest/gtest.h>
TEST(WilsonGaugeForce, GivesCorrectForceInStrongCouplingLimit)
{
  WilsonGaugeForce<Su3, 3> s;
  // unit field of size {2,3,4}:
  LinkLattice<Su3, 3> testgauge({ 2, 3, 4 }, Su3(1));

  // force term is zero, because P - P^+ is zero for all hermitian matrices.
  const auto result = Su3Algebra(0);

  auto tmp = s(testgauge);
  for (auto const &dir : tmp) {
    for (auto const &a : dir) {
      EXPECT_TRUE(a == result);
    }
  }
}

TEST(WilsonGaugeForce, TreeLevelForceForExpansionTypes)
{
  WilsonGaugeForceExpansion<Su3, 3, 4> s;
  // unit field of size {2,3,4}:
  LinkLattice<Expansion<Su3, 4>, 3> testgauge({ 2, 3, 4 },
                                              Expansion<Su3, 4>(1));

  // force term is zero, because P - P^+ is zero for all hermitian matrices.
  const auto result = Expansion<Su3Algebra, 4>(0);

  auto tmp = s(testgauge);
  for (auto const &dir : tmp) {
    for (auto const &a : dir) {
      EXPECT_TRUE(a == result);
    }
  }
}
