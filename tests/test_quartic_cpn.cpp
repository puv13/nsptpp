#include "actions/quartic_cpn.h"
#include <array>
#include <complex>
#include <gtest/gtest.h>


TEST(Quartic_CPN, CanInit)
{
  const std::size_t dim = 2;
  const std::size_t NN = 1;
  std::array<size_t, dim> dimsar{ 2, 2 };
  CP<std::complex<double>, NN> cp_init(1.0);


  SiteLattice<CP<std::complex<double>, NN>, dim> cplat(dimsar, cp_init);
  QuarticCPNAction<std::complex<double>, 2, 2> x;
}


TEST(Quartic_CPN, CanGetCorrectResultsForAction)
{
  const std::size_t dim = 2;
  const std::size_t NN = 2;
  double beta = 0.123456;
  std::array<size_t, dim> dimsar{ 2, 2 };
  CP<std::complex<double>, NN> cp_init(1.0 / std::sqrt(NN));

  SiteLattice<CP<std::complex<double>, NN>, dim> cplat(dimsar, cp_init);

  QuarticCPNAction<std::complex<double>, NN, dim> cp_act(cplat, beta);

  EXPECT_DOUBLE_EQ(-1.0 * beta * dim * static_cast<double>(NN),
                   cp_act.atSite(0));
  EXPECT_DOUBLE_EQ(-1.0 * beta * dim * static_cast<double>(NN * cplat.volume()),
                   cp_act.total());
}


TEST(Quartic_CPN, CanGetCorrectResultsForEngergyDensity)
{
  const std::size_t dim = 4;
  const std::size_t NN = 1;
  double beta = 0.123456;
  double testbeta = beta;
  std::array<size_t, dim> dimsar{ 2, 3, 5, 11 };
  CP<std::complex<double>, NN> cp_init(1.0 / std::sqrt(NN));

  SiteLattice<CP<std::complex<double>, NN>, dim> cplat(dimsar, cp_init);
  QuarticCPNAction<std::complex<double>, NN, dim> cp_act(cplat, beta);

  EXPECT_DOUBLE_EQ(1.0 * dim * static_cast<double>(NN), cp_act.energyDensity());
  EXPECT_DOUBLE_EQ(-1.0 * testbeta * dim * static_cast<double>(NN),
                   cp_act.average());
}
