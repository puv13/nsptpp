#include "actions/canonical_cpn.h"
#include <array>
#include <complex>
#include <gtest/gtest.h>


TEST(Canonical_CPN, CanInit)
{
  const std::size_t dim = 2;
  const std::size_t NN = 1;
  std::array<size_t, dim> dimsar{ 2, 2 };
  CP<std::complex<double>, NN> cp_init(1.0);
  U1 u1_init(1.0);

  FullLattice<CP<std::complex<double>, NN>, U1, dim> cplat(dimsar, cp_init,
                                                           u1_init);
  CanonicalCPNAction<std::complex<double>, 2, 2> x;
}


TEST(Canonical_CPN, CanGetCorrectResultsForAction)
{
  const std::size_t dim = 2;
  const std::size_t NN = 2;
  double beta = 0.123456;
  std::array<size_t, dim> dimsar{ 2, 2 };
  CP<std::complex<double>, NN> cp_init(1.0 / std::sqrt(NN));
  U1 u1_init(std::complex<double>(1.0, 1.0), true);
  double u1_real = std::real(u1_init.value());

  FullLattice<CP<std::complex<double>, NN>, U1, dim> cplat(dimsar, cp_init,
                                                           u1_init);

  CanonicalCPNAction<std::complex<double>, NN, dim> cp_act(cplat, beta);

  EXPECT_DOUBLE_EQ(-2.0 * beta * dim * static_cast<double>(NN) * u1_real,
                   cp_act.atSite(0));
  EXPECT_DOUBLE_EQ(-2.0 * beta * dim *
                     static_cast<double>(NN * cplat.volume()) * u1_real,
                   cp_act.total());
}


TEST(Canonical_CPN, CanGetCorrectResultsForEngergyDensity)
{
  const std::size_t dim = 4;
  const std::size_t NN = 1;
  double beta = 0.123456;
  double testbeta = beta;
  std::array<size_t, dim> dimsar{ 2, 3, 5, 11 };
  CP<std::complex<double>, NN> cp_init(1.0 / std::sqrt(NN));
  U1 u1_init(1.0);

  FullLattice<CP<std::complex<double>, NN>, U1, dim> cplat(dimsar, cp_init,
                                                           u1_init);

  CanonicalCPNAction<std::complex<double>, NN, dim> cp_act(cplat, beta);

  EXPECT_DOUBLE_EQ(2.0 * dim * static_cast<double>(NN), cp_act.energyDensity());
  EXPECT_DOUBLE_EQ(-2.0 * testbeta * dim * static_cast<double>(NN),
                   cp_act.average());
}

TEST(Canonical_CPN, CanGetCorrectResultsForPlaquette)
{
  const std::size_t dim = 3;
  const std::size_t NN = 2;
  double beta = 0.123456;
  std::array<size_t, dim> dimsar{ 2, 2, 2 };
  CP<std::complex<double>, NN> cp_init(1.0 / std::sqrt(NN));
  U1 u1_init(1.0);

  FullLattice<CP<std::complex<double>, NN>, U1, dim> cplat(
    dimsar, cp_init, { u1_init, u1_init, u1_init * 2. });

  CanonicalCPNObservables<std::complex<double>, NN, dim> cp_obs(cplat, beta);

  EXPECT_DOUBLE_EQ(1.0, cp_obs.meanPlaquette(0lu, 1lu));
  EXPECT_DOUBLE_EQ(4.0, cp_obs.meanPlaquette(0lu, 2lu));
  EXPECT_DOUBLE_EQ(4.0, cp_obs.meanPlaquette(1lu, 2lu));
  EXPECT_DOUBLE_EQ((4. + 4. + 1) / 3., cp_obs.meanPlaquette());

  // EXPECT_DOUBLE_EQ(-2.0*beta*dim*static_cast<double>(N*cplat.volume()),cp_act.total());
}

TEST(Canonical_CPN, CanTestNormality)
{
  const std::size_t dim = 2;
  const std::size_t NN = 3;
  double beta = 0.123456;
  std::array<size_t, dim> dimsar{ 2, 3 };
  CP<std::complex<double>, NN> cp_init(1.0);
  U1 u1_init(1.0);

  FullLattice<CP<std::complex<double>, NN>, U1, dim> cplat(dimsar, cp_init,
                                                           u1_init);


  CanonicalCPNObservables<std::complex<double>, NN, dim> cp_obs(cplat, beta);

  EXPECT_EQ(false, cp_obs.check_normalisation());
}


TEST(Canonical_CPN, CanGetCorrectResultsForLandauGaugeFunctional)
{
  const std::size_t dim = 2;
  const std::size_t NN = 3;
  double beta = 0.123456;
  std::array<size_t, dim> dimsar{ 11, 13 };
  CP<std::complex<double>, NN> cp_init(1.0);
  U1 u1_init(1.0);

  FullLattice<CP<std::complex<double>, NN>, U1, dim> cplat(dimsar, cp_init,
                                                           u1_init);


  CanonicalCPNGaugeFix<std::complex<double>, NN, dim> cp_gf(cplat, beta);

  EXPECT_DOUBLE_EQ(static_cast<double>(cplat.volume() * dim),
                   cp_gf.LandauGaugeFunctional());

  // std::cout << cp_gf.LandauGaugeFunctional() << std::endl;
  for (auto i = 0; i < 20; ++i) {
    cp_gf.LandauGaugeSweep();
    // std::cout << cp_gf.LandauGaugeFunctional() << std::endl;
  }
}
