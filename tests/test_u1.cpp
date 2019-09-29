/// file
///
/// Tests related to the U1 gauge group implementation


#include "gaugegroups/u1.h"
#include "latticefields/cp.h"
#include "test_helpers.h"
#include <gtest/gtest.h>

TEST(U1, ConstructsOne)
{
  std::complex<double> init(2., 2.);
  U1 u(1);

  EXPECT_EQ(u.norm(), 1);
}

TEST(U1, ArithmeticOperators)
{
  std::array<std::complex<double>, 2> init = randomVector<2>();
  U1 u(init[0]);
  U1 v(init[1]);

  EXPECT_DOUBLE_EQ(std::real((v + u).value()), std::real(init[0] + init[1]));
  EXPECT_DOUBLE_EQ(std::imag((v + u).value()), std::imag(init[0] + init[1]));
  EXPECT_DOUBLE_EQ(std::real((u - v).value()), std::real(init[0] - init[1]));
  EXPECT_DOUBLE_EQ(std::imag((u - v).value()), std::imag(init[0] - init[1]));
  EXPECT_DOUBLE_EQ(std::real((v * u).value()), std::real(init[0] * init[1]));
  EXPECT_DOUBLE_EQ(std::imag((v * u).value()), std::imag(init[0] * init[1]));
  EXPECT_DOUBLE_EQ(std::real((u / v).value()), std::real(init[0] / init[1]));
  EXPECT_DOUBLE_EQ(std::imag((u / v).value()), std::imag(init[0] / init[1]));

  v += u;
  EXPECT_EQ(v.value(), init[0] + init[1]);
  u -= v;
  EXPECT_EQ(u.value(), -init[1]);
  v /= v;
  EXPECT_EQ(v.value(), std::complex<double>(1., 0.));
  v *= u;
  EXPECT_EQ(v.value(), u.value());

  // Arithmetic with complex numbers
  std::complex<double> cplx(0., 1.);
  v = u * cplx;
  v = v * cplx;
  EXPECT_EQ(v.value(), -1. * u.value());

  v = cplx * v;
  v = cplx * v;
  EXPECT_EQ(v.value(), u.value());
}

TEST(U1, MultiplicationWithCP)
{
  std::complex<double> init(2., 2.);
  U1 u(init, true);

  const CP<std::complex<double>, 4> cp1(std::complex<double>(1.0));

  CP<std::complex<double>, 4> multl = u * cp1;
  CP<std::complex<double>, 4> multr = cp1 * u;

  for (auto el : multl) {
    EXPECT_EQ(el, u.value());
  }

  for (auto el : multr) {
    EXPECT_EQ(el, u.value());
  }
}

TEST(U1, RandomIsGood)
{
  double res = 0.;
  size_t numtests = 100000lu;
  double pi = std::atan(1.) * 4.;


  for (size_t i = 0; i < numtests; ++i) {
    auto ru1 = randomU1();
    res += ru1.phase();
  }

  double mean = res / (static_cast<double>(numtests));

  // The mean of the phases better be 0 (within 3*sigma)
  EXPECT_NEAR(std::abs(mean), 0.,
              3 * std::sqrt(pi * pi / (3. * static_cast<double>(numtests))));
}
