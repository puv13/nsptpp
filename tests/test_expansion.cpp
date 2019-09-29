#include "expansion.h"
#include "gaugegroups/su3.h"
#include "latticefields/cp.h"
#include "test_helpers.h"
#include <cmath>
#include <complex>
#include <gtest/gtest.h>

#define NORD 5ul
TEST(Expansion, GivesNumberOfOrders)
{
  Expansion<int, NORD> a;
  EXPECT_EQ(a.numberOfOrders(),
            static_cast<decltype(a.numberOfOrders())>(NORD));
}
TEST(Expansion, CanCheckIfTypeIsAnExpansion)
{
  auto isExpansion = is_expansion<Expansion<int, 1>>::value;
  auto isNotExpansion = is_expansion<int>::value;
  EXPECT_TRUE(isExpansion);
  EXPECT_TRUE(not isNotExpansion);
}
TEST(Expansion, CanCheckIfTwoTypesCanBeAdded)
{
  using namespace CHECK;
  bool res = canBeAdded<int, Su3>::value;
  EXPECT_FALSE(res);
  res = canBeAdded<Su3, Su3>::value;
  EXPECT_FALSE(res);
  res = canBeAdded<int, int>::value;
  EXPECT_TRUE(res);
  res = canBeAdded<Su3Algebra, Su3Algebra>::value;
  EXPECT_TRUE(res);
}
TEST(Expansion, ConstexprInitializationOfASingleElementArray)
{
  constexpr Expansion<double, 1> a;
  EXPECT_EQ(a.numberOfOrders(), 1ul);
  for (auto const &e : a.data) {
    EXPECT_DOUBLE_EQ(e, 0.0);
  }
}
TEST(Expansion, ConstexprInitializationToZero)
{
  constexpr Expansion<double, NORD> a;
  EXPECT_EQ(a.numberOfOrders(), NORD);
  for (auto const &e : a.data) {
    EXPECT_DOUBLE_EQ(e, 0.0);
  }
}
TEST(Expansion, ConstexprInitializationToOne)
{
  constexpr Expansion<double, NORD> a{
    Expansion<double, NORD>::one<double, NORD>
  };
  EXPECT_EQ(a.numberOfOrders(),
            static_cast<decltype(a.numberOfOrders())>(NORD));
  EXPECT_DOUBLE_EQ(a[0], 1.0) << " at order 0";
  int i = 1;
  for (auto it = std::next(std::begin(a.data)); it != std::end(a.data); ++it)
    EXPECT_DOUBLE_EQ(*it, 0.0) << " at order " << i++;
}
TEST(Expansion, ConstexprInitializationToSomeValue)
{
  constexpr Expansion<double, NORD> a(5.67);
  EXPECT_EQ(a.numberOfOrders(), NORD);
  EXPECT_DOUBLE_EQ(a[0], 5.67) << " at order 0";
  int i = 1;
  for (auto it = std::next(std::begin(a.data)); it != std::end(a.data); ++it)
    EXPECT_DOUBLE_EQ(*it, 0.0) << " at order " << i++;
}
TEST(Expansion, ComparesEqual)
{
  Expansion<int, NORD> a, b;
  for (std::size_t i = 0; i < NORD; ++i) {
    a[i] = static_cast<int>(i);
    b[i] = static_cast<int>(i);
  }
  EXPECT_TRUE(a == b);
}
TEST(Expansion, ComparesNonEqual)
{
  Expansion<int, NORD> a, b;
  for (std::size_t i = 0; i < NORD; ++i) {
    a[i] = static_cast<int>(i);
    b[i] = static_cast<int>(i + 1);
  }
  EXPECT_TRUE(a != b);
}
TEST(Expansion, SumOfTwoExpansionsIsCorrect)
{
  Expansion<std::size_t, NORD> a, b, c;
  for (std::size_t i = 0; i < NORD; ++i) {
    a[i] = i;
    b[i] = i + 1;
    c[i] = 2 * i + 1;
  }
  EXPECT_TRUE(a + b == c);
}
TEST(Expansion, CanConvertElementTypes)
{
  Expansion<double, NORD> a, c;
  Expansion<int, NORD> b;
  for (std::size_t i = 0; i < NORD; ++i) {
    a[i] = static_cast<double>(i);
    b[i] = static_cast<int>(i);
  }
  c = expansion_cast<double>(a);
  EXPECT_TRUE(a == c);
}
TEST(Expansion, ScalarMultiplicationFromLeft)
{
  Expansion<double, NORD> a, b;
  for (std::size_t i = 0; i < NORD; ++i) {
    a[i] = 0.5 * static_cast<double>(i);
    b[i] = static_cast<double>(i);
  }
  EXPECT_TRUE(2.0 * a == b);
}
TEST(Expansion, ScalarMultiplicationFromRight)
{
  Expansion<double, NORD> a, b;
  for (std::size_t i = 0; i < NORD; ++i) {
    a[i] = 0.5 * static_cast<double>(i);
    b[i] = static_cast<double>(i);
  }
  EXPECT_TRUE(a * 2.0 == b);
}
TEST(Expansion, OrderByOrderMultiplication)
{
  Expansion<int, 3> a, b, c;
  a[0] = 28;
  b[0] = 91;
  a[1] = 17;
  b[1] = 37;
  a[2] = 72;
  b[2] = 12;
  a[3] = 89;
  b[3] = 92;

  c[0] = a[0] * b[0];
  c[1] = a[0] * b[1] + a[1] * b[0];
  c[2] = a[0] * b[2] + a[1] * b[1] + a[2] * b[0];
  c[3] = a[3] * b[0] + a[2] * b[1] + a[1] * b[2] + a[0] * b[3];

  EXPECT_TRUE(a * b == c);
}
TEST(Expansion, OrderByOrderSquare)
{
  // the reasoning behind this test is that if both operands are given as
  // reference,
  // they must not change during the multiplication, because they would destroy
  // the result.
  Expansion<double, 3> a, c;
  a[0] = 1.2;
  a[1] = 4.5;
  a[2] = 7.8;

  c[0] = a[0] * a[0];
  c[1] = a[0] * a[1] + a[1] * a[0];
  c[2] = a[0] * a[2] + a[1] * a[1] + a[2] * a[0];

  EXPECT_TRUE(a * a == c);
}
TEST(Expansion, InitializationWithSingleElement)
{
  Expansion<double, NORD> a(123.), b;
  b[0] = 123.;
  // not necessary
  for (std::size_t i = 1; i < NORD; ++i) {
    b[i] = 0.0;
  }
  EXPECT_TRUE(a == b);
}
TEST(Expansion, AssignmentWithElement)
{
  Expansion<double, NORD> a, b(123);
  // not necessary, just to set it to sth different than zero.
  for (std::size_t i = 1; i < NORD; ++i) {
    a[i] = 2.0;
  }
  a = 123.;
  EXPECT_TRUE(a == b);
}
TEST(Expansion, CompareToElement)
{
  Expansion<double, NORD> a(123.);
  EXPECT_TRUE(a == 123.);
}
TEST(Expansion, InvertExpansion)
{
  Expansion<double, NORD> a, ainv, one(1);
  for (std::size_t i = 0; i < NORD; ++i) {
    a[i] = static_cast<double>(i + 1);
  }
  ainv = inverse(a);
  EXPECT_TRUE(one == a * ainv);
}
TEST(Expansion, InversionFailsForZeroZerothOrder)
{
  Expansion<double, NORD> a(0);
  EXPECT_THROW(inverse(a), ExpansionError);
}
TEST(Expansion, ExponentialCoefficients)
{
  constexpr auto res = ::expcoefficients<NORD>();
  double prev = 1.0;
  EXPECT_DOUBLE_EQ(res[0], 1.0);
  for (std::size_t i = 1; i < NORD; ++i) {
    EXPECT_DOUBLE_EQ(res[i], prev / static_cast<double>(i));
    prev = res[i];
  }
}
TEST(Expansion, ExpExpansionResemblesExp)
{
  Expansion<double, 10> test(0.1);
  const auto exptest = exp(test);
  EXPECT_NEAR(exptest[0], std::exp(test[0]), 1e-9);
  for (std::size_t ord = 1; ord < 10; ++ord) {
    EXPECT_DOUBLE_EQ(exptest[ord], 0.0);
  }
}
TEST(Expansion, LogCoefficients)
{
  constexpr auto res = ::logcoefficients<NORD>();
  EXPECT_DOUBLE_EQ(res[0], 0.0);
  for (std::size_t i = 1; i < NORD; ++i) {
    EXPECT_DOUBLE_EQ(res[i],
                     std::pow((-1.0), (i + 1)) / static_cast<double>(i));
  }
}
TEST(Expansion, LogExpansionResemblesLog)
{
  Expansion<double, 10> test(0.1);
  const auto logtest = log(test);
  EXPECT_NEAR(logtest[0], std::log(test[0]), 1e-9);
  for (std::size_t ord = 1; ord < 10; ++ord) {
    EXPECT_DOUBLE_EQ(logtest[ord], 0.0);
  }
}
TEST(Expansion, CanTakeExpOfDoubleExpansionZerothOrder)
{
  Expansion<double, 10> test(0.5); // this is an expansion (0.5, 0, 0, ...)

  auto exptest = exp(test);

  EXPECT_NEAR(exptest[0], std::exp(0.5), 1e-8);
  for (std::size_t i = 1; i < 10; ++i) {
    EXPECT_DOUBLE_EQ(exptest[i], 0.0);
  }
}
TEST(Expansion, GetMathematicaResultForExp)
{
  Expansion<double, 4> a;
  a[0] = 0.123; // very pseudo-random numbers.
  a[1] = 0.378;
  a[2] = 0.739;
  a[3] = 1.202;

  // got this from
  // Series[Exp[A0+A1*b+A2*b*b+A3*b*b*b], {b, 0, 3}]
  Expansion<double, 4> should;
  should[0] = exp(a[0]);
  should[1] = exp(a[0]) * a[1];
  should[2] = exp(a[0]) * (a[1] * a[1] + 2 * a[2]) * 0.5;
  should[3] =
    exp(a[0]) * (a[1] * a[1] * a[1] + 6 * a[1] * a[2] + 6 * a[3]) / 6.0;

  auto res = exp(a);
  for (std::size_t i = 0u; i < 4; ++i) {
    EXPECT_DOUBLE_EQ(res[i], should[i]) << "at order " << i;
  }
}
TEST(Expansion, CanTakeLogAndExpOfExpansion)
{
  Expansion<std::complex<double>, NORD> test(randomVector<NORD>());
  test *= 0.01;
  auto reconstr1 = log(exp(test));
  auto reconstr2 = exp(log(test));
  for (std::size_t i = 0; i < NORD; ++i) {
    EXPECT_NEAR(std::real(reconstr1[i]), std::real(test[i]), 1.e-9);
    EXPECT_NEAR(std::imag(reconstr1[i]), std::imag(test[i]), 1.e-9);
    EXPECT_NEAR(std::real(reconstr2[i]), std::real(test[i]), 1.e-9);
    EXPECT_NEAR(std::imag(reconstr2[i]), std::imag(test[i]), 1.e-9);
  }
}
TEST(Expansion, ExpComparesToOldFortranCode)
{
  Expansion<double, 7> testfield;
  testfield[0] = 0.0;
  testfield[1] = 0.1;
  testfield[2] = 0.2;
  testfield[3] = 0.3;
  testfield[4] = 0.4;
  testfield[5] = 0.5;
  testfield[6] = 0.6;
  auto res = exp(testfield);
  EXPECT_DOUBLE_EQ(res[0], 1.0);
  EXPECT_NEAR(res[1], testfield[1], 1e-7);
  EXPECT_NEAR(res[2], 0.20500000312924385, 1e-7);
  EXPECT_NEAR(res[3], 0.32016667919109265, 1e-7);
  EXPECT_NEAR(res[4], 0.45100417490725714, 1e-7);
  EXPECT_NEAR(res[5], 0.60353342133272314, 1e-7);
  EXPECT_NEAR(res[6], 0.78448419917942502, 1e-7);
}
TEST(Expansion, LogComparesToOldFortranCode)
{
  Expansion<double, 7> testfield;
  testfield[0] = 1.0;
  testfield[1] = 0.1;
  testfield[2] = 0.2;
  testfield[3] = 0.3;
  testfield[4] = 0.4;
  testfield[5] = 0.5;
  testfield[6] = 0.6;
  auto res = log(testfield);
  EXPECT_DOUBLE_EQ(res[0], 0.0);
  EXPECT_NEAR(res[1], testfield[1], 1e-7);
  EXPECT_NEAR(res[2], 0.19500000283122063, 1e-7);
  EXPECT_NEAR(res[3], 0.28033334467311699, 1e-7);
  EXPECT_NEAR(res[4], 0.35197500381320707, 1e-7);
  EXPECT_NEAR(res[5], 0.40680199590530985, 1e-7);
  EXPECT_NEAR(res[6], 0.44278651820920401, 1e-7);
}

TEST(Expansion, CanExpandCustomTypes)
{
  // Can we expand SU3?
  Expansion<Su3, NORD> s(Su3Consts::one);
  EXPECT_EQ(s[0], Su3Consts::one);
  for (std::size_t i = 1; i < NORD; ++i) {
    EXPECT_EQ(s[i], Su3Consts::zero);
  }

  // Can we expand CP
  // Expansion< CP<double,2>, NORD> c(1.);

  // std::cout << s[0] << std::endl;
  // std::cout << s[1] << std::endl;
  // std::cout << s[2] << std::endl;
  // std::cout << s[3] << std::endl;
  // std::cout << s[4] << std::endl;
}

TEST(Expansion, CanExpandAndDefaultConstructCustomTypes)
{
  // Can we use the default constructor?
  Expansion<Su3, NORD> s;
  for (std::size_t i = 0; i < NORD; ++i) {
    EXPECT_EQ(s[i], Su3Consts::zero);
  }
}
TEST(Expansion, CanExpandCustomTypesWithTypeConversionConstructor)
{
  /// this takes an int as a argument but initializes the expansion of Su3
  /// with static_cast<Su3>(0)
  Expansion<Su3, NORD> s(0);
}


TEST(Expansion, CanCastToLowerOrderExpansion)
{
  Expansion<double, NORD> b;
  Expansion<double, NORD - 1> c, a;

  // not necessary, just to set it to sth different than zero.
  for (std::size_t i = 1; i < NORD - 1; ++i) {
    b[i] = 2.0;
    c[i] = 2.0;
  }
  b[NORD - 1] = 4.;

  a = b.cast_lower<NORD - 1>();

  EXPECT_TRUE(a == c);
}

TEST(Expansion, CanCastToHigherOrderExpansion)
{
  Expansion<double, NORD> b;
  Expansion<double, NORD + 2> c, a;

  // not necessary, just to set it to sth different than zero.
  for (std::size_t i = 1; i < NORD - 1; ++i) {
    b[i] = 2.0;
    c[i] = 2.0;
  }
  c[NORD] = 0.0;
  c[NORD + 1] = 0.0;

  a = b.cast_higher<NORD + 2>();

  EXPECT_TRUE(a == c);
}


TEST(Expansion, CanUseRealAndImagForU1)
{

  Expansion<U1, NORD> u;
  Expansion<double, NORD> re, im;

  for (std::size_t i = 0; i < NORD - 1; ++i) {
    u[i] = std::complex<double>(1.0 * static_cast<double>(i),
                                2.0 * static_cast<double>(i));
    re[i] = 1.0 * static_cast<double>(i);
    im[i] = 2.0 * static_cast<double>(i);
  }

  auto ur = real(u);
  auto ui = imag(u);
  for (std::size_t i = 0; i < NORD - 1; ++i) {

    EXPECT_EQ(ur[i], re[i]);
    EXPECT_EQ(ui[i], im[i]);
  }
}
