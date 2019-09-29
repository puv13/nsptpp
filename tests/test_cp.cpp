///! \file
/// Automated tests for the CP class and related code

#include "lattice.h"
#include "latticefields/cp.h"
#include "nummat.h"
#include <array>
#include <complex>
#include <fstream>
#include <gtest/gtest.h>
#include <iostream>
#include <string>

TEST(CP, CanUseStdCtor)
{
  const CP<std::complex<double>, 1> cp1;
  EXPECT_EQ(cp1[0], 0.0);
  CP<int, 1> cp2;
  EXPECT_EQ(cp2[0], 0);
}


TEST(CP, CanUseInitCtor)
{
  const CP<std::complex<double>, 1> cp1(std::complex<double>(5.0));
  EXPECT_EQ(cp1[0], 5.0);
}

TEST(CP, CanUseCopyCtor)
{
  const CP<std::complex<double>, 2> cp1(std::complex<double>(5.0));
  CP<std::complex<double>, 2> cp2(cp1);
  EXPECT_EQ(cp1[1], cp2[1]);
}

TEST(CP, CanCopy)
{
  const CP<std::complex<double>, 2> cp1(std::complex<double>(5.0));
  CP<std::complex<double>, 2> cp2, cp3(std::complex<double>(.0, 1.));
}

TEST(CP, CanComputeNorm)
{
  const CP<std::complex<double>, 25> cp(std::complex<double>(1.0));
  EXPECT_DOUBLE_EQ(cp.norm(), 5.0);
  // std::cout << decltype(cp.norm()) << std::endl;
}

TEST(CP, CanUseArithmeticOps)
{
  const CP<std::complex<double>, 25> cp1(std::complex<double>(1.0));
  const CP<std::complex<double>, 25> cp2(std::complex<double>(2.0));

  CP<std::complex<double>, 25> diff(cp1);
  CP<std::complex<double>, 25> sum(cp1);
  diff -= cp2;
  sum += diff;

  for (std::size_t i = 0; i < 25; i++) {

    EXPECT_EQ(diff[i], std::complex<double>(-1., 0.));
    EXPECT_EQ(sum[i], std::complex<double>(0., 0.));
  }


  diff *= std::complex<double>(0., 1.);
  for (std::size_t i = 0; i < 25; i++) {

    EXPECT_EQ(diff[i], std::complex<double>(0., -1.));
  }

  const CP<std::complex<double>, 25> cp3 = cp1 * std::complex<double>(0., 2.);
  for (std::size_t i = 0; i < 25; i++) {

    EXPECT_EQ(cp3[i], std::complex<double>(0., 2.));
  }


  // Left multiplication
  const CP<std::complex<double>, 25> cp4 = std::complex<double>(0., 2.) * cp1;
  for (std::size_t i = 0; i < 25; i++) {

    EXPECT_EQ(cp4[i], std::complex<double>(0., 2.));
  }

  const CP<std::complex<double>, 25> cp5 = 2L * cp1;
  for (std::size_t i = 0; i < 25; i++) {

    EXPECT_EQ(cp5[i], std::complex<double>(2., 0.));
  }
}

TEST(CP, Equality)
{
  CP<std::complex<double>, 25> cp1(20.);
  CP<std::complex<double>, 25> cp2(cp1);

  EXPECT_EQ(cp1, cp2);

  cp2[12] = 20 - 1.e-14;
  EXPECT_EQ(cp1 == cp2, false);

  cp2[12] = 20;
  EXPECT_EQ(cp1, cp2);
}

TEST(CP, CanIterateOverField)
{
  CP<std::complex<double>, 25> cp;
  for (std::size_t i = 0; i < 25; i++) {
    cp[i] = static_cast<double>(i + 1);
  }


  size_t i = 0;
  for (auto el : cp) {
    EXPECT_EQ(el, std::complex<double>(static_cast<double>(++i)));
  }
}

TEST(CP, CanAccessFieldComponentsWithIterator)
{
  CP<std::complex<double>, 25> cp(0.);

  // Initialise
  size_t i = 0;
  for (auto it = cp.begin(); it < cp.end(); ++it) {
    *it = static_cast<double>(++i);
  }

  i = 0;
  for (auto el : cp) {
    EXPECT_EQ(el, std::complex<double>(static_cast<double>(++i)));
  }

  i = 5;
  auto it = CPiterator<CP<std::complex<double>, 25>>(i, cp);

  auto cconst = std::complex<double>(3, 3);
  *it = cconst;

  EXPECT_EQ(cp[i], cconst);
}


TEST(CP, ThrowsIfIndexOutOfBounds)
{
  // Check if runtime error is thrown
  CP<std::complex<double>, 5> cp(0.);
  // ASSERT_DEATH({cp.at(6);}, "Index out of range at \\s*\\S*:\\s*\\S*" );
  // This death test uses a compound statement.
  ASSERT_THROW(
    {
      size_t n = 6;
      cp.at(n);
    },
    std::runtime_error);
}


// TEST(CP, CanComputeTracePP)
// {

//     CP<std::complex<double>,5> cp(1);

//     auto tr = tracePP(cp,cp);

//     EXPECT_DOUBLE_EQ(tr,1.0);

// }


TEST(CP, RandomIsGood)
{
  int ntests = 1000;
  // Check for N=1
  std::ofstream CPfile;
  // Debug output to check uniformity
  // CPfile.open("./CPrand.dat",std::ios::out);
  for (auto i = 0; i < ntests; ++i) {
    auto CP1 = randomCP<1lu>();
    auto CP20 = randomCP<20lu>();
    EXPECT_DOUBLE_EQ(CP1.norm(), 1.);
    EXPECT_DOUBLE_EQ(CP20.norm(), 1.);

    // CPfile.precision(5);
    // CPfile << std::scientific << std::real(CP1[0]) << "\t" <<
    // std::imag(CP1[0]) << "\n";
  }
  // CPfile.close();


  int seed = 13;
  std::default_random_engine drng(seed);
  auto CPfirst = randomCP<100lu>(drng);

  for (auto i = 0; i < ntests; ++i) {
    auto CP10 = randomCP<10lu>(drng);
    auto CP40 = randomCP<40lu>(drng);
    EXPECT_DOUBLE_EQ(CP10.norm(), 1.);
    EXPECT_DOUBLE_EQ(CP40.norm(), 1.);
  }

  std::default_random_engine drng2(seed);
  auto CPseed = randomCP<100lu>(drng2);

  EXPECT_EQ(CPfirst, CPseed);
}

TEST(CP, ScalarProd)
{

  CP<std::complex<double>, 25> cp(0.);
  CP<std::complex<double>, 25> right(0.);

  double N = 25;
  double spexp = 1 / 6. * N * (N + 1) * (2 * N + 1);
  // Initialise
  size_t i = 0;
  for (auto it = cp.begin(); it < cp.end(); ++it) {
    *it = static_cast<double>(++i);
  }

  right = cp;

  double res = std::real(scalar_prod(cp, right));

  EXPECT_DOUBLE_EQ(res, spexp);

  std::complex<double> II(0, 1.);

  auto cres = scalar_prod(cp, right * II);

  EXPECT_NEAR(std::abs(cres - II * spexp), 0.0, 1.e-12);

  cres = scalar_prod(cp * II, right);

  EXPECT_NEAR(std::abs(cres + II * spexp), 0.0, 1.e-12);

  U1 link(II);
  cres = scalar_prod(cp * link, right);

  EXPECT_NEAR(std::abs(cres + II * spexp), 0.0, 1.e-12);
}


TEST(CP, TwistCanInit)
{
  std::array<std::complex<double>, 25> initar;
  for (auto i = 0LU; i < 25; ++i) {
    initar[i] = static_cast<double>(i);
  }


  CP<std::complex<double>, 25> cp(0.);
  CPtwistPhase<std::complex<double>, 25> cptp1;
  CPtwistPhase<std::complex<double>, 25> cptp2(2.0);
  CPtwistPhase<std::complex<double>, 25> cptp3(initar);

  for (auto i = 0LU; i < 25; ++i) {
    EXPECT_EQ(cptp1[i], std::complex<double>(1.0));
    EXPECT_EQ(cptp2[i], std::complex<double>(2.0));
    EXPECT_EQ(cptp3[i], initar[i]);
  }
}

TEST(CP, TwistArithemtics)
{
  std::array<std::complex<double>, 25> initar;
  for (auto i = 0LU; i < 25; ++i) {
    initar[i] = static_cast<double>(i);
  }

  CP<std::complex<double>, 25> cp(1.);
  std::complex<double> cplx(1., .2);
  U1 u1(-1.0);


  CPtwistPhase<std::complex<double>, 25> cptp1;

  auto cptp2 = cptp1;

  // // Exponent of twist phase
  // for (auto i=0LU; i<25; ++i)
  // {
  // 	EXPECT_EQ(cptp2[i],std::complex<double>(1.0));
  // }

  // Multiplication with complex number

  cptp2 = cptp2 * cplx;
  cp = cp * cptp2; // Left and Right mult. with CPtwistPhase are different!!!
  cp = cptp2 * cp;
  auto cptp3 = cplx * cptp2;

  for (auto i = 0LU; i < 25; ++i) {
    EXPECT_EQ(cptp2[i], cplx);
    EXPECT_EQ(cptp3[i], cplx * cplx);
    EXPECT_EQ(cp[i], cplx * std::conj(cplx));
  }

  // Multiplication with u1
  cptp2 = cptp1;
  cptp2 = cptp2 * u1;
  cp = cp * cptp2; // Left and Right mult. with CPtwistPhase are different!!!
  cp = cptp2 * cp;
  cptp3 = u1 * cptp2;

  for (auto i = 0LU; i < 25; ++i) {
    EXPECT_EQ(cptp2[i], u1.value());
    EXPECT_EQ(cptp3[i], u1.value() * u1.value());
    EXPECT_EQ(cp[i], u1.value() * cplx * std::conj(u1.value() * cplx));
  }
}

#ifdef BLAS_AVAIL
TEST(CP, CanMultiplyByNummat)
{
  CP<std::complex<double>, 3> rhs(1.);
  std::array<std::complex<double>, 3 * 3> matarr{ 1, 1, 1, 2, 2, 2, 3, 3, 3 };
  nummat<3> mat(matarr);

  // std::cout << rhs << std::endl << mat;

  auto res = mat * rhs;

  for (auto i = 0LU; i < 3LU; ++i) {
    EXPECT_EQ(static_cast<double>(i + 1) * 3, res[i]);
  }


  // std::cout << res ;
}
#endif // BLAS_AVAIL


// TEST(CP, TwistedBC)
// {
//     CP<std::complex<double>,4> cp(1.);
//     double pi2 = std::atan(1.)*2.;
//     CPtwistPhase<std::complex<double>,4> tp(pi2);

//     std::array<CPtwistPhase<std::complex<double>,4>, 2> tp_ar;
//     tp_ar.fill(tp);

//     BoundaryCondition<double, 2> pbc({1.,1.},"PERIODIC");
//     BoundaryCondition<CPtwistPhase<std::complex<double>,4>, 2>
//     abc(tp_ar,"AP_TWIST");

//     U1 u1_init(std::complex<double>(1.,0.));
//     FullLattice<CP<std::complex<double>,4>,U1, 2> cplat ({2,2},cp,u1_init);
//     auto  idx = cplat.coordToLinearIndex({0,0});
//     FullLatticeIterator<FullLattice<CP<std::complex<double>,4>,U1, 2> >
//     fli(idx,cplat);


//     auto orig = (cplat.at({0,1})).site();
//     auto back = fli.neighborSite(1,Direction::BACKWARD,abc);
//     auto fwd  = (cplat.at({1,1})).neighborSite(0,Direction::FORWARD,abc);

//     //std::cout << orig << "\n" << back << "\n" << fwd  << std::endl;


// }


// TEST(CP, CanNormalise)
// {
//     CP<std::complex<double>,25> cp(std::complex<double>(1.0));
//     cp.normalise();
//     EXPECT_DOUBLE_EQ(std::abs(cp.norm(),1.0);
// }
