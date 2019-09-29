///! \file
/// Automated tests for the CP class and related code
/// Needs a CBLAS compatible blas library!

#include "nummat.h"
#include <gtest/gtest.h>
#include <iostream>


TEST(NUMMAT, CanUseStdCtor)
{
  const nummat<2> nm1;

  EXPECT_EQ(nm1(0LU, 0LU), std::complex<double>(0.0));
}


TEST(NUMMAT, CanComputeFNorm)
{
  const nummat<2> nm1(1);

  EXPECT_DOUBLE_EQ(sqrt(2), normF(nm1));
}


TEST(NUMMAT, CanMultiplyMatrix)
{

  std::complex<double> II(0., 1.);

  std::array<std::complex<double>, 4> asx = { 0., 1, 1, 0 };
  std::array<std::complex<double>, 4> asy = { 0., -II, II, 0 };
  std::array<std::complex<double>, 4> asz = { 1, 0, 0, -1 };

  nummat<2> sx(asx), sy(asy), sz(asz), res;

  res = sx * sy - sy * sx;

  EXPECT_TRUE((res / (2. * II)) == sz);
  EXPECT_FALSE(sy == sz);

  // std::cout << ((res/(2.*II)) == sz) << std::endl;
}

TEST(NUMMAT, CanMultiplyScalar)
{

  std::complex<double> SC(12.35, 67.89);

  std::array<std::complex<double>, 4> Aar = { 1., 0., 0., 1. };
  std::array<std::complex<double>, 4> Sar = { SC, 0., 0., SC };

  nummat<2> A(Aar), S(Sar);

  // right scalar multiplication
  auto res = A * SC;
  EXPECT_TRUE(res == S);

  // left scalar multiplication
  res = SC * A;
  EXPECT_TRUE(res == S);

  // *= operator
  A *= SC;
  EXPECT_TRUE(A == S);
}

TEST(NUMMAT, MatrixMatrixAlgebra)
{

  std::array<std::complex<double>, 9> Aar = { 468, 279, 934, 588, 595,
                                              869, 210, 987, 769 };


  std::array<std::complex<double>, 9> Bar = { 336, 489, 735, 573, 708,
                                              815, 589, 926, 179 };


  std::array<std::complex<double>, 9> Car = { 132, -210, 199, 15, -113,
                                              54,  -379, 61,  590 };


  nummat<3> A(Aar), B(Bar), C(Car), D(C);


  // Subtraction
  EXPECT_TRUE(A - B == C);

  // Addition
  EXPECT_TRUE(A == C + B);


  // Operator +=
  D += B;
  EXPECT_TRUE(A == D);

  // Operator -=
  A -= B;
  EXPECT_TRUE(A == C);
}


TEST(NUMMAT, Trace)
{

  nummat<20> nm20(1.0);
  nummat<40> nm40(1.0);
  nummat<20> nm20_10(10.0);
  const std::array<std::complex<double>, 10> ar(
    { 1., 2., 3., 4., 5., 6., 7., 8., 9., 10 });
  auto nm10 = diag(ar);


  EXPECT_EQ(trace(nm20), 20.);
  EXPECT_EQ(trace(nm40), 40.);
  EXPECT_EQ(trace(nm20_10), 200.);
  EXPECT_EQ(trace(nm10), 55.);
}


TEST(NUMMAT, Exp)
{
  nummat<20> nm20 = diag(std::array<std::complex<double>, 20>{
    1.,  2.,  3.,  4.,  5.,  6.,  7.,  8.,  9.,  10.,
    11., 12., 13., 14., 15., 16., 17., 18., 19., 20 });


  const double error = 1.e-12;
  auto res = exp(nm20, error, 200);
  for (auto i = 0LU; i < 2LU; ++i) {
    // std::cout << std::real(res(i,i)) << "\t" << std::exp(i+1) <<
    // std::endl;

    EXPECT_NEAR(std::real(res(i, i)), std::exp(i + 1), error);
  }
}


TEST(NUMMAT, MultiplywithHermitianconjugate)
{
  std::complex<double> II(0., 1.);
  std::array<std::complex<double>, 4> ar2{ II, 0., 0., II };
  std::array<std::complex<double>, 4> res{ 1., 0., 0., 1. };
  nummat<2> nm2(ar2);
  nummat<2> nmres(res);
  nummat<2> nm3 = nm2.mult_conj(nm2);

  EXPECT_TRUE(nmres == nm3);
}


TEST(NUMMAT, CanHermitianConjugate)
{
  std::complex<double> II(0., 1.);
  std::array<std::complex<double>, 4> ar2{ II, 2. * II, -3. * II, II };
  std::array<std::complex<double>, 4> solar{ -1. * II, 3. * II, -2. * II,
                                             -1. * II };
  nummat<2> nm2(ar2);
  nummat<2> sol(solar);
  nummat<2> res;

  res = nm2.dagger();

  EXPECT_TRUE(res == sol);
}
