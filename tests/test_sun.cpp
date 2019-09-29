#include "gaugegroups/sun.h"
#include "test_helpers.h"
#include <gtest/gtest.h>

#ifdef BLAS_AVAIL

TEST(SuN, ConstructsOne)
{
  const std::size_t N = 5;
  SuN<N> one(1);
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) {
      if (i == j)
        EXPECT_EQ(one(i, j), std::complex<double>(1.0, 0.0));
      else
        EXPECT_EQ(one(i, j), std::complex<double>(0.0, 0.0));
    }
  }
}


TEST(SuN, GroupMultiplication)
{

  std::complex<double> II(0., 1.);

  std::array<std::complex<double>, 4> asx = { 0., 1, 1, 0 };
  std::array<std::complex<double>, 4> asy = { 0., -II, II, 0 };
  std::array<std::complex<double>, 4> asz = { 1, 0, 0, -1 };

  SuN<2> sx(asx), sy(asy), sz(asz), res;

  res = sx * sy - sy * sx;

  EXPECT_TRUE((res / (2. * II)) == sz);
  EXPECT_FALSE(sy == sz);
}

TEST(SuN, CanConstructGenerators)
{

  const std::size_t N = 3;

  auto res = Generators<N>();

  SuNAlgebra<N> Casimir(0.0);

  for (auto mat : res) {
    // std::cout << mat << std::endl;
    Casimir += (mat * mat);
  }


  auto tcm = std::abs(trace(Casimir));
  EXPECT_DOUBLE_EQ(tcm, 0.5 * static_cast<double>(N * N - 1LU));

  // std::cout << "Trace: " <<  tcm << std::endl;
}


#endif // CBLAS available
