#include "gaugegroups/su3.h"
#include "test_helpers.h"
#include <gtest/gtest.h>

TEST(Su3, ConstructsOne)
{
  Su3 one(1);
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      if (i == j)
        EXPECT_EQ(one(i, j), std::complex<double>(1.0, 0.0));
      else
        EXPECT_EQ(one(i, j), std::complex<double>(0.0, 0.0));
    }
  }
}
TEST(Su3, ComplexConstAddition)
{
  constexpr std::complex<double> a =
    std::complex<double>(1.2, 3.0) + std::complex<double>(4.3, 2.9);

  // make sure that we use the standard addition:
  std::complex<double> b = std::operator+(std::complex<double>(1.2, 3.0),
                                          std::complex<double>(4.3, 2.9));
  EXPECT_EQ(a, b);
}
TEST(Su3, ComplexConstSubtraction)
{
  constexpr std::complex<double> a =
    std::complex<double>(1.2, 3.0) - std::complex<double>(4.3, 2.9);

  // make sure that we use the standard subtraction:
  std::complex<double> b = std::operator-(std::complex<double>(1.2, 3.0),
                                          std::complex<double>(4.3, 2.9));
  EXPECT_EQ(a, b);
}
TEST(Su3, ComplexConstMultiplication)
{
  constexpr std::complex<double> a =
    std::complex<double>(1.2, 3.0) * std::complex<double>(4.3, 2.9);

  // make sure that we use the standard multiplication:
  std::complex<double> b = std::operator*(std::complex<double>(1.2, 3.0),
                                          std::complex<double>(4.3, 2.9));
  EXPECT_EQ(a, b);
}
TEST(Su3, HermitianConjugate)
{
  std::ranlux24 rd;
  Su3Algebra u = random(rd);
  ThreeByThreeMatrix udag = dagger(u);
  for (std::size_t a = 0; a < 3; ++a) {
    for (std::size_t b = 0; b < 3; ++b) {
      EXPECT_DOUBLE_EQ(std::real(udag(a, b)), std::real(std::conj(u(b, a))));
      EXPECT_DOUBLE_EQ(std::imag(udag(a, b)), std::imag(std::conj(u(b, a))));
    }
  }
}
TEST(Su3, CanTakeExpOfDiagonalMatrix)
{
  ThreeByThreeMatrix two(2.0);
  auto res = exp(two);
  EXPECT_DOUBLE_EQ(std::real(trace(res)) / 3.0, std::exp(2.0));
}
TEST(Su3, CanTakeExpOfMatrix)
{
  // compare to octave:
  ThreeByThreeMatrix a;
  a(0, 0) = 0.12;
  a(0, 1) = 0.23;
  a(0, 2) = 0.34;
  a(1, 0) = 0.45;
  a(1, 1) = 0.56;
  a(1, 2) = 0.67;
  a(2, 0) = 0.78;
  a(2, 1) = 0.89;
  a(2, 2) = 0.90;
  auto res = exp(a);
  EXPECT_NEAR(std::real(res(0, 0)), 1.48832, 1e-4);
  EXPECT_NEAR(std::imag(res(0, 0)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(0, 1)), 0.67668, 1e-4);
  EXPECT_NEAR(std::imag(res(0, 1)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(0, 2)), 0.83420, 1e-4);
  EXPECT_NEAR(std::imag(res(0, 2)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(1, 0)), 1.25900, 1e-4);
  EXPECT_NEAR(std::imag(res(1, 0)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(1, 1)), 2.55112, 1e-4);
  EXPECT_NEAR(std::imag(res(1, 1)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(1, 2)), 1.77866, 1e-4);
  EXPECT_NEAR(std::imag(res(1, 2)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(2, 0)), 1.95769, 1e-4);
  EXPECT_NEAR(std::imag(res(2, 0)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(2, 1)), 2.34040, 1e-4);
  EXPECT_NEAR(std::imag(res(2, 1)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(2, 2)), 3.53233, 1e-4);
  EXPECT_NEAR(std::imag(res(2, 2)), 0., 1e-15);
}
TEST(Su3, CanTakeLogOfDiagonalMatrix)
{
  ThreeByThreeMatrix a(0.0);
  a(0, 0) = 0.1;
  a(1, 1) = 0.2;
  a(2, 2) = 0.3;
  auto logA = log(a);
  ThreeByThreeMatrix shouldBeLogA(0.0);
  shouldBeLogA(0, 0) = std::log(0.1);
  shouldBeLogA(1, 1) = std::log(0.2);
  shouldBeLogA(2, 2) = std::log(0.3);
  for (std::size_t c = 0; c < 3; ++c) {
    for (std::size_t b = 0; b < 3; ++b) {
      EXPECT_NEAR(std::real(logA(c, b)), std::real(shouldBeLogA(c, b)), 1e-11);
      EXPECT_NEAR(std::imag(logA(c, b)), std::imag(shouldBeLogA(c, b)), 1e-11);
    }
  }
}
TEST(Su3, CanTakeSqrtOfMatrix)
{
  ThreeByThreeMatrix a;
  a(0, 0) = 0.12;
  a(0, 1) = 0.23;
  a(0, 2) = 0.34;
  a(1, 0) = 0.45;
  a(1, 1) = 0.56;
  a(1, 2) = 0.67;
  a(2, 0) = 0.78;
  a(2, 1) = 0.89;
  a(2, 2) = 0.90;
  auto b = a * a;
  ThreeByThreeMatrix sqrtB = sqrt(b);
  auto c = sqrtB * sqrtB;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      EXPECT_NEAR(std::real(c(i, j)), std::real(b(i, j)), 1e-11);
      EXPECT_NEAR(std::imag(c(i, j)), std::imag(b(i, j)), 1e-11);
    }
  }
}
TEST(Su3, CanTakeLogOfMatrix)
{
  ThreeByThreeMatrix a;
  a(0, 0) = 0.12;
  a(0, 1) = 0.23;
  a(0, 2) = 0.34;
  a(1, 0) = 0.45;
  a(1, 1) = 0.56;
  a(1, 2) = 0.67;
  a(2, 0) = 0.78;
  a(2, 1) = 0.89;
  a(2, 2) = 0.90;
  auto b = a * dagger(a);
  ASSERT_NEAR(std::real(b(0, 0)), 0.18290, 1e-5);
  ASSERT_NEAR(std::imag(b(0, 0)), 0., 1e-15);
  ASSERT_NEAR(std::real(b(0, 1)), 0.41060, 1e-5);
  ASSERT_NEAR(std::imag(b(0, 1)), 0., 1e-15);
  ASSERT_NEAR(std::real(b(0, 2)), 0.60430, 1e-5);
  ASSERT_NEAR(std::imag(b(0, 2)), 0., 1e-15);
  ASSERT_NEAR(std::real(b(1, 0)), 0.41060, 1e-5);
  ASSERT_NEAR(std::imag(b(1, 0)), 0., 1e-15);
  ASSERT_NEAR(std::real(b(1, 1)), 0.96500, 1e-5);
  ASSERT_NEAR(std::imag(b(1, 1)), 0., 1e-15);
  ASSERT_NEAR(std::real(b(1, 2)), 1.45240, 1e-5);
  ASSERT_NEAR(std::imag(b(1, 2)), 0., 1e-15);
  ASSERT_NEAR(std::real(b(2, 0)), 0.60430, 1e-5);
  ASSERT_NEAR(std::imag(b(2, 0)), 0., 1e-15);
  ASSERT_NEAR(std::real(b(2, 1)), 1.45240, 1e-5);
  ASSERT_NEAR(std::imag(b(2, 1)), 0., 1e-15);
  ASSERT_NEAR(std::real(b(2, 2)), 2.21050, 1e-5);
  ASSERT_NEAR(std::imag(b(2, 2)), 0., 1e-15);
  auto c = log(b);
  EXPECT_NEAR(std::real(c(0, 0)), -5.172148, 1e-5);
  EXPECT_NEAR(std::imag(c(0, 0)), 0., 1e-15);
  EXPECT_NEAR(std::real(c(0, 1)), 2.642320, 1e-5);
  EXPECT_NEAR(std::imag(c(0, 1)), 0., 1e-15);
  EXPECT_NEAR(std::real(c(0, 2)), 0.025184, 1e-5);
  EXPECT_NEAR(std::imag(c(0, 2)), 0., 1e-15);
  EXPECT_NEAR(std::real(c(1, 0)), 2.642320, 1e-5);
  EXPECT_NEAR(std::imag(c(1, 0)), 0., 1e-15);
  EXPECT_NEAR(std::real(c(1, 1)), -5.029185, 1e-5);
  EXPECT_NEAR(std::imag(c(1, 1)), 0., 1e-15);
  EXPECT_NEAR(std::real(c(1, 2)), 3.383145, 1e-5);
  EXPECT_NEAR(std::imag(c(1, 2)), 0., 1e-15);
  EXPECT_NEAR(std::real(c(2, 0)), 0.025184, 1e-5);
  EXPECT_NEAR(std::imag(c(2, 0)), 0., 1e-15);
  EXPECT_NEAR(std::real(c(2, 1)), 3.383145, 1e-5);
  EXPECT_NEAR(std::imag(c(2, 1)), 0., 1e-15);
  EXPECT_NEAR(std::real(c(2, 2)), -1.035712, 1e-5);
  EXPECT_NEAR(std::imag(c(2, 2)), 0., 1e-15);
}
TEST(Su3, CanInvertDiagonalMatrix)
{
  ThreeByThreeMatrix two(2.0);
  auto res = inverse(two);
  EXPECT_NEAR(std::real(trace(res)) / 3.0, 0.5, 1.e-15);
}
TEST(Su3, CanInvertMatrix)
{
  // compare to octave:
  ThreeByThreeMatrix a;
  a(0, 0) = 0.12;
  a(0, 1) = 0.23;
  a(0, 2) = 0.34;
  a(1, 0) = 0.45;
  a(1, 1) = 0.56;
  a(1, 2) = 0.67;
  a(2, 0) = 0.78;
  a(2, 1) = 0.89;
  a(2, 2) = 0.90;
  // Cond numver of a is ~ 10^2
  // So we can expect a precision of ~ 10^(-15+2) = 10^(-13) in the result
  auto res = inverse(a);
  EXPECT_NEAR(std::real(res(0, 0)), -25.427, 1e-3);
  EXPECT_NEAR(std::imag(res(0, 0)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(0, 1)), 26.336, 1e-3);
  EXPECT_NEAR(std::imag(res(0, 1)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(0, 2)), -10.000, 1e-3);
  EXPECT_NEAR(std::imag(res(0, 2)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(1, 0)), 32.397, 1e-3);
  EXPECT_NEAR(std::imag(res(1, 0)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(1, 1)), -43.306, 1e-3);
  EXPECT_NEAR(std::imag(res(1, 1)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(1, 2)), 20.000, 1e-3);
  EXPECT_NEAR(std::imag(res(1, 2)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(2, 0)), -10.000, 1e-3);
  EXPECT_NEAR(std::imag(res(2, 0)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(2, 1)), 20.000, 1e-3);
  EXPECT_NEAR(std::imag(res(2, 1)), 0., 1e-15);
  EXPECT_NEAR(std::real(res(2, 2)), -10.000, 1e-3);
  EXPECT_NEAR(std::imag(res(2, 2)), 0., 1e-15);
  auto prod = res * a;
  EXPECT_NEAR(1.0, std::real(prod(0, 0)), 1e-13);
  EXPECT_NEAR(0.0, std::imag(prod(0, 0)), 1e-15);
  EXPECT_NEAR(0.0, std::real(prod(0, 1)), 1e-13);
  EXPECT_NEAR(0.0, std::imag(prod(0, 1)), 1e-15);
  EXPECT_NEAR(0.0, std::real(prod(0, 2)), 1e-13);
  EXPECT_NEAR(0.0, std::imag(prod(0, 2)), 1e-15);
  EXPECT_NEAR(0.0, std::real(prod(1, 0)), 1e-13);
  EXPECT_NEAR(0.0, std::imag(prod(1, 0)), 1e-15);
  EXPECT_NEAR(1.0, std::real(prod(1, 1)), 1e-13);
  EXPECT_NEAR(0.0, std::imag(prod(1, 1)), 1e-15);
  EXPECT_NEAR(0.0, std::real(prod(1, 2)), 1e-13);
  EXPECT_NEAR(0.0, std::imag(prod(1, 2)), 1e-15);
  EXPECT_NEAR(0.0, std::real(prod(2, 0)), 1e-13);
  EXPECT_NEAR(0.0, std::imag(prod(2, 0)), 1e-15);
  EXPECT_NEAR(0.0, std::real(prod(2, 1)), 1e-13);
  EXPECT_NEAR(0.0, std::imag(prod(2, 1)), 1e-15);
  EXPECT_NEAR(1.0, std::real(prod(2, 2)), 1e-13);
  EXPECT_NEAR(0.0, std::imag(prod(2, 2)), 1e-15);
}
TEST(Su3, CanNotAssignToGroupFrom3x3Matrix)
{
  Su3Algebra a;
  Su3 g;
  ThreeByThreeMatrix m;

  m = a; // this works.
  m = g; // also this.
  /*
  a = m; // this should not.
  a = g; // this should not.
  g = a; // this should not.
  g = m; // this should not.
  */
  // but it does when explicitly casted:
  a = static_cast<Su3Algebra>(m);
  a = static_cast<Su3Algebra>(g);
  g = static_cast<Su3>(a);
  g = static_cast<Su3>(m);
}
TEST(Su3, TypeChangingOperators)
{
  Su3Algebra algebraelem;
  Su3 groupelem;
  ThreeByThreeMatrix matrixelem;

  // algebra <> algebra
  const bool algTimesAlgIsMatrix =
    std::is_same<decltype(algebraelem * algebraelem),
                 ThreeByThreeMatrix>::value;
  EXPECT_TRUE(algTimesAlgIsMatrix);
  const bool algPlusAlgIsAlg =
    std::is_same<decltype(algebraelem + algebraelem), Su3Algebra>::value;
  EXPECT_TRUE(algPlusAlgIsAlg);
  const bool algMinusAlgIsAlg =
    std::is_same<decltype(algebraelem - algebraelem), Su3Algebra>::value;
  EXPECT_TRUE(algMinusAlgIsAlg);

  // group <> group:
  const bool groupTimesGroupIsGroup =
    std::is_same<decltype(groupelem * groupelem), Su3>::value;
  EXPECT_TRUE(groupTimesGroupIsGroup);
  const bool groupPlusGroupIsMatrix =
    std::is_same<decltype(groupelem + groupelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(groupPlusGroupIsMatrix);
  const bool groupMinusGroupIsMatrix =
    std::is_same<decltype(groupelem - groupelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(groupMinusGroupIsMatrix);

  // group <> algebra:
  const bool algTimesGroupIsMatrix =
    std::is_same<decltype(algebraelem * groupelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(algTimesGroupIsMatrix);
  const bool algPlusGroupIsMatrix =
    std::is_same<decltype(algebraelem + groupelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(algPlusGroupIsMatrix);
  const bool algMinusGroupIsMatrix =
    std::is_same<decltype(algebraelem - groupelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(algMinusGroupIsMatrix);

  // algebra <> group:
  const bool groupTimesAlgIsMatrix =
    std::is_same<decltype(groupelem * algebraelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(groupTimesAlgIsMatrix);
  const bool groupPlusAlgIsMatrix =
    std::is_same<decltype(groupelem + algebraelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(groupPlusAlgIsMatrix);
  const bool groupMinusAlgIsMatrix =
    std::is_same<decltype(groupelem - algebraelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(groupMinusAlgIsMatrix);

  // matrix <> matrix:
  const bool matrixTimesMatrixIsMatrix =
    std::is_same<decltype(matrixelem * matrixelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(matrixTimesMatrixIsMatrix);
  const bool matrixPlusMatrixIsMatrix =
    std::is_same<decltype(matrixelem + matrixelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(matrixPlusMatrixIsMatrix);
  const bool matrixMinusMatrixIsMatrix =
    std::is_same<decltype(matrixelem - matrixelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(matrixMinusMatrixIsMatrix);

  // alg <> matrix:
  const bool algTimesMatrixIsMatrix =
    std::is_same<decltype(algebraelem * matrixelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(algTimesMatrixIsMatrix);
  const bool algPlusMatrixIsMatrix =
    std::is_same<decltype(algebraelem + matrixelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(algPlusMatrixIsMatrix);
  const bool algMinusMatrixIsMatrix =
    std::is_same<decltype(algebraelem - matrixelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(algMinusMatrixIsMatrix);

  // group <> matrix:
  const bool groupTimesMatrixIsMatrix =
    std::is_same<decltype(groupelem * matrixelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(groupTimesMatrixIsMatrix);
  const bool groupPlusMatrixIsMatrix =
    std::is_same<decltype(groupelem + matrixelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(groupPlusMatrixIsMatrix);
  const bool groupMinusMatrixIsMatrix =
    std::is_same<decltype(groupelem - matrixelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(groupMinusMatrixIsMatrix);

  // matrix <> group:
  const bool matrixTimesGroupIsMatrix =
    std::is_same<decltype(matrixelem * groupelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(matrixTimesGroupIsMatrix);
  const bool matrixPlusGroupIsMatrix =
    std::is_same<decltype(matrixelem + groupelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(matrixPlusGroupIsMatrix);
  const bool matrixMinusGroupIsMatrix =
    std::is_same<decltype(matrixelem - groupelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(matrixMinusGroupIsMatrix);

  // matrix <> algebra:
  const bool matrixTimesAlgebraIsMatrix =
    std::is_same<decltype(matrixelem * algebraelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(matrixTimesAlgebraIsMatrix);
  const bool matrixPlusAlgebraIsMatrix =
    std::is_same<decltype(matrixelem + algebraelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(matrixPlusAlgebraIsMatrix);
  const bool matrixMinusAlgebraIsMatrix =
    std::is_same<decltype(matrixelem - algebraelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(matrixMinusAlgebraIsMatrix);

  // scalar <> *:
  const bool scalarTimesAlgebraIsAlgebra =
    std::is_same<decltype(1.2 * algebraelem), Su3Algebra>::value;
  EXPECT_TRUE(scalarTimesAlgebraIsAlgebra);
  const bool scalarTimesGroupIsMatrix =
    std::is_same<decltype(1.2 * groupelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(scalarTimesGroupIsMatrix);
  const bool scalarTimesMatrixIsMatrix =
    std::is_same<decltype(1.2 * matrixelem), ThreeByThreeMatrix>::value;
  EXPECT_TRUE(scalarTimesMatrixIsMatrix);
}
TEST(Su3, Su3AlgebraPlusEqualOperator)
{
  std::ranlux24 rd;
  Su3Algebra a = random(rd);
  Su3Algebra b = random(rd);
  Su3Algebra c(a);

  c += b;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t k = 0; k < 3; ++k) {
      EXPECT_EQ(c(i, k), a(i, k) + b(i, k)) << " at i = " << i << ", k = " << k;
    }
  }
}
TEST(Su3, Su3AlgebraAdditionOperator)
{
  std::ranlux24 rd;
  Su3Algebra a = random(rd);
  Su3Algebra b = random(rd);

  auto c = a + b;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t k = 0; k < 3; ++k) {
      EXPECT_EQ(c(i, k), a(i, k) + b(i, k)) << " at i = " << i << ", k = " << k;
    }
  }
}
TEST(Su3, Su3AlgebraMinusEqualOperator)
{
  std::ranlux24 rd;
  Su3Algebra a = random(rd);
  Su3Algebra b = random(rd);
  Su3Algebra c(a);
  c -= b;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t k = 0; k < 3; ++k) {
      EXPECT_EQ(c(i, k), a(i, k) - b(i, k)) << " at i = " << i << ", k = " << k;
    }
  }
}
TEST(Su3, Su3AlgebraSubtractionPartOne)
{
  Su3Algebra a(std::complex<double>(0.3, 0.8)),
    b(std::complex<double>(8.1, 8.4)), c;
  c = a - b;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t k = 0; k < 3; ++k) {
      auto tmp = a(i, k) - b(i, k);
      EXPECT_DOUBLE_EQ(real(c(i, k)), real(tmp))
        << " at i = " << i << ", k = " << k << "a(i,k)=" << a(i, k)
        << ", b(i,k)=" << b(i, k);
      EXPECT_DOUBLE_EQ(imag(c(i, k)), imag(tmp))
        << " at i = " << i << ", k = " << k << "a(i,k)=" << a(i, k)
        << ", b(i,k)=" << b(i, k);
    }
  }
}
TEST(Su3, Su3AlgebraSubtraction)
{
  std::ranlux24 rd;
  Su3Algebra a = random(rd);
  Su3Algebra b = random(rd);
  auto c = a - b;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t k = 0; k < 3; ++k) {
      auto tmp = a(i, k) - b(i, k);
      EXPECT_DOUBLE_EQ(real(c(i, k)), real(tmp))
        << " at i = " << i << ", k = " << k << "a(i,k)=" << a(i, k)
        << ", b(i,k)=" << b(i, k);
      EXPECT_DOUBLE_EQ(imag(c(i, k)), imag(tmp))
        << " at i = " << i << ", k = " << k << "a(i,k)=" << a(i, k)
        << ", b(i,k)=" << b(i, k);
    }
  }
}
TEST(Su3, determinantOfAThreeByThreeMatrix)
{
  ThreeByThreeMatrix a;
  a(0, 0) = { 5., 0. };
  a(0, 1) = { -1., 0. };
  a(0, 2) = { 9., 0 };
  a(1, 0) = { -1., 0. };
  a(1, 1) = { 6., 0. };
  a(1, 2) = { -1., 0. };
  a(2, 0) = { 9., 0. };
  a(2, 1) = { -1., 0 };
  a(2, 2) = { 7., 0 };
  double d = real(det(a));
  EXPECT_DOUBLE_EQ(d, -270);
}
TEST(Su3, randomAlgebraElementIsInAlgebra)
{
  std::ranlux24 rd;
  Su3Algebra alg = random(rd);
  std::complex<double> tr = trace(alg);
  EXPECT_NEAR(real(tr), 0.0, 1e-15);
  EXPECT_NEAR(imag(tr), 0.0, 1e-15);
  bool hermitian = isHermitian(alg);
  EXPECT_TRUE(hermitian);
  std::complex<double> tr2 = trace(alg * alg);
  std::cout << tr2 << std::endl;
}
TEST(Su3, Su3Multiplication)
{
  std::ranlux24 rd;
  Su3 a = static_cast<Su3>(exp(random(rd)));
  Su3 b = static_cast<Su3>(exp(random(rd)));
  auto c = a * b;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t k = 0; k < 3; ++k) {
      std::complex<double> tmp(0.0, 0.0);
      for (std::size_t j = 0; j < 3; ++j) {
        tmp += a(i, j) * b(j, k);
      }
      EXPECT_NEAR(std::real(c(i, k)), std::real(tmp), 1.e-14);
      EXPECT_NEAR(std::imag(c(i, k)), std::imag(tmp), 1.e-14);
    }
  }
}
TEST(Su3, Su3Trace)
{
  constexpr Su3 one(1.0);
  Su3 testmat(std::array<std::array<std::complex<double>, 3>, 3>{
    { std::array<std::complex<double>, 3>{ { 1.2, 3.3, 5.3 } },
      std::array<std::complex<double>, 3>{ { 2.6, 8.3, 1.9 } },
      std::array<std::complex<double>, 3>{ { 8.3, 1.3, 8.5 } } } });

  constexpr std::complex<double> res = 1.2 + 8.3 + 8.5;
  EXPECT_TRUE(trace(one) == std::complex<double>(3.0, 0.0));
  EXPECT_TRUE(trace(testmat) == res);
}
TEST(Su3, Su3GeneratorsObeyAlgebra)
{
  for (std::size_t a = 0; a < 8; ++a) {
    for (std::size_t b = 0; b < 8; ++b) {
      double delta_ab = (a == b ? 1.0 : 0.0);
      EXPECT_DOUBLE_EQ(
        std::real(trace(Su3Consts::lambda[a] * Su3Consts::lambda[b])),
        2.0 * delta_ab)
        << " at a == " << a << " and b == " << b;
      EXPECT_DOUBLE_EQ(
        std::imag(trace(Su3Consts::lambda[a] * Su3Consts::lambda[b])), 0.0)
        << " at a == " << a << " and b == " << b;
      ;
    }
  }
}
TEST(Su3, Su3GeneratorsAreTraceless)
{
  for (std::size_t a = 0; a < 8; ++a) {
    EXPECT_DOUBLE_EQ(std::real(trace(Su3Consts::lambda[a])), 0.0);
    EXPECT_DOUBLE_EQ(std::imag(trace(Su3Consts::lambda[a])), 0.0);
  }
}
