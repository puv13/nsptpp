/// \file
///
/// SU(3)
#pragma once

#include "complex_helpers.h"
#include <array>
#include <complex>
#include <functional>
#include <iostream>
#include <limits>
#include <random>

// common base of Su3 and Su3Algebra
class ThreeByThreeMatrix {
 protected:
  std::array<std::array<std::complex<double>, 3>, 3> data;

 public:
  constexpr ThreeByThreeMatrix(ThreeByThreeMatrix const &other) = default;
  constexpr ThreeByThreeMatrix(
    std::array<std::array<std::complex<double>, 3>, 3> const &arr) noexcept
    : data(arr)
  {
  }
  constexpr ThreeByThreeMatrix(std::complex<double> const &c) noexcept
    : ThreeByThreeMatrix(std::array<std::array<std::complex<double>, 3>, 3>{
        { std::array<std::complex<double>, 3>{ { c, 0, 0 } },
          std::array<std::complex<double>, 3>{ { 0, c, 0 } },
          std::array<std::complex<double>, 3>{ { 0, 0, c } } } })
  {
  }
  constexpr ThreeByThreeMatrix() noexcept : ThreeByThreeMatrix(0){};
  explicit constexpr ThreeByThreeMatrix(
    std::array<std::complex<double>, 9> const &arr)
    : ThreeByThreeMatrix(std::array<std::array<std::complex<double>, 3>, 3>{
        { std::array<std::complex<double>, 3>{ { arr[0], arr[1], arr[2] } },
          std::array<std::complex<double>, 3>{ { arr[3], arr[4], arr[5] } },
          std::array<std::complex<double>, 3>{ { arr[6], arr[7], arr[8] } } } })
  {
  }

  // access to the elements:
  std::complex<double> &operator()(std::size_t i, std::size_t j)
  {
    return data[i][j];
  }
  constexpr std::complex<double> const &operator()(std::size_t i,
                                                   std::size_t j) const
  {
    return data[i][j];
  }
  ThreeByThreeMatrix &operator*=(std::complex<double> const &a)
  {
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        data[i][j] *= a;
      }
    }
    return *this;
  }
  ThreeByThreeMatrix &operator-=(ThreeByThreeMatrix const &other)
  {
    for (std::size_t i = 0; i < 3; ++i)
      for (std::size_t j = 0; j < 3; ++j)
        data[i][j] -= other(i, j);
    return *this;
  }
  ThreeByThreeMatrix &operator+=(ThreeByThreeMatrix const &other)
  {
    for (std::size_t i = 0; i < 3; ++i)
      for (std::size_t j = 0; j < 3; ++j)
        data[i][j] += other(i, j);
    return *this;
  }
};
// non-member functions
constexpr ThreeByThreeMatrix operator+(ThreeByThreeMatrix const &a,
                                       ThreeByThreeMatrix const &b)
{
  return ThreeByThreeMatrix(std::array<std::array<std::complex<double>, 3>, 3>{
    { std::array<std::complex<double>, 3>{
        { a(0, 0) + b(0, 0), a(0, 1) + b(0, 1), a(0, 2) + b(0, 2) } },
      std::array<std::complex<double>, 3>{
        { a(1, 0) + b(1, 0), a(1, 1) + b(1, 1), a(1, 2) + b(1, 2) } },
      std::array<std::complex<double>, 3>{
        { a(2, 0) + b(2, 0), a(2, 1) + b(2, 1), a(2, 2) + b(2, 2) } } } });
}
constexpr ThreeByThreeMatrix operator-(ThreeByThreeMatrix const &a,
                                       ThreeByThreeMatrix const &b)
{
  return ThreeByThreeMatrix(std::array<std::array<std::complex<double>, 3>, 3>{
    { std::array<std::complex<double>, 3>{
        { a(0, 0) - b(0, 0), a(0, 1) - b(0, 1), a(0, 2) - b(0, 2) } },
      std::array<std::complex<double>, 3>{
        { a(1, 0) - b(1, 0), a(1, 1) - b(1, 1), a(1, 2) - b(1, 2) } },
      std::array<std::complex<double>, 3>{
        { a(2, 0) - b(2, 0), a(2, 1) - b(2, 1), a(2, 2) - b(2, 2) } } } });
}
constexpr bool operator==(ThreeByThreeMatrix const &a,
                          ThreeByThreeMatrix const &b)
{
  return (a(0, 0) == b(0, 0) and a(0, 1) == b(0, 1) and a(0, 2) == b(0, 2) and
          a(1, 0) == b(1, 0) and a(1, 1) == b(1, 1) and a(1, 2) == b(1, 2) and
          a(2, 0) == b(2, 0) and a(2, 1) == b(2, 1) and a(2, 2) == b(2, 2));
}
constexpr bool operator!=(ThreeByThreeMatrix const &a,
                          ThreeByThreeMatrix const &b)
{
  return not(a == b);
}
constexpr ThreeByThreeMatrix dagger(ThreeByThreeMatrix const &u)
{
  return ThreeByThreeMatrix(std::array<std::array<std::complex<double>, 3>, 3>{
    { std::array<std::complex<double>, 3>{
        { conj(u(0, 0)), conj(u(1, 0)), conj(u(2, 0)) } },
      std::array<std::complex<double>, 3>{
        { conj(u(0, 1)), conj(u(1, 1)), conj(u(2, 1)) } },
      std::array<std::complex<double>, 3>{
        { conj(u(0, 2)), conj(u(1, 2)), conj(u(2, 2)) } } } });
}
constexpr std::complex<double> trace(ThreeByThreeMatrix const &a)
{
  return a(0, 0) + a(1, 1) + a(2, 2);
}
constexpr ThreeByThreeMatrix operator*(ThreeByThreeMatrix const &a,
                                       std::complex<double> const &b)
{
  return ThreeByThreeMatrix(std::array<std::array<std::complex<double>, 3>, 3>{
    { std::array<std::complex<double>, 3>{
        { b * a(0, 0), b * a(0, 1), b * a(0, 2) } },
      std::array<std::complex<double>, 3>{
        { b * a(1, 0), b * a(1, 1), b * a(1, 2) } },
      std::array<std::complex<double>, 3>{
        { b * a(2, 0), b * a(2, 1), b * a(2, 2) } } } });
}
constexpr ThreeByThreeMatrix operator*(std::complex<double> const &a,
                                       ThreeByThreeMatrix const &b)
{
  return b * a;
}

constexpr ThreeByThreeMatrix operator*(ThreeByThreeMatrix const &a,
                                       ThreeByThreeMatrix const &b)
{
  return ThreeByThreeMatrix(std::array<std::array<std::complex<double>, 3>, 3>{
    { std::array<std::complex<double>, 3>{
        { a(0, 0) * b(0, 0) + a(0, 1) * b(1, 0) + a(0, 2) * b(2, 0),
          a(0, 0) * b(0, 1) + a(0, 1) * b(1, 1) + a(0, 2) * b(2, 1),
          a(0, 0) * b(0, 2) + a(0, 1) * b(1, 2) + a(0, 2) * b(2, 2) } },
      std::array<std::complex<double>, 3>{
        { a(1, 0) * b(0, 0) + a(1, 1) * b(1, 0) + a(1, 2) * b(2, 0),
          a(1, 0) * b(0, 1) + a(1, 1) * b(1, 1) + a(1, 2) * b(2, 1),
          a(1, 0) * b(0, 2) + a(1, 1) * b(1, 2) + a(1, 2) * b(2, 2) } },
      std::array<std::complex<double>, 3>{
        { a(2, 0) * b(0, 0) + a(2, 1) * b(1, 0) + a(2, 2) * b(2, 0),
          a(2, 0) * b(0, 1) + a(2, 1) * b(1, 1) + a(2, 2) * b(2, 1),
          a(2, 0) * b(0, 2) + a(2, 1) * b(1, 2) + a(2, 2) * b(2, 2) } } } });
}
// Algebra SU3:
class Su3Algebra : public ThreeByThreeMatrix {
 public:
  explicit constexpr Su3Algebra(ThreeByThreeMatrix const &a)
    : ThreeByThreeMatrix(a)
  {
  }
  constexpr Su3Algebra(
    std::array<std::array<std::complex<double>, 3>, 3> const &a)
    : ThreeByThreeMatrix(a)
  {
  }
  constexpr Su3Algebra(std::complex<double> const &a) : ThreeByThreeMatrix(a) {}
  constexpr Su3Algebra() : ThreeByThreeMatrix(0) {}

  static const std::size_t Nc = 3;
};
// enable additions and subtractions within the algebra:
constexpr Su3Algebra operator+(Su3Algebra const &a, Su3Algebra const &b)
{
  return static_cast<Su3Algebra>(static_cast<ThreeByThreeMatrix>(a) +
                                 static_cast<ThreeByThreeMatrix>(b));
}
constexpr Su3Algebra operator-(Su3Algebra const &a, Su3Algebra const &b)
{
  return static_cast<Su3Algebra>(static_cast<ThreeByThreeMatrix>(a) -
                                 static_cast<ThreeByThreeMatrix>(b));
}
constexpr Su3Algebra operator*(Su3Algebra const &a,
                               std::complex<double> const &b)
{
  return static_cast<Su3Algebra>(static_cast<ThreeByThreeMatrix>(a) * b);
}
constexpr Su3Algebra operator*(std::complex<double> const &a,
                               Su3Algebra const &b)
{
  return static_cast<Su3Algebra>(a * static_cast<ThreeByThreeMatrix>(b));
}

// GROUP SU3:
// we provide a naive implementation here, that can then be optimized.
class Su3 : public ThreeByThreeMatrix {
 public:
  explicit constexpr Su3(ThreeByThreeMatrix const &a) : ThreeByThreeMatrix(a) {}
  constexpr Su3(std::array<std::array<std::complex<double>, 3>, 3> const &a)
    : ThreeByThreeMatrix(a)
  {
  }
  constexpr Su3(std::complex<double> const &a) : ThreeByThreeMatrix(a) {}
  constexpr Su3() : ThreeByThreeMatrix(1) {}
  typedef Su3Algebra algebra;
};
// enable multiplications within the group:
constexpr Su3 operator*(Su3 const &a, Su3 const &b)
{
  return static_cast<Su3>(static_cast<ThreeByThreeMatrix>(a) *
                          static_cast<ThreeByThreeMatrix>(b));
}

// since also constexpr is missing for sqrt, we add 1./sqrt(3) as a magic
// number:
static constexpr double invsqrt3 = 0.57735026918962576451;

namespace Su3Consts {
using cmplx = std::complex<double>;
using arr = std::array<cmplx, 3>;
using mat = std::array<arr, 3>;
static constexpr Su3 one = Su3(1);
static constexpr Su3Algebra zero = Su3Algebra(0);
constexpr std::array<Su3Algebra, 8> lambda{
  { // lambda_1
    Su3Algebra(
      mat{ { arr{ { cmplx(0., 0.), cmplx(1., 0.), cmplx(0., 0.) } },
             arr{ { cmplx(1., 0.), cmplx(0., 0.), cmplx(0., 0.) } },
             arr{ { cmplx(0., 0.), cmplx(0., 0.), cmplx(0., 0.) } } } }),
    // lambda_2
    Su3Algebra(
      mat{ { arr{ { cmplx(0., 0.), cmplx(0., -1.), cmplx(0., 0.) } },
             arr{ { cmplx(0., 1.), cmplx(0., 0.), cmplx(0., 0.) } },
             arr{ { cmplx(0., 0.), cmplx(0., 0.), cmplx(0., 0.) } } } }),
    // lambda_3
    Su3Algebra(
      mat{ { arr{ { cmplx(1., 0.), cmplx(0., 0.), cmplx(0., 0.) } },
             arr{ { cmplx(0., 0.), cmplx(-1., 0.), cmplx(0., 0.) } },
             arr{ { cmplx(0., 0.), cmplx(0., 0.), cmplx(0., 0.) } } } }),
    // lambda_4
    Su3Algebra(
      mat{ { arr{ { cmplx(0., 0.), cmplx(0., 0.), cmplx(1., 0.) } },
             arr{ { cmplx(0., 0.), cmplx(0., 0.), cmplx(0., 0.) } },
             arr{ { cmplx(1., 0.), cmplx(0., 0.), cmplx(0., 0.) } } } }),
    // lambda_5
    Su3Algebra(
      mat{ { arr{ { cmplx(0., 0.), cmplx(0., 0.), cmplx(0., -1.) } },
             arr{ { cmplx(0., 0.), cmplx(0., 0.), cmplx(0., 0.) } },
             arr{ { cmplx(0., 1.), cmplx(0., 0.), cmplx(0., 0.) } } } }),
    // lambda_6
    Su3Algebra(
      mat{ { arr{ { cmplx(0., 0.), cmplx(0., 0.), cmplx(0., 0.) } },
             arr{ { cmplx(0., 0.), cmplx(0., 0.), cmplx(1., 0.) } },
             arr{ { cmplx(0., 0.), cmplx(1., 0.), cmplx(0., 0.) } } } }),
    // lambda_7
    Su3Algebra(
      mat{ { arr{ { cmplx(0., 0.), cmplx(0., 0.), cmplx(0., 0.) } },
             arr{ { cmplx(0., 0.), cmplx(0., 0.), cmplx(0., -1.) } },
             arr{ { cmplx(0., 0.), cmplx(0., 1.), cmplx(0., 0.) } } } }),
    // lambda_8
    Su3Algebra(
      mat{ { arr{ { cmplx(invsqrt3, 0.), cmplx(0., 0.), cmplx(0., 0.) } },
             arr{ { cmplx(0., 0.), cmplx(invsqrt3, 0.), cmplx(0., 0.) } },
             arr{ { cmplx(0., 0.), cmplx(0., 0.),
                    cmplx(-2. * invsqrt3, 0.) } } } }) }
};
}
template <typename Rng>
inline Su3Algebra random(Rng &rng)
{
  Su3Algebra tmp;
  std::normal_distribution<double> dist(0.0, 1.0);
  // std::normal_distribution generates two values at once and caches one of
  // them. If the cache is empty we get a compiler warning. To avoid this call
  // std::normal_distribution once before we actually need it.
  // See also:
  // https://stackoverflow.com/questions/16008271/why-are-c11-random-distributions-mutable
  dist(rng);

  for (std::size_t i = 0; i < 8; ++i) {
    // the original fortran code samples with a variance of 2, but takes the T_a
    // = 0.5 * lambda_a.
    // we instead take the lambda_a and sample with a variance of 1:
    tmp += dist(rng) * Su3Consts::lambda[i];
  }
  return tmp;
}
inline ThreeByThreeMatrix exp(ThreeByThreeMatrix const &a)
{
  // computes exp(A) = \sum_{i=0} A^i / (n!)
  ThreeByThreeMatrix m(1.0), res(1.0);
  double residual = 1.0;
  double fac = 1.0;
  std::size_t iter = 1;
  while (residual > std::numeric_limits<double>::epsilon()) {
    fac /= static_cast<double>(iter);
    m = m * a;
    residual = fac * std::abs(trace(m)) / 3.0;

    res += fac * m;
    iter++;
  }
  return res;
}
inline std::complex<double> det(ThreeByThreeMatrix const &a)
{
  // computes det(A)
  return a(0, 0) * a(1, 1) * a(2, 2) + a(0, 1) * a(1, 2) * a(2, 0) +
         a(1, 0) * a(2, 1) * a(0, 2) -
         (a(2, 0) * a(1, 1) * a(0, 2) + a(1, 0) * a(0, 1) * a(2, 2) +
          a(2, 1) * a(1, 2) * a(0, 0));
}
inline std::complex<double> det(std::array<std::complex<double>, 4> const &a)
{
  // interprets the 4 numbers as 2x2 matrix and computes the determinant.
  // needed later
  return a[0] * a[3] - a[1] * a[2];
}
inline ThreeByThreeMatrix inverse(ThreeByThreeMatrix const &a)
{
  // computes A^{-1}
  const std::complex<double> detA = det(a);
  if (std::abs(detA) < std::numeric_limits<double>::epsilon()) {
    throw std::runtime_error("inverse: matrix is singular.");
  }
  const std::complex<double> oneOverDet = 1.0 / detA;
  return oneOverDet * ThreeByThreeMatrix(std::array<std::complex<double>, 9>{
                        det({ a(1, 1), a(1, 2), a(2, 1), a(2, 2) }),
                        det({ a(0, 2), a(0, 1), a(2, 2), a(2, 1) }),
                        det({ a(0, 1), a(0, 2), a(1, 1), a(1, 2) }),
                        det({ a(1, 2), a(1, 0), a(2, 2), a(2, 0) }),
                        det({ a(0, 0), a(0, 2), a(2, 0), a(2, 2) }),
                        det({ a(0, 2), a(0, 0), a(1, 2), a(1, 0) }),
                        det({ a(1, 0), a(1, 1), a(2, 0), a(2, 1) }),
                        det({ a(0, 1), a(0, 0), a(2, 1), a(2, 0) }),
                        det({ a(0, 0), a(0, 1), a(1, 0), a(1, 1) }) });
}
inline bool isHermitian(ThreeByThreeMatrix const &a)
{
  if (a(0, 1) == std::conj(a(1, 0)) and a(0, 2) == std::conj(a(2, 0)) and
      a(1, 2) == std::conj(a(2, 1)))
    return true;
  else
    return false;
}
inline void sqrtIteration(ThreeByThreeMatrix &y, ThreeByThreeMatrix &z,
                          ThreeByThreeMatrix &tmp)
{
  tmp = y;
  y = 0.5 * (y + inverse(z));
  z = 0.5 * (z + inverse(tmp));
}
inline ThreeByThreeMatrix sqrt(ThreeByThreeMatrix const &a)
{
  // no checks are done here..
  // for the algorithm, see Denman-Beavers iteration,
  // https://en.wikipedia.org/wiki/Square_root_of_a_matrix
  ThreeByThreeMatrix y(a), z(1.0), tmp;
  double residual = 1.0;
  std::size_t iteration = 0;
  constexpr std::size_t MAXITER = 100000;
  constexpr double TARGETRESIDUAL = 1e-10;

  while (residual > TARGETRESIDUAL and iteration < MAXITER) {
    sqrtIteration(y, z, tmp);
    iteration += 1;
    residual = std::abs(1.0 - det(a * z * z));
  }
  if (iteration >= MAXITER) {
    throw std::runtime_error("sqrt did not converge.");
  }
  // do a couple of extra iterations so that we get a better result (bad
  // residual definition)
  for (std::size_t i = 0u; i < 5; ++i) {
    sqrtIteration(y, z, tmp);
  }
  return y;
}
inline ThreeByThreeMatrix logMercator(ThreeByThreeMatrix const &a)
{
  // this does the log via the simple mercator series,
  // logA = log(1+K) = K - K^2/2 + K^3/3 - K^4/4 + ...
  // this works well if A is close to 1.
  ThreeByThreeMatrix k = a - ThreeByThreeMatrix(1.);
  double residual = 1.;
  ThreeByThreeMatrix kP(k), result(k), tmp;
  double iter(2.0), sign(-1.0);
  while (residual > std::numeric_limits<double>::epsilon() and iter < 1e5) {
    kP = kP * k;
    tmp = (sign / iter) * kP;
    result = result + tmp;
    residual = std::abs(trace(tmp)) / 3.0; // probably a bad estimator
    sign *= -1.0;
    iter += 1.0;
  }
  if (residual > std::numeric_limits<double>::epsilon()) {
    throw std::runtime_error("log did not converge.");
  }
  return result;
}
inline ThreeByThreeMatrix log(ThreeByThreeMatrix const &a)
{
  // this does some "preconditioning" for the log:
  // it uses the identity
  // log(A) = 2^k log(A^(1/(2^k)))
  // A^(1/(2^k)) is, for large enough k, close to 1 and
  // the log can be computed via the logMercator
  constexpr double MAXRES = 1e-2;
  ThreeByThreeMatrix s(sqrt(a));
  double residual = std::abs(3.0 - trace(s));
  std::size_t iter = 1;
  while (residual > MAXRES) {
    s = sqrt(s);
    iter += 1;
    residual = std::abs(3.0 - trace(s));
  }
  return static_cast<double>(std::pow(2, iter)) * logMercator(s);
}

constexpr Su3 dagger(Su3 const &u)
{
  return static_cast<Su3>(dagger(static_cast<ThreeByThreeMatrix>(u)));
}

// should be hermitian, but...
constexpr Su3Algebra dagger(Su3Algebra const &u)
{
  return static_cast<Su3Algebra>(dagger(static_cast<ThreeByThreeMatrix>(u)));
}
