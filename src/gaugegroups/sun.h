/// \file
///
/// SU(N) implementation (following the implementation of SU(3) in su3.h
/// Needs a CBLAS compatible blas library!

#pragma once

#include "../src/nummat.h"
#include "complex_helpers.h"
#include <random>

/// \brief SU(N) Algebra
template <std::size_t N>
class SuNAlgebra : public nummat<N> {
 public:
  constexpr SuNAlgebra() : nummat<N>(0.0) {}
  explicit constexpr SuNAlgebra(nummat<N> const &M) : nummat<N>(M) {}
  constexpr SuNAlgebra(std::complex<double> const &c) : nummat<N>(c) {}
  constexpr SuNAlgebra(std::array<std::complex<double>, N * N> const &ar)
    : nummat<N>(ar)
  {
  }
};

/// SU(N) Algebra Algebra

/// SU(N) Matrix Multiplication
template <std::size_t N>
constexpr SuNAlgebra<N> operator*(SuNAlgebra<N> const &left,
                                  SuNAlgebra<N> const &right)
{

  return SuNAlgebra<N>(static_cast<nummat<N>>(left) *
                       static_cast<nummat<N>>(right));
}

/// SU(N) Matrix Subtraction
template <std::size_t N>
constexpr SuNAlgebra<N> operator-(SuNAlgebra<N> const &left,
                                  SuNAlgebra<N> const &right)
{

  return SuNAlgebra<N>(static_cast<nummat<N>>(left) -
                       static_cast<nummat<N>>(right));
}

/// SU(N) Matrix Addition
template <std::size_t N>
constexpr SuNAlgebra<N> operator+(SuNAlgebra<N> const &left,
                                  SuNAlgebra<N> const &right)
{

  return SuNAlgebra<N>(static_cast<nummat<N>>(left) +
                       static_cast<nummat<N>>(right));
}


///  \brief SU(N) Group
template <std::size_t N>
class SuN : public nummat<N> {

 public:
  constexpr SuN() : nummat<N>(0.0) {}
  explicit constexpr SuN(nummat<N> const &M) : nummat<N>(M) {}
  constexpr SuN(std::complex<double> const &c) : nummat<N>(c) {}
  constexpr SuN(std::array<std::complex<double>, N * N> const &ar)
    : nummat<N>(ar)
  {
  }
};

/// SU(N) Group Algebra

/// SU(N) Matrix Multiplication
template <std::size_t N>
constexpr SuN<N> operator*(SuN<N> const &left, SuN<N> const &right)
{

  return SuN<N>(static_cast<nummat<N>>(left) * static_cast<nummat<N>>(right));
}

/// SU(N) Matrix Subtraction
template <std::size_t N>
constexpr SuN<N> operator-(SuN<N> const &left, SuN<N> const &right)
{

  return SuN<N>(static_cast<nummat<N>>(left) - static_cast<nummat<N>>(right));
}

/// SU(N) Matrix Addition
template <std::size_t N>
constexpr SuN<N> operator+(SuN<N> const &left, SuN<N> const &right)
{

  return SuN<N>(static_cast<nummat<N>>(left) + static_cast<nummat<N>>(right));
}


template <std::size_t N>
constexpr std::array<SuNAlgebra<N>, N * N - 1> Generators()
{
  static_assert(N > 1, "SU(N) is trivial for N<2.");

  std::array<SuNAlgebra<N>, N * N - 1> gen;
  gen.fill(SuNAlgebra<N>(0.0));

  auto num_sym = (N * (N - 1)) / (2LU);
  auto num_asm = (N * (N - 1)) / (2LU);


  for (auto k = 1LU; k < N; ++k) {
    for (auto j = 0LU; j < k; ++j) {
      auto s = k + j - 1;
      auto a = s + num_sym;
      // std::cout << "k:" << k << "\tj:" << j << "\ti:" << i << std::endl;

      // N(N-1)/2 symmetric generators
      gen[s](j, k) = 0.5;
      gen[s](k, j) = 0.5;

      // N(N-1)/2 anti-symmetric generators
      gen[a](j, k) = std::complex<double>(0.0, -0.5);
      gen[a](k, j) = std::complex<double>(0.0, 0.5);
    }
  }

  // (N-1) diagonal generators

  for (auto l = 0LU; l < N - 1; ++l) {
    auto d = l + num_sym + num_asm;

    auto rt = 0.5 * sqrt(2. / (static_cast<double>((l + 2) * (l + 1))));

    for (auto j = 0LU; j < l + 1; ++j) {
      gen[d](j, j) = rt;
    }
    gen[d](l + 1, l + 1) = -1.0 * rt * static_cast<double>(l + 1);
  }

  return gen;
}


/// This generates Gaussian Langevin noise for the SUN group.
/// Note that the noise itself is not group valued!
/// NormalGenerator has to have an overloaded operator()
/// that returns a normal (0,1) distributed double.
/// By default the noise has a standard deviation of sqrt(2) (variance of 2).
template <std::size_t N, typename NormalGenerator>
SuNAlgebra<N> noiseSUN(NormalGenerator &rng, double stddev = std::sqrt(2.))
{

  SuNAlgebra<N> res;
  static auto gens = Generators<N>();

  for (auto G : gens) {
    double noise = rng() * stddev;
    res += G * noise;
  }

  return res;
}


/// This generates a random SUN group element.
/// NormalGenerator has to have an overloaded operator()
/// that returns a normal (0,1) distributed double.
template <std::size_t N, typename Generator>
SuN<N> randomSUN(Generator &rng)
{

  nummat<N> res;
  static auto gens = Generators<N>();
  std::normal_distribution<double> ndist(0.0, 1.0);

  for (auto G : gens) {
    double noise = ndist(rng);
    res += G * noise;
  }

  res *= std::complex<double>(0., 1.);

  return SuN<N>(exp(res));
}
