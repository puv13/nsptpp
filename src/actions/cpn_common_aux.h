/// \file
///
/// \brief Stuff needed by both CP(N-1) action implementations


#pragma once

#include "../src/gaugegroups/sun.h"
#include "../src/lattice.h"
#include "../src/latticefields/cp.h"
#include "action.h"
#include "expansion.h"
#include "nummat.h"
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <limits>
#include <numeric>

// Global constants known at compile time
// ----------------------------------------------------------------------
#ifdef SET_N_
constexpr std::size_t N = SET_N_;
#else
constexpr std::size_t N = 3;
#endif

// We need pi a lot here, so ...
#ifndef M_PI
#define M_PI 3.14159265358979323846 // Double precision pi
#endif

// Constant for the Langevin noise generation
constexpr double CP_Langevin_noise_scale = std::sqrt(2.);


// Type names
// ----------------------------------------------------------------------
using size_t = std::size_t;
using cplx = std::complex<double>;
using CPN = CP<cplx, N>;


/// Computes the Hermitian matrix
/// scalar*left*dagger(right)+right*dagger(scalar*left) and adds it to M.
/// Warning! This only fills the upper part of the Hermitian matrix. M should be
/// a compressed Hermitian matrix in CblasUpper form!
template <size_t N>
void external_product(std::complex<double> const &scalar,
                      CP<std::complex<double>, N> const &left,
                      CP<std::complex<double>, N> const &right, nummat<N> &M)
{

  cblas_zher2(CBLAS_ORDER::CblasRowMajor, CBLAS_UPLO::CblasUpper, N,
              reinterpret_cast<const double *>(&scalar),
              reinterpret_cast<const double *>(&left[0]), 1,
              reinterpret_cast<const double *>(&right[0]), 1,
              reinterpret_cast<double *>(M.data()), N);
}

/// operator * for CP fields
template <size_t N>
auto operator*(CP<std::complex<double>, N> const &l,
               CP<std::complex<double>, N> const &r)
{
  return scalar_prod(l, r);
}


namespace CommonCPNExpHelpers {
}


/// += operator for complex double expansion and double expansion
template <size_t order>
auto operator+=(Expansion<std::complex<double>, order> &one,
                Expansion<double, order> const &other)
  -> Expansion<std::complex<double>, order>
{

  for (size_t i = 0; i < order; ++i) {
    one[i] += other[i];
  }
  return one;
}


/// Scalar Product for expansions
template <typename T, size_t order>
auto scalar_prod(Expansion<T, order> const &left,
                 Expansion<T, order> const &right)
  -> Expansion<decltype(scalar_prod(left[0], right[0])), order>
{

  Expansion<decltype(scalar_prod(left[0], right[0])), order> res;

  for (auto n = 0LU; n < order; ++n) {
    for (auto i = 0LU; i <= n; ++i) {
      res[n] += scalar_prod(left[n - i], right[i]);
    }
  }

  return res;
}


template <size_t N, size_t order>
Expansion<nummat<N>, order> operator+(Expansion<nummat<N>, order> L,
                                      nummat<N> const &R)
{

  if (order > 1) {
    L[1] = L[1] + R;
  }

  return L;
}


template <size_t N, size_t order>
auto operator*(Expansion<nummat<N>, order> const &left,
               Expansion<CPN, order> const &right) -> Expansion<CPN, order>
{
  Expansion<CPN, order> res;

  for (size_t n = 0; n < order; ++n) {
    for (size_t i = 0; i <= n; ++i) {
      res[n] += left[n - i] * right[i];
    }
  }

  return res;
}


template <size_t N, size_t order>
auto operator*(Expansion<std::complex<double>, order> const &left,
               Expansion<nummat<N>, order> const &right)
  -> Expansion<nummat<N>, order>
{
  Expansion<nummat<N>, order> res;

  for (size_t n = 0; n < order; ++n) {
    for (size_t i = 0; i <= n; ++i) {
      res[n] += left[n - i] * right[i];
    }
  }
  return res;
}


template <size_t N, size_t order>
auto operator*(Expansion<nummat<N>, order> const &left,
               Expansion<std::complex<double>, order> const &right)
  -> Expansion<nummat<N>, order>
{
  return right * left;
}


template <size_t N, size_t order>
auto operator*(Expansion<std::complex<double>, order> const &left,
               Expansion<CP<cplx, N>, order> const &right)
  -> Expansion<CP<cplx, N>, order>
{
  Expansion<CP<cplx, N>, order> res;

  for (size_t n = 0; n < order; ++n) {
    for (size_t i = 0; i <= n; ++i) {
      res[n] += left[n - i] * right[i];
    }
  }

  return res;
}


template <size_t N, size_t order>
auto operator*(Expansion<std::complex<double>, order> const &left,
               nummat<N> const &right) -> Expansion<nummat<N>, order>
{
  Expansion<nummat<N>, order> res;

  for (size_t n = 0; n < order; ++n) {
    res[n] = left[n] * right;
  }

  return res;
}
