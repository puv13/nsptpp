/// \file
///
/// \brief Implementation of the CP(N-1) field


#pragma once

#include "../gaugegroups/complex_helpers.h"
#include "../gaugegroups/u1.h"
#include "nummat.h"
#include <array>
#include <complex>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

using size_t = std::size_t;
template <typename T, size_t N>
class CP;
template <typename T, size_t N>
class CPtwistPhase;


/// Scalar Product for CP
template <typename T, size_t N>
T scalar_prod(CP<T, N> const &left, CP<T, N> const &right)
{
  auto res = T(0.);

  for (size_t i = 0; i < N; ++i) {
    res += std::conj(left[i]) * right[i];
  }
  return res;
}


/// Scalar Product specialisation for complex<double>
template <size_t N>
std::complex<double> scalar_prod(CP<std::complex<double>, N> const &left,
                                 CP<std::complex<double>, N> const &right)
{
  // Debug stuff
  // std::cout << "COMPLEX SPECIALISATION OF SCALAR PROD" << std::endl;
  std::complex<double> res = 0.;

  cblas_zdotc_sub(N, reinterpret_cast<double const *>(&left[0]), 1,
                  reinterpret_cast<double const *>(&right[0]), 1,
                  reinterpret_cast<__complex__ double *>(&res));
  return res;
}


/// Matrix-Vector multiplication
template <size_t N>
constexpr CP<std::complex<double>, N>
operator*(nummat<N> const &lhs, CP<std::complex<double>, N> const &rhs)
{

  CP<std::complex<double>, N> res(0.0);
  constexpr std::complex<double> alpha(1.0), beta(0.0);

  cblas_zgemv(CBLAS_ORDER::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans, N, N,
              reinterpret_cast<const double *>(&alpha),
              reinterpret_cast<const double *>(lhs.data()), N,
              reinterpret_cast<const double *>(&rhs[0]), 1,
              reinterpret_cast<const double *>(&beta),
              reinterpret_cast<double *>(&res[0]), 1);

  return res;
}


/// Generate a (uniformly sampled) random CP field
/// \details
/// The basic idea is to uniformly pick a vector on the 2N-unit-sphere
/// For this we need 2N Gaussian random numbers.
/// A normalised vector of this 2N random numbers lies on the 2N-unit-sphere
/// and constructing several such vectors samples the sphere uniformly.
/// (See  mathworld.wolfram.com/HyperspherePointPicking.html)
/// We then take the entries of the 2N-dim vector as real and imaginary parts
/// of a N-dim complex vector.
template <size_t N>
inline CP<std::complex<double>, N> randomCP()
{

  std::random_device rd;
  std::seed_seq sseq({ rd(), rd(), rd() });
  std::ranlux48 generator(sseq);

  std::normal_distribution<double> ndist(0.0, 1.0);


  CP<std::complex<double>, N> res;

  for (auto &el : res) {
    el = std::complex<double>{ ndist(generator), ndist(generator) };
  }

  res.normalise();
  return res;
}


/// \brief Generate a (uniformly sampled) random CP field using generator gen
/// For details see the description of the function randomCP()
template <size_t N, typename RNG>
inline CP<std::complex<double>, N> randomCP(RNG &gen)
{

  std::normal_distribution<double> ndist(0.0, 1.0);

  CP<std::complex<double>, N> res;

  for (auto &el : res) {
    el = std::complex<double>{ ndist(gen), ndist(gen) };
  }

  res.normalise();
  return res;
}


/// Iterator for CP class
///
template <typename T>
class CPiterator : public std::iterator<std::bidirectional_iterator_tag, T> {

 private:
  size_t pos;
  T *parent;

 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************

  explicit CPiterator(size_t pos_, T &parent_) : pos(pos_), parent(&parent_) {}

  explicit CPiterator(T &parent_) : pos(0ul), parent(&parent_) {}

  // **********************************************************************
  // Deconstructor
  // **********************************************************************
  ~CPiterator() = default;

  // **********************************************************************
  // Member functions
  // **********************************************************************
  constexpr size_t index() const { return pos; }

  /// Prefix increment
  CPiterator &operator++()
  {
    ++pos;
    return *this;
  }

  /// Prefix decrement
  CPiterator &operator--()
  {
    --pos;
    return *this;
  }

  /// Postfix increment
  CPiterator operator++(int)
  {
    auto copy = *this;
    ++pos;
    return copy;
  }

  /// Postfix decrement
  CPiterator operator--(int)
  {
    auto copy = *this;
    --pos;
    return copy;
  }

  bool operator==(CPiterator const &other) const
  {
    return (pos == other.pos and parent == other.parent);
  }

  bool operator!=(CPiterator const &other) const
  {
    return not((*this) == other);
  }

  bool operator<(CPiterator const &other) const
  {
    if (parent != other.parent) {
      throw std::runtime_error("Operator < only makes sense for same parent");
    }
    return (pos < other.pos);
  }


  CPiterator operator[](int offset) const
  {
    return CPiterator(pos + offset, *parent);
  }

  auto operator*() const -> decltype(parent->operator[](pos))
  {
    return parent->operator[](pos);
  }
};


/// \brief Class for CP(N-1) Fields
template <typename T, size_t N>
class CP {

 private:
  /// Dimension. Keep in mind that N is the same as in CP(N-1)
  /// i.e., CP(2) <-> N=3 !!!
  std::array<T, N> data;

 public:
  // **********************************************************************
  // Typedefs
  // **********************************************************************
  typedef CPiterator<const CP<T, N>> const_iterator;
  typedef CPiterator<CP<T, N>> iterator;

  // **********************************************************************
  // Constructors
  // **********************************************************************

  /// Standard Ctor
  CP() { data.fill(T(0)); }

  /// \brief Init Ctor
  /// \details Set all field components to the given value
  explicit CP(T const &init, const bool &normalise = false)
  {
    data.fill(init);
    if (normalise) {
      this->normalise();
    }
  }

  explicit CP(std::array<T, N> const &init_ar, const bool &normalise = false)
  {
    for (size_t i = 0; i < N; ++i) {
      this->data[i] = init_ar[i];
      if (normalise) {
        this->normalise();
      }
    }
  }

  /// Copy Ctor
  // CP(const CP<T,N>  &orig) {
  //     data = orig.data;
  // }


  // **********************************************************************
  // Member Functions
  // **********************************************************************

  /// Access elements of cp field
  T &operator[](std::size_t const idx) { return data[idx]; }

  /// Access elemenst of constant cp field
  constexpr T const &operator[](const std::size_t idx) const
  {
    return data[idx];
  }

  /// \brief Compute norm of cp field
  /// \details  Use C++ magic to infer return type. Note that for
  /// std::complex<T> the return type of std::abs is T, i.e.,
  /// std::abs returns a double if invoked for a complex<double>.
  /// \todo Exception if data is empty
  constexpr auto norm() const -> decltype(std::abs(data.front()))
  {
    auto res = std::abs(data[0]);
    res = 0;
    for (const auto &x : data) {
      auto abs = std::abs(x);
      res += abs * abs;
    }
    return std::sqrt(res);
  }

  void normalise()
  {
    auto nrm = this->norm();
    if (nrm < std::abs(T(1.)) * 1e-16) {
      throw std::runtime_error("Norm is numerically zero. Cannot normalise.");
    }

    for (auto idx = data.begin(); idx != data.end(); ++idx) {
      *idx /= nrm;
    }
  }

  constexpr static size_t size() { return N; }

  // **********************************************************************
  // Operators
  // **********************************************************************

  // Assignment
  CP<T, N> &operator=(const CP<T, N> &orig)
  {
    this->data = orig.data;
    return *this;
  }

  // Arithmetic

  /// Minus assignment operator
  CP<T, N> &operator-=(CP<T, N> const &rhs)
  {
    for (size_t i = 0; i < N; ++i) {
      data[i] -= rhs[i];
    }

    return *this;
  }

  /// Minus
  CP<T, N> operator-(CP<T, N> const &rhs)
  {
    auto res = *this;

    res -= rhs;

    return res;
  }

  /// Equality
  bool operator==(CP<T, N> const &rhs) const
  {

    for (auto i = 0ul; i < N; i++) {

      auto eq(this->data[i] == rhs[i]);

      if (not eq) {
        return false;
      }
    }
    return true;
  }

  /// Plus assignment operator
  CP<T, N> &operator+=(CP<T, N> const &rhs)
  {
    for (size_t i = 0; i < N; ++i) {
      data[i] += rhs[i];
    }

    return *this;
  }

  /// Multiplication assignment (scalar multiplication)
  CP<T, N> &operator*=(T const &rhs)
  {
    for (size_t i = 0; i < N; ++i) {
      data[i] *= rhs;
    }
    return *this;
  }

  /// Multiplication with scalar
  CP<T, N> operator*(T const &rhs) const
  {
    auto res = *this;
    for (size_t i = 0; i < N; ++i) {
      res[i] *= rhs;
    }
    return res;
  }

  /// Multiplication with U1
  CP<T, N> operator*(U1 const &rhs) const
  {
    auto res = *this;
    for (size_t i = 0; i < N; ++i) {
      res[i] *= rhs.value();
    }
    return res;
  }

  /// Multiplication with CPtwistPhase
  // CP<T,N>  operator*( CPtwistPhase<T,N> const& rhs) const
  // {
  //     auto res = *this;
  //     for (size_t i=0; i<N; ++i)
  //     {
  // 	res[i] *= rhs[i];
  //     }
  //     return res;
  // }

  // **********************************************************************
  // Friend functions
  // **********************************************************************

  /// Multiplication with U1 from the left
  friend CP<T, N> operator*(U1 const &lhs, CP<T, N> const &rhs)
  {
    auto res = rhs;
    return res * lhs;
  }


  // friend CP<T,N> operator *(CP<T,N> const & lhs, double  const &rhs);
  friend CP<T, N> operator*(std::complex<double> lhs, CP<T, N> const &rhs)
  {
    auto res = rhs;
    for (size_t i = 0; i < N; ++i) {
      res.data[i] *= lhs;
    }

    return res;
  }

  /// scalar product
  friend T scalar_prod<T>(CP<T, N> const &left, CP<T, N> const &right);

  /// MC update
  template <typename Generator>
  friend CP<T, N> update(CP<T, N> const &orig, Generator &gen,
                         double const &eps = 5.e-2)
  {

    std::uniform_real_distribution<double> dist(0, eps);

    // Get random CP(N-1) field
    CP<T, N> diff = randomCP<N>(gen);
    // Rescale by random epsilon
    diff *= dist(gen);

    // add original field
    diff += orig;
    // In general we will have to normalise the new field
    diff.normalise();

    return diff;
  }

  // **********************************************************************
  // Iterators over CP field components
  // **********************************************************************
  iterator begin() { return iterator(0, *this); }
  iterator end() { return iterator(data.size(), *this); }
  /// Element access with boundary check via iterator
  iterator at(const size_t &idx)
  {
    if (idx >= this->data.size()) {
      throw std::runtime_error("Index out of range in CPiterator");
    }
    return iterator(idx, *this);
  }


  const_iterator begin() const { return const_iterator(0, *this); }
  const_iterator end() const
  {
    return const_iterator(data.size(), static_cast<CP<T, N> const &>(*this));
  }
  /// Element access with boundary check via const iterator
  const_iterator at(const size_t &idx) const
  {
    if (idx >= this->data.size()) {
      throw std::runtime_error("Index out of range in const CPiterator ");
    }
    return const_iterator(idx, *this);
  }

  // **********************************************************************
  // Friends
  // **********************************************************************

  /// Simple printing of CP(n-1) field
  friend std::ostream &operator<<(std::ostream &stream, CP<T, N> const &cp)
  {

    stream.precision(6);
    stream << std::scientific << "[ ";

    for (auto i = 0ul; i < cp.size() - 1; ++i) {
      stream << "(" << std::setw(13) << std::real(cp[i]) << "," << std::setw(13)
             << std::imag(cp[i]) << "), ";
    }
    stream << "(" << std::setw(13) << std::real(cp[cp.size() - 1]) << ","
           << std::setw(13) << std::imag(cp[cp.size() - 1]) << ")"
           << " ]";
    return stream;
  }
};


/// \brief Helper class for twisted CP(N-1) boundary conditions
template <typename T, size_t N>
class CPtwistPhase {

 private:
  /// Dimension. Keep in mind that N is the same as in CP(N-1)
  /// i.e., CP(2) <-> N=3 !!!
  std::array<T, N> phases;

 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************

  ///  \brief Standard Ctor
  ///  \details Set all phases to 1
  explicit CPtwistPhase() { phases.fill(T(1.0)); }

  /// \brief Init Ctor
  /// \details Set all phases to the given value
  explicit CPtwistPhase(T const &init) { phases.fill(init); }

  /// \brief Init Ctor
  /// \details Set all phases according to given values in array
  explicit CPtwistPhase(std::array<T, N> const &init_ar)
  {
    for (size_t i = 0; i < N; ++i) {
      this->phases[i] = init_ar[i];
    }
  }


  // **********************************************************************
  // Member Functions
  // **********************************************************************


  /// Access phases of CPtwistPhase
  T &operator[](std::size_t const idx) { return phases[idx]; }

  /// Access phases of constant CPtwistPhase
  constexpr T const &operator[](const std::size_t idx) const
  {
    return phases[idx];
  }

  /// Get size
  constexpr static size_t size() { return N; }

  /// Multiplication with CP(N-1) field
  CP<T, N> operator*(CP<T, N> const &field) const
  {

    auto res = field;
    for (size_t i = 0; i < N; ++i) {
      res[i] *= this->phases[i];
    }
    return res;
  }

  /// Multiplication with U1
  CPtwistPhase<T, N> operator*(U1 const &rhs) const
  {
    auto res = *this;
    for (size_t i = 0; i < N; ++i) {
      res[i] *= rhs.value();
    }
    return res;
  }

  /// Multiplication with a scalar
  CPtwistPhase<T, N> operator*(T const &rhs) const
  {
    auto res = *this;
    for (size_t i = 0; i < N; ++i) {
      res[i] *= rhs;
    }
    return res;
  }


  // **********************************************************************
  // Friend functions
  // **********************************************************************

  /// Multiplication with U1 from the left
  friend CPtwistPhase<T, N> operator*(U1 const &lhs,
                                      CPtwistPhase<T, N> const &rhs)
  {
    auto res = rhs;
    return res * lhs;
  }

  /// Multiply CP(N-1) field from left
  friend CP<T, N> operator*(CP<T, N> const &field,
                            CPtwistPhase<T, N> const &rhs)
  {
    auto res = field;
    for (size_t i = 0; i < N; ++i) {
      res[i] *= std::conj(rhs[i]);
    }
    return res;
  }


  /// Multiplication with complex<double> from the left
  friend CPtwistPhase<T, N> operator*(std::complex<double> lhs,
                                      CPtwistPhase<T, N> const &rhs)
  {
    auto res = rhs;
    for (size_t i = 0; i < N; ++i) {
      res.phases[i] *= lhs;
    }

    return res;
  }

  /// Exp defined for compatibility with current implementation of
  /// boundary conditions.
  // friend CPtwistPhase<T,N> exp(CPtwistPhase<T,N> const &orig){
  //     auto res = orig;

  //     for (size_t i=0; i<N; ++i)
  //     {
  // 	res[i] = std::exp(res[i]);
  //     }

  //     return res;

  // }
};
