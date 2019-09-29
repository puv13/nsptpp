/// \file
///
/// Implement U(1) gauge field with all the features required by README.md


#pragma once
#include "complex_helpers.h"
#include <complex>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <random>
//#include "../latticefields/cp.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846 // Double precision pi
#endif

using cmplx = std::complex<double>;
using size_t = std::size_t;
// const auto II  =  std::complex<double>(0.,1.);


class U1 {
 private:
  /// This should be a complex number with norm 1
  cmplx data;

 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************

  U1() : data(1.) {}
  /// \todo Automatic enforcement of group property in gauge group
  /// constructions ?
  U1(const cmplx &init, const bool &normalise = false) : data(init)
  {
    if (normalise) {
      this->normalise();
    }
  }


  // **********************************************************************
  // Destructor
  // **********************************************************************
  ~U1() = default;


  // **********************************************************************
  // Member functions
  // **********************************************************************

  double norm() const { return std::abs(data); }

  cmplx value() const
  {
    auto res = data;
    return res;
  }

  U1 conj() const
  {
    auto res = *this;
    res.data = std::conj(res.data);
    return res;
  }

  void normalise()
  {

    cmplx nrm = this->norm();
    this->data /= nrm;
  }

  double phase() const { return std::arg(data); }


  /// Assignment
  U1 &operator=(U1 const &orig)
  {
    this->data = orig.data;
    return *this;
  }

  // Arithmetic

  /// Minus assignment operator
  U1 &operator-=(U1 const &rhs)
  {
    this->data -= rhs.data;

    return *this;
  }

  /// Minus
  U1 operator-(U1 const &rhs)
  {
    auto res = *this;

    res -= rhs;

    return res;
  }

  /// Plus assignment operator
  U1 &operator+=(U1 const &rhs)
  {
    this->data += rhs.data;

    return *this;
  }

  /// Plus
  U1 operator+(U1 const &rhs)
  {
    auto res = *this;

    res += rhs;

    return res;
  }

  /// Multiplication assignment operator
  U1 &operator*=(U1 const &rhs)
  {
    this->data *= rhs.data;

    return *this;
  }

  /// Multiplication
  U1 operator*(U1 const &rhs) const
  {
    auto res = *this;

    res *= rhs;

    return res;
  }

  /// Multiplication with scalar
  U1 operator*(std::complex<double> const &rhs)
  {
    auto res = *this;
    res.data *= rhs;
    return res;
  }

  /// Division assignment operator
  U1 &operator/=(U1 const &rhs)
  {
    this->data /= rhs.data;

    return *this;
  }

  /// Division
  U1 operator/(U1 const &rhs)
  {
    auto res = *this;

    res /= rhs;

    return res;
  }


  /// Equality
  bool operator==(U1 const &rhs) const { return (this->data == rhs.data); }

  // **********************************************************************
  // Friend functions
  // **********************************************************************

  friend cmplx log(U1 const &orig) { return std::log(orig.data); }

  friend cmplx exp(U1 const &orig) { return std::exp(orig.data); }

  friend cmplx real(U1 const &orig) { return std::real(orig.data); }

  friend cmplx imag(U1 const &orig) { return std::imag(orig.data); }

  friend U1 inverse(U1 const &orig) { return 1. / orig.data; }

  friend U1 dagger(U1 const &orig)
  {
    auto res = orig.conj();
    return res;
  }

  friend std::complex<double> trace(U1 const &orig) { return orig.value(); }

  friend U1 operator*(U1 const &lhs, std::complex<double> const &rhs)
  {
    auto res = lhs;
    return res * rhs;
  }

  friend U1 operator*(std::complex<double> const &lhs, U1 const &rhs)
  {
    return rhs * lhs;
  }

  /// Simple printing of U1 field
  friend std::ostream &operator<<(std::ostream &stream, U1 const &u1)
  {
    stream.precision(6);
    auto out = u1.value();
    stream << std::scientific << "(" << std::setw(13) << std::real(out) << ","
           << std::setw(13) << std::imag(out) << ")";
    return stream;
  }


  template <typename Generator>
  friend U1 update(U1 const &orig, Generator &gen, double const &eps = 5.e-2)
  {
    auto II = std::complex<double>(0., 1.);
    std::uniform_real_distribution<double> dist(-eps, eps);

    U1 newU1(orig);

    double y = dist(gen);

    auto factor = (1 + II * y) / std::sqrt(1 + y * y);

    newU1.data *= factor;

    // newU1.normalise();

    return newU1;
  }
};

// \todo How should U1 random(Distribution &d, Generator &gen) be implemented?
template <typename Distribution, typename Generator>
U1 randomU1(Distribution &d, Generator &gen)
{

  double randnum = d(gen);

  return U1(std::polar(1., randnum));
}

/// Returns a group element of U1
/// with uniform distribution in U1.
template <typename Generator>
U1 randomU1(Generator &gen)
{

  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  double phase = M_PI * dist(gen);

  return U1(std::polar(1., phase));
}


/// This generates Gaussian Langevin noise for the U1 group.
/// Note that the noise itself is not group valued!
/// NormalGenerator has to have an overloaded operator()
/// that returns a normal (0,1) distributed double.
/// By default the noise has a standard deviation of sqrt(2) (variance of 2).
template <typename NormalGenerator>
double noiseU1(NormalGenerator &gen, double stddev = std::sqrt(2.))
{

  double res = stddev * gen();

  return res;
}


/// Returns a (uniformly sampled) random U1 field (
inline U1 randomU1()
{

  // Init ranlux with true random number (Well, at least we try. The C++
  // standard does not guaranty std::random_device() gives a true random
  // number)
  std::random_device rd;
  std::seed_seq sseq({ rd(), rd(), rd() });
  // std::seed_seq sseq ({1,2,3,4});
  std::ranlux48 generator(sseq);

  // std::ostream_iterator<unsigned> out (std::cout," ");
  // sseq.param(out); std::cout << "\n" << generator << "\n" << std::endl;
  // Uniformly pick a phase from [-pi,pi)
  // Note that this is consistent with the std::arg function from std::complex
  std::uniform_real_distribution<double> dist(-1.0, 1.0);
  double phase = M_PI * dist(generator);

  return U1(std::polar(1., phase));
}
