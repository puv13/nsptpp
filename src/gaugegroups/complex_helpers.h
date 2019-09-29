/// \file
///
/// Fhis file fixes the missing constexpr qualifiers in <complex>
/// see also http://open-std.org/JTC1/SC22/WG21/docs/papers/2016/p0415r0.html


#pragma once

#include <complex>

constexpr std::complex<double> operator*(std::complex<double> const &a,
                                         std::complex<double> const &b)
{
  return std::complex<double>(a.real() * b.real() - a.imag() * b.imag(),
                              a.real() * b.imag() + a.imag() * b.real());
}
constexpr std::complex<double> operator+(std::complex<double> const &a,
                                         std::complex<double> const &b)
{
  return std::complex<double>(a.real() + b.real(), a.imag() + b.imag());
}
constexpr std::complex<double> operator-(std::complex<double> const &a,
                                         std::complex<double> const &b)
{
  return std::complex<double>((a.real() - b.real()), (a.imag() - b.imag()));
}
constexpr std::complex<double> conj(std::complex<double> const &a)
{
  return std::complex<double>(a.real(), -a.imag());
}
