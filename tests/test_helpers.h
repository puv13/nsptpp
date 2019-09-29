#pragma once
#include <array>
#include <complex>
#include <random>

template <std::size_t len>
std::array<std::complex<double>, len> randomVector()
{
  static std::ranlux24 rd;
  std::normal_distribution<> d(0, 2);

  std::array<std::complex<double>, len> e;
  for (auto &el : e) {
    el = std::complex<double>(d(rd), d(rd));
  }
  return e;
}
