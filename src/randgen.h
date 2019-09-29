/// \file
/// Definitions for a simple random class to be used to conveniently pass a
/// random number generator and  distribution to a function


#pragma once

#include <array>
#include <iostream>
#include <random>
#include <sstream>
#include <string>
#include <type_traits>
#include <vector>

using std::string;

/// Base class for random number generation
template <typename Generator, typename Distribution,
          std::size_t InitLength = 312>
class rnd_base {
 private:
  Generator gen;
  Distribution dist;

 public:
  /// Standard constructor
  rnd_base() { hardware_init(); }

  /// SeedSeq constructor
  explicit rnd_base(std::seed_seq &seedseq) { gen = Generator{ seedseq }; }


  /// Init Random number generator with real random numbers.
  ///
  /// (At least we try. The standard does not enforce random_device to
  /// be a true random number generator.)
  void hardware_init()
  {
    std::random_device rd;
    std::uniform_int_distribution<size_t> udint;
    std::vector<std::size_t> seeds;
    for (auto i = 0ul; i < InitLength; ++i) {
      seeds.push_back(udint(rd));
    }
    std::seed_seq sseq(seeds.begin(), seeds.end());

    gen = Generator{ sseq };
  }

  /// Returns a random number sampled from the distribution
  /// 'dist' using the generator 'gen'.
  auto operator()() { return dist(gen); }

  /// Returns a string encoding the internal state of the GENERATOR
  string get_generator_state() const
  {

    std::stringstream state;
    state << gen;

    return state.str();
  }

  /// Returns a string encoding the internal state of the DISTRIBUTION
  string get_distribution_state() const
  {

    std::stringstream state;
    state << dist;

    return state.str();
  }

  /// Returns an array of strings encoding the internal state of
  /// the generator and the distribution.
  auto get_state() const
  {

    std::array<std::string, 2> statear;

    statear[0] = this->get_generator_state();
    statear[1] = this->get_distribution_state();

    return statear;
  }


  /// Set the internal state of the GENERATOR with a string
  void set_generator_state(string statestr)
  {

    std::stringstream state;
    state << statestr;

    state >> gen;
  }

  /// Set the internal state of the DISTRIBUTION with a string
  void set_distribution_state(string statestr)
  {

    std::stringstream state;
    state << statestr;

    state >> dist;
  }

  /// Set the internal state of the generator and the distribution
  /// with an array of strings.
  void set_state(std::array<std::string, 2> const &statesar)
  {

    std::stringstream diststate;
    std::stringstream genstate;

    genstate << statesar[0];
    genstate >> gen;

    diststate << statesar[1];
    diststate >> dist;
  }

  /// Equality operator
  template <typename G2, typename D2>
  bool operator==(rnd_base<G2, D2, InitLength> const &other) const
  {

    bool test = (this->get_state() == other.get_state());

    test = (test && std::is_same<Generator, G2>::value &&
            std::is_same<Distribution, D2>::value);

    return test;
  }
};


// #############################################################################
// Specialisations for often used cases
// #############################################################################


// -----------------------------------------------------------------------------
// RANLUX
// -----------------------------------------------------------------------------

/// This corresponds to ranlux with a luxury level of 3 (p=223)
/// of the original Fortran code from F. James
/// Computer Physics Communications, 79 (1994) 111–114
template <typename Distribution>
using ranlux24 = rnd_base<std::ranlux24, Distribution, 48LU>;


/// Standard normal distribution (mean=0,stddev=1) based on the ranlux generator
template <typename output_type>
using ranlux24_normal = ranlux24<std::normal_distribution<output_type>>;

/// Uniform distribution in [0,1) based on the ranlux generator
template <typename output_type>
using ranlux24_uniform = ranlux24<std::uniform_real_distribution<output_type>>;


// -----------------------------------------------------------------------------
// MERSENNE TWISTER
// -----------------------------------------------------------------------------

/// Mersenne Twister
template <typename Distribution>
using mt19937 =
  rnd_base<std::mt19937_64, Distribution, std::mt19937_64::state_size>;

/// Standard normal distribution (mean=0,stddev=1) based on the Mersenne
/// generator
template <typename output_type>
using mt19937_normal = mt19937<std::normal_distribution<output_type>>;


/// Uniform distribution in [0,1) based on the Mersenne generator
template <typename output_type>
using mt19937_uniform = mt19937<std::uniform_real_distribution<output_type>>;
