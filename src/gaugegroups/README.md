# Gaugegroups #

This folder contains specializations for various gauge groups.
In the following the typename will be abbreviated "T".

Each gauge group also must have an implementation of its algebra,
abbreviated A.

T and A do not need to derive from anything but should implement
at least:

T should implement:

  * a constructor that represents "one" in the sense of a
    group and can be called `T(1)`.
  * a member typedef pointing to the generating algebra, called `algebra`,
    i.e. `public typedef T::algebra = A`.
  * group multiplication: `constexpr operator*(T const &, T const &)`
  * complex scalar multiplication: `constexpr operator*(std::complex<double> const &, T const &)`
  * the group inverse of an element: `T inverse(T const &)`
  * (non-) equality: `constexpr operator==(T const &, T const &)`
    and `constexpr operator!=(T const &, T const &)`
  * `constexpr T dagger(T const & u)`

A should implement:
  
  * a constructor that represents the neutral element when called `A(0)`
  * addition and subtraction, `A operator+(A const &, A const &)` and
    `A operator-(A const &, A const &)`
  * a function
	```
	template< typename Generator >
	A random( Generator &gen );
	```
    that returns a random element from the algebra. 
    The template parameter `Generator` must be a Uniform Random Bit Generator
    http://en.cppreference.com/w/cpp/concept/UniformRandomBitGenerator
  * (non-) equality: `constexpr operator==(T const &, T const &)`
  * `constexpr std::complex<double> trace(T const &)`
