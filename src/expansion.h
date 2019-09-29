/// \file
///
/// Series expansion of lattice fields
#pragma once

#include <array>
#include <complex>
#include <exception>
#include <iostream>
#include <string>
#include <tuple>
#include <typeindex>

/// this allows a contexpr array initialization with the order (size) being a
/// template parameter.
template <typename Tuple, typename VTuple = std::remove_reference_t<Tuple>,
          std::size_t... Indices>

constexpr std::array<
  std::common_type_t<std::tuple_element_t<Indices, VTuple>...>,
  sizeof...(Indices)>
to_array(Tuple &&tuple, std::index_sequence<Indices...>)
{
  return { { std::get<Indices>(std::forward<Tuple>(tuple))... } };
}

template <typename Tuple, typename VTuple = std::remove_reference_t<Tuple>>
constexpr decltype(auto) to_array(Tuple &&tuple)
{
  return to_array(std::forward<Tuple>(tuple),
                  std::make_index_sequence<std::tuple_size<VTuple>::value>{});
}

template <typename T, class Function, std::size_t... Indices>
constexpr std::array<T, sizeof...(Indices)>
make_array_helper(Function f, std::index_sequence<Indices...>)
{
  return { { f(Indices)... } };
}

template <typename T, std::size_t N, class Function>
constexpr std::array<T, N> make_array(Function f)
{
  return make_array_helper<T>(f, std::make_index_sequence<N>{});
}

template <typename T>
struct arr_firstOrderSetter {
  const T value;
  constexpr arr_firstOrderSetter(T val) : value(val){};
  constexpr T operator()(std::size_t const &elem)
  {
    return (elem == 0 ? value : T(0));
  }
};


namespace CHECK {
struct DoesNotExist {
};

template <typename A, typename B>
DoesNotExist operator+(A const &a, B const &b);

template <typename A, typename B>
struct canBeAdded {
 private:
  // TODO: make this work!!
  static constexpr typename std::remove_reference<A>::type *a =
    static_cast<typename std::remove_reference<A>::type *>(nullptr);
  static constexpr typename std::remove_reference<B>::type *b =
    static_cast<typename std::remove_reference<B>::type *>(nullptr);

 public:
  static constexpr bool value =
    not std::is_same<decltype(*a + *b), DoesNotExist>::value;
};

} // namespace CHECK

struct ExpansionError : std::runtime_error {
  ExpansionError(std::string const &msg) : std::runtime_error(msg) {}
};

/// \brief Expansion class
///
/// A template for expansions with coefficients of type T.
template <typename T, std::size_t order>
class Expansion {
 public:
  std::array<T, order> data;
  typedef T scalartype;

  /// Init Ctor from array
  constexpr Expansion(std::array<T, order> const &arr) : data{ arr }
  {
    static_assert(sizeof(Expansion<T, order>) == order * sizeof(T),
                  "Cannot construct an Expansion of types that hold "
                  "more than just a collection of doubles."
                  "Expansion and members do not have the correct sizes. "
                  "Maybe member has additional fields?");
  }

  /// Standard Ctor
  constexpr Expansion()
    : Expansion(make_array<T, order>(arr_firstOrderSetter<T>(T(0))))
  {
  }


  /// Iterator to first expansion order
  auto begin() { return data.begin(); };

  /// Iterator to expansion order zero
  auto const begin() const { return data.begin(); };

  /// Iterator to expansion order zero
  auto end() { return data.end(); };

  /// Iterator to highest expansion order
  auto const end() const { return data.end(); };

  /*static_assert( order > 0,
                   "cannot construct an expansion without elements." );
    for( auto &elem : data )
        elem = T( 0 );
  */
  /// Copy Ctor
  constexpr Expansion(Expansion<T, order> const &other) = default;


  /// Copy Ctor from expansion of different type
  template <typename T2>
  Expansion(Expansion<T2, order> const &other) : Expansion()
  {
    for (std::size_t i = 0; i < order; ++i) {
      data[i] = static_cast<T>(other[i]);
    }
  }

  // template< typename T2 >
  //  constexpr Expansion( T2 const &elem ) :
  //  data((elem*(Expansion<T,order>::one<T,order>)).data) {}

  /// Init Ctor from single element
  template <typename T2>
  constexpr Expansion(T2 const &elem)
    : data(make_array<T, order>(arr_firstOrderSetter<T>(static_cast<T>(elem))))
  {
  }

  /// Ctor. for (multiplicative) identity expansion
  template <typename T2, std::size_t order2>
  static constexpr auto one =
    Expansion<T2, order2>(make_array<T2, order2>(arr_firstOrderSetter<T2>(1)));

  /// Ctor. for zero expansion
  template <typename T2, std::size_t order2>
  static constexpr auto zero = Expansion<T2, order2>();

  /// Get expansion order
  static constexpr std::size_t numberOfOrders() { return order; }

  /// Element (expansion coefficient) access
  constexpr T const &operator[](std::size_t index) const
  {
    return data[static_cast<std::size_t>(index)];
  }

  /// Element (expansion coefficient) access
  T &operator[](std::size_t index)
  {
    return data[static_cast<std::size_t>(index)];
  }

  /// Check for equality
  bool operator==(Expansion<T, order> const &other) const
  {
    auto it_this = this->data.begin();
    auto it_other = other.data.begin();
    for (; it_this != this->data.end() and it_other != other.data.end();
         ++it_this, ++it_other) {
      if (*it_this != *it_other)
        return false;
    }
    return true;
  }


  /// Expansions of different base type are not equal
  template <typename T2, std::size_t order2>
  bool operator==(Expansion<T2, order2> const &) const
  {
    return false;
  }

  /// Check if two expansion are different
  bool operator!=(Expansion<T, order> const &other) const
  {
    return not((*this) == other);
  }


  /// -= operator for Expansions
  auto operator-=(Expansion<T, order> const &other)
    -> Expansion<decltype(data[0] - other[0]), order>
  {
    auto it_this = this->data.begin();
    auto it_other = other.data.cbegin();
    for (; it_this != this->data.end() and it_other != other.data.cend();
         ++it_this, ++it_other)
      *it_this -= *it_other;
    return *this;
  }

  /// *= operator for multiplication of 'scalar' and expansion
  auto operator*=(T const &other) -> Expansion<decltype(data[0] * other), order>
  {
    auto it_this = this->data.begin();
    for (; it_this != this->data.end(); ++it_this) {
      *it_this *= other;
    }
    return *this;
  }

  /// Calculate the memory requirements of the Expansion in multiples of
  /// the size of a double
  static constexpr std::size_t numberOfDoubles()
  {
    return order * sizeof(T) / sizeof(double);
  }


  /// \brief Serialise an expansion
  ///
  /// Represent the expansion with an array of doubles
  std::array<double, numberOfDoubles()> serialize() const
  {
    std::array<double, numberOfDoubles()> res;
    auto it = res.begin();
    for (std::size_t i = 0; i < order; ++i) {
      for (std::size_t j = 0; j < numberOfDoubles() / order; ++j) {
        *it = reinterpret_cast<double const *>(&(data[i]))[j];
        ++it;
      }
    }
    return res;
  }

  /// Cast expansion to an expansion of lower order
  /// Higher order coefficients will be lost!
  template <std::size_t lower>
  Expansion<T, lower> cast_lower() const
  {
    // Make sure this can only  be used intentionally
    static_assert(lower < order,
                  "Can only cast to an expansion of lower order");

    Expansion<T, lower> res;
    for (std::size_t i = 0; i < lower; ++i) {
      res[i] = this->data[i];
    }

    return res;
  }


  /// Cast expansion to an expansion of higher order
  template <std::size_t higher>
  Expansion<T, higher> cast_higher() const
  {
    // Make sure this can only  be used intentionally
    static_assert(higher > order,
                  "Can only cast to an expansion of higher order");

    Expansion<T, higher> res;
    for (std::size_t i = 0; i < order; ++i) {
      res[i] = this->data[i];
    }

    return res;
  }
}; // Class Expansion

/// += operator for expansions
template <typename T, std::size_t order>
auto operator+=(Expansion<T, order> &one, Expansion<T, order> const &other) ->
  typename std::enable_if<CHECK::canBeAdded<T, T>::value,
                          Expansion<T, order> &>::type
{
  /*
    auto it_this  = this->data.begin();
    auto it_other = other.data.cbegin();
    for( ; it_this != this->data.end() and it_other != other.data.cend();
    ++it_this, ++it_other )
    *it_this += *it_other;
    return *this;
  */
  for (std::size_t i = 0; i < order; ++i) {
    one[i] += other[i];
  }
  return one;
}

/// -= operator for expansions
template <typename T, std::size_t order>
auto operator-=(Expansion<T, order> &one, Expansion<T, order> const &other) ->
  typename std::enable_if<CHECK::canBeAdded<T, T>::value,
                          Expansion<T, order> &>::type
{
  for (std::size_t i = 0; i < order; ++i) {
    one[i] -= other[i];
  }
  return one;
}

/// '+' operator for expansions with different coefficient types
template <typename T1, typename T2, std::size_t order>
auto operator+(Expansion<T1, order> const &one,
               Expansion<T2, order> const &other) ->
  typename std::enable_if<CHECK::canBeAdded<T1, T2>::value,
                          Expansion<decltype(one[0] + other[0]), order>>::type
{
  Expansion<decltype(one[0] + other[0]), order> res;
  for (std::size_t i = 0; i < order; ++i) {
    res[i] = one[i] + other[i];
  }
  return res;
}


/// \brief Deserialise and array of doubles
///
/// Construct expansion with coefficients of type T from
/// an array of doubles.
template <typename T, std::size_t order>
Expansion<T, order>
deserialize(std::array<double, order * sizeof(T) / sizeof(double)> arr)
{
  Expansion<T, order> res;
  constexpr std::size_t membersize = res.numberOfDoubles() / order;
  for (std::size_t i = 0; i < order; ++i) {
    for (std::size_t j = 0; j < membersize; ++j) {
      reinterpret_cast<double *>(&(res[i]))[j] = arr[i * membersize + j];
    }
  }
  return res;
}

/// Utility struct to check if something is an expansion
template <typename T>
struct is_expansion {
  static constexpr bool value = false;
  static constexpr std::size_t ord = 0;
};

/// Utility struct to check if something is an expansion
template <typename T, std::size_t order>
struct is_expansion<Expansion<T, order>> {
  static constexpr bool value = true;
  static constexpr std::size_t ord = order;
};

/// Multiplication of expansion with 'scalar' from the right.
///
/// Wraps the member function.
template <typename T, std::size_t order>
auto operator*(Expansion<T, order> left, T const &right)
  -> Expansion<decltype(left[0] * right), order>
{
  left *= right;
  return left;
}

/// Multiplication of expansion with a 'scalar' from the left.
///
/// This impl requires that T2 is not an expansion, but an elementary type that
/// can multiply the elements of the expansion.., therefore, we check if T2 is
/// NOT an
/// expansion. SFINAE ensures that this impl only exists if T2 fulfills this
/// requirement.
///
/// \todo: THIS DOES NOT WORK IF "left*right[0]" cannot be cast to decltype(T).
template <typename T2, typename T, std::size_t order>
auto operator*(T2 const &left, Expansion<T, order> right)
  -> std::enable_if_t<not is_expansion<T2>::value,
                      Expansion<decltype(left * right[0]), order>>
{
  // This does not work in general!
  // for( auto it = right.data.begin(); it != right.data.end(); ++it )
  //     *it = left * ( *it );
  // return right;
  Expansion<decltype(left * right[0]), order> res;

  for (auto i = 0LU; i < order; ++i) {
    res[i] = left * right[i];
  }

  return res;
}

/// Multiplication of expansion with a 'scalar' from the right.
///
/// Check the 'scalar' multiplication from the left for details.
template <typename T2, typename T, std::size_t order>
auto operator*(Expansion<T, order> left, T2 const &right)
  -> std::enable_if_t<not is_expansion<T2>::value,
                      Expansion<decltype(left[0] * right), order>>
{

  Expansion<decltype(left[0] * right), order> res;

  for (auto i = 0LU; i < order; ++i) {
    res[i] = left[i] * right;
  }

  return res;
}

/// Cast expansion to expansion with coefficients of different type
template <typename Tto, typename Tfrom, std::size_t order>
Expansion<Tto, order> expansion_cast(Expansion<Tfrom, order> const &from)
{
  Expansion<Tto, order> res;
  for (std::size_t i = 0; i < order; ++i) {
    res[i] = static_cast<Tto>(from[i]);
  }
  return res;
}

/// Multiplication of two expansions
///
/// T*T is not necessarily of type T (the scalar product of vectors is a
/// counterexample)
/// Is there a nicer alternative to writing decltype(T()*T())?
template <typename T, std::size_t order>
Expansion<decltype(T() * T()), order> operator*(Expansion<T, order> const &a,
                                                Expansion<T, order> const &b)
{

  Expansion<decltype(T() * T()), order> res;
  for (std::size_t n = 0; n < order; ++n) {
    for (std::size_t i = 0; i <= n; ++i) {
      res[n] += a[n - i] * b[i];
    }
  }
  return res;
}

/// \brief Inverse for simple floating point types.
///
/// Is needed for inverse(Expansion).
/// Must be implemented elsewhere for more complicated types.
template <typename T>
inline std::enable_if_t<std::is_floating_point<T>::value, T> inverse(T const &a)
{
  if (a == 0.0)
    throw ExpansionError("Cannot invert element.");
  return 1. / a;
}

/// Inverse for complex double.
inline std::complex<double> inverse(std::complex<double> const &a)
{
  if (real(a) == 0.0 and imag(a) == 0.0)
    throw ExpansionError("Cannot invert element.");
  return 1. / a;
}

/// Sum of expansions
template <typename T, std::size_t order>
Expansion<T, order> operator+(Expansion<T, order> left,
                              Expansion<T, order> const &right)
{
  left += right;
  return left;
}

/// Difference of expansions
template <typename T, std::size_t order>
Expansion<T, order> operator-(Expansion<T, order> left,
                              Expansion<T, order> const &right)
{
  left -= right;
  return left;
}

/// \brief Utility function to compute the inverse of an expansion.
///
/// The inverse of the zeroth order coefficient has to be supplied.
/// (Only the inverse of the zeroth order coefficient is needed to calculate
/// the inverse of the expansion!)
template <typename T, std::size_t order>
Expansion<T, order> inverse(Expansion<T, order> const &a,
                            T const &zerothOrderInverse)
{
  // start with the zeroth order inverse:
  Expansion<T, order> res(zerothOrderInverse);
  T tmp;
  for (std::size_t n = 1; n < order; ++n) {
    tmp = T(0);
    for (std::size_t i = 0; i < n; ++i)
      tmp += a[n - i] * res[i];
    res[n] = -res[0] * tmp;
  }
  return res;
}

/// \brief Inverse of an expansion
///
/// See inverse(Expansion<T, order> const &a,T const &z) for details.
template <typename T, std::size_t order>
Expansion<T, order> inverse(Expansion<T, order> const &a)
{
  return inverse(a, inverse(a[0]));
}

namespace {

/// \brief Utility function to compute expansion of the exponential function
///
/// This gives the ith coefficient c_i of
/// exp(x) = \sum_{i=0} c_i x^i
/// where c_i = 1 / (i!).
/// Uses template magic. See also 'constexpr std::array< double, 1 >
/// expcoefficients()'
template <std::size_t order>
constexpr std::array<double, order> expcoefficients()
{
  return std::array<double, order>{ to_array(std::tuple_cat(
    expcoefficients<order - 1>(),
    std::array<double, 1>{
      { const_cast<double &>(static_cast<std::array<double, order - 1> const &>(
          expcoefficients<order - 1>())[order - 2]) /
        static_cast<double>(order - 1) } })) };
}

/// \brief Utility function for exponential of expansion
template <>
constexpr std::array<double, 1> expcoefficients()
{
  return std::array<double, 1>{ { 1.0 } };
}

/// \brief Utility function to compute expansion of the natural log of (1+x)
///
/// This should give the ith coefficient c_i of
/// log(1+x) = \sum_{i=0} c_i x^i
/// with c_i = \partial^i log(1+x)|x=0  / (i!)
/// derivative part only:
/// log(1+x) ->  0
/// 1/(1+x) = (1+x)^{-1}  -> 1
/// 1/(1+x)^2 * (-1) = -(1+x)^{-2} -> -1
/// 2/(1+x)^3 = 2*(1+x)^{-3} -> 2
/// -6/(1+x)^4 = -6*(1+x)^{-4} -> -6
/// complete with /i! --> c_i = (-1)^(i+1) * (i-1)! / (i!) = (-1)^(i+1)/i
/// for i > 0
constexpr double logcoefficient(std::size_t const i)
{
  return i == 0 ? 0.0 : (i % 2 == 0 ? -1. / static_cast<double>(i) :
                                    /*i odd*/ 1. / static_cast<double>(i));
}

/// \brief Utility function for log of expansion
template <std::size_t order>
constexpr std::array<double, order> logcoefficients()
{
  return make_array<double, order>(logcoefficient);
}

/// \brief Utility function to insert expansion in a power series
template <typename T, std::size_t order>
/* constexpr ?? */ Expansion<T, order>
taylor(std::array<double, order> const &coefficients,
       Expansion<T, order> const &expansion)
{
  /// res = \sum_i c_i A^i
  /// powers = (A^i)^(j) = 1 + A + A^2 + A^3 + ...
  /// i.e. powers contains:
  /// powers[0] = 1, 0, 0, ...
  /// powers[1] = A^(0), A^(1), A^(2), ...
  /// powers[2] = A^(0)*A^(0), A^(1)*A^(0)+A^(0)*A^(1), ...
  /// ...
  std::array<Expansion<T, order>, order> powers;
  powers[0] = Expansion<T, order>{ T(1) };
  for (std::size_t i = 1; i < order; ++i) {
    powers[i] = powers[i - 1] * expansion;
  }
  /// resum, i.e. sort by powers of the final expansion parameter:
  /// res[0] = c_0 + c_1*A^(0) + c_2*A(0)*A(0) + ...
  /// res[1] =       c_1*A^(1) + c_2*(A^(1)*A^(0)+A^(0)*A^(1)) + ...
  /// res[2] =                   c_2*(A^(2)*
  Expansion<T, order> res;
  for (std::size_t j = 0; j < order; ++j) {
    for (std::size_t i = j; i < order; ++i) {
      res[i] += coefficients[j] * powers[j][i];
    }
  }
  return res;
}
} // anonymous namespace

/// \brief Exponential function of an expansion
///
/// \warning This function makes use of the identity \f$e^{a+b} = e^a
///  \cdot e^b\f$.  The reason is that we calculate the exponential
///  for the zeroth order separately. This is done to ensure that the
///  exponential is again *exact* up to the given order \f$m\f$. Let
///  \f$X\f$ be an expansion with expansion parameter \f$z\f$. If it
///  starts at zeroth order \f$(z^0)\f$, any term \f$X^n\f$ in the
///  series expansion of the exponential will also start with order
///  zero.  Consider for example the Taylor expansion of \f$e^S\f$
///  where \f$S = 1 + \epsilon x + \epsilon^2 y\f$.  The Taylor series
///  of which will have infinitely many terms of order 0 (in
///  \f$\epsilon\f$).  (This is important when your expansion is not
///  around the vacuum but, for example, around an instanton
///  background.)  (Convince yourself that this is true!) Therefore,
///  infinitely many terms in the exponential series would have to be
///  taken into account.  If, however, \f$X\f$ starts at order
///  \f$z^k\f$ with \f$k>0\f$ the power \f$X^n\f$ will be (at least)
///  of order \f$z^n\f$ in the expansion parameter. Therefore, all
///  terms after \f$X^{(m)}\f$ in the exponential series can be
///  ignored at order \f$z^{(m)}\f$.  Be extremely careful when doing
///  this for non-commuting types, where in general \f$e^{a+b} \neq
///  e^a \cdot e^b\f$.
template <typename T, std::size_t order>
Expansion<T, order> exp(Expansion<T, order> const &x)
{
  auto tmp = x;
  tmp[0] = T(0); // x -> x - x[0]
  return exp(x[0]) * taylor(expcoefficients<order>(), tmp);
}


///  \brief Logarithm of an expansion
///
///  To calculate the logarithm the series expansion of \f$\ln(1+x)\f$ is
///  used. Therefore, we first have to bring the expansion in a form that is
///  suitable for this series expansion.
///  \warning This function makes use of the identity \f$\ln(A\cdot
///  B)=\ln(A)+\ln(B)\f$. Be
///  careful if the expansion type is non-commuting!
///  See also  exp(Expansion< T, order > const &x ).
template <typename T, std::size_t order>
Expansion<T, order> log(Expansion<T, order> const &x)
{
  /*
  std::cout << "log("; for( auto const & elem : x ) std::cout << elem << " ";
  std::cout << ") = ";
  auto logZerothOrder = log(x[0]);
  std::cout << ") = " << logZerothOrder << "+log(";
  auto l = inverse(x[0]) * x;
  for( auto const & elem : l ) std::cout << elem << " ";
  std::cout << ") = ";
  std::cout << ") = " << logZerothOrder << "+(";
  l = inverse(x[0]) * x;
  l[0] = T(0);
  l = taylor(logcoefficients<order>(), l);
  for( auto const & elem : l ) std::cout << elem << " ";
  l = Expansion<T, order>(logZerothOrder)+l;
  std::cout << " = ";
  for( auto const & elem : l ) std::cout << elem << " ";
  std::cout << std::endl;
  */

  // x -> (x-x0)/x0
  auto tmp = x;
  T zerothOrderInverse = inverse(x[0]);
  tmp[0] = T(0);                  // x -> x - x[0] \sim x[0] = 0
  tmp = zerothOrderInverse * tmp; // x -> x / x[0]
  return Expansion<T, order>(log(x[0])) + taylor(logcoefficients<order>(), tmp);
}

/// \brief Compute the Hermitian conjugate of the expansion.
/// The function dagger() has to be defined for the expansion type.
template <typename T, std::size_t order>
Expansion<T, order> dagger(Expansion<T, order> const &x)
{
  auto res = x;
  for (auto &i : res)
    i = dagger(i);
  return res;
}

/// \brief Compute the trace of an expansion.
/// Trace has to be defined for the expansion type.
template <typename T, std::size_t order>
decltype(auto) trace(Expansion<T, order> const &x)
{
  Expansion<decltype(trace(x[0])), order> res;
  for (std::size_t i = 0; i < order; ++i) {
    res[i] = trace(x[i]);
  }
  return res;
}

/// \brief Calculate real part of the expansion
template <typename T, std::size_t order>
decltype(auto) real(Expansion<T, order> const &x)
{
  using namespace std;
  Expansion<decltype(real(x[0])), order> res;
  for (std::size_t i = 0; i < order; ++i) {
    res[i] = real(x[i]);
  }
  return res;
}

/// \brief Calculate imaginary  part of the expansion
template <typename T, std::size_t order>
decltype(auto) imag(Expansion<T, order> const &x)
{
  using namespace std;
  Expansion<decltype(imag(x[0])), order> res;
  for (std::size_t i = 0; i < order; ++i) {
    res[i] = imag(x[i]);
  }
  return res;
}

/// \brief Complex conjugate of the expansion
template <typename T, std::size_t order>
decltype(auto) conj(Expansion<T, order> const &x)
{
  using namespace std;
  Expansion<decltype(std::conj(x[0])), order> res;
  for (std::size_t i = 0; i < order; ++i) {
    res[i] = std::conj(x[i]);
  }
  return res;
}
