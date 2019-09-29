/// \file
/// A helper class defining numerical complex matrices for use with BLAS.
/// Needs a CBLAS compatible BLAS library!

#pragma once
#include <array>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>

#ifndef USE_MKL
#include <cblas.h>
#else
#include <mkl.h>
#endif

// **********************************************************************
/// Numerical complex matrix class
/// \brief This class defines a complex matrix with n rows and n columns
///
/// This class implements a complex matrix. The entries are stored in an
/// std::array in row-major order, so that it is easy to implement linear
/// algebra functions using the standard CBLAS interface.
// **********************************************************************
template <std::size_t n>
class nummat {

 private:
  std::array<std::complex<double>, n * n> mat;

 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************

  constexpr nummat() { mat.fill(std::complex<double>(0.0)); }
  /// Construct a matrix from an array
  constexpr nummat(std::array<std::complex<double>, n * n> const &init)
    : mat(init)
  {
  }
  constexpr nummat(std::complex<double> const &c)
  {
    mat.fill(std::complex<double>(0.0));
    for (auto i = 0LU; i < n; ++i) {
      mat[i + n * i] = c;
    }
  }

  ///! Copy constructor
  constexpr nummat(nummat<n> const &orig) = default;

  // **********************************************************************
  // Methods
  // **********************************************************************

  /// Access elements of matrix
  std::complex<double> &operator()(std::size_t const i, std::size_t const j)
  {
    return mat[i * n + j];
  }

  constexpr std::complex<double> const &operator()(std::size_t const i,
                                                   std::size_t const j) const
  {
    return mat[i * n + j];
  }

  std::complex<double> *data() noexcept { return mat.data(); }

  std::complex<double> const *data() const noexcept { return mat.data(); }


  /// Matrix multiplication
  constexpr nummat<n> operator*(nummat<n> const &rhs) const
  {

    constexpr std::complex<double> alpha(1.0);
    constexpr std::complex<double> beta(0.0);

    nummat<n> res;

    cblas_zgemm(CBLAS_ORDER::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans,
                CBLAS_TRANSPOSE::CblasNoTrans, n, n, n,
                reinterpret_cast<const double *>(&alpha),
                reinterpret_cast<const double *>(this->data()), n,
                reinterpret_cast<const double *>(rhs.data()), n,
                reinterpret_cast<const double *>(&beta),
                reinterpret_cast<double *>(res.data()), n);
    return res;
  }

  /// Matrix subtraction
  constexpr nummat<n> operator-(nummat<n> const &rhs) const
  {

    constexpr std::complex<double> alpha(-1.0);

    nummat<n> res(this->mat);

    cblas_zaxpy(n * n, reinterpret_cast<const double *>(&alpha),
                reinterpret_cast<const double *>(rhs.data()), 1,
                reinterpret_cast<double *>(res.data()), 1);
    return res;
  }


  /// Matrix addition
  constexpr nummat<n> operator+(nummat<n> const &rhs) const
  {

    std::complex<double> alpha(1.0);

    nummat<n> res(this->mat);

    cblas_zaxpy(n * n, reinterpret_cast<const double *>(&alpha),
                reinterpret_cast<const double *>(rhs.data()), 1,
                reinterpret_cast<double *>(res.data()), 1);
    return res;
  }


  /// Operator "+=" for matrix addition
  void operator+=(nummat<n> const &rhs)
  {

    constexpr std::complex<double> alpha(1.0);

    cblas_zaxpy(n * n, reinterpret_cast<const double *>(&alpha),
                reinterpret_cast<const double *>(rhs.data()), 1,
                reinterpret_cast<double *>(this->data()), 1);
  }

  /// Operator "-=" for matrix subtraction
  void operator-=(nummat<n> const &rhs)
  {

    constexpr std::complex<double> alpha(-1.0);

    cblas_zaxpy(n * n, reinterpret_cast<const double *>(&alpha),
                reinterpret_cast<const double *>(rhs.data()), 1,
                reinterpret_cast<double *>(this->data()), 1);
  }


  /// Scalar multiplication (from right)
  nummat<n> operator*(std::complex<double> const &rhs) const
  {

    nummat<n> res(this->mat);

    cblas_zscal(n * n, reinterpret_cast<const double *>(&rhs),
                reinterpret_cast<double *>(res.data()), 1);
    return res;
  }

  /// Operator "*=" for scalar multiplication. (Rescaling)
  void operator*=(std::complex<double> const &a)
  {
    cblas_zscal(n * n, reinterpret_cast<const double *>(&a),
                reinterpret_cast<double *>(this->data()), 1);
  }

  /// Scalar division
  constexpr nummat<n> operator/(std::complex<double> const &rhs) const
  {
    return (*this) * (1. / rhs);
  }


  /// Equality check \todo Better float equal?
  constexpr bool operator==(nummat<n> const &right) const
  {
    bool res = true;
    for (auto i = 0LU; i < n; ++i) {
      for (auto j = 0LU; j < n; ++j) {
        res = ((*this)(i, j) == right(i, j));
        if (not res)
          return res;
      }
    }
    return res;
  }


  // Multiplication with Hermitian conjugate
  constexpr nummat<n> mult_conj(nummat<n> const &rhs) const
  {

    constexpr std::complex<double> alpha(1.0);
    constexpr std::complex<double> beta(0.0);

    nummat<n> res;

    cblas_zgemm(CBLAS_ORDER::CblasRowMajor, CBLAS_TRANSPOSE::CblasNoTrans,
                CBLAS_TRANSPOSE::CblasConjTrans, n, n, n,
                reinterpret_cast<const double *>(&alpha),
                reinterpret_cast<const double *>(this->data()), n,
                reinterpret_cast<const double *>(rhs.data()), n,
                reinterpret_cast<const double *>(&beta),
                reinterpret_cast<double *>(res.data()), n);
    return res;
  }


  // Real part
  // constexpr nummat<n> real() const {

  //     nummat<n> res;

  //     for (auto i=0LU; i<n;++i){
  // 	for(auto j=0LU; j<n;++j){
  // 	    res(i,j)= this->mat[i*n+j];
  // 	}
  //     }

  //     return res;
  // }

  /// Hermitian conjugate
  constexpr nummat<n> dagger() const
  {
    nummat<n> res;

    for (auto i = 0LU; i < n; ++i) {
      for (auto j = 0LU; j < n; ++j) {
        res(i, j) = std::conj(this->mat[j * n + i]);
      }
    }

    return res;
  }


  /// Matrix inverse \todo Implement for general matrix
  constexpr nummat<n> inverse() const
  {

    if (*this == nummat<n>(1.)) {
      return nummat<n>(this->mat);
    }
    else {
      throw std::runtime_error("Inverse only implemented for unit matrix!");
    }
  }


  /// Matrix log \todo Implement for general matrix
  constexpr nummat<n> log() const
  {

    nummat<n> res(0.0);

    if (*this == nummat<n>(1.)) {
      return res;
    }
    else {
      throw std::runtime_error("Log only implemented for unit matrix!");
    }
  }


  // **********************************************************************
  // Friends
  // **********************************************************************

  /// Multiplication with scalar from the left
  friend constexpr nummat<n> operator*(std::complex<double> const &lhs,
                                       nummat<n> const &rhs)
  {
    auto res = rhs;
    return res * lhs;
  }

  /// Simple printing of nummat
  friend std::ostream &operator<<(std::ostream &stream, nummat<n> const &nm)
  {

    stream.precision(6);

    for (auto i = 0ul; i < n * n; ++i) {
      stream << std::fixed << "(" << std::setw(13) << std::real(nm.mat[i])
             << "," << std::setw(13) << std::imag(nm.mat[i]) << ")\t";

      if (0 == (i + 1) % n) {
        stream << "\n";
      }
    }
    return stream;
  }
};


// Non-member functions

/// Trace
template <std::size_t n>
constexpr std::complex<double> trace(nummat<n> const &M)
{
  std::complex<double> res(0.0);
  for (auto i = 0LU; i < n; ++i) {
    res += M(i, i);
  }

  return res;
}

/// Dagger (needed as a non-member function for expansions)
template <std::size_t n>
constexpr nummat<n> dagger(nummat<n> const &M)
{
  return M.dagger();
}

/// Inverse (needed as a non-member function for expansions)
///
/// Wraps member function
template <std::size_t n>
constexpr nummat<n> inverse(nummat<n> const &M)
{
  return M.inverse();
}

/// Matrix Log (needed as a non-member function for expansions)
///
/// Wraps member function
template <std::size_t n>
constexpr nummat<n> log(nummat<n> const &M)
{
  return M.log();
}


/// Construct diagonal matrix from array
template <std::size_t n>
constexpr nummat<n> diag(std::array<std::complex<double>, n> const &init)
{
  nummat<n> M;
  for (auto i = 0LU; i < n; ++i) {
    M(i, i) = init[i];
  }

  return M;
}

/// Frobenius norm of a matrix
template <std::size_t n>
double normF(nummat<n> const &M)
{
  return cblas_dznrm2(n * n, reinterpret_cast<const double *>(M.data()), 1);
}

/// Matrix exponential
/// (Via Taylor series. This is probably neither fast nor stable.)
template <std::size_t n>
nummat<n> exp(nummat<n> const &M,
              const double &error = std::numeric_limits<double>::epsilon(),
              const size_t &maxiter = 100LU)
{

  double prec_goal = error;
  // Sanity check
  if (prec_goal < std::numeric_limits<double>::epsilon()) {
    prec_goal = std::numeric_limits<double>::epsilon();
    std::cout
      << "Warning! Precision goal was reset to machine epsilon in matrix exp"
      << prec_goal << std::endl;
  }

  nummat<n> res{ 1.0 }, old{ M }, aux{ 1.0 };

  double eps = 1.0;
  double fac = 1.0;

  size_t iter = 0LU;
  while (eps > error && iter < maxiter) {

    ++iter;
    fac *= static_cast<double>(iter);
    aux = aux * M;

    res += aux / fac;
    // Error estimate. (Quite pessimistic)
    eps = normF(old - res);
    // std::cout << "Error: " << eps << std::endl;
    old = res;
  }

  // std::cout << iter << std::endl;
  return res;
}

/// Inflates compressed Hermitian matrix (CblasUpper format) to a dense matrix
template <size_t NN>
void inflate(nummat<NN> &M)
{
  for (auto i = 1LU; i < NN; ++i) {
    for (auto j = 0LU; j < i; ++j) {
      M(i, j) = std::conj(M(j, i));
    }
  }
}
