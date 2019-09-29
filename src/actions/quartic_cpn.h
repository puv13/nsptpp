/// \file
///
/// \brief Quartic CP(N) action


#pragma once

#include "./cpn_common_aux.h"

template <typename T, size_t N, size_t dim>
class QuarticCPNAction {

 private:
  SiteLattice<CP<T, N>, dim> const *lattice;
  double beta;
  typedef CP<T, N> CPtype;

 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************

  /// Default constructor
  QuarticCPNAction() : beta(1.0) { lattice = nullptr; }

  /// Construct with given lattice
  QuarticCPNAction(SiteLattice<CP<T, N>, dim> const &orig,
                   const double &coupling = 1.0)
    : lattice(&orig), beta(coupling)
  {
  }


  // **********************************************************************
  // Deconstructor
  // **********************************************************************
  ~QuarticCPNAction() = default;


  // **********************************************************************
  // Member Functions
  // **********************************************************************

  /// Setter for beta
  void setbeta(const double &coupling) { beta = coupling; }
  /// Getter for beta
  double getbeta() { return this->beta; }

  template <typename P = double>
  auto atSite(size_t const &idx, const BoundaryCondition<P, dim> &bc =
                                   BoundaryCondition<P, dim>()) const
    -> decltype(std::abs(T(1.)))
  {
    auto lat = *lattice;
    SiteLatticeIterator<const SiteLattice<CPtype, dim>> fli(idx, lat);

    auto sum = std::abs(T(0.0));

    auto site = fli.site();

    for (auto i = 0lu; i < dim; i++) {
      auto sp = std::abs(
        scalar_prod(site, fli.neighborSite(i, Direction::FORWARD, bc)));
      sum += sp * sp;
    }

    return -1.0 * (this->beta) * static_cast<double>(N) * sum;
  }


  template <typename P = double>
  auto
  total(const BoundaryCondition<P, dim> &bc = BoundaryCondition<P, dim>()) const
    -> decltype(std::abs(T(1.)))
  {
    // const auto lat = *lattice;
    // SiteLatticeIterator<const SiteLattice<CPtype, dim>> fli(lat);

    auto res = std::abs(T(0.0));

    // Loop over all sites
    for (auto fli = lattice->begin(); fli != lattice->end(); ++fli) {
      auto site = fli.site();

      for (auto i = 0lu; i < dim; i++) {
        auto sp = std::abs(
          scalar_prod(site, fli.neighborSite(i, Direction::FORWARD, bc)));
        res += sp * sp;
      }
    }

    return -1. * (this->beta) * static_cast<double>(N) * res;
  }


  template <typename P = double>
  auto average(const BoundaryCondition<P, dim> &bc =
                 BoundaryCondition<P, dim>()) const -> decltype(std::abs(T(1.)))
  {
    auto vol = static_cast<double>(lattice->volume());
    return total(bc) / vol;
  }

  template <typename P = double>
  auto energyDensity(
    const BoundaryCondition<P, dim> &bc = BoundaryCondition<P, dim>()) const
    -> decltype(std::abs(T(1.)))
  {
    auto vol = static_cast<double>(lattice->volume());
    return -total(bc) / (vol * this->beta);
  }
};


// ####################################################################################################
// ####################################################################################################
// Langevin related code
// ####################################################################################################
// ####################################################################################################


// --------------------------------------------------------------------------------
// Helper functions
// --------------------------------------------------------------------------------
namespace QuarticCPNHelpers {

template <size_t N, size_t dim>
struct QuarticCPForceStruct {
  nummat<N> site{ 0.0 };
};

/// Simple printing of QuarticCPForce struct
template <size_t N, size_t dim>
std::ostream &operator<<(std::ostream &stream,
                         QuarticCPForceStruct<N, dim> const &cfs)
{
  stream.precision(6);
  stream << std::scientific << "Force site: " << cfs.site;
  stream << std::endl;
  return stream;
}


template <size_t N>
struct QuarticCPCasimirStruct {
  const double site = static_cast<double>(N);
};

/// Simple printing of Casimir struct
template <size_t N>
std::ostream &operator<<(std::ostream &stream,
                         QuarticCPCasimirStruct<N> const &cps)
{
  stream.precision(6);
  stream << std::scientific << "Casimir site: " << cps.site << std::endl;
  return stream;
}

template <size_t N>
nummat<N> operator*(nummat<N> const &left, QuarticCPCasimirStruct<N> const &CA)
{

  nummat<N> res;

  res = left * CA.site;

  return res;
}


template <size_t N>
nummat<N> operator*(QuarticCPCasimirStruct<N> const &CA, nummat<N> const &right)
{
  nummat<N> res;
  res = CA.site * right;
  return res;
}
}


namespace QuarticCPNExpHelpers {

/// Computes the Hermitian matrix
/// scalar*left*dagger(right)+right*dagger(scalar*left) and adds it to M.
/// Warning! This only fills the upper part of the Hermitian matrix. M should be
/// a compressed Hermitian matrix in CblasUpper form!
template <size_t N, size_t order>
void external_product(Expansion<std::complex<double>, order> const &scalar,
                      Expansion<CP<std::complex<double>, N>, order> const &x,
                      Expansion<CP<std::complex<double>, N>, order> const &y,
                      Expansion<nummat<N>, order> &M)
{
  // Multiply  x by scalar
  auto xbar = scalar * x;
  constexpr std::complex<double> dummy{ 1., 0. };

  for (auto n = 0LU; n < order; ++n) {
    for (auto l = 0LU; l <= n; ++l) {
      cblas_zher2(CBLAS_ORDER::CblasRowMajor, CBLAS_UPLO::CblasUpper, N,
                  reinterpret_cast<const double *>(&dummy),
                  reinterpret_cast<const double *>(&xbar[l][0]), 1,
                  reinterpret_cast<const double *>(&y[n - l][0]), 1,
                  reinterpret_cast<double *>(M[n].data()), N);
    }
  }
}


template <size_t N>
using QuarticCPForceStruct = nummat<N>;


template <size_t N, size_t dim, typename Generator>
class QuarticrandomCPdf {
 public:
  QuarticCPForceStruct<N> operator()(Generator &rng) const
  {

    QuarticCPForceStruct<N> rd;
    rd = noiseSUN<N, Generator>(rng, CP_Langevin_noise_scale);

    return rd;
  }
};
}


// --------------------------------------------------------------------------------
// Drift Force
// --------------------------------------------------------------------------------
namespace NoExpansion {

using namespace QuarticCPNHelpers;

/// CP(N-1) Drift Force for quartic action
template <size_t N, size_t dim, typename P = double>
class QuarticCPForce
  : public Force<SiteLattice<CP<std::complex<double>, N>, dim>,
                 std::vector<nummat<N>>> {
 private:
  BoundaryCondition<P, dim> bc;

 public:
  typedef nummat<N> QuarticCPdf;
  typedef std::vector<QuarticCPdf> QuarticCPForceField;
  typedef SiteLattice<CP<std::complex<double>, N>, dim> QuarticCPlattice;

  /// Returns force field using standard boundary cond.
  QuarticCPForceField operator()(QuarticCPlattice const &lat) const
  {

    // Pre-allocate memory
    QuarticCPForceField FF;
    FF.reserve(lat.volume());


    for (auto it = lat.begin(); it != lat.end(); ++it) {

      QuarticCPdf res;

      const auto II = std::complex<double>(0.0, 1.0);

      auto unit_matrix = cmplx(0.0);
      for (auto i = 0lu; i < dim; i++) {
        auto field = it.site();

        auto field_fwd = it.neighborSite(i, Direction::FORWARD, bc);
        auto field_bwd = it.neighborSite(i, Direction::BACKWARD, bc);

        // CP field
        // ----------------------------------------------------------------------

        // First
        // -----------------------
        // External product part
        auto factor = -1. * II * 0.5 * scalar_prod(field, field_fwd);
        // std::cout << "nu: " << i << "\tfactor:" << factor << std::endl;
        external_product(factor * N, field, field_fwd, res);
        // Diagonal part
        unit_matrix += 2 * std::real(scalar_prod(field_fwd, field) * factor);
        // std::cout << "res:\n" << res.site <<std::endl;

        // Second (relative minus sign!)
        // -----------------------
        factor = -1.0 * II * 0.5 * scalar_prod(field, field_bwd);
        // std::cout << "nu: " << i << "\tfactor:" << factor << std::endl;
        // External product part
        external_product(factor * N, field, field_bwd, res);
        // Diagonal part
        auto s_f_ffwd = 2 * std::real(scalar_prod(field_bwd, field) * factor);
        // std::cout << "res:\n" << res.site <<std::endl;
        unit_matrix += s_f_ffwd;
      }


      inflate(res);
      // std::cout << "Unit Matrix: " << unit_matrix << std::endl;
      res += (unit_matrix * nummat<N>(-1.0));

      FF.push_back(res);
      // ################################################################################
    }

    return FF;
  }

  /// Constructors
  explicit QuarticCPForce() : bc(BoundaryCondition<P, dim>()) {}
  explicit QuarticCPForce(const BoundaryCondition<P, dim> &_bc) : bc(_bc) {}
  ~QuarticCPForce(){};
  /// Compiler complains about constexpr here ???
  bool isGroupValued() const { return false; }
};
}

namespace Expan {

using namespace QuarticCPNExpHelpers;
using namespace CommonCPNExpHelpers;

/// CP(N-1) Drift Force for quartic action
template <size_t N, size_t dim, size_t order, typename P = double>
class QuarticCPForce
  : public Force<
      SiteLattice<Expansion<CP<std::complex<double>, N>, order>, dim>,
      std::vector<Expansion<nummat<N>, order>>> {
 private:
  BoundaryCondition<P, dim> bc;

 public:
  typedef Expansion<nummat<N>, order> QuarticCPdf;
  typedef std::vector<QuarticCPdf> QuarticCPForceField;
  typedef SiteLattice<Expansion<CP<std::complex<double>, N>, order>, dim>
    QuarticCPlattice;

  /// Returns force field using bc boundary cond.
  QuarticCPForceField operator()(QuarticCPlattice const &lat) const
  {

    // Pre-allocate memory
    QuarticCPForceField FF;
    FF.reserve(lat.volume());


    for (auto it = lat.begin(); it != lat.end(); ++it) {

      QuarticCPdf res;

      const auto II = std::complex<double>(0.0, 1.0);

      Expansion<std::complex<double>, order> unit_matrix(
        std::complex<double>(0.0));
      for (auto i = 0lu; i < dim; i++) {
        auto field = it.site();

        auto field_fwd = it.neighborSite(i, Direction::FORWARD, bc);
        auto field_bwd = it.neighborSite(i, Direction::BACKWARD, bc);

        // CP field
        // ----------------------------------------------------------------------

        // First
        // -----------------------
        // External product part
        auto factor = scalar_prod(field, field_fwd) * (-1. * II * 0.5);
        // std::cout << "nu: " << i << "\tfactor:" << factor << std::endl;
        external_product(factor * static_cast<double>(N), field, field_fwd,
                         res);
        // Diagonal part
        unit_matrix += 2 * real(scalar_prod(field_fwd, field) * factor);
        // std::cout << "res:\n" << res.site <<std::endl;

        // Second (relative minus sign!)
        // -----------------------
        factor = scalar_prod(field, field_bwd) * (-1.0 * II * 0.5);
        // std::cout << "nu: " << i << "\tfactor:" << factor << std::endl;
        // External product part
        external_product(factor * static_cast<double>(N), field, field_bwd,
                         res);
        // Diagonal part
        auto s_f_ffwd = 2 * real(scalar_prod(field_bwd, field) * factor);
        // std::cout << "res:\n" << res.site <<std::endl;
        unit_matrix += s_f_ffwd;
      }


      for (auto i = 0LU; i < order; ++i) {
        inflate(res[i]);
      }
      // std::cout << "Unit Matrix: " << unit_matrix << std::endl;
      res += (unit_matrix) * (nummat<N>(-1.0));

      FF.push_back(res);
      // ################################################################################
    }

    return FF;
  }

  /// Constructors
  explicit QuarticCPForce() : bc(BoundaryCondition<P, dim>()) {}
  explicit QuarticCPForce(const BoundaryCondition<P, dim> &_bc) : bc(_bc) {}
  ~QuarticCPForce(){};
  /// Compiler complains about constexpr here ???
  bool isGroupValued() const { return false; }
};
}
