/// \file
///
/// \brief Canonical CP(N) action


#pragma once

#include "../src/gaugegroups/u1.h"
#include "./cpn_common_aux.h"


template <typename T, size_t N, size_t dim>
class CanonicalCPNAction {

 private:
  FullLattice<CP<T, N>, U1, dim> const *lattice;
  double beta;
  typedef CP<T, N> CPtype;

 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************

  /// Default constructor
  CanonicalCPNAction() : beta(1.0) { lattice = nullptr; }

  /// Construct with given lattice
  CanonicalCPNAction(FullLattice<CP<T, N>, U1, dim> const &orig,
                     const double &coupling = 1.0)
    : lattice(&orig), beta(coupling)
  {
  }


  // **********************************************************************
  // Deconstructor
  // **********************************************************************
  ~CanonicalCPNAction() = default;


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
    FullLatticeIterator<const FullLattice<CPtype, U1, dim>> fli(idx, lat);

    auto sum = T(0.0);

    auto site = fli.site();

    for (auto i = 0lu; i < dim; i++) {
      sum += scalar_prod(site, fli.neighborSite(i, Direction::FORWARD, bc)) *
             fli.link(i).value();
    }

    return -2. * (this->beta) * N * std::real(sum);
  }


  template <typename P = double>
  auto
  total(const BoundaryCondition<P, dim> &bc = BoundaryCondition<P, dim>()) const
    -> decltype(std::abs(T(1.)))
  {
    const auto lat = *lattice;
    FullLatticeIterator<const FullLattice<CPtype, U1, dim>> fli(lat);

    auto res = T(0.0);

    // Loop over all sites
    for (fli = lat.begin(); fli != lat.end(); ++fli) {
      auto site = fli.site();

      for (auto i = 0lu; i < dim; i++) {
        res += scalar_prod(site, fli.neighborSite(i, Direction::FORWARD, bc)) *
               fli.link(i).value();
      }
    }

    return -2. * (this->beta) * N * (std::real(res));
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


template <typename T, size_t N, size_t dim>
class CanonicalCPNObservables {

 private:
  FullLattice<CP<T, N>, U1, dim> const *lattice;
  double beta;
  typedef CP<T, N> CPtype;

 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************

  /// Default constructor
  CanonicalCPNObservables() : beta(1.0) { lattice = nullptr; }

  /// Construct with given lattice
  CanonicalCPNObservables(FullLattice<CP<T, N>, U1, dim> const &orig,
                          const double &coupling = 1.0)
    : lattice(&orig), beta(coupling)
  {
  }


  // **********************************************************************
  // Deconstructor
  // **********************************************************************

  ~CanonicalCPNObservables() = default;


  // **********************************************************************
  // Member Functions
  // **********************************************************************

  bool check_normalisation(double norm_eps = 1.e-9)
  {
    const auto lat = *lattice;
    FullLatticeIterator<const FullLattice<CPtype, U1, dim>> fli(lat);
    /// Loop over all sites
    for (fli = lat.begin(); fli != lat.end(); ++fli) {

      auto site = fli.site();

      auto snorm = site.norm();

      if (std::abs(snorm - 1) > norm_eps) {
        // std::cout << "Site norm at index " << fli.index() << " is " << snorm
        // 	      << std::endl;
        return false;
      }

      for (auto i = 0lu; i < dim; i++) {
        auto lnorm = fli.link(i).norm();
        if (std::abs(lnorm - 1) > norm_eps) {
          std::cout << "Link norm at index " << fli.index() << " and site " << i
                    << " is " << snorm << std::endl;

          return false;
        }
      }
    }

    return true;
  }

  /// Returns (the real part of the trace of ) the plaquette
  template <typename P = double>
  auto plaquette(
    size_t site, size_t mu, size_t nu,
    BoundaryCondition<P, dim> const &bc = BoundaryCondition<P, dim>()) const
    -> decltype(std::real(trace(U1(1.))))
  {
    if (site >= lattice->volume()) {
      throw std::runtime_error("Site out of range");
    }
    if (mu >= dim) {
      throw std::runtime_error("mu out of range");
    }
    if (nu >= dim) {
      throw std::runtime_error("nu out of range");
    }
    if (mu == nu) {
      throw std::runtime_error("mu must not equal nu");
    }

    auto it = lattice->begin();
    it = it[site];

    // U_µν(x) =U_µ(x)U_ν(x+µ)dagger(U_µ(n+ν))dagger(U_ν(x))

    auto plaq = it.link(mu);
    plaq *= it.neighborLink(mu, nu, Direction::FORWARD, bc);
    plaq *= dagger(it.neighborLink(nu, mu, Direction::FORWARD, bc));
    plaq *= dagger(it.link(nu));

    // auto coord = lattice->linearIndexToCoord(site);
    // std::cout << "( ";
    // for (auto && el : coord)
    // {
    // 	std::cout << el << "\t";
    // }
    // std::cout << ") ";
    // std::cout << "Site: " << site << "\tIt: " << it.index()
    // 	                  << "\tIt_nu: " << it_nu.index()
    // 	      << "\tIt_mu: " << it_mu.index() << "\n";

    return std::real(trace(plaq));
  }


  /// Compute mean plaquette in the mu-nu plane
  template <typename P = double>
  auto meanPlaquette(size_t mu, size_t nu, BoundaryCondition<P, dim> const &bc =
                                             BoundaryCondition<P, dim>()) const
    -> decltype(std::real(trace(U1(1.))))
  {

    auto res = plaquette(0lu, mu, nu, bc);
    for (auto i = 1ul; i < lattice->volume(); ++i) {
      res += plaquette(i, mu, nu);
    }

    return res / static_cast<double>(lattice->volume());
  }


  /// Compute global plaquette average
  template <typename P = double>
  auto meanPlaquette(
    BoundaryCondition<P, dim> const &bc = BoundaryCondition<P, dim>()) const
    -> decltype(std::real(trace(U1(1.))))
  {

    decltype(std::real(trace(U1(1.)))) res = 0.0;

    for (auto mu = 0lu; mu < dim; ++mu) {
      for (auto nu = mu + 1lu; nu < dim; ++nu) {
        auto res_mu_nu = plaquette(0ul, mu, nu);
        for (auto i = 1ul; i < lattice->volume(); ++i) {
          res_mu_nu += plaquette(i, mu, nu, bc);
        }

        res_mu_nu *= 1. / static_cast<double>(lattice->volume());

        res += res_mu_nu;
      }
    }

    return 2.0 * res / (static_cast<double>(dim * (dim - 1)));
  }


  /// Compute trace of P(x)*P(y), where P is the tensor product of a CP(N-1)
  /// field z, i.e
  /// P_ij=conj(z_i)*z_J
  auto tracePP(size_t x, size_t y) -> decltype(std::abs(T(1.)))
  {
    auto V = lattice->volume();
    if (x >= V || y >= V) {
      throw std::runtime_error("Site out of range");
    }

    auto sp = scalar_prod(lattice->getSite(x), lattice->getSite(y));

    auto res = std::abs(sp);
    res *= res;

    res = res - 1. / static_cast<double>(N);

    return res;
  }


  /// Magnetic susceptibility chi_m
  auto chi_m() -> decltype(std::abs(T(1.)))
  {
    if (dim != 2) {
      throw std::runtime_error("chi_m not implemented for dim != 2.");
    }


    auto V = lattice->volume();

    auto res = std::abs(T(0.0));
    for (auto x = 0ul; x < V; ++x) {
      for (auto y = 0ul; y < V; ++y) {
        auto trace = tracePP(x, y);
        res += trace;
      }
    }

    return res / static_cast<double>(V);
  }

  /// Correlation length xi_G
  /// For a defintion see for example Campostrini,Rossi and Vicari PRD 46,
  /// 6, 1992, p 2647
  /// The time component has index zero!
  /// Not implemented for dim != 2
  auto xi_G() -> decltype(std::abs(T(1.)))
  {
    if (dim != 2) {
      throw std::runtime_error("xi_G not implemented for dim != 2.");
    }

    std::complex<double> II(0., 1.);
    constexpr double pi_2 = 2. * M_PI;

    T G_00(0.0);
    T G_10(0.0);
    auto L = lattice->dimensionsArray()[0];
    auto factor = II * pi_2 / static_cast<double>(L);

    for (auto x = 0ul; x < lattice->volume(); ++x) {
      auto x_cord = lattice->linearIndexToCoord(x);
      auto t_x = x_cord[0];

      for (auto y = 0ul; y < lattice->volume(); ++y) {
        auto y_cord = lattice->linearIndexToCoord(y);
        auto t_y = y_cord[0];

        auto trace = tracePP(x, y);

        // time difference
        auto t_diff = static_cast<double>((t_x - t_y + L) % L);

        G_00 += trace;
        G_10 += trace * std::exp(factor * t_diff);
      }
    }

    auto ratio = G_00 / G_10 - 1.;
    auto sin_factor = 2.0 * std::sin(pi_2 / (2. * static_cast<double>(L)));
    sin_factor *= sin_factor;

    // This should be real ... shouldn't it ? ...
    auto res = std::sqrt(ratio / sin_factor);

    // Well ...
    return std::abs(res);
  }

  /// Compute all Wall-Wall correlators
  /// Not implemented for dim != 2
  auto WW_correlators() -> std::vector<decltype(std::abs(T(1.)))>
  {
    if (dim != 2) {
      throw std::runtime_error(
        "Wall-Wall correlators not implemented for dim != 2.");
    }


    auto dims = lattice->dimensionsArray();
    auto Lt = dims[0];
    auto Ls = dims[1];

    std::vector<decltype(std::abs(T(1.)))> WW_res(Lt, std::abs(T(0.0)));

    for (auto x = 0ul; x < lattice->volume(); ++x) {
      auto x_cord = lattice->linearIndexToCoord(x);
      auto t_x = x_cord[0];

      for (auto y = 0ul; y < lattice->volume(); ++y) {
        auto y_cord = lattice->linearIndexToCoord(y);
        auto t_y = y_cord[0];

        auto trace = tracePP(x, y);

        // time difference
        auto t_diff = (t_x - t_y + Lt) % Lt;

        WW_res[t_diff] += trace;
      }
    }


    // Normalise
    for (auto i = 0LU; i < WW_res.size(); ++i) {
      WW_res[i] /= static_cast<double>(Ls);
    }


    return WW_res;
  }

  /// Top. Plaquette, an observable needed for our top. charge definition
  /// on the lattice.
  template <typename P = double>
  auto
  top_plaq(size_t site,
           const BoundaryCondition<P, dim> &bc = BoundaryCondition<P, dim>())
    -> decltype(std::abs(T(1.0)))
  {
    if (site >= lattice->volume()) {
      throw std::runtime_error("Site out of range");
    }

    if (dim != 2) {
      throw std::runtime_error(
        "'Topological plaquette' not implemented for dim != 2.");
    }

    auto res = std::abs(T(0.0));

    auto it = lattice->begin();
    it = it[site];
    auto neib = it.neighbor(0, 1);

    auto z_x = it.site();
    auto z_x_mu = it.neighborSite(0, Direction::FORWARD, bc);
    auto z_x_nu = it.neighborSite(1, Direction::FORWARD, bc);
    auto z_x_mu_nu = neib.neighborSite(1, Direction::FORWARD, bc);

    // + θx,µ
    res += std::arg(scalar_prod(z_x, z_x_mu));

    // + θx+µ,ν
    res += std::arg(scalar_prod(z_x_mu, z_x_mu_nu));

    // - θx+ν,µ
    res -= std::arg(scalar_prod(z_x_nu, z_x_mu_nu));

    // - θx,ν
    res -= std::arg(scalar_prod(z_x, z_x_nu));


    // Bring to (-pi, pi] interval

    while (res <= -1. * M_PI) {
      res += 2 * M_PI;
    }

    while (res > M_PI) {
      res -= 2 * M_PI;
    }

    return res;
  }

  /// Top. Plaquette, an observable needed for our top. charge definition
  /// on the lattice.
  template <typename P = double>
  auto Q_top(const BoundaryCondition<P, dim> &bc = BoundaryCondition<P, dim>())
    -> decltype(std::abs(T(1.0)))
  {

    auto res = std::abs(T(0.0));

    for (auto x = 0ul; x < lattice->volume(); ++x) {
      res += top_plaq(x, bc);
    }

    return res / (2. * M_PI);
  }

  auto GaugeFieldSpatialSum()
    -> std::vector<std::array<decltype(std::abs(T(1.))), dim>>
  {

    auto dims = lattice->dimensionsArray();
    std::array<decltype(std::abs(T(1.))), dim> A;
    std::vector<std::array<decltype(std::abs(T(1.))), dim>> Avec;

    A.fill(std::abs(T(0.)));

    // Output vector has size of temporal lattice extend
    for (auto l = 0LU; l < dims[0]; ++l) {
      Avec.push_back(A); // Init with "zero"
    }


    // Loop over all lattice sites
    // std::vector<size_t> tvec;
    // tvec.resize(Avec.size());
    for (auto fli = lattice->begin(); fli != lattice->end(); ++fli) {
      auto cord = lattice->linearIndexToCoord(fli.index());
      auto t = cord[0];


      // tvec[t]++;
      for (auto mu = 0LU; mu < dim; ++mu) {
        (Avec[t])[mu] += std::sin((fli.link(mu)).phase());
      }
    }

    // std::cout << "TVEC: ";
    // std::for_each(tvec.begin(),tvec.end(),[](auto a){std::cout << a <<
    // "\t";});
    // No normalisation !!!
    return Avec;
  }
};


/// Class that implements gauge fixing and related functions
/// \brief Gauge fixing class
///
/// \note The functions defined here only make sense for PERIODIC boundary
/// conditions.
template <typename T, size_t N, size_t dim>
class CanonicalCPNGaugeFix {

 private:
  FullLattice<CP<T, N>, U1, dim> *const lattice;
  double beta;
  typedef CP<T, N> CPtype;

 public:
  // **********************************************************************
  // Constructors
  // **********************************************************************

  /// Default constructor
  CanonicalCPNGaugeFix() : beta(1.0) { lattice = nullptr; }

  /// Construct with given lattice
  CanonicalCPNGaugeFix(FullLattice<CP<T, N>, U1, dim> &orig,
                       const double &coupling = 1.0)
    : lattice(&orig), beta(coupling)
  {
  }


  // **********************************************************************
  // Deconstructor
  // **********************************************************************

  ~CanonicalCPNGaugeFix() = default;


  // **********************************************************************
  // Member Functions
  // **********************************************************************
  void randomGaugeTransformation()
  {
    auto &lat = *lattice;
    FullLatticeIterator<FullLattice<CPtype, U1, dim>> fli(0, lat);

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


    for (fli = lat.begin(); fli != lat.end(); ++fli) {
      double rphase = M_PI * dist(generator);

      // Apply gauge transformation to site ...
      std::complex<double> ii(0, 1);
      (*fli).site *= std::exp(-ii * rphase);

      // ... and to all adjacent links
      for (auto i = 0UL; i < dim; ++i) {
        auto n_dir = fli.neighbor(i, -1);
        (*fli).links[i] *= std::exp(-ii * rphase);
        (*n_dir).links[i] *= std::exp(ii * rphase);
      }
    }
  }

  double LandauGaugeFunctional() const
  {

    double res = 0;
    for (auto it = lattice->begin(); it != lattice->end(); ++it) {
      auto links = it.linkarray();

      for (size_t i = 0; i < lattice->dimensions(); ++i) {
        res += std::cos(links[i].phase());
      }
    }
    return res;
  }

  void localLandauGauge(size_t const &idx, const double or_param = 1.0)
  {
    auto &lat = *lattice;
    FullLatticeIterator<FullLattice<CPtype, U1, dim>> fli(idx, lat);

    // Calculate the phase to locally fix the config to Landau gauge
    double num = 0.0;
    double den = 0.0;

    for (auto i = 0UL; i < dim; ++i) {
      auto n_dir = fli.neighbor(i, -1);
      auto here = (fli.link(i)).phase();
      auto neib = (n_dir.link(i)).phase();

      num += (std::sin(here) - std::sin(neib));
      den += (std::cos(here) + std::cos(neib));
    }

    auto phase = or_param * std::atan(num / den);

    // Apply gauge transformation to site ...
    std::complex<double> ii(0, 1);
    (*fli).site *= std::exp(-ii * phase);

    // ... and to all adjacent links
    for (auto i = 0UL; i < dim; ++i) {
      auto n_dir = fli.neighbor(i, -1);
      (*fli).links[i] *= std::exp(-ii * phase);
      (*n_dir).links[i] *= std::exp(ii * phase);
    }
  }

  double localLandauGaugeQuality(size_t const &idx)
  {
    auto &lat = *lattice;
    FullLatticeIterator<FullLattice<CPtype, U1, dim>> fli(idx, lat);

    double res = 0.0;
    for (auto i = 0UL; i < dim; ++i) {
      auto n_dir = fli.neighbor(i, -1);
      auto here = (fli.link(i)).phase();
      auto neib = (n_dir.link(i)).phase();

      res += (std::sin(here) - std::sin(neib));
    }

    return res * res;
  }


  /// \Todo Implement boundary condition
  void LandauGaugeSweep(const double or_param = 1.0)
  {

    auto lat = *lattice;
    // FullLatticeIterator<FullLattice<CPtype,U1, dim>> fli(lat);

    // size_t count[] = {0,0};
    // Loop over all even and odd sites separately
    for (auto eo = 0; eo < 2; ++eo) {
// for (fli=lat.begin(); fli != lat.end(); ++fli)
#pragma omp parallel for
      for (auto idx = 0LU; idx < lat.volume(); ++idx) {
        // auto idx = fli.index();
        auto coords = lat.linearIndexToCoord(idx);

        auto csum = std::accumulate(coords.begin(), coords.end(), 0);
        if (eo == csum % 2) {
          localLandauGauge(idx, or_param);
          // count[eo]++;
        }
      }
    }

    // std::cout << "Sites: " << count[0] << "/" << count[1] << std::endl;
    // std::cout << "Volume: " << lat.volume() << std::endl;
  }

  double LandauGaugeQuality()
  {

    auto lat = *lattice;

    double res = 0.0;
    for (auto idx = 0LU; idx < lat.volume(); ++idx) {
      res += localLandauGaugeQuality(idx);
    }

    return res / static_cast<double>(lat.volume());
  }

  size_t LandauGaugeDriver(const size_t gc_num, const size_t sw_num,
                           const double or_param = 1.0)
  {

    double minus_inf = std::numeric_limits<double>::lowest();
    FullLattice<CPtype, U1, dim> gauge_copy(*lattice), best_copy(*lattice);

    double best = minus_inf; // Init with large negative number
    double last = 0;

    size_t iter = 0LU;

    for (auto i = 0UL; i < gc_num; ++i) {

      CanonicalCPNGaugeFix<T, N, dim> gc_gf(gauge_copy, this->beta);
      if (0LU != i) // Use the original config once
      {
        gc_gf.randomGaugeTransformation();
      }

      // double last_LF=gc_gf.LandauGaugeFunctional();
      for (auto s = 0UL; s < sw_num; ++s) {
        gc_gf.LandauGaugeSweep(or_param);
        ++iter;
        // double new_LF=gc_gf.LandauGaugeFunctional();

        // if (std::abs(new_LF-last_LF) <
        // 1.e-6/static_cast<double>((*lattice).volume()))
        if (gc_gf.LandauGaugeQuality() < 1.e-9) {
          // last_LF = new_LF;
          last = gc_gf.LandauGaugeFunctional();
          // std::cout << "\t\t Best local LF: " << last_LF << "\t (best global:
          // " << best << ")"
          // 	  << "(" << ++s << " sweeps)" <<  std::endl;
          break;
        }
        // last_LF = new_LF;
      }

      // last = last_LF;
      if (last > best) {
        best_copy = gauge_copy;
        best = last;
      }
    }

    *lattice = best_copy;
    return iter;
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
namespace CanonicalCPNHelpers {

template <size_t N, size_t dim>
struct CPForceStruct {
  nummat<N> site{ 0.0 };
  std::array<std::complex<double>, dim> links;
};

/// Simple printing of CPForce struct
template <size_t N, size_t dim>
std::ostream &operator<<(std::ostream &stream, CPForceStruct<N, dim> const &cfs)
{

  stream.precision(6);
  stream << std::scientific << "Force site: " << cfs.site;
  for (auto i = 0LU; i < dim; i++) {
    stream << "\n"
           << "Force links:" << cfs.links[i];
  }
  stream << std::endl;
  return stream;
}

template <size_t N>
struct CPCasimirStruct {
  const double site = static_cast<double>(N);
  const double links = 0.;
};

/// Simple printing of Casimir struct
template <size_t N>
std::ostream &operator<<(std::ostream &stream, CPCasimirStruct<N> const &cps)
{

  stream.precision(6);
  stream << std::scientific << "Casimir site: " << cps.site << "\n"
         << "Casimir links:" << cps.links;
  return stream;
}


template <size_t N, size_t dim>
SiteAndLinks<CP<std::complex<double>, N>, U1, dim>
operator*(CPForceStruct<N, dim> const &cpdf,
          SiteAndLinks<CP<std::complex<double>, N>, U1, dim> const &sl)
{

  SiteAndLinks<CP<std::complex<double>, N>, U1, dim> res;

  res.site = cpdf.site * sl.site;
  for (auto i = 0LU; i < dim; ++i) {
    res.links[i] = cpdf.links[i] * sl.links[i];
  }
  return res;
}


template <size_t N, size_t dim>
CPForceStruct<N, dim> exp(CPForceStruct<N, dim> const &cpdf)
{
  CPForceStruct<N, dim> res;

  res.site = exp(cpdf.site);
  for (auto i = 0LU; i < dim; ++i) {
    res.links[i] = exp(cpdf.links[i]);
  }
  return res;
}


template <size_t N, size_t dim>
CPForceStruct<N, dim> operator*(CPForceStruct<N, dim> const &cpdf,
                                std::complex<double> const &scalar)
{

  CPForceStruct<N, dim> res;

  res.site = scalar * cpdf.site;
  for (auto i = 0LU; i < dim; ++i) {
    res.links[i] = cpdf.links[i] * scalar;
  }
  return res;
}


template <size_t N, size_t dim>
CPForceStruct<N, dim> operator*(CPForceStruct<N, dim> const &lcpdf,
                                CPCasimirStruct<N> const &CA)
{

  CPForceStruct<N, dim> res;

  res.site = lcpdf.site * CA.site;
  for (auto i = 0LU; i < dim; ++i) {
    res.links[i] = lcpdf.links[i] * CA.links;
  }
  return res;
}


template <size_t N, size_t dim>
CPForceStruct<N, dim> operator*(CPCasimirStruct<N> const &CA,
                                CPForceStruct<N, dim> const &rcpdf)
{

  CPForceStruct<N, dim> res;

  res.site = CA.site * rcpdf.site;
  for (auto i = 0LU; i < dim; ++i) {
    res.links[i] = CA.links * rcpdf.links[i];
  }
  return res;
}


template <size_t N, size_t dim>
CPForceStruct<N, dim> operator+(CPForceStruct<N, dim> const &L,
                                CPForceStruct<N, dim> const &R)
{

  CPForceStruct<N, dim> res;

  res.site = L.site + R.site;
  for (auto i = 0LU; i < dim; ++i) {
    res.links[i] = L.links[i] + R.links[i];
  }
  return res;
}


template <size_t N, size_t dim>
CPForceStruct<N, dim> operator-(CPForceStruct<N, dim> const &L,
                                CPForceStruct<N, dim> const &R)
{

  CPForceStruct<N, dim> res;

  res.site = L.site - R.site;
  for (auto i = 0LU; i < dim; ++i) {
    res.links[i] = L.links[i] - R.links[i];
  }
  return res;
}

template <size_t N, size_t dim>
CPForceStruct<N, dim> operator*(std::complex<double> const &scalar,
                                CPForceStruct<N, dim> const &cpdf)
{
  return cpdf * scalar;
}

} // End of namespace "CanonicalCPNHelpers "


namespace CanonicalCPNExpHelpers {


/// Multiplication of U1 and CP expansions
template <size_t N, size_t order>
auto operator*(Expansion<U1, order> const &left,
               Expansion<CP<std::complex<double>, N>, order> const &right)
  -> Expansion<CP<std::complex<double>, N>, order>
{

  Expansion<CP<std::complex<double>, N>, order> res;

  for (auto n = 0LU; n < order; ++n) {
    for (auto i = 0LU; i <= n; ++i) {
      res[n] += left[n - i] * right[i];
    }
  }


  return res;
}

template <size_t N, size_t order>
auto operator*(Expansion<CP<std::complex<double>, N>, order> const &left,
               Expansion<U1, order> const &right)
  -> Expansion<CP<std::complex<double>, N>, order>
{
  return right * left;
}


/// Multiplication of U1 and std::complex<double> expansions
template <size_t order>
auto operator*(Expansion<U1, order> const &left,
               Expansion<std::complex<double>, order> const &right)
  -> Expansion<std::complex<double>, order>
{

  Expansion<std::complex<double>, order> res;

  for (auto n = 0LU; n < order; ++n) {
    for (auto i = 0LU; i <= n; ++i) {
      res[n] += left[n - i].value() * right[i];
    }
  }

  return res;
}


template <size_t order>
auto operator*(Expansion<std::complex<double>, order> const &left,
               Expansion<U1, order> const &right)
  -> Expansion<std::complex<double>, order>
{
  return right * left;
}


template <typename T, size_t order>
std::ostream &operator<<(std::ostream &stream, Expansion<T, order> const &e)
{


  for (auto i = 0ul; i < order; ++i) {
    stream << "order: " << i << " --> " << e[i] << std::endl;
  }
  return stream;
}

/// Computes the Hermitian matrix
/// scalar*left*dagger(right)+right*dagger(scalar*left) and adds it to M.
/// Warning! This only fills the upper part of the Hermitian matrix. M should be
/// a compressed Hermitian matrix in CblasUpper form!
template <size_t N, size_t order>
void external_product(Expansion<U1, order> const &scalar,
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


template <size_t N, size_t dim>
struct CPForceStruct {
  nummat<N> site{ 0.0 };
  std::array<std::complex<double>, dim> links;
};


template <typename T, size_t N, size_t dim>
CPForceStruct<N, dim> operator*(T const &lhs, CPForceStruct<N, dim> rhs)
{
  rhs.site = lhs * rhs.site;
  for (auto i = 0LU; i < dim; ++i) {
    rhs.links[i] = lhs * rhs.links[i];
  }

  return rhs;
}


template <size_t N, size_t dim, typename Generator>
class randomCPdf {
 public:
  CPForceStruct<N, dim> operator()(Generator &rng) const
  {

    CPForceStruct<N, dim> rd;
    rd.site = noiseSUN<N, Generator>(rng, CP_Langevin_noise_scale);


    for (auto i = 0LU; i < rd.links.size(); ++i) {
      rd.links[i] = noiseU1(rng, CP_Langevin_noise_scale);
    }

    return rd;
  }
};

template <size_t N, size_t dim, size_t order>
struct CPForceStructExp {
  Expansion<nummat<N>, order> site{ 0.0 };
  std::array<Expansion<std::complex<double>, order>, dim> links;
};

template <size_t N, size_t dim, size_t order>
CPForceStructExp<N, dim, order> operator+(CPForceStruct<N, dim> const &lhs,
                                          CPForceStructExp<N, dim, order> rhs)
{
  if (order > 1) {
    rhs.site[1] = lhs.site + rhs.site[1];
    for (auto i = 0LU; i < dim; ++i) {
      rhs.links[i][1] = lhs.links[i] + rhs.links[i][1];
    }
  }

  return rhs;
}

template <size_t N, size_t dim, size_t order>
CPForceStructExp<N, dim, order> operator+(CPForceStructExp<N, dim, order> lhs,
                                          CPForceStruct<N, dim> const &rhs)
{
  return (rhs + lhs);
}


template <size_t N>
struct CPCasimirStruct {
  const double site = static_cast<double>(N);
  const double links = 0.;
};


template <size_t N, size_t dim, size_t order>
using SiteAndLinksExp =
  SiteAndLinks<Expansion<CP<std::complex<double>, N>, order>,
               Expansion<U1, order>, dim>;

template <size_t N, size_t dim, size_t order>
SiteAndLinksExp<N, dim, order>
operator*(CPForceStructExp<N, dim, order> const &cpdf,
          SiteAndLinksExp<N, dim, order> const &sl)
{

  SiteAndLinksExp<N, dim, order> res;

  res.site = cpdf.site * sl.site;
  for (auto i = 0LU; i < dim; ++i) {
    res.links[i] = cpdf.links[i] * sl.links[i];
  }
  return res;
}


template <size_t N, size_t dim, size_t order>
CPForceStructExp<N, dim, order> exp(CPForceStructExp<N, dim, order> const &cpdf)
{
  CPForceStructExp<N, dim, order> res;

  res.site = exp(cpdf.site);
  for (auto i = 0LU; i < dim; ++i) {
    res.links[i] = exp(cpdf.links[i]);
  }
  return res;
}


template <size_t N, size_t dim, size_t order>
CPForceStructExp<N, dim, order>
operator*(CPForceStructExp<N, dim, order> const &cpdf,
          std::complex<double> const &scalar)
{

  CPForceStructExp<N, dim, order> res;

  res.site = scalar * cpdf.site;
  for (auto i = 0LU; i < dim; ++i) {
    res.links[i] = cpdf.links[i] * scalar;
  }
  return res;
}


template <size_t N, size_t dim, size_t order>
CPForceStructExp<N, dim, order>
operator*(CPForceStructExp<N, dim, order> const &lcpdf,
          CPCasimirStruct<N> const &CA)
{

  CPForceStructExp<N, dim, order> res;

  res.site = lcpdf.site * CA.site;
  for (auto i = 0LU; i < dim; ++i) {
    res.links[i] = lcpdf.links[i] * CA.links;
  }
  return res;
}


template <size_t N, size_t dim, size_t order>
CPForceStructExp<N, dim, order>
operator*(CPCasimirStruct<N> const &CA,
          CPForceStructExp<N, dim, order> const &rcpdf)
{
  return rcpdf * CA;
}


template <size_t N, size_t dim, size_t order>
CPForceStructExp<N, dim, order>
operator+(CPForceStructExp<N, dim, order> const &L,
          CPForceStructExp<N, dim, order> const &R)
{

  CPForceStructExp<N, dim, order> res;

  res.site = L.site + R.site;
  for (auto i = 0LU; i < dim; ++i) {
    res.links[i] = L.links[i] + R.links[i];
  }
  return res;
}


template <size_t N, size_t dim, size_t order>
CPForceStructExp<N, dim, order>
operator-(CPForceStructExp<N, dim, order> const &L,
          CPForceStructExp<N, dim, order> const &R)
{

  CPForceStructExp<N, dim, order> res;

  res.site = L.site - R.site;
  for (auto i = 0LU; i < dim; ++i) {
    res.links[i] = L.links[i] - R.links[i];
  }
  return res;
}

template <size_t N, size_t dim, size_t order>
CPForceStructExp<N, dim, order>
operator*(std::complex<double> const &scalar,
          CPForceStructExp<N, dim, order> const &cpdf)
{
  return cpdf * scalar;
}

// template<size_t order>
// auto operator+=( Expansion< std::complex<double>, order > &one,
// 		     Expansion< double, order > const &other)
// 	->  Expansion< std::complex<double>, order >
// {

// 	for( size_t i=0; i < order; ++i) {
// 	    one[i] += other[i];
// 	}
// 	return one;
// }


// template<size_t N, size_t order>
// auto operator*( Expansion<std::complex<double>,order> const &left, nummat<N>
// const & right)
// 	-> Expansion< nummat<N>, order >
// {
// 	Expansion< nummat<N>, order > res;

// 	for (auto i=0LU; i<order; ++i){
// 	    res[i]=left[i]*right;
// 	}

// 	return res;
// }

// template<size_t N, size_t order>
// auto operator*( Expansion<std::complex<double>,order> const &left,
// Expansion<nummat<N>,order> const & right)
// 	-> Expansion< nummat<N>, order >
// {
// 	Expansion< nummat<N>, order > res;

// 	for( size_t n = 0; n < order; ++n )
// 	{
// 	    for( size_t i = 0; i <= n; ++i )
// 	    {
// 		res[ n ] += left[ n - i ] * right[ i ];
// 	    }
// 	}

// 	return res;
// }

// template<size_t N, size_t order>
// auto operator*( Expansion<nummat<N>,order> const &left,
// Expansion<std::complex<double>,order> const & right)
// 	-> Expansion< nummat<N>, order >
// {
// 	return right*left;
// }

// template<size_t N, size_t order>
// auto operator*( Expansion<nummat<N>,order> const &left,
// Expansion<CP<std::complex<double>,N>,order> const & right)
// 	-> Expansion< CP<std::complex<double>,N>, order >
// {
// 	Expansion< CP<std::complex<double>,N>, order > res;

// 	for( size_t n = 0; n < order; ++n )
// 	{
// 	    for( size_t i = 0; i <= n; ++i )
// 	    {
// 		res[ n ] += left[ n - i ] * right[ i ];
// 	    }
// 	}

// 	return res;
// }


} // End of namespace "CanonicalCPNExpHelpers"


// --------------------------------------------------------------------------------
// Drift Force
// --------------------------------------------------------------------------------
namespace NoExpansion {

using namespace CanonicalCPNHelpers;

/// CP(N-1) Drift Force
template <size_t N, size_t dim, typename P = double, typename L = double>
class CPForce : public Force<FullLattice<CP<std::complex<double>, N>, U1, dim>,
                             std::vector<CPForceStruct<N, dim>>> {

 private:
  BoundaryCondition<P, dim> sbc;
  BoundaryCondition<L, dim> lbc;

 public:
  typedef CPForceStruct<N, dim> CPdf;
  typedef std::vector<CPdf> CPForceField;
  typedef FullLattice<CP<std::complex<double>, N>, U1, dim> CPlattice;


  /// Returns force field using standard boundary cond.
  CPForceField operator()(CPlattice const &lat) const
  {

    // Pre-allocate memory
    CPForceField FF;
    FF.reserve(lat.volume());


    for (auto it = lat.begin(); it != lat.end(); ++it) {


      // ################################################################################

      CPdf res;

      const auto II = std::complex<double>(0.0, 1.0);

      auto unit_matrix = cmplx(0.0);
      for (auto i = 0lu; i < dim; i++) {
        auto link = it.link(i).value();
        auto link_bwd = it.neighborLink(i, i, Direction::BACKWARD, lbc).value();

        auto field = it.site();

        auto field_fwd = it.neighborSite(i, Direction::FORWARD, sbc);
        auto field_bwd = it.neighborSite(i, Direction::BACKWARD, sbc);

        // CP field
        // ----------------------------------------------------------------------

        // First
        // -----------------------
        // External product part
        auto factor = -1. * link_bwd * II * 0.5;
        // std::cout << "nu: " << i << "\tfactor:" << factor << std::endl;
        external_product(factor * N, field, field_bwd, res.site);
        // Diagonal part
        unit_matrix += 2 * std::real(scalar_prod(field_bwd, field) * factor);
        // std::cout << "res:\n" << res.site <<std::endl;

        // Second (relative minus sign!)
        // -----------------------
        factor = link * II * 0.5;
        // std::cout << "nu: " << i << "\tfactor:" << factor << std::endl;
        // External product part
        external_product(factor * N, field_fwd, field, res.site);
        // Diagonal part
        auto s_f_ffwd = 2 * std::real(scalar_prod(field, field_fwd) * factor);
        // std::cout << "res:\n" << res.site <<std::endl;
        unit_matrix += s_f_ffwd;

        // Links
        // ----------------------------------------------------------------------
        // res.links[i] = -2*s_f_ffwd*NN;
        res.links[i] = 2 * N * std::imag(scalar_prod(field, field_fwd) * link);
        /* std::cout << 2*std::imag(link*s_f_ffwd) << std::endl; */
      }

      inflate(res.site);
      // std::cout << "Unit Matrix: " << unit_matrix << std::endl;
      res.site += (unit_matrix * nummat<N>(-1.0));

      FF.push_back(res);
      // ################################################################################
    }

    return FF;
  }


  /// Constructor
  explicit CPForce()
    : sbc(BoundaryCondition<P, dim>()), lbc(BoundaryCondition<L, dim>())
  {
  }
  explicit CPForce(const BoundaryCondition<P, dim> &_sbc,
                   const BoundaryCondition<L, dim> &_lbc)
    : sbc(_sbc), lbc(_lbc)
  {
  }
  ~CPForce() {}
  /// Compiler complains about constexpr here ???
  bool isGroupValued() const { return false; }
};
} // End of namespace "NoExpansion"


namespace Expan {

using namespace CanonicalCPNExpHelpers;
using namespace CommonCPNExpHelpers;
const double CP_Langevin_noise_scale = sqrt(2.);

/// CP(N-1) Drift Force
template <size_t N, size_t dim, size_t order, typename P = double,
          typename L = double>
class CPForce
  : public Force<FullLattice<Expansion<CP<std::complex<double>, N>, order>,
                             Expansion<U1, order>, dim>,
                 std::vector<CPForceStructExp<N, dim, order>>> {

 private:
  BoundaryCondition<P, dim> sbc;
  BoundaryCondition<L, dim> lbc;

 public:
  typedef CPForceStructExp<N, dim, order> CPdf;
  typedef std::vector<CPdf> CPForceField;
  typedef FullLattice<Expansion<CP<std::complex<double>, N>, order>,
                      Expansion<U1, order>, dim>
    CPlattice;


  /// Returns force field using standard boundary cond.
  CPForceField operator()(CPlattice const &lat) const
  {

    // Pre-allocate memory
    CPForceField FF;
    FF.reserve(lat.volume());


    for (auto it = lat.begin(); it != lat.end(); ++it) {


      // ################################################################################

      CPdf res;

      const auto II = std::complex<double>(0.0, 1.0);

      Expansion<std::complex<double>, order> unit_matrix(
        std::complex<double>(0.0));
      for (auto i = 0lu; i < dim; i++) {
        auto link = it.link(i);
        auto link_bwd = it.neighborLink(i, i, Direction::BACKWARD, lbc);

        auto field = it.site();

        auto field_fwd = it.neighborSite(i, Direction::FORWARD, sbc);
        auto field_bwd = it.neighborSite(i, Direction::BACKWARD, sbc);

        // CP field
        // ----------------------------------------------------------------------

        // First
        // -----------------------
        // External product part
        auto factor = (-1. * II * 0.5) * link_bwd;
        // std::cout << "nu: " << i << "\tfactor:" << factor << std::endl;
        external_product(factor * static_cast<double>(N), field, field_bwd,
                         res.site);
        // Diagonal part
        unit_matrix += 2. * real(scalar_prod(field_bwd, field) * factor);
        // std::cout << "res:\n" << res.site <<std::endl;

        // Second (relative minus sign!)
        // -----------------------
        factor = link * (II * 0.5);
        // std::cout << "nu: " << i << "\tfactor:" << factor << std::endl;
        // External product part
        external_product(factor * static_cast<double>(N), field_fwd, field,
                         res.site);
        // Diagonal part
        auto s_f_ffwd = 2 * real(scalar_prod(field, field_fwd) * factor);
        // std::cout << "res:\n" << res.site <<std::endl;
        unit_matrix += s_f_ffwd;

        // Links
        // ----------------------------------------------------------------------
        // res.links[i] = -2*s_f_ffwd*NN;
        res.links[i] = 2. * static_cast<double>(N) *
                       imag(scalar_prod(field, field_fwd) * link);
        /* std::cout << 2*std::imag(link*s_f_ffwd) << std::endl; */
      }

      for (auto i = 0LU; i < order; ++i) {
        inflate(res.site[i]);
      }

      // //std::cout << "Unit Matrix: " << unit_matrix << std::endl;
      res.site += (unit_matrix) * (nummat<N>(-1.0));

      FF.push_back(res);
      // ################################################################################
    }

    return FF;
  }

  ~CPForce(){};
  explicit CPForce()
    : sbc(BoundaryCondition<P, dim>()), lbc(BoundaryCondition<L, dim>())
  {
  }
  explicit CPForce(const BoundaryCondition<P, dim> &_sbc,
                   const BoundaryCondition<L, dim> &_lbc)
    : sbc(_sbc), lbc(_lbc){};

  /// Compiler complains about constexpr here ???
  bool isGroupValued() const { return false; }
};
} // End of namespace "Expan"
