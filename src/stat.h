/// \file
/// Definitions for a few useful statistical functions


#pragma once
#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <valarray>
#include <vector>

/// Computes the mean over a vector
template <typename T>
T vec_mean(std::vector<T> const &vec)
{
  T sum = std::accumulate(vec.begin(), vec.end(), T(0));
  return sum / static_cast<T>(vec.size());
}

/// Computes the standard error of the mean over a vector
template <typename T>
T vec_std_err(std::vector<T> const &vec)
{
  T mean = vec_mean(vec);
  T stdev = T(0);
  std::for_each(vec.begin(), vec.end(), [&stdev, &mean](T const &elem) {
    stdev += (elem - mean) * (elem - mean);
  });

  T NN = static_cast<T>(vec.size());
  return std::sqrt(stdev / ((NN - 1) * NN));
}

/// Compute expectation value of x^2 for a vector of x values
template <typename T>
T vec_exp2(std::vector<T> const &vec)
{
  T exp2 = T(0);
  std::for_each(vec.begin(), vec.end(),
                [&exp2](T const &elem) { exp2 += elem * elem; });

  T NN = static_cast<T>(vec.size());
  return exp2 / NN;
}

/// Compute standard error for expectation value of x^2 for a vector of x values
template <typename T>
T vec_exp2_std_err(std::vector<T> const &vec)
{
  T exp2 = vec_exp2(vec);
  T stdev = T(0);

  std::for_each(vec.begin(), vec.end(), [&stdev, &exp2](T const &elem) {
    stdev += (elem * elem - exp2 * exp2) * (elem * elem - exp2 * exp2);
  });

  T NN = static_cast<T>(vec.size());
  return std::sqrt(stdev / ((NN - 1) * NN));
}


/// Compute mean of "vector observables"
/// Input is given by a vector of vectors (which should all have the same
/// size S). Output is a vector of size S which has as i-th entry the mean over
/// the i-th entries of the provided vectors.
template <typename T>
std::vector<T> vec_vec_mean(std::vector<std::vector<T>> const &vec)
{

  auto length = vec[0].size();
  auto mean_vec = vec[0];


  for (auto i = 1LU; i < vec.size(); ++i) {
    if (length != vec[i].size()) {
      throw std::runtime_error("Wrong input for vec_vec_mean."
                               " Vectors have different length.");
    }

    for (auto j = 0LU; j < length; ++j) {
      mean_vec[j] += (vec[i])[j];
    }
  }

  // Normalise
  for (auto j = 0LU; j < length; ++j) {
    mean_vec[j] /= static_cast<T>(vec.size());
  }

  return mean_vec;
}


/// Computes standard error for vector of vectors. See vec_vec_mean for details
template <typename T>
std::vector<T> vec_vec_std_err(std::vector<std::vector<T>> const &vec)
{

  auto mean_vec = vec_vec_mean(vec);
  auto length = vec[0].size();
  std::vector<T> std_err_vec(length, T(0.0));

  for (auto i = 0LU; i < vec.size(); ++i) {
    if (length != vec[i].size()) {
      throw std::runtime_error("Wrong input for vec_vec_mean."
                               " Vectors have different length.");
    }

    for (auto j = 0LU; j < length; ++j) {
      auto fact = (vec[i])[j] - mean_vec[j];
      fact *= fact;

      std_err_vec[j] += fact;
    }
  }

  auto NN = static_cast<T>(vec.size());

  // Normalise
  for (auto j = 0LU; j < length; ++j) {
    std_err_vec[j] = std::sqrt(std_err_vec[j] / ((NN - 1) * NN));
  }

  return std_err_vec;
}


/// Bootstrap mean and error
template <typename T, typename Generator>
std::vector<T> bootstrap(std::vector<T> const &data, std::size_t bs_samples,
                         Generator &gen)
{

  std::uniform_int_distribution<std::size_t> dist(0, data.size() - 1);

  const double percentile = (1. - 0.682689492137086) * 0.5;

  std::vector<T> res_vec;
  std::vector<T> bs_vec;
  bs_vec.resize(data.size());

  for (auto i = 0LU; i < bs_samples; ++i) {

    // Get BS sample
    for (auto j = 0LU; j < data.size(); ++j) {
      bs_vec[j] = data[dist(gen)];
    }

    // Compute mean
    res_vec.push_back(vec_mean(bs_vec));
    // std::cout << res_vec.back() << "\n";
  }

  // Sorting res_vec to compute percentiles  (modeled by a normal distribution
  // where ~ 68% of the values lie within one standard deviation from the mean)
  std::sort(res_vec.begin(), res_vec.end(), std::less<T>());

  auto mean = vec_mean(res_vec);
  auto idx = static_cast<std::size_t>(
    std::trunc(percentile * (static_cast<double>(res_vec.size()))));

  auto std_err = res_vec[res_vec.size() - idx] - res_vec[idx];
  std_err *= 0.5;


  res_vec.clear();
  res_vec.push_back(mean);
  res_vec.push_back(std_err);

  return res_vec;
}


template <typename T>
std::vector<T> binning(std::vector<T> const &data, std::size_t bin_size)
{

  auto orig_size = data.size();

  if (bin_size > orig_size) {
    throw std::runtime_error("Bin size larger than original vector.");
  }

  std::vector<T> res_vec;

  // Take into account that bin_size might not divide the vector size
  auto last = orig_size / bin_size;
  last *= bin_size;
  for (auto i = 0LU; i < orig_size - bin_size; i += bin_size) {
    T res = T(0.0);
    for (auto j = 0LU; j < bin_size; ++j) {
      res += data[i + j];
    }
    res_vec.push_back(res / static_cast<double>(bin_size));
  }

  // Simply put the remaining values in the last bin
  T res = T(0.0);
  auto cnt = 1LU;
  for (auto j = last; j < orig_size; ++j) {
    res += data[j];
    ++cnt;
  }
  res_vec.push_back(res / static_cast<double>(cnt));


  return res_vec;
}


/// \brief A class for running statistics calculations
///
/// Computes the cumulative running mean, variance and standard deviation of a
/// stream of data. The data is not saved and the statistical information is
/// updated 'on the fly'.
///
/// See:
/// Welford, B. P., Technometrics 3 (1962), pp 419--412
/// doi:10.1080/00401706.1962.10490022
/// https://amstat.tandfonline.com/doi/pdf/10.1080/00401706.1962.10490022
template <typename T>
class RunningStatistics {

 private:
  T m;               // Running mean
  T var;             // Running variance
  std::size_t count; // Current sample size

 public:
  RunningStatistics() : m(T(0)), var(T(0)), count(0LU) {}

  void push(T const &x)
  {

    count++;
    double NN = static_cast<double>(count);

    auto new_mean = (NN > 1) ? ((NN - 1.) * m + x) / NN : x;
    // Variance is zero for sample size one
    auto new_var = (NN > 1) ? var + (x - m) * (x - new_mean) : T(0);

    m = new_mean;
    var = new_var;
  }

  T mean() const { return m; }

  T variance() const
  {
    return count > 1 ? var / (static_cast<double>(count) - 1.) : 0.0;
  }

  T std() const { return std::sqrt(this->variance()); }
};


/// \brief A class for running statistics calculations for "array observables"
template <typename T, size_t len>
class RunningStatisticsArray {
 private:
  std::valarray<T> m;   // Running mean
  std::valarray<T> var; // Running variance
  std::size_t count;    // Current sample size

 public:
  RunningStatisticsArray() : m(T(0), len), var(T(0), len), count(0LU) {}

  void push(std::valarray<T> const &x)
  {

    count++;
    double NN = static_cast<double>(count);

    auto new_mean = (NN > 1) ? ((NN - 1.) * m + x) / NN : x;
    // Variance is zero for sample size one
    auto new_var =
      (NN > 1) ? var + (x - m) * (x - new_mean) : std::valarray<T>(0.0, len);

    m = new_mean;
    var = new_var;
  }

  std::valarray<T> mean() const { return m; }

  std::valarray<T> variance() const
  {
    return count > 1 ? var / (static_cast<double>(count) - 1.)
                     : std::valarray<T>(0.0, len);
  }

  std::valarray<T> std() const { return std::sqrt(this->variance()); }
};
