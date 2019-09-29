/// file
///
/// Tests for the statistical functions in stat.h

#include "stat.h"
#include <gtest/gtest.h>
#include <valarray>
#include <vector>


TEST(STAT, CanComputeRunningStatistics)
{

  RunningStatistics<double> rs;


  double max = 1000.;
  for (auto x = 1.0; x <= max; ++x) {

    rs.push(x);
  }

  EXPECT_NEAR(rs.mean(), (max + 1.) / 2., 1.e-13 * max);
  EXPECT_NEAR(rs.variance(), (max + 1.) * max / 12., 1.e-13 * max);
  EXPECT_NEAR(rs.std() * rs.std(), (max + 1.) * max / 12., 1.e-13 * max);
}


TEST(STAT, CanComputeRunningStatisticsArray)
{
  RunningStatisticsArray<double, 3L> rsa;

  double max = 1000.;
  for (auto x = 1.0; x <= max; ++x) {

    std::valarray<double> X(x, 3);
    rsa.push(X);
  }

  for (auto i = 0LU; i < 3LU; ++i) {
    EXPECT_NEAR(rsa.mean()[i], (max + 1.) / 2., 1.e-13 * max);
    EXPECT_NEAR(rsa.variance()[i], (max + 1.) * max / 12., 1.e-13 * max);
    EXPECT_NEAR(rsa.std()[i] * rsa.std()[i], (max + 1.) * max / 12.,
                1.e-13 * max);
  }
}
