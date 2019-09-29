/// file
///
/// Tests related to random number generation

#include "randgen.h"
#include <gtest/gtest.h>
#include <vector>


TEST(RANDGEN, Construction)
{
  std::normal_distribution<double> nd;

  /// Standard Ctor should use hardware random number generator
  {
    rnd_base<std::ranlux24, decltype(nd), 48LU> rnd1, rnd2;
    EXPECT_FALSE(rnd1 == rnd2);
  }


  /// Two generators with same initialisation should be the same
  {
    std::vector<std::size_t> seeds;
    for (auto i = 0ul; i < 48LU; ++i) {
      seeds.push_back(i + 1);
    }
    std::seed_seq sseq(seeds.begin(), seeds.end());

    rnd_base<std::ranlux24, decltype(nd), 48LU> rnd1(sseq), rnd2(sseq);
    EXPECT_TRUE(rnd1 == rnd2);
    EXPECT_DOUBLE_EQ(rnd1(), rnd2());
  }
}

TEST(RANDGEN, Equality)
{
  /// Two generators with same state but different engine or dist should
  //  NOT be the same
  std::vector<std::size_t> seeds;
  for (auto i = 0ul; i < 48LU; ++i) {
    seeds.push_back(i + 1);
  }
  std::seed_seq sseq(seeds.begin(), seeds.end());

  std::uniform_real_distribution<double> ud;
  std::cauchy_distribution<double> cd;
  rnd_base<std::ranlux24, decltype(cd), 48LU> rnd1(sseq);
  rnd_base<std::ranlux24, decltype(ud), 48LU> rnd2(sseq);
  rnd2.set_state(rnd1.get_state());

  // std::cout << rnd1.get_state()[0] << "\n"  << rnd1.get_state()[1] << "\n"
  //  	      << rnd2.get_state()[0] << "\n"  << rnd2.get_state()[1] << "\n";

  EXPECT_FALSE(rnd1 == rnd2);
  EXPECT_TRUE(rnd1.get_state() == rnd2.get_state());

  /// Two generators with same engine and dist should but different state
  //  should NOT be the same
  rnd_base<std::ranlux24, decltype(cd), 48LU> rnd3(sseq);

  EXPECT_TRUE(rnd1 == rnd3);
  rnd1();
  EXPECT_FALSE(rnd1 == rnd3);
  rnd3();
  EXPECT_TRUE(rnd1 == rnd3);
}


TEST(RANDGEN, CanRestoreState)
{

  std::knuth_b engine;
  std::gamma_distribution<double> gd;


  rnd_base<decltype(engine), decltype(gd)> rnd1, rnd2;


  for (auto i = 0ul; i < 1000LU; ++i) {
    rnd1();
  }

  rnd2.set_state(rnd1.get_state());

  bool test = true;
  for (auto i = 0ul; i < 1000LU; ++i) {
    test = (test && rnd1() == rnd2());
    if (not test) {
      break;
    }
  }

  EXPECT_TRUE(test);
}
