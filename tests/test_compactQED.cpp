/// file
///
/// Tests related to the compact QED implementation

#include "actions/compact_qed.h"
#include "gaugegroups/u1.h"
#include "test_helpers.h"
#include <gtest/gtest.h>


TEST(compact_QED, CanInit)
{
  const std::size_t dim = 2;
  double beta = 0.987;
  std::array<size_t, dim> dimsar{ 2, 2 };
  U1 u1_init(1.0);

  LinkLattice<U1, dim> qplat(dimsar, u1_init);
  QEDAction<2> x;
  QEDAction<2> y(qplat, beta);
}

TEST(compact_QED, CanGetCorrectResultsForAction)
{
  const std::size_t dim = 3;
  double beta = 0.123456;
  std::array<size_t, dim> dimsar;
  dimsar.fill(2);
  U1 u1_init(std::complex<double>(1.0, 3.0), true);

  LinkLattice<U1, dim> qplat(dimsar, u1_init);
  QEDAction<dim> QED_act(qplat, beta);

  double mul = static_cast<double>((dim - 1) * 2);
  double plaq = std::real(std::pow(u1_init.value(), dim - 1) *
                          std::pow(std::conj(u1_init.value()), dim - 1));

  // EXPECT_DOUBLE_EQ(beta*mul*(1.-plaq),QED_act.atLink(0,0));
  EXPECT_NEAR(beta * mul * (1. - plaq), QED_act.atLink(0, 0), 1.e-15);
}


TEST(compact_QED, CanMC)
{
  const std::size_t dim = 4;
  const std::size_t nhit = 1;
  double beta = 2;
  std::array<size_t, dim> dimsar{ 4, 4, 4, 4 };
  U1 u1_init(1.0);

  LinkLattice<U1, dim> qplat(dimsar, u1_init);

  Compact_QED_update<dim> mc(&qplat, beta);
  QEDAction<dim> QED_act(qplat, beta);

  std::random_device rd;
  std::uniform_int_distribution<size_t> udint;
  // std::size_t udint_max=udint.max();
  std::vector<size_t> seeds;
  for (auto i = 0ul; i < std::mt19937_64::state_size; ++i) {
    seeds.push_back(udint(rd));
  }

  std::seed_seq sseq(seeds.begin(), seeds.end());
  std::mt19937_64 generator(sseq);

  // Make sure overrelaxation does not change action
  mc.multihit_MC(generator, nhit, 0.6);
  auto init = QED_act.energyDensity();
  mc.overrelaxation();
  auto final = QED_act.energyDensity();

  EXPECT_NEAR(init, final, 1.e-14);
  // for (auto i=0; i<100000; ++i){

  // 	auto acc = mc.multihit_MC(generator,nhit,0.6);

  // 	if (i%1000 == 0){
  // 	    std::cout << "\t"
  // 	    	      <<
  // static_cast<double>(acc)/(static_cast<double>(nhit*qplat.volume()*dim))
  // 	    	      << "\t" << QED_act.energyDensity();

  // 	    mc.overrelaxation();
  // 	    std::cout << "\t" << QED_act.energyDensity() << std::endl;
  // 	}
  // }
}
