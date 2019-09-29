#include "actions/principal_chiral_model.h"
#include <array>
#include <complex>
#include <gtest/gtest.h>


TEST(PCM, CanInit)
{
  const std::size_t dim = 2;
  const std::size_t NN = 1;
  std::array<size_t, dim> dimsar{ 2, 2 };
  nummat<NN> init(1.0);
  SiteLattice<nummat<NN>, dim> pcmplat(dimsar, init);
  PCMAction<std::complex<double>, 2, 2> x;
}


TEST(PCM, CanGetCorrectResultsForAction)
{
  const std::size_t dim = 2;
  const std::size_t NN = 6;
  double beta = 0.123456;
  std::array<size_t, dim> dimsar{ 2, 2 };

  nummat<NN> init(1.0);
  SiteLattice<nummat<NN>, dim> pcmlat(dimsar, init);
  PCMAction<std::complex<double>, NN, dim> pcmact(pcmlat, beta);

  EXPECT_DOUBLE_EQ(
    static_cast<double>(pcmlat.volume() * pcmlat.dimensions() * NN),
    pcmact.total());
}


TEST(PCM, CanLeftMultiplyWithPCMtwistPhase)
{
  // const std::size_t dim=2;
  const std::size_t NN = 3;
  auto II = std::complex<double>(0., 1.);

  std::array<std::complex<double>, NN *NN> init_ar = { 1.0, 1.0, 1.0, 1.0, 1.0,
                                                       1.0, 1.0, 1.0, 1.0 };
  std::array<std::complex<double>, NN *NN> final_ar = {
    4., 6., 10. * II, 6., 9., 15 * II, -10. * II, -15 * II, 25.
  };

  std::array<std::complex<double>, NN> phase_ar = { 2 * II, 3 * II, 5 };
  nummat<NN> init(init_ar);
  nummat<NN> final(final_ar);


  PCMtwistPhase<std::complex<double>, NN> phase(phase_ar);

  auto tst = phase * init;

  EXPECT_EQ(tst, final);
}


TEST(PCM, CanRightMultiplyWithPCMtwistPhase)
{
  // const std::size_t dim=2;
  const std::size_t NN = 3;
  auto II = std::complex<double>(0., 1.);

  std::array<std::complex<double>, NN *NN> init_ar = { 1.0, 1.0, 1.0, 1.0, 1.0,
                                                       1.0, 1.0, 1.0, 1.0 };
  std::array<std::complex<double>, NN *NN> final_ar = {
    4., 6., -10. * II, 6., 9., -15 * II, 10. * II, 15 * II, 25.
  };

  std::array<std::complex<double>, NN> phase_ar = { 2 * II, 3 * II, 5 };
  nummat<NN> init(init_ar);
  nummat<NN> final(final_ar);


  PCMtwistPhase<std::complex<double>, NN> phase(phase_ar);

  auto tst = init * phase;

  EXPECT_EQ(tst, final);
}


TEST(PCM, CanLeftMultiplyWithPCMtwistMatrix)
{
  // const std::size_t dim=2;
  const std::size_t NN = 3;
  auto II = std::complex<double>(0., 1.);

  std::array<std::complex<double>, NN *NN> init_ar = { 1.0, 1.0, 1.0, 1.0, 1.0,
                                                       1.0, 1.0, 1.0, 1.0 };
  std::array<std::complex<double>, NN *NN> final_ar = {
    4., 6., 10. * II, 6., 9., 15 * II, -10. * II, -15 * II, 25.
  };

  std::array<std::complex<double>, NN> phase_ar = { 2 * II, 3 * II, 5 };
  nummat<NN> init(init_ar);
  nummat<NN> final(final_ar);


  PCMtwistMatrix<NN> mat(phase_ar);

  auto tst = mat * init;

  EXPECT_EQ(tst, final);
}


TEST(PCM, CanRightMultiplyWithPCMtwistMatrix)
{
  // const std::size_t dim=2;
  const std::size_t NN = 3;
  auto II = std::complex<double>(0., 1.);

  std::array<std::complex<double>, NN *NN> init_ar = { 1.0, 1.0, 1.0, 1.0, 1.0,
                                                       1.0, 1.0, 1.0, 1.0 };
  std::array<std::complex<double>, NN *NN> final_ar = {
    4., 6., -10. * II, 6., 9., -15 * II, 10. * II, 15 * II, 25.
  };

  std::array<std::complex<double>, NN> phase_ar = { 2 * II, 3 * II, 5 };
  nummat<NN> init(init_ar);
  nummat<NN> final(final_ar);


  PCMtwistMatrix<NN> mat(phase_ar);

  auto tst = init * mat;

  EXPECT_EQ(tst, final);
}
