/// \file
///
/// \brief Tests for the lattice code


#include "expansion.h"
#include "gaugegroups/su3.h"
#include "lattice.h"
#include <array>
#include <gtest/gtest.h>
#include <iostream>

/// Define constant for EXPECT_NEAR tests
const double near_eps = 1.e-15;


TEST(SiteLattice, GivesCorrectVolume)
{
  SiteLattice<int, 4> lat({ 3, 4, 5, 6 });
  EXPECT_EQ(lat.volume(), 3 * 4 * 5 * 6ul);
}

TEST(SiteLattice, GivesCorrectDimension)
{
  SiteLattice<int, 3> lat({ 3, 5, 6 });
  EXPECT_EQ(lat.dimensions(), 3ul);
}

TEST(SiteLattice, ConstructionIsInitialization)
{
  SiteLattice<int, 4> lat({ 3, 3, 3, 3 }, 2);
  EXPECT_EQ(lat.volume(), 3 * 3 * 3 * 3ul);
  // expect that there is a two at every lattice point:
  for (std::size_t i = 0; i < lat.volume(); ++i) {
    EXPECT_EQ(lat[i], 2);
  }
}

TEST(SiteLattice, DefaultConstructionWorks)
{
  SiteLattice<double, 4> lat;
  EXPECT_EQ(lat.volume(), 0ul);
}

TEST(SiteLattice, CopyConstructionWorks)
{
  SiteLattice<double, 4> lat({ 5, 2, 4, 3 }, 2);
  SiteLattice<double, 4> copy(lat);

  for (std::size_t i = 0; i < lat.volume(); ++i) {
    EXPECT_EQ(lat[i], copy[i]);
  }
}

TEST(SiteLattice, AssignmentWorks)
{
  SiteLattice<double, 4> lat({ 5, 2, 4, 3 }, 2);
  SiteLattice<double, 4> copy;

  copy = lat;

  for (std::size_t i = 0; i < lat.volume(); ++i) {
    EXPECT_EQ(lat[i], copy[i]);
  }
}


TEST(SiteLattice, CanUseIteratorToRunLinearlyOverLattice)
{
  SiteLattice<int, 3> lat({ 2, 3, 4 });
  int i = 0;

  // Initialise lattice
  for (auto it = lat.begin(); it != lat.end(); ++it) {
    *it = (++i);
  }

  // Test indexing
  i = 0;
  for (auto it = lat.begin(); it != lat.end(); ++it) {
    EXPECT_EQ(*it, ++i);
  }

  // Test for correct values at given site
  for (std::size_t j = 0; j < lat.volume(); ++j) {
    EXPECT_EQ(static_cast<size_t>(lat[j]), (j + 1));
  }
}

TEST(SiteLattice, CanUseIteratorToRunOverCartesianLattice)
{
  constexpr std::size_t Nx{ 5 }, Ny{ 7 };
  SiteLattice<int, 2> lat({ Nx, Ny });
  std::size_t x(0), y(0);
  auto it = lat.begin();
  // this should run linearly over the lattice:
  auto ind = 0;
  for (; y < Ny; ++y) {
    x = 0;
    for (; x < Nx; ++x) {
      *it = static_cast<int>(x * 1000 + y);
      EXPECT_EQ(it.index(), static_cast<size_t>(ind));
      it.advanceInDirection(0, 1);
      ind++;
    }
    it.advanceInDirection(1, 1);
  }
  // after running over the lattice, ind  should now be the volume:
  EXPECT_EQ(static_cast<size_t>(ind), lat.volume());
  x = 0;
  y = 0;
  for (auto const &elem : lat) {
    EXPECT_EQ(x * 1000 + y, static_cast<size_t>(elem));
    x++;
    if (x == Nx) {
      x = 0;
      y++;
    }
  }
}
TEST(SiteLattice, CanUseIteratorToRunOverSlicesOfLattice)
{
  constexpr std::size_t Nx{ 5 }, Ny{ 7 };
  SiteLattice<int, 2> lat({ Nx, Ny }, 0);
  auto it = lat.begin();
  // this should wrap around the x direction and inc the value by one.
  for (std::size_t x = 0; x < 2 * Nx + 3; ++x) {
    *it += 1;
    it.advanceInDirection(0, 1);
  }
  // after running over the x slices 2times + 3, the values of the xslice at y=0
  // should be
  // 3 for x <= 3, 2 for x = 4;
  for (it = lat.begin(); it != lat.end(); ++it) {
    auto pos = it.coord();
    if (pos[0] < 3 and pos[1] == 0)
      EXPECT_EQ(*it, 3);
    else if (pos[0] >= 3 and pos[1] == 0)
      EXPECT_EQ(*it, 2);
    else
      EXPECT_EQ(*it, 0);
  }
}
TEST(SiteLattice, MultiSiteHoppingWorks)
{
  constexpr std::size_t Nx{ 5 }, Ny{ 7 };
  SiteLattice<int, 2> lat({ Nx, Ny }, 0);
  auto it = lat.begin();
  std::advance(it, 10); // go somewhere else.

  auto otherit = it; // copy

  // now, advance the iterator in y-dir for 2 hops:
  it.advanceInDirection(1, 2);

  // is it the same as advancing twice?
  otherit.advanceInDirection(1, 1);
  otherit.advanceInDirection(1, 1);

  EXPECT_TRUE(it == otherit);
}
TEST(SiteLattice, LatticeSum)
{
  SiteLattice<int, 4> lat({ 2, 3, 4, 5 }, 2.);
  int s = sum(lat);
  EXPECT_EQ(s, 2 * 3 * 4 * 5 * 2);
}
TEST(SiteLattice, ManipulateOnePointInLattice)
{
  SiteLattice<int, 4> lat({ 2, 3, 4, 5 }, 2.);
  *(lat.at({ 1, 2, 3, 4 })) = 123;
  // volume-1 points are 2, one point is 123:
  EXPECT_EQ(static_cast<size_t>(sum(lat)), 123 + (lat.volume() - 1) * 2);
  // check that the first and last point is correct:
  EXPECT_EQ(lat.at({ 0, 0, 0, 0 }).index(), 0ul);
  EXPECT_EQ(lat.at({ 1, 2, 3, 4 }).index(), lat.volume() - 1);
  // check that this is indeed the last point, i.e. that the next point is
  // end():
  EXPECT_TRUE((++(lat.at({ 1, 2, 3, 4 }))) == lat.end());
}

TEST(SiteLattice, CanObeyPeriodicBoundaryCond)
{
  BoundaryCondition<double, 2> pbc({ 1., 1. }, "PERIODIC");
  SiteLattice<std::complex<double>, 2> lat({ 3, 2 }, 1.);
  *lat.at({ 1, 0 }) = 2.;
  *lat.at({ 1, 1 }) = 2.;
  *lat.at({ 2, 0 }) = 3.;
  *lat.at({ 2, 1 }) = 3.;

  SiteLatticeIterator<SiteLattice<std::complex<double>, 2>> sli(0, lat);

  auto idx = lat.coordToLinearIndex({ 0, 0 });
  auto nidx = lat.coordToLinearIndex({ 1, 0 });

  EXPECT_EQ(lat[nidx], sli.neighborSite(0, Direction::FORWARD, pbc));
  EXPECT_EQ((*lat.at({ 2, 0 })), sli.neighborSite(0, Direction::BACKWARD));

  idx = lat.coordToLinearIndex({ 2, 1 });
  sli = SiteLatticeIterator<SiteLattice<std::complex<double>, 2>>(idx, lat);

  EXPECT_EQ(*lat.at({ 0, 1 }), sli.neighborSite(0, Direction::FORWARD, pbc));
  EXPECT_EQ(*lat.at({ 2, 0 }), sli.neighborSite(1, Direction::FORWARD, pbc));
}


TEST(SiteLattice, CanObeyAntiPeriodicBoundaryCond)
{

  BoundaryCondition<double, 2> abc({ -1, -1 }, "ANTI-PERIODIC");
  SiteLattice<std::complex<double>, 2> lat({ 3, 2 }, 1.);
  *lat.at({ 1, 0 }) = 2.;
  *lat.at({ 1, 1 }) = 2.;
  *lat.at({ 2, 0 }) = 3.;
  *lat.at({ 2, 1 }) = 3.;

  SiteLatticeIterator<SiteLattice<std::complex<double>, 2>> sli(0, lat);

  auto idx = lat.coordToLinearIndex({ 0, 0 });
  auto nidx = lat.coordToLinearIndex({ 1, 0 });

  EXPECT_EQ(lat[nidx], sli.neighborSite(0, Direction::FORWARD, abc));
  EXPECT_NEAR(std::abs((*lat.at({ 2, 0 })) +
                       sli.neighborSite(0, Direction::BACKWARD, abc)),
              0., near_eps);

  idx = lat.coordToLinearIndex({ 2, 1 });
  sli = SiteLatticeIterator<SiteLattice<std::complex<double>, 2>>(idx, lat);
  EXPECT_NEAR(
    std::abs(*lat.at({ 0, 1 }) + sli.neighborSite(0, Direction::FORWARD, abc)),
    0.0, near_eps);


  EXPECT_NEAR(
    std::abs(*lat.at({ 2, 0 }) + sli.neighborSite(1, Direction::FORWARD, abc)),
    0.0, near_eps);
}


// TEST(SiteLattice, CanObeyPhaseBoundaryCond){

//     double phase = std::atan(1.)*0.3;
//     auto II      = std::complex<double>(0.,1.);
//     auto bwd     = std::exp(II*phase);
//     auto fwd     = std::exp(-II*phase);

//     BoundaryCondition<double, 2> pbc({phase,phase},"PHASE");
//     SiteLattice<std::complex<double>, 2> lat({3,2}, 1.);
//      *lat.at({1,0})=2.;
//      *lat.at({1,1})=2.;
//      *lat.at({2,0})=3.;
//      *lat.at({2,1})=3.;

//      auto  idx = lat.coordToLinearIndex({0,0});
//      auto nidx = lat.coordToLinearIndex({1,0});

//      SiteLatticeIterator< SiteLattice<std::complex<double>, 2>> sli(idx,lat);

//      EXPECT_EQ(lat[nidx],sli.neighborSite(0,Direction::FORWARD,pbc));
//      EXPECT_EQ(bwd*(*lat.at({2,0})),sli.neighborSite(0,Direction::BACKWARD,pbc));

//      idx = lat.coordToLinearIndex({2,1});
//      sli = SiteLatticeIterator< SiteLattice<std::complex<double>,
//      2>>(idx,lat);

//      EXPECT_EQ(fwd*(*lat.at({0,1})),sli.neighborSite(0,Direction::FORWARD,pbc));
//      EXPECT_EQ(*lat.at({2,0})*fwd,
//      sli.neighborSite(1,Direction::FORWARD,pbc));
// }


TEST(LinkLattice, GivesCorrectVolume)
{
  LinkLattice<int, 4> lat({ 3, 4, 5, 6 });
  EXPECT_EQ(lat.volume(), 3 * 4 * 5 * 6ul);
}

TEST(LinkLattice, GivesCorrectDimension)
{
  LinkLattice<int, 3> lat({ 3, 5, 6 });
  EXPECT_EQ(lat.dimensions(), 3ul);
}

TEST(LinkLattice, ConstructionWithLinkIsInitialization)
{
  LinkLattice<int, 4> lat({ 3, 3, 3, 3 }, 2);
  EXPECT_EQ(lat.volume(), 3 * 3 * 3 * 3ul);
  std::array<int, 4> ar;
  ar.fill(2);
  // expect that all links are set to 2:
  for (std::size_t i = 0; i < lat.volume(); ++i) {
    EXPECT_EQ(lat[i], ar);
  }
}

TEST(LinkLattice, ConstructionWithLinkFieldIsInitialization)
{

  const std::array<int, 4> ar = { 1, 2, 3, 4 };
  LinkLattice<int, 4> lat({ 3, 3, 3, 3 }, ar);
  EXPECT_EQ(lat.volume(), 3 * 3 * 3 * 3ul);
  // expect that all links are set to ar:
  for (std::size_t i = 0; i < lat.volume(); ++i) {
    EXPECT_EQ(lat[i], ar);
  }
}


TEST(LinkLattice, DefaultConstructionWorks)
{
  LinkLattice<double, 4> lat;
  EXPECT_EQ(lat.volume(), 0ul);
}

TEST(LinkLattice, CopyConstructionWorks)
{
  LinkLattice<double, 4> lat({ 5, 2, 4, 3 }, 2.);
  LinkLattice<double, 4> copy(lat);

  for (std::size_t i = 0; i < lat.volume(); ++i) {
    EXPECT_EQ(lat[i], copy[i]);
  }
}

TEST(LinkLattice, AssignmentWorks)
{
  LinkLattice<double, 4> lat({ 5, 2, 4, 3 }, { 2., 7., 9., 23. });
  LinkLattice<double, 4> copy;

  copy = lat;

  for (std::size_t i = 0; i < lat.volume(); ++i) {
    EXPECT_EQ(lat[i], copy[i]);
  }
}


TEST(LinkLattice, CanAccessLinkArray)
{
  const std::array<int, 4> a = { 1, 2, 3, 4 }, b = { 2, 3, 4, 5 },
                           c = { -1, -1, -1, -1 };
  LinkLattice<int, 4> lat({ 2, 2, 2, 2 }, 2.);

  lat[0] = a;
  lat[1] = b;
  EXPECT_EQ(lat[0], a);
  EXPECT_EQ(lat[1], b);

  // lat[0]-b should equal c
  std::transform(lat[0].begin(), lat[0].end(), b.begin(), lat[0].begin(),
                 std::minus<int>());
  EXPECT_EQ(lat[0], c);
}

TEST(LinkLattice, CanAccessLink)
{

  const std::array<int, 4> a = { 0, 0, 0, 0 }, b = { 2, 4, 6, 8 };
  LinkLattice<int, 4> lat({ 2, 2, 2, 2 }, 2.);

  lat[0] = a;
  for (std::size_t i = 0; i < 4; i++) {
    lat(0, i) = 2 * (static_cast<int>(i) + 1);

    EXPECT_EQ(lat(0, i), b[i]);
  }

  EXPECT_EQ(lat[0], b);
}

TEST(LinkLattice, CanUseIteratorToRunLinearlyOverLinkLattice)
{
  LinkLattice<int, 3> lat({ 2, 3, 4 });
  int i = 0;
  std::array<int, 3> a = { 0, 0, 0 };

  // Initialise lattice
  for (auto it = lat.begin(); it != lat.end(); ++it) {
    a.fill(++i);
    *it = a;
  }

  // Test indexing
  i = 0;
  for (auto it = lat.begin(); it != lat.end(); ++it) {
    std::array<int, 3> b;
    b.fill(++i);
    EXPECT_EQ(*it, b);
  }

  // Test for correct values at given site
  for (std::size_t j = 0; j < lat.volume(); ++j) {
    std::array<int, 3> b;
    b.fill(static_cast<int>(j) + 1);
    EXPECT_EQ(lat[j], b);
  }
}


TEST(LinkLattice, CanObeyPeriodicBoundaryCond)
{
  BoundaryCondition<double, 2> pbc({ 1., 1. }, "PERIODIC");
  std::array<std::complex<double>, 2> init_ar;
  LinkLattice<std::complex<double>, 2> lat({ 3, 2 }, 1.);
  init_ar.fill(2.);
  *lat.at({ 1, 0 }) = init_ar;
  init_ar.fill(2.);
  *lat.at({ 1, 1 }) = init_ar;
  init_ar.fill(3.);
  *lat.at({ 2, 0 }) = init_ar;
  init_ar.fill(3.);
  *lat.at({ 2, 1 }) = init_ar;

  LinkLatticeIterator<LinkLattice<std::complex<double>, 2>> lli(lat);

  // No Boundary crossed
  auto nidx = lat.coordToLinearIndex({ 1, 0 });

  EXPECT_EQ(lat(nidx, 0), lli.neighborLink(0, 0, Direction::FORWARD, pbc));
  EXPECT_EQ(lat(nidx, 1), lli.neighborLink(0, 1, Direction::FORWARD, pbc));


  // 0-dir. boundary crossed in forward dir
  auto idx = lat.coordToLinearIndex({ 2, 0 });
  nidx = lat.coordToLinearIndex({ 0, 0 });
  lli = LinkLatticeIterator<LinkLattice<std::complex<double>, 2>>(idx, lat);

  EXPECT_EQ(lat(nidx, 0), lli.neighborLink(0, 0, Direction::FORWARD, pbc));

  // 1-dir. boundary crossed in forward dir
  idx = lat.coordToLinearIndex({ 1, 1 });
  nidx = lat.coordToLinearIndex({ 1, 0 });
  lli = LinkLatticeIterator<LinkLattice<std::complex<double>, 2>>(idx, lat);

  EXPECT_EQ(lat(nidx, 1), lli.neighborLink(1, 1, Direction::FORWARD, pbc));
  // Cross check with no crossing
  EXPECT_EQ(lat(nidx, 1), lli.neighborLink(1, 1, Direction::BACKWARD, pbc));

  // 1-dir. boundary crossed in backward dir
  idx = lat.coordToLinearIndex({ 0, 0 });
  nidx = lat.coordToLinearIndex({ 2, 0 });
  lli = LinkLatticeIterator<LinkLattice<std::complex<double>, 2>>(idx, lat);
  EXPECT_EQ(lat(nidx, 1), lli.neighborLink(0, 1, Direction::BACKWARD, pbc));
}

TEST(LinkLattice, CanObeyAntiPeriodicBoundaryCond)
{

  BoundaryCondition<double, 2> abc({ -1, -1 }, "ANTIPERIODIC");
  std::array<std::complex<double>, 2> init_ar;
  LinkLattice<std::complex<double>, 2> lat({ 3, 2 }, 1.);
  init_ar.fill(2.);
  *lat.at({ 1, 0 }) = init_ar;
  init_ar.fill(2.);
  *lat.at({ 1, 1 }) = init_ar;
  init_ar.fill(3.);
  *lat.at({ 2, 0 }) = init_ar;
  init_ar.fill(3.);
  *lat.at({ 2, 1 }) = init_ar;

  LinkLatticeIterator<LinkLattice<std::complex<double>, 2>> lli(lat);


  // No Boundary crossed
  auto nidx = lat.coordToLinearIndex({ 1, 0 });

  EXPECT_EQ(lat(nidx, 0), lli.neighborLink(0, 0, Direction::FORWARD, abc));
  EXPECT_EQ(lat(nidx, 1), lli.neighborLink(0, 1, Direction::FORWARD, abc));


  // 0-dir. boundary crossed in forward dir
  auto idx = lat.coordToLinearIndex({ 2, 0 });
  nidx = lat.coordToLinearIndex({ 0, 0 });
  lli = LinkLatticeIterator<LinkLattice<std::complex<double>, 2>>(idx, lat);

  EXPECT_NEAR(
    std::abs(lat(nidx, 0) + lli.neighborLink(0, 0, Direction::FORWARD, abc)),
    0.0, near_eps);
  EXPECT_NEAR(
    std::abs(lat(nidx, 1) + lli.neighborLink(0, 1, Direction::FORWARD, abc)),
    0.0, near_eps);

  // 0-dir. boundary crossed in backward dir
  idx = lat.coordToLinearIndex({ 0, 0 });
  nidx = lat.coordToLinearIndex({ 2, 0 });
  lli = LinkLatticeIterator<LinkLattice<std::complex<double>, 2>>(idx, lat);

  EXPECT_NEAR(
    std::abs(lat(nidx, 0) + lli.neighborLink(0, 0, Direction::BACKWARD, abc)),
    0.0, near_eps);
  EXPECT_NEAR(
    std::abs(lat(nidx, 1) + lli.neighborLink(0, 1, Direction::BACKWARD, abc)),
    0.0, near_eps);

  // 1-dir. boundary crossed in forward dir
  idx = lat.coordToLinearIndex({ 1, 1 });
  nidx = lat.coordToLinearIndex({ 1, 0 });
  lli = LinkLatticeIterator<LinkLattice<std::complex<double>, 2>>(idx, lat);

  EXPECT_NEAR(
    std::abs(lat(nidx, 1) + lli.neighborLink(1, 1, Direction::FORWARD, abc)),
    0.0, near_eps);
  EXPECT_NEAR(
    std::abs(lat(nidx, 0) + lli.neighborLink(1, 0, Direction::FORWARD, abc)),
    0.0, near_eps);
  // Cross check with no crossing
  EXPECT_EQ(lat(nidx, 1), lli.neighborLink(1, 1, Direction::BACKWARD, abc));

  // 1-dir. boundary crossed in backward dir
  idx = lat.coordToLinearIndex({ 1, 0 });
  nidx = lat.coordToLinearIndex({ 1, 1 });
  lli = LinkLatticeIterator<LinkLattice<std::complex<double>, 2>>(idx, lat);
  EXPECT_NEAR(
    std::abs(lat(nidx, 1) + lli.neighborLink(1, 1, Direction::BACKWARD, abc)),
    0.0, near_eps);
  EXPECT_NEAR(
    std::abs(lat(nidx, 0) + lli.neighborLink(1, 0, Direction::BACKWARD, abc)),
    0.0, near_eps);
}


// TEST(LinkLattice, CanObeyPhaseBoundaryCond){

//     double phase  = std::atan(1.)*0.12045; //"Random" numbers for phase
//     double phase1 = std::atan(1.)*0.269858;
//     auto II       = std::complex<double>(0.,1.);
//     auto bwd      = std::exp(II*phase);
//     auto fwd      = std::exp(-II*phase);
//     auto bwd1     = std::exp(II*phase1);
//     auto fwd1     = std::exp(-II*phase1);
//     BoundaryCondition<double, 2> abc({phase,phase1},"PHASE");
//     std::array<std::complex<double>, 2> init_ar;
//     LinkLattice<std::complex<double>, 2> lat({3,2}, 1. );
//     init_ar.fill(2.);  *lat.at({1,0})=init_ar;
//     init_ar.fill(2.);  *lat.at({1,1})=init_ar;
//     init_ar.fill(3.);  *lat.at({2,0})=init_ar;
//     init_ar.fill(3.);  *lat.at({2,1})=init_ar;

//     // 0-dir. boundary crossed in forward dir
//     auto idx  = lat.coordToLinearIndex({2,0});
//     auto nidx = lat.coordToLinearIndex({0,0});
//     auto lli  = LinkLatticeIterator< LinkLattice<std::complex<double>,
//     2>>(idx,lat);

//     EXPECT_NEAR(std::abs(fwd*lat(nidx,0)-lli.neighborLink(0,0,Direction::FORWARD,abc)),0.0,near_eps);
//     EXPECT_NEAR(std::abs(fwd*lat(nidx,1)-lli.neighborLink(0,1,Direction::FORWARD,abc)),0.0,near_eps);

//     // 0-dir. boundary crossed in backward dir
//     idx  = lat.coordToLinearIndex({0,0});
//     nidx = lat.coordToLinearIndex({2,0});
//     lli  = LinkLatticeIterator< LinkLattice<std::complex<double>,
//     2>>(idx,lat);

//     EXPECT_NEAR(std::abs(bwd*lat(nidx,0)-lli.neighborLink(0,0,Direction::BACKWARD,abc)),0.0,near_eps);
//     EXPECT_NEAR(std::abs(bwd*lat(nidx,1)-lli.neighborLink(0,1,Direction::BACKWARD,abc)),0.0,near_eps);

//     // 1-dir. boundary crossed in forward dir
//     idx  = lat.coordToLinearIndex({1,1});
//     nidx = lat.coordToLinearIndex({1,0});
//     lli  = LinkLatticeIterator< LinkLattice<std::complex<double>,
//     2>>(idx,lat);

//     EXPECT_NEAR(std::abs(fwd1*lat(nidx,1)-lli.neighborLink(1,1,Direction::FORWARD,abc)),0.0,near_eps);
//     EXPECT_NEAR(std::abs(fwd1*lat(nidx,0)-lli.neighborLink(1,0,Direction::FORWARD,abc)),0.0,near_eps);
//     // Cross check with no crossing
//     EXPECT_EQ(lat(nidx,1),lli.neighborLink(1,1,Direction::BACKWARD,abc));

//     // 1-dir. boundary crossed in backward dir
//     idx  = lat.coordToLinearIndex({1,0});
//     nidx = lat.coordToLinearIndex({1,1});
//     lli  = LinkLatticeIterator< LinkLattice<std::complex<double>,
//     2>>(idx,lat);
//     EXPECT_NEAR(std::abs(bwd1*lat(nidx,0)-lli.neighborLink(1,0,Direction::BACKWARD,abc)),0.0,near_eps);
//     EXPECT_NEAR(std::abs(bwd1*lat(nidx,1)-lli.neighborLink(1,1,Direction::BACKWARD,abc)),0.0,near_eps);

// }


TEST(FullLattice, ConstructionIsInitialisation)
{
  FullLattice<double, double, 4> lat({ 3, 4, 5, 6 });
  FullLattice<double, double, 5> lat_init({ 3, 4, 5, 6, 7 }, 0.0,
                                          { 1., 2., 3., 4., 5. });


  for (size_t i = 0; i < lat.volume(); ++i) {
    for (size_t nu = 0; nu < lat.dimensions(); ++nu) {
      EXPECT_DOUBLE_EQ(lat.getLink(i, nu), 0.);
    }
    EXPECT_DOUBLE_EQ(lat.getSite(i), 0.);
  }


  for (size_t i = 0; i < lat_init.volume(); ++i) {
    for (size_t nu = 0; nu < lat_init.dimensions(); ++nu) {
      EXPECT_DOUBLE_EQ(lat_init.getLink(i, nu), 1. + static_cast<double>(nu));
    }
    EXPECT_DOUBLE_EQ(lat_init.getSite(i), 0.);
  }
}

TEST(FullLattice, CanUseIteratorToRunLinearlyOverLattice)
{
  FullLattice<int, int, 3> lat({ 2, 3, 4 });
  int i = 0;

  // Initialise lattice
  for (auto it = lat.begin(); it != lat.end(); ++it) {
    (*it).site = (++i);
    (*it).links = std::array<int, 3>{ i * 2, i * 3, i * 4 };
  }

  // Test indexing
  i = 0;
  for (auto it = lat.begin(); it != lat.end(); ++it) {
    EXPECT_EQ((*it).site, ++i);
  }

  // Test for correct values at given site
  for (std::size_t j = 0; j < lat.volume(); ++j) {
    EXPECT_EQ(static_cast<size_t>(lat[j].site), (j + 1));
    i = static_cast<int>(j + 1);
    EXPECT_EQ(lat[j].links, (std::array<int, 3>{ i * 2, i * 3, i * 4 }));
  }
}


TEST(FullLattice, CanObeyPeriodicBoundaryCond)
{
  BoundaryCondition<double, 2> pbc({ 1., 1. }, "PERIODIC");
  // std::array<std::complex<double>, 2> init_ar;
  FullLattice<std::complex<double>, std::complex<double>, 2> lat({ 3, 2 }, 1.,
                                                                 0.);
  SiteAndLinks<std::complex<double>, std::complex<double>, 2> slstruct;

  slstruct.links.fill(2.);
  slstruct.site = 2.;
  (*lat.at({ 1, 0 })) = slstruct;
  slstruct.links.fill(2.);
  slstruct.site = 2.;
  (*lat.at({ 1, 1 })) = slstruct;
  slstruct.links.fill(3.);
  slstruct.site = 3.;
  (*lat.at({ 2, 0 })) = slstruct;
  slstruct.links.fill(3.);
  slstruct.site = 3.;
  (*lat.at({ 2, 1 })) = slstruct;

  FullLatticeIterator<
    FullLattice<std::complex<double>, std::complex<double>, 2>>
    fli(lat);

  auto idx = lat.coordToLinearIndex({ 0, 0 });
  auto nidx = lat.coordToLinearIndex({ 1, 0 });


  // Links
  EXPECT_EQ(lat.getLink(nidx, 0),
            fli.neighborLink(0, 0, Direction::FORWARD, pbc));
  EXPECT_EQ(lat.getLink(nidx, 1),
            fli.neighborLink(0, 1, Direction::FORWARD, pbc));
  // Sites
  EXPECT_EQ(lat.getSite(nidx), fli.neighborSite(0, Direction::FORWARD, pbc));
  EXPECT_EQ((*lat.at({ 2, 0 })).site, fli.neighborSite(0, Direction::BACKWARD));


  fli = FullLatticeIterator<
    FullLattice<std::complex<double>, std::complex<double>, 2>>(idx, lat);

  EXPECT_EQ(lat.getLink(nidx, 0), fli.neighborLink(0, 0, Direction::FORWARD));

  idx = lat.coordToLinearIndex({ 1, 1 });
  nidx = lat.coordToLinearIndex({ 1, 0 });

  fli = FullLatticeIterator<
    FullLattice<std::complex<double>, std::complex<double>, 2>>(idx, lat);

  // Links
  EXPECT_EQ(lat.getLink(nidx, 1), fli.neighborLink(1, 1, Direction::FORWARD));
  EXPECT_EQ(lat.getLink(nidx, 1), fli.neighborLink(1, 1, Direction::BACKWARD));
  // Site
  EXPECT_EQ(lat.getSite(nidx), fli.neighborLink(1, 1, Direction::FORWARD));
  EXPECT_EQ(lat.getSite(nidx), fli.neighborLink(1, 1, Direction::BACKWARD));
}


TEST(FullLattice, CanObeyMixedBoundaryCond)
{

  BoundaryCondition<double, 2> pbc({ 1., 1. }, "PERIODIC");
  BoundaryCondition<double, 2> abc({ -1, -1 }, "PERIODIC");
  // std::array<std::complex<double>, 2> init_ar;
  FullLattice<std::complex<double>, std::complex<double>, 2> lat({ 3, 2 }, 1.,
                                                                 0.);
  SiteAndLinks<std::complex<double>, std::complex<double>, 2> slstruct;

  slstruct.links.fill(2.);
  slstruct.site = 2.;
  (*lat.at({ 1, 0 })) = slstruct;
  slstruct.links.fill(2.);
  slstruct.site = 2.;
  (*lat.at({ 1, 1 })) = slstruct;
  slstruct.links.fill(3.);
  slstruct.site = 3.;
  (*lat.at({ 2, 0 })) = slstruct;
  slstruct.links.fill(3.);
  slstruct.site = 3.;
  (*lat.at({ 2, 1 })) = slstruct;

  auto idx = lat.coordToLinearIndex({ 2, 0 });
  auto nidx = lat.coordToLinearIndex({ 0, 0 });
  FullLatticeIterator<
    FullLattice<std::complex<double>, std::complex<double>, 2>>
    fli(idx, lat);

  // Periodic Links anti-periodic sites
  // Links
  EXPECT_EQ(lat.getLink(nidx, 0),
            fli.neighborLink(0, 0, Direction::FORWARD, pbc));
  EXPECT_EQ(lat.getLink(nidx, 1),
            fli.neighborLink(0, 1, Direction::FORWARD, pbc));
  // Sites
  EXPECT_NEAR(std::abs((*lat.at({ 0, 0 })).site +
                       fli.neighborSite(0, Direction::FORWARD, abc)),
              0., near_eps);
  EXPECT_NEAR(std::abs((*lat.at({ 2, 1 })).site -
                       fli.neighborSite(1, Direction::FORWARD, abc)),
              0., near_eps);
  EXPECT_NEAR(std::abs((*lat.at({ 2, 1 })).site +
                       fli.neighborSite(1, Direction::BACKWARD, abc)),
              0., near_eps);


  // Anti-Periodic Links and periodic sites
  // Links
  idx = lat.coordToLinearIndex({ 0, 0 });
  nidx = lat.coordToLinearIndex({ 2, 0 });
  // Test bc in 0-Direction
  fli = FullLatticeIterator<
    FullLattice<std::complex<double>, std::complex<double>, 2>>(idx, lat);
  EXPECT_NEAR(std::abs(lat.getLink(nidx, 0) +
                       fli.neighborLink(0, 0, Direction::BACKWARD, abc)),
              0., near_eps);
  EXPECT_NEAR(std::abs(lat.getLink(nidx, 1) +
                       fli.neighborLink(0, 1, Direction::BACKWARD, abc)),
              0., near_eps);
  // Test bc in 1-Direction
  idx = lat.coordToLinearIndex({ 0, 1 });
  nidx = lat.coordToLinearIndex({ 0, 0 });
  EXPECT_NEAR(std::abs(lat.getLink(nidx, 0) +
                       fli.neighborLink(1, 0, Direction::FORWARD, abc)),
              0., near_eps);
  EXPECT_NEAR(std::abs(lat.getLink(nidx, 1) +
                       fli.neighborLink(1, 1, Direction::FORWARD, abc)),
              0., near_eps);

  // Sites
  EXPECT_EQ(lat.getSite(nidx), fli.neighborSite(1, Direction::FORWARD, pbc));
  EXPECT_EQ(lat.getSite(nidx), fli.neighborSite(1, Direction::BACKWARD, pbc));
  idx = lat.coordToLinearIndex({ 2, 0 });
}
