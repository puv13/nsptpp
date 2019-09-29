#include "./gaugegroups/u1.h"
#include "./latticefields/cp.h"
#include "expansion.h"
#include "io.h"
#include "lattice.h"
#include "test_helpers.h"
#include <cstdio>
#include <gtest/gtest.h>

TEST(IO, CanInstantiateWriter)
{
  std::remove("test_file.h5");
  Hdf5File wtr("test_file.h5");
}
TEST(IO, CanPushIntoNonexistingGroups)
{
  Hdf5File wtr("test_file.h5");
  wtr.push("testgroup");
}
TEST(IO, CanPushIntoExistingGroup)
{
  // exists from last test
  Hdf5File wtr("test_file.h5");
  wtr.push("testgroup");
}
TEST(IO, CanPopFromGroup)
{
  Hdf5File wtr("test_file.h5");
  wtr.push("testgroup");
  wtr.pop();
}
TEST(IO, CanSerializeUnderlyingTypeObjectComplexDouble)
{
  Expansion<std::complex<double>, 5> expansion = randomVector<5>();
  auto serialized = expansion.serialize();
  EXPECT_EQ(serialized.size(), 5UL * 2UL);
  for (std::size_t i = 0; i < 5; ++i) {
    EXPECT_EQ(serialized[2 * i + 0], std::real(expansion[i]));
    EXPECT_EQ(serialized[2 * i + 1], std::imag(expansion[i]));
  }
}
TEST(IO, CanDeserializeUnderlyingTypeObjectComplexDouble)
{
  Expansion<std::complex<double>, 5> expansion = randomVector<5>();
  auto serialized = expansion.serialize();

  Expansion<std::complex<double>, 5> other; // is set to 0
  other = deserialize<std::complex<double>, 5>(serialized);
  EXPECT_TRUE(expansion == other);
}

// TEST(IO, CanWriteAndReadComplexExpansionSiteLatticeToAndFromFile)
// {
//   SiteLattice<Expansion<std::complex<double>, 5>, 3> lat({3,4,5});
//   //pseudo
//   for( auto & e : lat )
//     e = randomVector<5>();

//   std::remove("test_SL_expansion.h5");
//   Hdf5File f("./test_SL_expansion.h5");
//   f.write_lattice(lat);
//   SiteLattice<Expansion<std::complex<double>, 5>, 3> read_lat =
//   f.read_lattice<Expansion<std::complex<double>, 5>, 3>();

//   ASSERT_EQ(lat.volume(), read_lat.volume());
//   for( std::size_t i = 0; i < lat.volume(); ++i ) {
//     ASSERT_TRUE(lat[i] == read_lat[i]);
//   }
// }

// TEST(IO, CanWriteAndReadComplexExpansionLinkLatticeToAndFromFile)
// {
//   LinkLattice<Expansion<std::complex<double>, 5>, 3> lat({3,4,5});
//   //pseudo
//   for( auto & e : lat )
//     e = randomVector<5>();

//   std::remove("test_LL_expansion.h5");
//   Hdf5File f("./test_LL_expansion.h5");
//   f.write_lattice(lat);
//   LinkLattice<Expansion<std::complex<double>, 5>, 3> read_lat =
//   f.read_lattice<Expansion<std::complex<double>, 5>, 3>();

//   ASSERT_EQ(lat.volume(), read_lat.volume());
//   for( std::size_t i = 0; i < lat.volume(); ++i ) {
//     ASSERT_TRUE(lat[i] == read_lat[i]);
//   }
// }


// TEST(IO, CanCreateGroupHirachy)
// {
//     // Hdf5File f("./test_groups.h5");

//     f.push("A");
//     f.push("a1");
//     f.push("aa1");
//     f.pop();
//     f.push("aa2");
//     f.pop();
//     f.pop();
//     f.pop();
//     f.push("B");

//     FullLattice<std::complex<double>,std::complex<double>,2> lat({2,2});

//     f.write_lattice(lat);
// }

TEST(IO, CanSerializeAndDeserializeArbitraryTypes)
{

  // CP(N)
  CP<std::complex<double>, 4> cp1(randomVector<1>()[0], true), cp2;
  std::array<double, sizeof(cp1) / sizeof(double)> serial_cp;

  serial_cp = serialize_type(cp1);
  // for (auto el : serial_cp){
  //  	std::cout << el << "\t";
  // }

  cp2 = deserialize_type<decltype(cp2)>(serial_cp);

  EXPECT_EQ(cp1, cp2);

  // U1
  U1 u(randomVector<1>()[0], true), u2;
  std::array<double, sizeof(u) / sizeof(double)> serial_u1;
  serial_u1 = serialize_type(u);

  u2 = deserialize_type<decltype(u2)>(serial_u1);

  EXPECT_EQ(u, u2);
}

TEST(IO, CanWriteAndReadBoundaryConditions)
{
  BoundaryCondition<double, 3> const bc({ 1.2, 2.3, 3.4 },
                                        "testBoundaryConditions");
  std::remove("./test_bc.h5");
  Hdf5File f("./test_bc.h5");
  f.write_boundaryCondition(bc);
  auto read_bc = f.read_boundaryCondition<double, 3>();

  ASSERT_EQ(read_bc.getName(), bc.getName());
  for (std::size_t i = 0; i < 3; ++i) {
    ASSERT_EQ(read_bc[i], bc[i]);
  }
}
// TEST(IO, ThrowsIfIncompatibleNumberOfDimensionsIsRead)
// {
//   SiteLattice<Expansion<std::complex<double>, 5>, 3>
//   lat(std::array<std::size_t, 3>({{3,4,5}}));
//   //pseudo
//   for( auto & e : lat )
//     e = randomVector<5>();

//   std::remove("test_expansion.h5");
//   Hdf5File file("./test_expansion.h5");
//   file.write_lattice(lat);

//   auto fun = [&file](SiteLattice<Expansion<std::complex<double>, 5>, 1>
//   &rl){file.read_lattice<Expansion<std::complex<double>, 5>, 1>(rl);};
//   // 1 != 3:
//   EXPECT_THROW(fun(), std::runtime_error);
// }
// TEST(IO, ThrowsIfIncompatibleNumberOfOrdersIsRead)
// {
//   SiteLattice<Expansion<std::complex<double>, 3>, 3>
//   lat(std::array<std::size_t, 3>({{3,4,5}}));
//   //pseudo
//   for( auto & e : lat )
//     e = randomVector<3>();

//   std::remove("test_expansion.h5");
//   Hdf5File file("./test_expansion.h5");
//   file.write_lattice(lat);

//   // 5 != 3:
//   auto fun = [&file](){file.read_lattice<Expansion<std::complex<double>, 5>,
//   3>();};
//   EXPECT_THROW(fun(), std::runtime_error);
// }


TEST(IO, CanWriteFullLatticeToFile)
{
  FullLattice<double, double, 2> lat_init({ 3, 7 }, 1.0, { 2., 3. }),
    lat_read({ 3, 7 }, 13.0, { 31., 42. });

  std::remove("./test_FullLattice.h5");
  Hdf5File f("./test_FullLattice.h5");
  // Hdf5File h("./test_FullLatticeRead.h5");

  f.write_lattice(lat_init, "testdataset");

  f.read_lattice(lat_read, "testdataset");

  auto init_it = lat_init.begin();
  auto read_it = lat_read.begin();

  for (; init_it != lat_init.end() and read_it != lat_read.end();
       ++init_it, ++read_it) {

    EXPECT_EQ((*init_it).site, (*read_it).site);

    for (auto i = 0UL; i < lat_init.dimensions(); ++i) {
      EXPECT_EQ((*init_it).links[i], (*read_it).links[i]);
    }
  }
}

TEST(IO, CanWriteLinkLatticeToFile)
{
  LinkLattice<double, 2> lat_init({ 3, 7 }, 1.0), lat_read({ 3, 7 }, 13.0);

  std::remove("./test_LinkLattice.h5");
  Hdf5File f("./test_LinkLattice.h5");
  f.write_lattice(lat_init, "testdataset");

  f.read_lattice(lat_read, "testdataset");

  auto init_it = lat_init.begin();
  auto read_it = lat_read.begin();

  for (; init_it != lat_init.end() and read_it != lat_read.end();
       ++init_it, ++read_it) {

    for (auto i = 0UL; i < lat_init.dimensions(); ++i) {
      EXPECT_EQ((*init_it)[i], (*read_it)[i]);
    }
  }
}

TEST(IO, CanWriteSiteLatticeToFile)
{
  SiteLattice<double, 2> lat_init({ 3, 7 }, 1.0), lat_read({ 3, 7 }, 13.0);

  std::remove("./test_SiteLattice.h5");
  Hdf5File f("./test_SiteLattice.h5");
  f.write_lattice(lat_init, "testdataset");

  f.read_lattice(lat_read, "testdataset");

  auto init_it = lat_init.begin();
  auto read_it = lat_read.begin();

  for (; init_it != lat_init.end() and read_it != lat_read.end();
       ++init_it, ++read_it) {

    EXPECT_EQ((*init_it), (*read_it));
  }
}


TEST(IO, CanWriteFullLatticeExpansionToFile)
{
  FullLattice<Expansion<std::complex<double>, 5>,
              Expansion<std::complex<double>, 5>, 2>
    lat_init({ 3, 7 }, 1.0, { 2., 3. }), lat_read({ 3, 7 }, 13.0, { 31., 42. });

  Hdf5File f("./test_FullLatticeExpansion.h5");
  // Hdf5File h("./test_FullLatticeRead.h5");

  f.write_lattice(lat_init, "testdataset");

  f.read_lattice(lat_read, "testdataset");

  auto init_it = lat_init.begin();
  auto read_it = lat_read.begin();

  for (; init_it != lat_init.end() and read_it != lat_read.end();
       ++init_it, ++read_it) {

    EXPECT_EQ((*init_it).site, (*read_it).site);

    for (auto i = 0UL; i < lat_init.dimensions(); ++i) {
      EXPECT_EQ((*init_it).links[i], (*read_it).links[i]);
    }
  }
}


TEST(IO, ThrowsIfIncompatibleDataSpaceExistsFullLattice)
{
  std::remove("./test_FullLatticeExistingDspace.h5");
  Hdf5File f("./test_FullLatticeExistingDspace.h5");
  FullLattice<double, double, 2> lat_1({ 3, 7 }), lat_2({ 4, 8 });

  f.write_lattice(lat_1);


  auto fun = [&f](FullLattice<double, double, 2> lat) { f.write_lattice(lat); };
  EXPECT_THROW(fun(lat_2), std::runtime_error);
}


TEST(IO, ThrowsIfIncompatibleSiteOrLinkDataTypeFullLattice)
{
  std::remove("./test_FullLatticeExistingDspace.h5");
  Hdf5File f("./test_FullLatticeExistingDspace.h5");
  FullLattice<double, double, 2> lat_1({ 3, 7 });
  FullLattice<std::complex<double>, std::complex<double>, 2> lat_2({ 3, 7 });

  f.write_lattice(lat_1);


  auto fun =
    [&f](FullLattice<std::complex<double>, std::complex<double>, 2> lat) {
      f.write_lattice(lat);
    };
  EXPECT_THROW(fun(lat_2), std::runtime_error);
}
