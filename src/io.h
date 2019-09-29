/// \file
///
/// \brief Code to simplify HDF5 input and output

#pragma once

#include "expansion.h"
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wold-style-cast"
#pragma GCC diagnostic ignored "-Wsign-conversion"
#pragma GCC diagnostic ignored "-Wreorder"
#include "highfive/H5Attribute.hpp"
#include "highfive/H5DataSet.hpp"
#include "highfive/H5DataSpace.hpp"
#include "highfive/H5File.hpp"
#pragma GCC diagnostic pop
#include "lattice.h"
#include "nsptpp_config.h"
#include <fstream>
#include <string>


// -----------------------------------------------------------------------------
// Helper functions
// -----------------------------------------------------------------------------

/// Check if file exists
inline bool fileExists(std::string const &fname)
{
  std::ifstream f(fname.c_str());
  return f.good();
}

/// Serialise a generic data type as an array of doubles
template <typename T>
inline std::array<double, sizeof(T) / sizeof(double)>
serialize_type(T const &orig)
{

  const std::size_t numberOfDoubles = sizeof(T) / sizeof(double);
  // Assert that size of T is multiple of the size of double
  static_assert(numberOfDoubles * sizeof(double) == sizeof(T),
                "Cannot serialize a data type that has a size that is not"
                " an integer multiple of sizeof(double)");


  std::array<double, numberOfDoubles> res;
  auto it = res.begin();
  for (std::size_t i = 0; i < numberOfDoubles; ++i) {
    // Cast orig as double pointer and fill res with doubles
    *it = reinterpret_cast<double const *>(&orig)[i];
    ++it;
  }

  return res;
}

/// Deserialize an array of doubles to get the type T (inverse of serialize)
template <typename T>
inline T
deserialize_type(std::array<double, sizeof(T) / sizeof(double)> const &serial)
{
  T res;
  const std::size_t numberOfDoubles = sizeof(T) / sizeof(double);
  // Assert that size of T is multiple of the size of double
  static_assert(numberOfDoubles * sizeof(double) == sizeof(T),
                "Cannot deserialize an array that has a size that is not"
                " equal to sizeof(T)");

  for (std::size_t i = 0; i < numberOfDoubles; ++i) {
    // Cast orig as double pointer and fill res with doubles
    reinterpret_cast<double *>(&res)[i] = serial[i];
  }

  return res;
}

// -----------------------------------------------------------------------------
// Main Class
// -----------------------------------------------------------------------------

/// Class for HDF5 file IO.
class Hdf5File {

 private:
  HighFive::Group curgroup;
  bool isWritable;
  HighFive::File file;
  std::vector<std::string> groups;
  /// Returns the current group
  std::string cwd() const
  {
    std::string res("/");
    for (auto const &elem : groups)
      res += elem + "/";
    return res;
  }

 public:
  /// Constructor
  Hdf5File(std::string const &filename, bool readOnly = false)
    : isWritable(not readOnly),
      file(HighFive::File(filename, fileExists(filename)
                                      ? (readOnly ? HighFive::File::ReadOnly
                                                  : HighFive::File::ReadWrite)
                                      : HighFive::File::ReadWrite |
                                          HighFive::File::Create |
                                          HighFive::File::Truncate))
  {
    curgroup = file.getGroup(".");
  }

  /// Sets current group to 'groupname'
  void push(std::string const &groupname)
  {
    groups.push_back(groupname);
    if (curgroup.exist(cwd()))
      curgroup = curgroup.getGroup(cwd());
    else if (isWritable)
      curgroup = curgroup.createGroup(cwd());
    else
      throw std::runtime_error("cannot create group in a read-only file.");
  }

  /// Sets current group to next higher level in group hierarchy
  void pop()
  {
    groups.pop_back();
    curgroup = curgroup.getGroup(cwd());
  }


  /// Convenience function to write Boundary Condition
  template <typename P, std::size_t dim>
  void write_boundaryCondition(BoundaryCondition<P, dim> const &bc)
  {

    // Serialise "Phases"
    auto bc_serial = serialize_type(bc.getPhases());
    auto bc_size = bc_serial.size();

    HighFive::DataSet dataset = curgroup.createDataSet<double>(
      "phases", HighFive::DataSpace({ bc_size }));

    std::vector<double> phases(bc_size);
    for (std::size_t i = 0u; i < bc_size; ++i) {
      phases[i] = bc_serial[i];
    }

    dataset.write(phases);

    // Write name
    std::vector<std::string> namelist;
    namelist.push_back(bc.getName());
    auto nameset = curgroup.createDataSet<std::string>(
      "name", HighFive::DataSpace::From(namelist));
    nameset.write(namelist);
  }


  /// Convenience function to read Boundary Condition
  template <typename P, std::size_t dim>
  BoundaryCondition<P, dim> read_boundaryCondition()
  {

    const std::size_t numberOfDoubles = (sizeof(P) * dim) / sizeof(double);
    // using namespace HighFive;

    auto phasesdset = curgroup.getDataSet("phases");

    // Get serialised phases
    std::vector<double> phases;
    phasesdset.read(phases);

    if (phases.size() != numberOfDoubles) {
      throw std::runtime_error(
        "number of phases (dimensions) does not match template argument");
    }

    // Convert to array
    std::array<double, numberOfDoubles> aux_arr;
    for (std::size_t i = 0u; i < dim; ++i) {
      aux_arr[i] = phases[i];
    }


    auto phasesarr = deserialize_type<std::array<P, dim>>(aux_arr);

    std::vector<std::string> namelist;
    auto namedset = curgroup.getDataSet("name");
    namedset.read(namelist);

    return BoundaryCondition<P, dim>(phasesarr, namelist.front());
  }


  /// Write an (scalar) attribute in current group and do some sanity checks
  template <typename T>
  void addAttribute(const std::string &attr_name, T const &attr_data,
                    const std::string &dset_name = "data")
  {
    HighFive::DataSet dataset = curgroup.getDataSet(dset_name);

    // Check for existence
    if (not dataset.hasAttribute(attr_name)) {
      HighFive::Attribute attr = dataset.createAttribute<T>(
        attr_name, HighFive::DataSpace::From(attr_data));

      attr.write(attr_data);
    }
    else {
      HighFive::Attribute attr = dataset.getAttribute(attr_name);
      // No way to check for data space of Attribute in HighFive?
      std::cout << "\tWarning: Overwriting existing data in Attribute"
                << " " << attr_name << ". \n";
      attr.write(attr_data);
    }
  }

  /// Read an (scalar) attribute in current group and do some sanity checks
  template <typename T>
  void readAttribute(T &attr_data, const std::string &attr_name,
                     const std::string &dset_name = "data")
  {


    HighFive::DataSet dataset = curgroup.getDataSet(dset_name);
    HighFive::Attribute attr = dataset.getAttribute(attr_name);

    attr.read(attr_data);
  }

  /// Read an (vector) attribute in current group and do some sanity checks
  template <typename T>
  void readAttributeVector(std::vector<T> &attr_data,
                           const std::string &attr_name,
                           const std::string &dset_name = "data")
  {


    HighFive::DataSet dataset = curgroup.getDataSet(dset_name);
    HighFive::Attribute attr = dataset.getAttribute(attr_name);

    attr.read(attr_data);
  }

  /// Write a vector attribute in current group and do some sanity checks
  template <typename T>
  void addAttributeVector(const std::string &attr_name,
                          const std::vector<T> &attr_data,
                          const std::string &dset_name = "data")
  {

    HighFive::DataSet dataset = curgroup.getDataSet(dset_name);

    // Check for existence
    if (not dataset.hasAttribute(attr_name)) {
      HighFive::Attribute attr = dataset.createAttribute<T>(
        attr_name, HighFive::DataSpace::From(attr_data));

      attr.write(attr_data);
    }
    else {
      HighFive::Attribute attr = dataset.getAttribute(attr_name);
      // No way to check for data space of Attribute in HighFive?
      std::cout << "\tWarning: Overwriting existing data in Attribute"
                << " " << attr_name << ". \n";
      attr.write(attr_data);
    }
  }

  /// Write Lattice params as attributes
  /// Should work with any lattice class derived from BareLattice
  template <typename L>
  void writeCommonLatticeAttributes(L const &lat,
                                    std::string &dset_name = "data")
  {

    // Git hash of commit used to build this code
    //
    // \note (If we don't explicitly specify std::string we get the linker
    // error "undefined reference to `HighFive::AtomicType<char
    // [41]>::AtomicType()'")
    addAttribute<std::string>("githash", NSPTPP_GIT_HASH, dset_name);

    // Build time of this code
    //
    // \note (If we don't explicitly specify std::string we get the linker
    // error "undefined reference to `HighFive::AtomicType<char
    // [41]>::AtomicType()'")
    addAttribute<std::string>("buildtime", NSPTPP_BUILD_TIME, dset_name);

    /// Lattice dimensions
    auto dimensionsarr = lat.dimensionsArray();
    std::vector<std::size_t> dimensions;
    dimensions.reserve(dimensionsarr.size());
    for (auto const &d : dimensionsarr) {
      dimensions.push_back(d);
    }
    std::size_t ndims(dimensions.size());
    addAttribute("number_of_dimensions", ndims, dset_name);
    addAttributeVector<std::size_t>("dimensions", dimensions, dset_name);
  }

  /// Return a list of Objects in the current group
  std::vector<std::string> ObjectNames()
  {
    return this->curgroup.listObjectNames();
  }


  // -----------------------------------------------------------------------------
  // FullLattice specific stuff
  // -----------------------------------------------------------------------------

  /// Specialisation for FullLattice with any types (incl. expansions)
  template <typename S, typename L, std::size_t dim>
  void write_lattice(FullLattice<S, L, dim> const &lat,
                     std::string dset_name = "data")
  {
    // Sanity check
    if (not isWritable) {
      throw std::runtime_error("Cannot write lattice to a read-only file.");
    }

    const std::size_t numberOfDoublesL = sizeof(L) / sizeof(double);
    const std::size_t numberOfDoublesS = sizeof(S) / sizeof(double);

    const std::size_t numberOfDoubles =
      numberOfDoublesS + dim * numberOfDoublesL;


    // Transform lattice to an std::vector<std::vector<double>> for output
    // Memory is cheap, so we generate a copy of the data.

    std::vector<std::vector<double>> latcopy(lat.volume());
    auto out_it = latcopy.begin();
    auto in_it = lat.begin();

    for (; out_it != latcopy.end() and in_it != lat.end(); ++in_it, ++out_it) {
      out_it->resize(numberOfDoubles); // Prepare output vector

      auto tmp = (*in_it); // Pointer to lattice site

      // Output vector
      auto output_vec = out_it->data();

      // Serialise current lattice site
      auto cur_site = serialize_type(tmp.site);
      for (auto i = 0UL; i < numberOfDoublesS; ++i) {
        output_vec[i] = cur_site[i];
      }

      // Serialise current lattice links
      for (auto i = 0UL; i < dim; ++i) {
        auto cur_link = serialize_type(tmp.links[i]);
        auto offset = numberOfDoublesS + i * numberOfDoublesL;

        for (auto j = 0LU; j < numberOfDoublesL; ++j) {
          output_vec[offset + j] = cur_link[j];
        }
      }
    }

    // Create HDF5 dataset
    if (not isWritable)
      throw std::runtime_error("Cannot write lattice to a read-only file");


    if (not curgroup.exist(dset_name)) // Data set does not exist
    {
      HighFive::DataSet dataset = curgroup.createDataSet<double>(
        dset_name, HighFive::DataSpace({ lat.volume(), numberOfDoubles }));

      // Write data
      dataset.write(latcopy);
    }
    // Data set with right dimensions does exist
    else {
      HighFive::DataSet dataset = curgroup.getDataSet(dset_name);
      HighFive::DataSpace dspace = dataset.getSpace();
      auto dspace_dims = dspace.getDimensions();
      auto dspace_ndim = dspace.getNumberDimensions();

      if (dspace_ndim == 2 && dspace_dims[0] == lat.volume() &&
          dspace_dims[1] == numberOfDoubles) {
        // Write data
        std::cout << "\tWarning: Overwriting existing data in dataset "
                  << dset_name << ". \n";
        dataset.write(latcopy);
      }
      else {
        throw std::runtime_error(
          "Existing data set with different dimensions.");
      }
    }

    // Write Lattice parameters as attributes
    writeCommonLatticeAttributes(lat, dset_name);

    // This assumes that both expansions have the same order (which
    // is the only sensible option anyway)
    if (is_expansion<S>::value) {
      constexpr std::size_t ord = is_expansion<S>::ord;
      addAttribute<std::size_t>("number_of_orders", ord, dset_name);
    }

    return;
  }

  /// read_lattice specialisation for FullLattice with any types (incl.
  /// expansions)
  template <typename S, typename L, std::size_t dim>
  void read_lattice(FullLattice<S, L, dim> &res, std::string dset_name = "data")
  {

    // Calculate required "storage space" per site
    const std::size_t numberOfDoublesL = sizeof(L) / sizeof(double);
    const std::size_t numberOfDoublesS = sizeof(S) / sizeof(double);

    // Get dataset
    HighFive::DataSet dataset = curgroup.getDataSet(dset_name);
    HighFive::DataSpace dataspace = dataset.getSpace();

    // Calculate dataspace storage size
    auto storage = 1LU;
    for (auto entry : dataspace.getDimensions()) {
      storage *= entry;
    }

    // Read number of dimensions
    std::size_t ndims(0);
    readAttribute(ndims, "number_of_dimensions", dset_name);
    if (ndims != dim) {
      throw std::runtime_error(
        "read configuration has different number of dimensions.");
    }

    // Read Lattice layout
    std::vector<std::size_t> dims(ndims);
    readAttributeVector<std::size_t>(dims, "dimensions", dset_name);
    std::array<std::size_t, dim> dims_arr;
    auto vol = 1LU;
    for (std::size_t i = 0u; i < dim; ++i) {
      dims_arr[i] = dims[i];
      vol *= dims[i];
    }

    // Consistency Check
    auto req_storage = vol * (numberOfDoublesL * dim + numberOfDoublesS);
    if (req_storage != storage) {
      throw std::runtime_error(
        "Cannot read lattice. Link and/or Site data type of "
        "stored lattice does not match expectations.");
    }

    // Now we can read the sites and links
    std::vector<std::vector<double>> latcopy(res.volume());
    dataset.read(latcopy);
    auto in_it = latcopy.begin();
    auto out_it = res.begin();


    for (; in_it != latcopy.end() and out_it != res.end(); ++in_it, ++out_it) {
      auto &tmp = *out_it; // & because otherwise we get a copy of
                           // *out_it


      // Deserialise current lattice site
      std::array<double, numberOfDoublesS> cur_site_array;
      for (auto i = 0UL; i < numberOfDoublesS; ++i) {
        // Copy data to cur_site_array
        cur_site_array[i] = (*in_it)[i];
        // std::cout << "in_it[" << i << "] : " << (*in_it)[i] << std::endl;
      }

      auto cur_site = deserialize_type<S>(cur_site_array);
      tmp.site = cur_site;

      // Deserialise current lattice links
      for (auto i = 0UL; i < dim; ++i) {
        //     auto cur_link = serialize_type(tmp.links[i]);
        std::array<double, numberOfDoublesL> cur_link_array;
        auto offset = numberOfDoublesS + i * numberOfDoublesL;

        for (auto j = 0LU; j < numberOfDoublesL; ++j) {
          cur_link_array[j] = (*in_it)[offset + j];
        }

        auto cur_link = deserialize_type<L>(cur_link_array);
        tmp.links[i] = cur_link;
      }
    }

    return;
  }


  // -----------------------------------------------------------------------------
  // LinkLattice specific stuff
  // -----------------------------------------------------------------------------

  /// read_lattice specialisation for LinkLattice with any type
  /// (incl. expansions)
  template <typename L, std::size_t dim>
  void read_lattice(LinkLattice<L, dim> &res, std::string dset_name = "data")
  {

    // Calculate required "storage space" per site
    const std::size_t numberOfDoublesL = sizeof(L) / sizeof(double);

    // Get dataset
    HighFive::DataSet dataset = curgroup.getDataSet(dset_name);
    HighFive::DataSpace dataspace = dataset.getSpace();

    // Calculate dataspace storage size
    auto storage = 1LU;
    for (auto entry : dataspace.getDimensions()) {
      storage *= entry;
    }

    // Read number of dimensions
    std::size_t ndims(0);
    readAttribute(ndims, "number_of_dimensions", dset_name);
    if (ndims != dim) {
      throw std::runtime_error(
        "read configuration has different number of dimensions.");
    }

    // Read Lattice layout
    std::vector<std::size_t> dims(ndims);
    readAttributeVector<std::size_t>(dims, "dimensions", dset_name);
    std::array<std::size_t, dim> dims_arr;
    auto vol = 1LU;
    for (std::size_t i = 0u; i < dim; ++i) {
      dims_arr[i] = dims[i];
      vol *= dims[i];
    }

    // Consistency Check
    auto req_storage = vol * (numberOfDoublesL * dim);
    if (req_storage != storage) {
      throw std::runtime_error("Cannot read lattice. Link data type of "
                               "stored lattice does not match expectations.");
    }

    // Now we can read the links
    std::vector<std::vector<double>> latcopy(res.volume());
    dataset.read(latcopy);
    auto in_it = latcopy.begin();
    auto out_it = res.begin();


    for (; in_it != latcopy.end() and out_it != res.end(); ++in_it, ++out_it) {
      auto &tmp = *out_it; // & because otherwise we get a copy of
                           // *out_it

      // Deserialise current lattice links
      for (auto i = 0UL; i < dim; ++i) {
        //     auto cur_link = serialize_type(tmp.links[i]);
        std::array<double, numberOfDoublesL> cur_link_array;
        auto offset = i * numberOfDoublesL;

        for (auto j = 0LU; j < numberOfDoublesL; ++j) {
          cur_link_array[j] = (*in_it)[offset + j];
        }

        auto cur_link = deserialize_type<L>(cur_link_array);
        tmp[i] = cur_link;
      }
    }
    return;
  }


  /// Specialisation for LinkLattice with any type (incl. expansions)
  template <typename L, std::size_t dim>
  void write_lattice(LinkLattice<L, dim> const &lat,
                     std::string dset_name = "data")
  {
    // Sanity check
    if (not isWritable) {
      throw std::runtime_error("Cannot write lattice to a read-only file.");
    }

    const std::size_t numberOfDoublesL = sizeof(L) / sizeof(double);
    const std::size_t numberOfDoubles = dim * numberOfDoublesL;


    // Transform lattice to an std::vector<std::vector<double>> for
    // output. Memory is cheap, so we generate a copy of the data.
    std::vector<std::vector<double>> latcopy(lat.volume());
    auto out_it = latcopy.begin();
    auto in_it = lat.begin();

    for (; out_it != latcopy.end() and in_it != lat.end(); ++in_it, ++out_it) {
      out_it->resize(numberOfDoubles); // Prepare output vector

      auto tmp = (*in_it); // Pointer to lattice site

      // Output vector
      auto output_vec = out_it->data();

      // Serialise current lattice links
      for (auto i = 0UL; i < dim; ++i) {
        auto cur_link = serialize_type(tmp[i]);
        auto offset = i * numberOfDoublesL;

        for (auto j = 0LU; j < numberOfDoublesL; ++j) {
          output_vec[offset + j] = cur_link[j];
        }
      }
    }

    // Create HDF5 dataset
    if (not isWritable)
      throw std::runtime_error("Cannot write lattice to a read-only file");


    if (not curgroup.exist(dset_name)) // Data set does not exist
    {
      HighFive::DataSet dataset = curgroup.createDataSet<double>(
        dset_name, HighFive::DataSpace({ lat.volume(), numberOfDoubles }));

      // Write data
      dataset.write(latcopy);
    }
    // Data set with right dimensions does exist
    else {
      HighFive::DataSet dataset = curgroup.getDataSet(dset_name);
      HighFive::DataSpace dspace = dataset.getSpace();
      auto dspace_dims = dspace.getDimensions();
      auto dspace_ndim = dspace.getNumberDimensions();

      if (dspace_ndim == 2 && dspace_dims[0] == lat.volume() &&
          dspace_dims[1] == numberOfDoubles) {
        // Write data
        std::cout << "\tWarning: Overwriting existing data "
                  << "in dataset " << dset_name << ". \n";
        dataset.write(latcopy);
      }
      else {
        throw std::runtime_error(
          "Existing data set with different dimensions.");
      }
    }

    // Write Lattice parameters as attributes
    writeCommonLatticeAttributes(lat, dset_name);

    if (is_expansion<L>::value) {
      constexpr std::size_t ord = is_expansion<L>::ord;
      addAttribute<std::size_t>("number_of_orders", ord, dset_name);
    }
    return;
  }


  // -----------------------------------------------------------------------------
  // SiteLattice specific stuff
  // -----------------------------------------------------------------------------

  /// read_lattice specialisation for SiteLattice with any type
  /// (incl. expansions)
  template <typename S, std::size_t dim>
  void read_lattice(SiteLattice<S, dim> &res, std::string dset_name = "data")
  {

    // Calculate required "storage space" per site
    const std::size_t numberOfDoubles = sizeof(S) / sizeof(double);

    // Get dataset
    HighFive::DataSet dataset = curgroup.getDataSet(dset_name);
    HighFive::DataSpace dataspace = dataset.getSpace();

    // Calculate dataspace storage size
    auto storage = 1LU;
    for (auto entry : dataspace.getDimensions()) {
      storage *= entry;
    }

    // Read number of dimensions
    std::size_t ndims(0);
    readAttribute(ndims, "number_of_dimensions", dset_name);
    if (ndims != dim) {
      throw std::runtime_error(
        "read configuration has different number of dimensions.");
    }

    // Read Lattice layout
    std::vector<std::size_t> dims(ndims);
    readAttributeVector<std::size_t>(dims, "dimensions", dset_name);
    std::array<std::size_t, dim> dims_arr;
    auto vol = 1LU;
    for (std::size_t i = 0u; i < dim; ++i) {
      dims_arr[i] = dims[i];
      vol *= dims[i];
    }

    // Consistency Check
    auto req_storage = vol * (numberOfDoubles);
    if (req_storage != storage) {
      throw std::runtime_error("Cannot read lattice. Site data type of "
                               "stored lattice does not match expectations.");
    }

    // Now we can read the sites
    std::vector<std::vector<double>> latcopy(res.volume());
    dataset.read(latcopy);
    auto in_it = latcopy.begin();
    auto out_it = res.begin();


    for (; in_it != latcopy.end() and out_it != res.end(); ++in_it, ++out_it) {
      auto &tmp = *out_it; // & because otherwise we get a copy of
                           // *out_it

      // Deserialise current lattice site
      std::array<double, numberOfDoubles> cur_site_array;
      for (auto i = 0UL; i < numberOfDoubles; ++i) {
        // Copy data to cur_site_array
        cur_site_array[i] = (*in_it)[i];
      }
      auto cur_site = deserialize_type<S>(cur_site_array);
      tmp = cur_site;
    }

    return;
  }


  /// Specialisation for SiteLattice with any type (incl. expansions)
  template <typename S, std::size_t dim>
  void write_lattice(SiteLattice<S, dim> const &lat,
                     std::string dset_name = "data")
  {
    // Sanity check
    if (not isWritable) {
      throw std::runtime_error("Cannot write lattice to a read-only file.");
    }

    const std::size_t numberOfDoubles = sizeof(S) / sizeof(double);


    // Transform lattice to an std::vector<std::vector<double>> for output
    // Memory is cheap, so we generate a copy of the data.

    std::vector<std::vector<double>> latcopy(lat.volume());
    auto out_it = latcopy.begin();
    auto in_it = lat.begin();

    for (; out_it != latcopy.end() and in_it != lat.end(); ++in_it, ++out_it) {
      out_it->resize(numberOfDoubles); // Prepare output vector

      auto tmp = (*in_it); // Pointer to lattice site

      // Output vector
      auto output_vec = out_it->data();

      // Serialise current lattice site
      auto cur_site = serialize_type(tmp);
      for (auto i = 0UL; i < numberOfDoubles; ++i) {
        output_vec[i] = cur_site[i];
      }
    }

    // Create HDF5 dataset
    if (not isWritable)
      throw std::runtime_error("Cannot write lattice to a read-only file");


    if (not curgroup.exist(dset_name)) // Data set does not exist
    {
      HighFive::DataSet dataset = curgroup.createDataSet<double>(
        dset_name, HighFive::DataSpace({ lat.volume(), numberOfDoubles }));

      // Write data
      dataset.write(latcopy);
    }
    // Data set with right dimensions does exist
    else {
      HighFive::DataSet dataset = curgroup.getDataSet(dset_name);
      HighFive::DataSpace dspace = dataset.getSpace();
      auto dspace_dims = dspace.getDimensions();
      auto dspace_ndim = dspace.getNumberDimensions();

      if (dspace_ndim == 2 && dspace_dims[0] == lat.volume() &&
          dspace_dims[1] == numberOfDoubles) {
        // Write data
        std::cout << "\tWarning: Overwriting existing data "
                  << "in dataset " << dset_name << ". \n";
        dataset.write(latcopy);
      }
      else {
        throw std::runtime_error(
          "Existing data set with different dimensions.");
      }
    }

    // Write Lattice parameters as attributes
    writeCommonLatticeAttributes(lat, dset_name);

    if (is_expansion<S>::value) {
      constexpr std::size_t ord = is_expansion<S>::ord;
      addAttribute<std::size_t>("number_of_orders", ord, dset_name);
    }

    return;
  }
};
