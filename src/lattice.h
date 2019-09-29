/// \file
///
/// \brief Generic lattice code
/// \note
/// The different lattice classes and lattice iterator classes are indentionally
/// not implemented using a class hierarchy. This is done to avoid possible
/// performance issues due to vtable look-ups.


#pragma once

#include <array>
#include <cassert>
#include <complex>
#include <iostream>
#include <iterator>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

// *****************************************************************************
// Enums and typedefs
// *****************************************************************************

/// Direction enum for more convenient handling of neighbours
enum Direction : int { FORWARD = 0, BACKWARD = 1 };

using size_t = std::size_t;
using std::exp;


// *****************************************************************************
// Helper Functions
// *****************************************************************************

/// Product of all members of an array
template <std::size_t dim>
constexpr std::size_t product(std::array<std::size_t, dim> const &arr)
{
  size_t res = 1;
  for (size_t i = 0; i < dim; ++i) {
    res *= arr[i];
  }
  return res;
}


// *****************************************************************************
/// Boundary condition class
/// \brief Helper for the lattice class.
///
/// This class implements simple boundary conditions (Periodic,
/// Anti-Periodic, ...), than can be implemented by "multiplying"
/// the fields at the boundary by a constant factor. "Multiplication" can be
/// any binary operation and has to be implemented by overloading *-Operator
/// for a custom data type.
/// In the simplest case we multiply the fields by a (complex) number. Periodic
/// boundary conditions are equivalent to multiplying by 1., for example.
// *****************************************************************************
template <typename T, size_t dim>
class BoundaryCondition {

 private:
  /// Array with the size of the number of lattice dimensions.
  /// \details For each lattice direction the "phase" that is
  /// acquired by crossing the border is stored. The phase factor
  /// is equal to 1 for periodic BC and -1 for anti-
  /// periodic BC.
  std::array<T, dim> bcPhase;
  /// Name for the BC.
  std::string name;

 public:
  // ***************************************************************************
  // Constructors
  // ***************************************************************************

  BoundaryCondition() : name("Default BC (Periodic)")
  {
    std::array<T, dim> tmp;
    tmp.fill(T(1.0));
    bcPhase = tmp;
  }
  BoundaryCondition(const std::array<T, dim> &_bcPhase, std::string str)
    : bcPhase(_bcPhase), name(str)
  {
  }

  BoundaryCondition(const BoundaryCondition<T, dim> &other)
  {
    bcPhase = other.getPhases();
    name = other.getName();
  }

  // ***************************************************************************
  // Destructor
  // ***************************************************************************
  ~BoundaryCondition() = default;

  // ***************************************************************************
  // Member Functions
  // ***************************************************************************

  /// Get boundary phase factor
  inline T const &operator[](size_t i) const { return bcPhase[i]; }
  /// Get boundary phase factor
  inline T &operator[](size_t i) { return bcPhase[i]; }
  std::array<T, dim> getPhases() const { return bcPhase; }
  std::string getName() const { return name; }
};


// *****************************************************************************
/// Lattice iterator class
/// \brief
// *****************************************************************************
template <typename T>
class raw_lattice_iterator
  : public std::iterator<std::bidirectional_iterator_tag, T> {
  size_t pos;
  T *parent;

 public:
  explicit raw_lattice_iterator(size_t pos_, T &parent_)
    : pos(pos_), parent(&parent_)
  {
  }

  virtual ~raw_lattice_iterator() = default;

  size_t index() const { return pos; }

  std::array<size_t, T::dimensions()> coord() const
  {
    return parent->linearIndexToCoord(pos);
  }

  /// Prefix increment
  raw_lattice_iterator &operator++()
  {
    ++pos;
    return *this;
  }

  /// Prefix decrement
  raw_lattice_iterator &operator--()
  {
    --pos;
    return *this;
  }

  /// Postfix increment
  raw_lattice_iterator operator++(int)
  {
    auto copy = *this;
    ++pos;
    return copy;
  }

  /// Postfix decrement
  raw_lattice_iterator operator--(int)
  {
    auto copy = *this;
    --pos;
    return copy;
  }


  bool operator==(raw_lattice_iterator const &other) const
  {
    return (pos == other.pos and parent == other.parent);
  }

  bool operator!=(raw_lattice_iterator const &other) const
  {
    return not((*this) == other);
  }

  raw_lattice_iterator operator[](int offset) const
  {
    return raw_lattice_iterator(pos + offset, *parent);
  }

  // this return const if the template argument is const and non-const
  // otherwise:
  typename std::conditional_t<std::is_const<T>::value,
                              std::add_const_t<typename T::scalartype> &,
                              typename T::scalartype &>
  operator*() const
  {
    return parent->operator[](pos);
  }

  raw_lattice_iterator &advanceInDirection(size_t dir, int n = 1)
  {
    // moves the iterator one unit in direction dir, i.e.
    // if dir=0, and the current position is (n1,n2,...), the position will be
    // (n1+1,n2,...)
    // after the call
    Direction fbwd;
    size_t step = 0;

    if (n > 0) {
      fbwd = Direction::FORWARD;
      step = static_cast<size_t>(n);
    }
    else if (n < 0) {
      fbwd = Direction::BACKWARD;
      step = static_cast<size_t>(-n);
    }
    else {
      throw std::runtime_error("advancing by 0 is meaningless.");
    }

    if (dir >= T::dimensions())
      throw std::runtime_error("direction out of range.");


    // step is positive, here
    for (size_t i = 0; i < step; ++i) {
      // reset position to new point:
      pos = parent->neighbor(pos, dir, fbwd);
    }

    return *this;
  }
  /// does NOT change *this, i.e. const version of advanceInDirection
  /// (makes a copy which is cheap in this case)
  raw_lattice_iterator neighbor(size_t dir, int n = 1) const
  {
    auto res = *this;
    res.advanceInDirection(dir, n);
    return res;
  }
};

// *****************************************************************************
/// SiteLattice iterator class
/// \brief Iterator for SiteLattice class
// *****************************************************************************
template <typename T>
class SiteLatticeIterator
  : public std::iterator<std::bidirectional_iterator_tag, T> {
 private:
  size_t pos;
  T *parent;

 public:
  // ***************************************************************************
  // Constructors
  // ***************************************************************************
  explicit SiteLatticeIterator(size_t pos_, T &parent_)
    : pos(pos_), parent(&parent_)
  {
  }

  explicit SiteLatticeIterator(T &parent_) : pos(0ul), parent(&parent_) {}

  /// Copy constructor
  // SiteLatticeIterator(const SiteLatticeIterator &orig)
  //     : pos(orig.pos), parent(orig.parent) {}


  // ***************************************************************************
  // Deconstructor
  // ***************************************************************************
  ~SiteLatticeIterator() = default;

  // ***************************************************************************
  // Member functions
  // ***************************************************************************
  size_t index() const { return pos; }

  std::array<size_t, T::dimensions()> coord() const
  {
    return parent->linearIndexToCoord(pos);
  }

  /// Prefix increment
  SiteLatticeIterator &operator++()
  {
    ++pos;
    return *this;
  }

  /// Prefix decrement
  SiteLatticeIterator &operator--()
  {
    --pos;
    return *this;
  }

  /// Postfix increment
  SiteLatticeIterator operator++(int)
  {
    auto copy = *this;
    ++pos;
    return copy;
  }

  /// Postfix decrement
  SiteLatticeIterator operator--(int)
  {
    auto copy = *this;
    --pos;
    return copy;
  }

  bool operator==(SiteLatticeIterator const &other) const
  {
    return (pos == other.pos and parent == other.parent);
  }

  bool operator!=(SiteLatticeIterator const &other) const
  {
    return not((*this) == other);
  }

  SiteLatticeIterator operator[](int offset) const
  {
    return SiteLatticeIterator(pos + offset, *parent);
  }

  // this return const if the template argument is const and non-const
  // otherwise:
  typename std::conditional_t<std::is_const<T>::value,
                              std::add_const_t<typename T::scalartype> &,
                              typename T::scalartype &>
  operator*() const
  {
    return parent->operator[](pos);
  }

  SiteLatticeIterator &advanceInDirection(size_t dir, int n = 1)
  {
    // moves the iterator one unit in direction dir, i.e.
    // if dir=0, and the current position is (n1,n2,...), the position will be
    // (n1+1,n2,...)
    // after the call
    Direction fbwd;
    size_t step = 0;

    if (n > 0) {
      fbwd = Direction::FORWARD;
      step = static_cast<size_t>(n);
    }
    else if (n < 0) {
      fbwd = Direction::BACKWARD;
      step = static_cast<size_t>(-n);
    }
    else {
      throw std::runtime_error("advancing by 0 is meaningless.");
    }

    if (dir >= T::dimensions()) {
      throw std::runtime_error("direction out of range.");
    }
    // step is positive, here
    for (size_t i = 0; i < step; ++i) {
      // reset position to new point:
      pos = parent->neighbor(pos, dir, fbwd);
    }
    return *this;
  }

  /// does NOT change *this, i.e. const version of advanceInDirection
  /// (makes a copy which is cheap in this case)
  SiteLatticeIterator neighbor(size_t dir, int n = 1) const
  {
    auto res = *this;
    res.advanceInDirection(dir, n);
    return res;
  }

  /// Get site field at current pos
  auto site() -> typename T::scalartype const
  {
    return parent->operator[](pos);
  }

  template <typename P = double>
  auto neighborSite(size_t dir, Direction fbwd,
                    const BoundaryCondition<P, T::dimensions()> &bc =
                      BoundaryCondition<P, T::dimensions()>()) ->
    typename T::scalartype const
  {
    if (dir >= T::dimensions()) {
      throw std::runtime_error("direction out of range.");
    }

    // Get neighbor index and and coordinates  "
    auto nidx = parent->neighbor(pos, dir, fbwd);
    auto ncoord = parent->linearIndexToCoord(nidx);

    auto field = parent->operator[](nidx);

    // Implement boundary cond.
    if (ncoord[dir] == 0 && Direction::FORWARD == fbwd) {
      // We have crossed the boundary in forward direction
      // CAUTION!
      // bc and link may not commute
      // Multiply link from left for forward crossing
      return bc[dir] * field;
    }
    else if (ncoord[dir] == (parent->dimensionsArray())[dir] - 1 &&
             Direction::BACKWARD == fbwd) {
      // We have crossed the boundary in backward direction
      // CAUTION!
      // bc and link may not commute
      // Multiply link from right for backward crossing
      return field * bc[dir];
    }
    else {
      // No boundary crossed
      return field;
    }
  }
};


// *****************************************************************************
/// LinkLattice iterator class
/// \brief  Iterator for LinkLattice class
/// \todo Common base class for lattice iterators?
// *****************************************************************************
template <typename T>
class LinkLatticeIterator
  : public std::iterator<std::bidirectional_iterator_tag, T> {
 private:
  size_t pos;
  T *parent;

 public:
  // ***************************************************************************
  // Constructors
  // ***************************************************************************
  explicit LinkLatticeIterator(size_t pos_, T &parent_)
    : pos(pos_), parent(&parent_)
  {
  }

  explicit LinkLatticeIterator(T &parent_) : pos(0ul), parent(&parent_) {}


  // ***************************************************************************
  // Deconstructor
  // ***************************************************************************
  ~LinkLatticeIterator() = default;

  // ***************************************************************************
  // Member functions
  // ***************************************************************************
  size_t index() const { return pos; }

  std::array<size_t, T::dimensions()> coord() const
  {
    return parent->linearIndexToCoord(pos);
  }

  /// Prefix increment
  LinkLatticeIterator &operator++()
  {
    ++pos;
    return *this;
  }

  /// Prefix decrement
  LinkLatticeIterator &operator--()
  {
    --pos;
    return *this;
  }

  /// Postfix increment
  LinkLatticeIterator operator++(int)
  {
    auto copy = *this;
    ++pos;
    return copy;
  }

  /// Postfix decrement
  LinkLatticeIterator operator--(int)
  {
    auto copy = *this;
    --pos;
    return copy;
  }

  bool operator==(LinkLatticeIterator const &other) const
  {
    return (pos == other.pos and parent == other.parent);
  }

  bool operator!=(LinkLatticeIterator const &other) const
  {
    return not((*this) == other);
  }

  LinkLatticeIterator operator[](size_t offset) const
  {
    return LinkLatticeIterator(pos + offset, *parent);
  }

  // this returns const if the template argument is const and non-const
  // otherwise:
  typename std::conditional_t<std::is_const<T>::value,
                              std::add_const_t<typename T::linkarraytype> &,
                              typename T::linkarraytype &>
  operator*() const
  {
    return parent->operator[](pos);
  }

  /// Moves along SITES !!
  LinkLatticeIterator &advanceInDirection(size_t dir, int n = 1)
  {
    // moves the iterator one unit in direction dir, i.e.
    // if dir=0, and the current position is (n1,n2,...), the position will be
    // (n1+1,n2,...)
    // after the call
    Direction fbwd;
    size_t step = 0;

    if (n > 0) {
      fbwd = Direction::FORWARD;
      step = static_cast<size_t>(n);
    }
    else if (n < 0) {
      fbwd = Direction::BACKWARD;
      step = static_cast<size_t>(-n);
    }
    else {
      throw std::runtime_error("advancing by 0 is meaningless.");
    }

    // step is positive, here
    for (size_t i = 0; i < step; ++i) {
      // reset position to new point:
      pos = parent->neighbor(pos, dir, fbwd);
    }
    return *this;
  }

  /// does NOT change *this, i.e. const version of advanceInDirection
  /// (makes a copy which is cheap in this case)
  LinkLatticeIterator neighbor(size_t dir, int n = 1) const
  {
    auto res = *this;
    res.advanceInDirection(dir, n);
    return res;
  }

  /// Get the link of the neighbor in dir direction that points in
  /// link_dir direction.
  /// NOTE that link_dir is always positive here.
  template <typename P = double>
  auto neighborLink(size_t dir, size_t link_dir, Direction fbwd,
                    const BoundaryCondition<P, T::dimensions()> &bc =
                      BoundaryCondition<P, T::dimensions()>()) ->
    typename T::linktype const
  {
    if (dir >= T::dimensions()) {
      throw std::runtime_error("direction out of range.");
    }
    if (link_dir >= T::dimensions()) {
      throw std::runtime_error("link direction out of range.");
    }


    // Get neighbor index and and coordinates
    auto nidx = parent->neighbor(pos, dir, fbwd);
    auto ncoord = parent->linearIndexToCoord(nidx);

    auto linkArray = parent->operator[](nidx);
    auto link = linkArray[link_dir];

    // Implement boundary cond
    if (ncoord[dir] == 0 && Direction::FORWARD == fbwd) {
      // We have crossed the boundary in forward direction
      // CAUTION!
      // bc and link may not commute
      // Multiply link from left for forward crossing
      return bc[dir] * link;
    }
    else if (ncoord[dir] == (parent->dimensionsArray())[dir] - 1 &&
             Direction::BACKWARD == fbwd) {
      // We have crossed the boundary in backward direction
      // CAUTION!
      // bc and link may not commute
      // Multiply link from right for backward crossing
      return link * bc[dir];
    }
    else {
      // No boundary crossed
      return link;
    }
  }


  template <typename P = double>
  auto neighborLinkArray(size_t dir, Direction fbwd,
                         const BoundaryCondition<P, T::dimensions()> &bc =
                           BoundaryCondition<P, T::dimensions()>()) ->
    typename T::linkarraytype const
  {

    if (dir >= T::dimensions()) {
      throw std::runtime_error("direction out of range.");
    }

    // Get neighbor index and and coordinates
    auto nidx = parent->neighbor(pos, dir, fbwd);
    auto ncoord = parent->linearIndexToCoord(nidx);

    // Get Links
    auto linkArray = parent->operator[](nidx);

    // Implement boundary cond.

    if (ncoord[dir] == 0 && Direction::FORWARD == fbwd) {
      // We have crossed the boundary in forward direction
      // CAUTION!
      // bc and link may not commute
      // Multiply link from left for forward crossing
      // Enforce boundary cond. for all links
      for (auto i = 0lu; i < T::dimensions(); ++i) {
        linkArray[i] = bc[dir] * linkArray[i];
      }
      return linkArray;
    }
    else if (ncoord[dir] == (parent->dimensionsArray())[dir] - 1 &&
             Direction::BACKWARD == fbwd) {
      // We have crossed the boundary in backward direction
      // CAUTION!
      // bc and link may not commute
      // Multiply link from right for backward crossing
      for (auto i = 0lu; i < T::dimensions(); ++i) {
        linkArray[i] = linkArray[i] * bc[dir];
      }
      return linkArray;
    }
    else {
      // No boundary crossed
      return linkArray;
    }
  }

  /// Get link at current pos in direction dir
  auto link(size_t dir) -> typename T::linktype const
  {
    return (parent->operator[](pos))[dir];
  }

  /// Get link array at current pos
  auto linkarray() -> typename T::linkarraytype const
  {
    return parent->operator[](pos);
  }
};


template <typename T>
class FullLatticeIterator
  : public std::iterator<std::bidirectional_iterator_tag, T> {

 private:
  size_t pos;
  T *parent;

 public:
  // ***************************************************************************
  // Constructors
  // ***************************************************************************
  explicit FullLatticeIterator(size_t pos_, T &parent_)
    : pos(pos_), parent(&parent_)
  {
  }

  explicit FullLatticeIterator(T &parent_) : pos(0ul), parent(&parent_) {}


  // ***************************************************************************
  // Member functions
  // ***************************************************************************
  size_t index() const { return pos; }

  std::array<size_t, T::dimensions()> coord() const
  {
    return parent->linearIndexToCoord(pos);
  }

  /// Prefix increment
  FullLatticeIterator &operator++()
  {
    ++pos;
    return *this;
  }

  /// Prefix decrement
  FullLatticeIterator &operator--()
  {
    --pos;
    return *this;
  }

  /// Postfix increment
  FullLatticeIterator operator++(int)
  {
    auto copy = *this;
    ++pos;
    return copy;
  }

  /// Postfix decrement
  FullLatticeIterator operator--(int)
  {
    auto copy = *this;
    --pos;
    return copy;
  }

  /// Equality check
  bool operator==(FullLatticeIterator const &other) const
  {
    return (pos == other.pos and parent == other.parent);
  }

  bool operator!=(FullLatticeIterator const &other) const
  {
    return not((*this) == other);
  }

  FullLatticeIterator operator[](size_t offset) const
  {
    return FullLatticeIterator(pos + offset, *parent);
  }

  // this returns const if the template argument is const and non-const
  // otherwise:
  typename std::conditional_t<std::is_const<T>::value,
                              std::add_const_t<typename T::fieldstructtype> &,
                              typename T::fieldstructtype &>
  operator*() const
  {
    return parent->operator[](pos);
  }

  /// Moves along lattice sites
  FullLatticeIterator &advanceInDirection(size_t dir, int n = 1)
  {
    // moves the iterator one unit in direction dir, i.e.
    // if dir=0, and the current position is (n1,n2,...), the position will be
    // (n1+1,n2,...)
    // after the call
    Direction fbwd;
    size_t step = 0;

    if (n > 0) {
      fbwd = Direction::FORWARD;
      step = static_cast<size_t>(n);
    }
    else if (n < 0) {
      fbwd = Direction::BACKWARD;
      step = static_cast<size_t>(-n);
    }
    else {
      throw std::runtime_error("advancing by 0 is meaningless.");
    }

    if (dir >= T::dimensions()) {
      throw std::runtime_error("direction out of range.");
    }

    // step is positive, here
    for (size_t i = 0; i < step; ++i) {
      // reset position to new point:
      pos = parent->neighbor(pos, dir, fbwd);
    }
    return *this;
  }

  /// does NOT change *this, i.e. const version of advanceInDirection
  /// (makes a copy which is cheap in this case)
  FullLatticeIterator neighbor(size_t dir, int n = 1) const
  {
    auto res = *this;
    res.advanceInDirection(dir, n);
    return res;
  }

  /// Get site field at current pos
  auto site() -> typename T::sitetype const
  {
    return parent->operator[](pos).site;
  }

  /// Get link at current pos in direction dir
  auto link(size_t dir) -> typename T::linktype const
  {
    return parent->operator[](pos).links[dir];
  }

  /// Get link array at current pos
  auto linkarray() -> typename T::linkarraytype const
  {
    return parent->operator[](pos).links;
  }


  template <typename P = double>
  auto neighborSite(size_t dir, Direction fbwd,
                    const BoundaryCondition<P, T::dimensions()> &bc =
                      BoundaryCondition<P, T::dimensions()>()) ->
    typename T::sitetype const
  {
    if (dir >= T::dimensions()) {
      throw std::runtime_error("direction out of range.");
    }

    // Get neighbor index and and coordinates
    auto nidx = parent->neighbor(pos, dir, fbwd);
    auto ncoord = parent->linearIndexToCoord(nidx);

    auto field = parent->operator[](nidx).site;


    // Implement boundary cond.
    if (ncoord[dir] == 0 && Direction::FORWARD == fbwd) {
      // We have crossed the boundary in forward direction
      // CAUTION!
      // bc and link may not commute
      // Multiply link from left for forward crossing
      return bc[dir] * field;
    }
    else if (ncoord[dir] == (parent->dimensionsArray())[dir] - 1 &&
             Direction::BACKWARD == fbwd) {
      // We have crossed the boundary in backward direction
      // CAUTION!
      // bc and link may not commute
      // Multiply link from right for backward crossing
      return field * bc[dir];
    }
    else {
      // No boundary crossed
      return field;
    }
  }


  /// Get the link of the neighbor in dir direction that points in
  /// link_dir direction.
  /// NOTE that link_dir is always positive here.
  template <typename P = double>
  auto neighborLink(size_t dir, size_t link_dir, Direction fbwd,
                    const BoundaryCondition<P, T::dimensions()> &bc =
                      BoundaryCondition<P, T::dimensions()>()) ->
    typename T::linktype const
  {
    if (dir >= T::dimensions()) {
      throw std::runtime_error("direction out of range.");
    }
    if (link_dir >= T::dimensions()) {
      throw std::runtime_error("link direction out of range.");
    }

    // Get neighbor index and and coordinates
    auto nidx = parent->neighbor(pos, dir, fbwd);
    auto ncoord = parent->linearIndexToCoord(nidx);

    auto linkArray = parent->operator[](nidx).links;
    auto link = linkArray[link_dir];

    if (ncoord[dir] == 0 && Direction::FORWARD == fbwd) {
      // We have crossed the boundary in forward direction
      // CAUTION!
      // bc and link may not commute
      // Multiply link from left for forward crossing
      return bc[dir] * link;
    }
    else if (ncoord[dir] == (parent->dimensionsArray())[dir] - 1 &&
             Direction::BACKWARD == fbwd) {
      // We have crossed the boundary in backward direction
      // CAUTION!
      // bc and link may not commute
      // Multiply link from right for backward crossing
      return link * bc[dir];
    }
    else {
      // No boundary crossed
      return link;
    }
  }


  template <typename P = double>
  auto neighborLinkArray(size_t dir, Direction fbwd,
                         const BoundaryCondition<P, T::dimensions()> &bc =
                           BoundaryCondition<P, T::dimensions()>()) ->
    typename T::linkarraytype const
  {

    if (dir >= T::dimensions()) {
      throw std::runtime_error("direction out of range.");
    }

    // Get neighbor index and and coordinates
    auto nidx = parent->neighbor(pos, dir, fbwd);
    auto ncoord = parent->linearIndexToCoord(nidx);

    auto linkArray = parent->operator[](nidx).links;

    // Implement boundary cond.

    if (ncoord[dir] == 0 && Direction::FORWARD == fbwd) {
      // We have crossed the boundary in forward direction
      // CAUTION!
      // bc and link may not commute
      // Multiply link from left for forward crossing
      // Enforce boundary cond. for all links
      for (auto i = 0lu; i < T::dimensions(); ++i) {
        linkArray[i] = bc[dir] * linkArray[i];
      }
      return linkArray;
    }
    else if (ncoord[dir] == (parent->dimensionsArray())[dir] - 1 &&
             Direction::BACKWARD == fbwd) {
      // We have crossed the boundary in backward direction
      // CAUTION!
      // bc and link may not commute
      // Multiply link from right for backward crossing
      for (auto i = 0lu; i < T::dimensions(); ++i) {
        linkArray[i] = linkArray[i] * bc[dir];
      }
      return linkArray;
    }
    else {
      // No boundary crossed
      return linkArray;
    }
  }
};

// *****************************************************************************
/// Lattice base class
///
/// This class takes care of indexing and projecting the lattice onto a
/// torus. It also provides methods to compute the lattice volume,
/// number of dimensions, etc. ...
// *****************************************************************************
template <size_t dim>
class BareLattice {

  // Private, except for subclasses
 protected:
  std::array<size_t, dim> dims; ///< Array holding the Lattice dimensions
  /// Look-up table for converting lattice index to coordinates
  std::vector<std::array<size_t, dim>> lut_indexToCoord;
  /// Look-up table for neighbours:  lut_neighbor[bwd/fwd][site][dir]
  std::array<std::vector<std::array<size_t, dim>>, 2> lut_neighbors;
  /// Wrap lattice around a torus
  void toTorus(std::array<size_t, dim> &coord) const
  {
    for (size_t i = 0; i < dim; ++i) {
      coord[i] = (coord[i] + dims[i]) % dims[i];
    }
  }

  /// Compute look-up table for Sites
  void fillSiteLut()
  {
    // Sanity check
    // assert(data.size() == product(dims));
    size_t lattice_size = product(dims);
    // Resize lut (no need to initialise)
    lut_indexToCoord.resize(lattice_size);
    std::array<size_t, dim> pos;
    // start at i = 0 <==> pos=(0,0,0,0...)
    for (auto &d : pos) {
      d = 0;
    }
    // Go over all indices
    for (size_t i = 0; i < lattice_size; ++i) {
      // For a given index construct the corresponding lattice
      // position

      // current position in cartesian coordinates is pos
      lut_indexToCoord[i] = pos;

      // Increase position index. The lowest index is running fastest
      pos[0] += 1;
      // Shift increase to appropriate coordinate
      for (size_t idim = 0; idim < dim - 1; ++idim) {

        if (pos[idim] == dims[idim]) {
          pos[idim] = 0;      // Wrap around and
          pos[idim + 1] += 1; // go to higher index
        }
      }
    }
  }

  /// Compute look-up table for neighbours
  void fillNeighborSiteLut()
  {
    // Sanity check
    // assert(lut_indexToCoord.size() == data.size());
    size_t lattice_size = product(dims);
    // Resize lut (no need to initialise)
    lut_neighbors[Direction::FORWARD].resize(lattice_size);
    lut_neighbors[Direction::BACKWARD].resize(lattice_size);

    for (size_t i = 0; i < lut_indexToCoord.size(); ++i) {
      for (size_t idim = 0; idim < dim; ++idim) {
        auto pos_cart = lut_indexToCoord[i];
        // FORWARD:
        pos_cart[idim] += 1; // Neighbour in idim direction
        toTorus(pos_cart);   // Make sure we stay on torus
        lut_neighbors[Direction::FORWARD][i][idim] =
          coordToLinearIndex(pos_cart); // Get index of neighbour
        // BACKWARD:
        // pos_cart[idim] = lut_indexToCoord[i][idim] - 1;
        pos_cart = lut_indexToCoord[i];
        pos_cart[idim] -= 1; // Neighbour in idim direction
        toTorus(pos_cart);   // Make sure we stay on torus
        lut_neighbors[Direction::BACKWARD][i][idim] =
          coordToLinearIndex(pos_cart); // Get index of neighbour
      }
    }
  }

 public:
  // ***************************************************************************
  // Constructors
  // ***************************************************************************

  /// Constructor
  BareLattice(std::array<size_t, dim> const &_Dim) : dims(_Dim)
  {
    // Calculate Look-Up Tables
    fillSiteLut();
    fillNeighborSiteLut();
  }

  /// Default constructor
  BareLattice()
  {
    dims = std::array<size_t, dim>();
    dims.fill(0);

    fillSiteLut();
    fillNeighborSiteLut();
  }

  /// Copy constructor
  BareLattice(BareLattice<dim> const &orig)
  {
    dims = orig.dims;
    lut_indexToCoord = orig.lut_indexToCoord;
    lut_neighbors = orig.lut_neighbors;
  }

  // ***************************************************************************
  // Destructor
  // ***************************************************************************
  virtual ~BareLattice() = default;


  // ***************************************************************************
  // Member functions
  // ***************************************************************************

  /// get lattice volume
  size_t volume() const { return product(dims); }

  /// get number of lattice dimensions
  static inline constexpr size_t dimensions() { return dim; }

  /// get array with dimensions
  inline std::array<size_t, dim> dimensionsArray() const { return dims; }

  /// Convert lattice index to lattice coordinate array
  inline std::array<size_t, dim> linearIndexToCoord(size_t coord) const
  {
    return lut_indexToCoord[coord];
  }

  /// Get index of neighbor
  inline size_t neighbor(size_t index, size_t mu, Direction dir) const
  {
    return lut_neighbors[dir][index][mu];
  }

  /// get lattice coordinates for a given lattice index
  inline size_t coordToLinearIndex(std::array<size_t, dim> coord) const
  {
    // first coord is fastest running, last coord is slowest running:
    // N[0]*N[1]*...*N[dim-2]*coord[dim-1]+
    // N[0]*N[1]*...*N[dim-3]*coord[dim-2]+
    // ..
    // N[0]*coord[1]+
    // coord[0] = index
    // Faster, but less clear  (save multiplications)
    // index = coord[0]+N[0]*(coord[1]+N[1]*(coord[2]+ N[2]*(coord[3]+ ...
    // N[idim-1]*coord[idim])
    size_t tmp = coord[static_cast<size_t>(dim - 1)]; // coord[idim]
    for (size_t idim = dim - 1; idim > 0; --idim) {
      tmp = (coord[static_cast<size_t>(idim - 1)] + dims[idim - 1] * tmp);
    }
    return tmp;
  }
};


// *****************************************************************************
/// Class for Lattice where the fields live only on the lattice sites
///
/// The class is implemented as a template and kept as generic as
/// possible.
// *****************************************************************************
template <typename T, size_t dim>
class SiteLattice : public BareLattice<dim> {
 private:
  std::vector<T> data; ///< Stores the fields

 public:
  // ***************************************************************************
  // Typedefs
  // ***************************************************************************

  typedef SiteLatticeIterator<const SiteLattice<T, dim>> const_iterator;
  typedef SiteLatticeIterator<SiteLattice<T, dim>> iterator;
  // typedef raw_lattice_iterator<const SiteLattice<T, dim>> const_iterator;
  // typedef raw_lattice_iterator<      SiteLattice<T, dim>> iterator;
  typedef T scalartype;

  // ***************************************************************************
  // Constructors
  // ***************************************************************************

  SiteLattice(std::array<size_t, dim> const &Dim_, T const &init)
    : BareLattice<dim>(Dim_), data(product(Dim_), init)
  {
    // Sanity check
    assert(data.size() == product(this->dims));

    // And another sanity check
    assert(this->lut_indexToCoord.size() == data.size());
  }

  /// Construct default  (T zero initialised) site lattice with given
  /// dimensions
  SiteLattice(std::array<size_t, dim> const &Dim_) : SiteLattice(Dim_, T(0)) {}

  /// Default constructor
  SiteLattice() : BareLattice<dim>()
  {

    std::vector<T> vtmp;
    vtmp.resize(this->volume());

    data = vtmp;
  }

  /// Copy constructor
  /// \Todo Explicit copy ctor breaks test_io, which seems fishy
  // SiteLattice(SiteLattice<T,dim> const &orig)
  //     : BareLattice<dim>(orig)
  // {
  //     data = orig.data;
  // }

  // ***************************************************************************
  // Destructor
  // ***************************************************************************
  ~SiteLattice() = default;

  // ***************************************************************************
  // Member Functions
  // ***************************************************************************

  inline T const &operator[](size_t i) const { return data[i]; }
  inline T &operator[](size_t i) { return data[i]; }


  // Iterators
  iterator begin() { return iterator(0, *this); }
  iterator end() { return iterator(data.size(), *this); }
  iterator at(std::array<size_t, dim> const &point)
  {
    return iterator(this->coordToLinearIndex(point), *this);
  }

  const_iterator begin() const { return const_iterator(0, *this); }
  const_iterator end() const
  {
    return const_iterator(data.size(),
                          static_cast<SiteLattice<T, dim> const &>(*this));
  }
  const_iterator at(std::array<size_t, dim> const &point) const
  {
    return const_iterator(this->coordToLinearIndex(point), *this);
  }
};

// *****************************************************************************
/// Class for Lattice where the fields live only on the lattice links
///
/// The class is implemented as a template and kept as generic as
/// possible.
// *****************************************************************************
template <typename T, size_t dim>
class LinkLattice : public BareLattice<dim> {
 private:
  /// Stores the fields as array of dim links per site
  std::vector<std::array<T, dim>> data;

 public:
  // ***************************************************************************
  // Typedefs
  // ***************************************************************************
  typedef LinkLatticeIterator<const LinkLattice<T, dim>> const_iterator;
  typedef LinkLatticeIterator<LinkLattice<T, dim>> iterator;
  typedef T linktype;
  typedef std::array<T, dim> linkarraytype;


  // ***************************************************************************
  // Constructors
  // ***************************************************************************

  /// Constructor to initialise all LINKS with init
  LinkLattice(std::array<size_t, dim> const &Dim_, T const &init)
    : BareLattice<dim>(Dim_)
  {

    // Make an array filled with init
    std::array<T, dim> init_ar;
    init_ar.fill(init);

    // Fill the vector with our array
    this->data =
      std::move(std::vector<std::array<T, dim>>(product(this->dims), init_ar));

    // Sanity check
    assert(this->lut_indexToCoord.size() == data.size());
  }

  /// Constructor to initialise all SITES with the LINK field
  LinkLattice(std::array<size_t, dim> const &Dim_,
              std::array<T, dim> const &init_ar)
    : BareLattice<dim>(Dim_)
  {
    // Fill the vector with our array
    data = std::vector<std::array<T, dim>>(product(this->dims), init_ar);

    // Sanity check
    assert(this->lut_indexToCoord.size() == data.size());
  }

  /// Constructor for default  (T zero initialised) link lattice with given
  /// dimenions
  LinkLattice(std::array<size_t, dim> const &Dim_) : LinkLattice(Dim_, T(0)) {}

  /// Default constructor
  LinkLattice() : BareLattice<dim>() {}

  /// Copy constructor
  // LinkLattice(LinkLattice<T,dim> const &orig)
  //     : BareLattice<dim>(orig)
  // {
  //     data = orig.data;
  // }

  // ***************************************************************************
  // Destructor
  // ***************************************************************************
  ~LinkLattice() = default;

  // ***************************************************************************
  // Member Functions
  // ***************************************************************************

  /// get array of all links at the site i
  inline std::array<T, dim> const &operator[](size_t i) const
  {
    return this->data[i];
  }
  /// get array of all links at the site i
  inline std::array<T, dim> &operator[](size_t i) { return this->data[i]; }

  /// get link in direction nu at site i
  inline T const &operator()(size_t i, size_t nu) const
  {
    if (nu >= this->dimensions()) {
      throw std::runtime_error("nu out of range.");
    }
    return this->data[i][nu];
  }
  /// get link in direction nu at site i
  inline T &operator()(size_t i, size_t nu)
  {
    if (nu >= this->dimensions()) {
      throw std::runtime_error("nu out of range.");
    }
    return this->data[i][nu];
  }

  // Iterators
  // Note that this iterates over sites!
  iterator begin() { return iterator(0, *this); }
  iterator end() { return iterator(data.size(), *this); }
  iterator at(std::array<size_t, dim> const &point)
  {
    return iterator(this->coordToLinearIndex(point), *this);
  }

  const_iterator begin() const { return const_iterator(0, *this); }
  const_iterator end() const
  {
    return const_iterator(data.size(),
                          static_cast<LinkLattice<T, dim> const &>(*this));
  }
  const_iterator at(std::array<size_t, dim> const &point) const
  {
    return const_iterator(this->coordToLinearIndex(point), *this);
  }
};


// *****************************************************************************
/// Struct to hold site and link fields. Helper for FullLattice class.
// *****************************************************************************
template <typename S, typename L, size_t dim>
struct SiteAndLinks {
  S site;
  std::array<L, dim> links;
};

/// Printing the SiteAndLinks Struct

template <typename S, typename L, size_t dim>
std::ostream &operator<<(std::ostream &stream,
                         SiteAndLinks<S, L, dim> const &sl)
{
  stream.precision(6);
  stream << std::scientific;

  stream << "Site:\n" << sl.site << std::endl;
  stream << "Links: \n";

  for (auto l : sl.links) {
    stream << l << std::endl;
  }

  return stream;
}


// *****************************************************************************
/// Class for Lattice with fields on the lattice sites AND links
///
/// The class is implemented as a template and kept as generic as
/// possible.
///
/// \TODO The name "FullLattice" is stupid. Rename this simply to
///       "Lattice" once the dust has settled?
// *****************************************************************************
template <typename S, typename L, size_t dim>
class FullLattice : public BareLattice<dim> {

 private:
  /// Store links and site in vector of  SiteAndLinks structs
  std::vector<SiteAndLinks<S, L, dim>> data;

 public:
  // ***************************************************************************
  // Typedefs
  // ***************************************************************************
  typedef FullLatticeIterator<const FullLattice<S, L, dim>> const_iterator;
  typedef FullLatticeIterator<FullLattice<S, L, dim>> iterator;
  typedef SiteAndLinks<S, L, dim> fieldstructtype;
  typedef S sitetype;
  typedef L linktype;
  typedef std::array<L, dim> linkarraytype;


  // ***************************************************************************
  // Constructors
  // ***************************************************************************

  /// Default constructor
  FullLattice() : BareLattice<dim>()
  {
    // Reserve vector memory
    data.resize(this->volume());
    // std::array has no well defined default constructor
    std::array<L, dim> ar;
    ar.fill(L()); // Fill array with default of L

    for (auto &&el : data) {

      el.site = S();
      el.links = ar;
    }


    // Finally init link vector
    // ldata = std::vector<std::array<L,dim>>();
    // ldata.resize(this->volume());

    // Sanity check
    assert(this->lut_indexToCoord.size() == data.size());
    // assert(this->lut_indexToCoord.size() == ldata.size());
  }

  /// Constructor to initialise Lattice with given site and link
  FullLattice(std::array<size_t, dim> const &Dim_, S const &sinit,
              L const &linit)
    : BareLattice<dim>(Dim_)
  {
    // Reserve vector memory
    data.resize(this->volume());
    // Make an array filled with init
    std::array<L, dim> linit_ar;
    linit_ar.fill(linit);

    // Fill data vector
    for (auto &&el : data) {

      el.site = sinit;
      el.links = linit_ar;
    }

    // Sanity check
    assert(this->lut_indexToCoord.size() == data.size());
  }

  /// Constructor to initialise Lattice with given site and link array
  FullLattice(std::array<size_t, dim> const &Dim_, S const &sinit,
              const std::array<L, dim> &linit_ar)
    : BareLattice<dim>(Dim_)
  {
    // Reserve vector memory
    data.resize(this->volume());

    // Fill data vector
    for (auto &&el : data) {
      el.site = sinit;
      el.links = linit_ar;
    }


    // Sanity check
    assert(this->lut_indexToCoord.size() == data.size());
  }

  /// Construct default (S and L zero initialised) lattice with given dimensions
  FullLattice(std::array<size_t, dim> const &Dim_)
    : FullLattice(Dim_, S(0), L(0))
  {
  }

  /// Copy constructor
  // 	FullLattice(const FullLattice &orig)
  // 	    : BareLattice<dim>(orig)
  // {
  // 	    data = orig.data;
  // }

  // ***************************************************************************
  /// Deconstructor
  // ***************************************************************************
  ~FullLattice() = default;

  // ***************************************************************************
  // Member Functions
  // ***************************************************************************

  /// get array of all links at the site i
  inline SiteAndLinks<S, L, dim> const &operator[](size_t i) const
  {
    return this->data[i];
  }
  /// get array of all links at the site i
  inline SiteAndLinks<S, L, dim> &operator[](size_t i) { return this->data[i]; }

  /// Get Link at site idx in direction dir
  inline L const &getLink(size_t idx, size_t dir) const
  {
    return this->data[idx].links[dir];
  }

  /// Get site at index  idx
  inline S const &getSite(size_t idx) const { return this->data[idx].site; }


  // Iterators
  iterator begin() { return iterator(0, *this); }
  iterator end() { return iterator(data.size(), *this); }
  iterator at(std::array<size_t, dim> const &point)
  {
    return iterator(this->coordToLinearIndex(point), *this);
  }

  const_iterator begin() const { return const_iterator(0, *this); }
  const_iterator end() const
  {
    return const_iterator(data.size(),
                          static_cast<FullLattice<S, L, dim> const &>(*this));
  }
  const_iterator at(std::array<size_t, dim> const &point) const
  {
    return const_iterator(this->coordToLinearIndex(point), *this);
  }
};

/// Lattice class
/// \brief Implements a generic lattice object
///
/// The class is implemented as a template and kept as generic as
/// possible.
// *****************************************************************************
// template <typename T, size_t dim> class Lattice: public BareLattice<dim> {

//     private:
// 	std::vector<T>                data;         ///< Stores the fields
// 	BoundaryCondition<dim>        bc;           ///< Boundary condition
// 	// Look-up table for boundary:  lut_boundary[dir][bwd/fwd][pos] ?

//     public:
// 	// Typedefs
// 	typedef raw_lattice_iterator<const Lattice<T, dim>> const_iterator;
// 	typedef raw_lattice_iterator<Lattice<T, dim>> iterator;
// 	typedef T scalartype;

// 	// **********************************************************************
// 	// Constructors
// 	// **********************************************************************

// 	Lattice(std::array<size_t, dim> const &Dim_, T const &init,
// 		BoundaryCondition<dim> _bc =  BoundaryCondition<dim>())
// 	    : BareLattice<dim>(Dim_),
// 	      data(product(Dim_), init),
// 	      bc(_bc)
// 	    {
// 		// Sanity check
// 		assert(data.size() == product(this->dims));

// 		// Calculate Look-Up Tables
// 		this->fillSiteLut();
// 		this->fillNeighborSiteLut();

// 		// And another sanity check
// 		assert(this->lut_indexToCoord.size() == data.size());
// 	    }

// 	/// \TODO: Check constructors
// 	Lattice(std::array<size_t, dim> const &Dim_)
// 	    : Lattice(Dim_, T(0)){}

//         /// Default constructor
// 	Lattice()
// 	    : BareLattice<dim>(){

// 	    std::vector<T> vtmp;


// 	    data = vtmp;
// 	    bc=BoundaryCondition<dim>();

// 	    this->fillSiteLut();
// 	    this->fillNeighborSiteLut();
// 	}

// 	/// Copy constructor
// 	/// \Todo Explicit copy ctor breaks test_io, which seems fishy
// 	Lattice(Lattice<T,dim> const &orig)
// 	    : BareLattice<dim>(orig)
// 	{
// 	    this->dims = orig.dims;
// 	    data = orig.data;
// 	    this->lut_indexToCoord  = orig.lut_indexToCoord;
// 	    this->lut_neighbors     = orig.lut_neighbors;
// 	    bc = orig.bc;

// 	}

// 	// **********************************************************************
// 	// Member Functions
// 	// **********************************************************************

// 	inline T const &operator[](size_t i) const { return data[i]; }
// 	inline T &operator[](size_t i) { return data[i]; }


// 	inline T get_neighbor_field(size_t index, size_t bearing, Direction step)
// const{

// 	    // sanity check
// 	    if (bearing >= this->dimensions()){
// 		throw std::runtime_error("direction out of range.");
// 	    }


// 	    if ( bc.type() == BCtype::SITE)
// 	    {


// 		// Make step in the direction given by 'bearing'
// 		// and make sure we stay on the torus

// 		auto idx = this->neighbor(index,bearing, step);
// 		auto pos = this->linearIndexToCoord(idx);

// 		if (pos[bearing]==0 && Direction::FORWARD==step)
// 		{
// 		    // We have crossed the boundary in forward direction
// 		    auto II  = std::complex<double>(0.,1.);
// 		    return data[idx]*std::exp(-II*bc[bearing]);
// 		}
// 		else if (pos[bearing]==this->dims[bearing]-1 &&
// Direction::BACKWARD==step)
// 		{
// 		    // We have crossed the boundary in backward direction
// 		    // NOTE THE SIGN CHANGE IN THE PHASE !!!
// 		    auto II  = std::complex<double>(0.,1.);
// 		    return data[idx]*exp(II*bc[bearing]);
// 		}
// 		else
// 		{
// 		    //No boundary crossed
// 		    return data[idx];
// 		}
// 	    }
// 	    else
// 	    {
// 		throw std::runtime_error("Only 'SITE' boundary conditions
// implemented!");
// 	    }
// 	}

// 	const_iterator begin() const { return const_iterator(0, *this); }
// 	const_iterator end() const {
// 	    return const_iterator(data.size(),
// 				  static_cast<Lattice<T, dim> const &>(*this));
// 	}
// 	iterator begin() { return iterator(0, *this); }
// 	iterator end() { return iterator(data.size(), *this); }
// 	iterator at(std::array<size_t, dim> const & point) {
// 	    return iterator(this->coordToLinearIndex(point), *this); }
// 	const_iterator at(std::array<size_t, dim> const & point) const {
// 	    return const_iterator(this->coordToLinearIndex(point), *this); }
// 	// Will implement iterators pointing to first element of shadow, boundary,
// 	// etc.
// };


template <typename T, size_t dim>
T sum(SiteLattice<T, dim> const &lat)
{
  return std::accumulate(lat.begin(), lat.end(), T(0));
}
