/// \file
/// A simple progress bar

#pragma once
using size_t = std::size_t;

/// A class to display a simple progress bar
class progress_bar {

 private:
  std::string indicator;
  std::string offset;
  double final;
  size_t width;

 public:
  // ---------------------------------------------------------------------
  // Constructors
  // ---------------------------------------------------------------------
  progress_bar()
    : indicator(std::string("=")), offset(std::string("\0")), final(100.),
      width(70LU)
  {
  }

  progress_bar(double max, size_t w = 70LU, std::string off = std::string("\0"),
               std::string ind = std::string("="))
    : indicator(ind), offset(off), final(max), width(w)
  {
  }

  // ---------------------------------------------------------------------
  /// Destructor
  // ---------------------------------------------------------------------
  ~progress_bar() = default;

  // ---------------------------------------------------------------------
  // Member functions
  // ---------------------------------------------------------------------

  void progress(const double &cur) const
  {
    double prog = cur / final;
    size_t pos = static_cast<size_t>(prog * static_cast<double>(width));

    // std::cout << cur/diff << std::endl;
    // std::cout << "diff: " << diff << "\tcur:" << cur << std::endl;
    std::cout << offset << "[";
    auto i = 0LU;
    while (i < pos) {
      std::cout << indicator;
      ++i;
    }

    if (i < width) {
      std::cout << ">";
    }
    ++i;


    while (i < width) {
      std::cout << " ";
      ++i;
    }

    std::cout << "]"
              << " " << std::fixed << std::setprecision(2) << std::setw(6)
              << std::right << prog * 100 << "%\r";
    std::cout.flush();
  }
};
