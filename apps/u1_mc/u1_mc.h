///! \file
/// Header for u1_mc.cpp

#pragma once

#include "../src/actions/compact_qed.h"
#include "../src/gaugegroups/u1.h"
#include "../src/lattice.h"
#include "../src/progress_bar.h"
#include "../src/stat.h"
#include "./u1_mc_params.h"
#include <algorithm>
#include <chrono>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>


using size_t = std::size_t;
using clck = std::chrono::high_resolution_clock;
