/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *  
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *  
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *  
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 *  
 *  Contact the Developers:
 *    E-Mail: xsli@uw.edu
 *  
 */
#pragma once

// System Headers (standard C++)
#include <string>        
#include <cassert>
#include <iterator>
#include <array>         
#include <unordered_map> 
#include <map> 
#include <vector>        
#include <numeric>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <exception>
#include <sstream>
#include <cctype>
#include <locale>
#include <algorithm>
#include <unistd.h>
#include <complex>
#include <memory>
#include <valarray>
#include <random>
#include <chrono>
#include <typeindex>
#include <system_error>
#include <unordered_set>
#include <set>
#include <bitset>

#include <chronusq_config.hpp> // Configuration header


#ifdef _OPENMP
  #include <omp.h>
#endif


// Standard typedefs
typedef std::complex<double> dcomplex;

