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

#include <coupledcluster/TAManager.hpp>
#include <utility>

//#define DEBUG_TAMANAGER

namespace ChronusQ {

  TA::TiledRange TAManager::toRange(std::string ranges) const {
    std::vector<TA::TiledRange1> range1s;
    range1s.reserve(ranges.length());
    for (char c : ranges)
      range1s.push_back(rangeTypes_.at(c));
    return TA::TiledRange(range1s.begin(), range1s.end());
  }

  std::vector<std::pair<size_t, size_t>> TAManager::toBlockRange(std::string ranges) const {
    std::vector<std::pair<size_t,size_t>> blockRanges;
    blockRanges.reserve(ranges.length());
    for (char c : ranges)
      blockRanges.push_back(blockRangeTypes_.at(c));
    return blockRanges;
  }

  size_t TAManager::elem_per_TA(std::string ranges) const {
    size_t elem = 1;
    for (char c : ranges)
      elem *= rangeTypes_.at(c).extent();
    return elem;
  }

  template<>
  TA::TArray<double> TAManager::malloc_fresh(std::string ranges) {
    if (rTAstat_.count(ranges) == 0) {
      rTAstat_[ranges] = {0, 0, 0, sizeof(double) * elem_per_TA(ranges)};
      rTAs_[ranges] = std::vector<TA::TArray<double>>();
    }
    rTAstat_[ranges].cur_total++;
    rTAstat_[ranges].hist_total++;
    rTAstat_[ranges].peak = std::max(rTAstat_[ranges].cur_total, rTAstat_[ranges].peak);

    cur_mem_ += rTAstat_[ranges].mem_each_;
    peak_mem_ = std::max(cur_mem_, peak_mem_);

    return TA::TArray<double>(TA::get_default_world(), toRange(ranges));
  }

  template<>
  TA::TArray<dcomplex> TAManager::malloc_fresh(std::string ranges) {
#ifdef DEBUG_TAMANAGER
    std::cout << "TAManager::malloc_fresh(" << ranges << ")" << std::endl;
#endif
    if (cTAstat_.count(ranges) == 0) {
      cTAstat_[ranges] = {0, 0, 0, sizeof(dcomplex) * elem_per_TA(ranges)};
      cTAs_[ranges] = std::vector<TA::TArray<dcomplex>>();
    }
    cTAstat_[ranges].cur_total++;
    cTAstat_[ranges].hist_total++;
    cTAstat_[ranges].peak = std::max(cTAstat_[ranges].cur_total, cTAstat_[ranges].peak);

    cur_mem_ += cTAstat_[ranges].mem_each_;
    peak_mem_ = std::max(cur_mem_, peak_mem_);

    return TA::TArray<dcomplex>(TA::get_default_world(), toRange(ranges));
  }

  template<>
  TA::TArray<double> TAManager::malloc(std::string ranges) {
    if (rTAs_.count(ranges) == 0 or rTAs_[ranges].size() == 0) {
      TA::TArray<double> tmp = malloc_fresh<double>(ranges);
      tmp.fill(0.0);
      return tmp;
    }

    TA::TArray<double> tmp(std::move(rTAs_[ranges].back()));

    rTAs_[ranges].pop_back();
    return tmp;
  }

  template<>
  TA::TArray<dcomplex> TAManager::malloc(std::string ranges) {
#ifdef DEBUG_TAMANAGER
    std::cout << "TAManager::malloc(" << ranges << ")" << std::endl;
#endif
    if (cTAs_.count(ranges) == 0 or cTAs_[ranges].size() == 0) {
#ifdef DEBUG_TAMANAGER
      std::cout << "|-";
#endif
      TA::TArray<dcomplex> tmp = malloc_fresh<dcomplex>(ranges);
      tmp.fill(0.0);
      return tmp;
    }

    TA::TArray<dcomplex> tmp(std::move(cTAs_[ranges].back()));

    cTAs_[ranges].pop_back();
    return tmp;
  }

  template<>
  void TAManager::free(std::string ranges, TA::TArray<double> &&ta, bool discard) {
    if (rTAs_.count(ranges) == 0)
      CErr("Requested range label does not exist in TAManager!");
    if (discard) {
      rTAstat_[ranges].cur_total--;
      cur_mem_ -= rTAstat_[ranges].mem_each_;
    } else
      rTAs_[ranges].push_back(std::forward<TA::TArray<double>>(ta));
    ta = TA::TArray<double>();
  }

  template<>
  void TAManager::free(std::string ranges, TA::TArray<dcomplex> &&ta, bool discard) {
#ifdef DEBUG_TAMANAGER
    std::cout << "TAManager::free(" << ranges << (discard ? ", discard" : ", to cache") << ")" << std::endl;
#endif
    if (cTAs_.count(ranges) == 0)
      CErr("Requested range label does not exist in TAManager!");
    if (discard) {
      cTAstat_[ranges].cur_total--;
      cur_mem_ -= cTAstat_[ranges].mem_each_;
    } else
      cTAs_[ranges].push_back(std::forward<TA::TArray<dcomplex>>(ta));
    ta = TA::TArray<dcomplex>();
  }

  void TAManager::discard_cache() {
#ifdef DEBUG_TAMANAGER
    std::cout << "TAManager::discard_cache()" << std::endl;
#endif
    for (auto &sTA : rTAs_) {
      rTAstat_[sTA.first].cur_total -= sTA.second.size();
      cur_mem_ -= sTA.second.size() * rTAstat_[sTA.first].mem_each_;
      sTA.second.clear();
    }
    for (auto &sTA : cTAs_) {
      cTAstat_[sTA.first].cur_total -= sTA.second.size();
      cur_mem_ -= sTA.second.size() * cTAstat_[sTA.first].mem_each_;
#ifdef DEBUG_TAMANAGER
      std::cout << "|-  " << sTA.first << " : -" << sTA.second.size() << std::endl;
#endif
      sTA.second.clear();
    }
    TA::get_default_world().gop.fence();
  }

  void TAManager::reset(bool reset_range_types) {
    peak_mem_ = 0;
    cur_mem_ = 0;
    if (reset_range_types) {
      rangeTypes_.clear();
      blockRangeTypes_.clear();
    }
    rTAstat_.clear();
    cTAstat_.clear();
    rTAs_.clear();
    cTAs_.clear();
    TA::get_default_world().gop.fence();
  }

  std::pair<double, char> memSize(size_t mem);

  std::ostream& operator<<(std::ostream& out, const TAManager& manager) {
    out << "--- TiledArray manager status ---" << std::endl;
    out << std::endl << "  - Types of ranges:" << std::endl;
    for (const auto &cr : manager.rangeTypes_) {
      out << "    " << cr.first << " : " << cr.second << std::endl;
    }

    if (not manager.rTAstat_.empty()) {
      out << std::endl << "  - Real valued TA objects:" << std::endl;
      out << "    " << "----------------------------------------------" << std::endl;
      out << "    " << " Dim.    size   Total   Avail.  Alloc.  Peak" << std::endl;
      out << "    " << "------  ------  ------  ------  ------  ------" << std::endl;
      for (const auto &sTA : manager.rTAstat_) {
        out << "    " << std::fixed << std::setw(6) << std::left << sTA.first;
        std::pair<double, char> mem_postfix = memSize(sTA.second.mem_each_);
        out << "  " << std::fixed << std::right << std::setw(5) << std::setprecision(1)
            << mem_postfix.first << mem_postfix.second;
        out << "  " << std::setw(4) << sTA.second.cur_total;
        out << "    " << std::setw(4) << manager.rTAs_.at(sTA.first).size();
        out << "    " << std::setw(4) << sTA.second.hist_total;
        out << "    " << std::setw(4) << sTA.second.peak << std::endl;
      }
      out << "    " << "----------------------------------------------" << std::endl;
    }

    if (not manager.cTAstat_.empty()) {
      out << std::endl << "  - Complex valued TA objects:" << std::endl;
      out << "    " << "----------------------------------------------" << std::endl;
      out << "    " << " Dim.    size   Total   Avail.  Alloc.  Peak" << std::endl;
      out << "    " << "------  ------  ------  ------  ------  ------" << std::endl;
      for (const auto &sTA : manager.cTAstat_) {
        out << "    " << std::fixed << std::setw(6) << std::left << sTA.first;
        std::pair<double, char> mem_postfix = memSize(sTA.second.mem_each_);
        out << "  " << std::fixed << std::right << std::setw(5) << std::setprecision(1)
        << mem_postfix.first << mem_postfix.second;
        out << "  " << std::setw(4) << sTA.second.cur_total;
        out << "    " << std::setw(4) << manager.cTAs_.at(sTA.first).size();
        out << "    " << std::setw(4) << sTA.second.hist_total;
        out << "    " << std::setw(4) << sTA.second.peak << std::endl;
      }
      out << "    " << "----------------------------------------------" << std::endl;
    }

    std::pair<double, char> mem_postfix = memSize(manager.cur_mem_);
    out << "  - TA objects total memory:" << std::endl
        << "    * Current allocated by manager "
        << std::fixed << std::right << std::setw(5) << std::setprecision(1)
        << mem_postfix.first << mem_postfix.second << "B" << std::endl;

    mem_postfix = memSize(manager.peak_mem_);
    out << "    * Manager high-water mark      "
        << std::fixed << std::right << std::setw(5) << std::setprecision(1)
        << mem_postfix.first << mem_postfix.second << "B" << std::endl;

    size_t ta_high_water_mark = TA::hostEnv::instance()->host_allocator().getHighWatermark();
    mem_postfix = memSize(ta_high_water_mark);
    out << "    * TiledArray high-water mark" << std::endl
        << "                      |- this rank "
        << std::fixed << std::right << std::setw(5) << std::setprecision(1)
        << mem_postfix.first << mem_postfix.second << "B" << std::endl;

    TA::get_default_world().gop.template reduce(&ta_high_water_mark, 1, std::plus<size_t>());
    mem_postfix = memSize(ta_high_water_mark);
    out << "                      |- overall   "
        << std::fixed << std::right << std::setw(5) << std::setprecision(1)
        << mem_postfix.first << mem_postfix.second << "B" << std::endl;

    return out;
  }

}; // namespace ChronusQ
