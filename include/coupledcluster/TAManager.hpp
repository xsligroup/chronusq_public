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

#include <chronusq_sys.hpp>

#ifdef CQ_HAS_TA
#include <tiledarray.h>
#include <cerr.hpp>

namespace ChronusQ {

  class TAManager;

  std::ostream& operator<<(std::ostream& out, const TAManager& manager);

  struct TAStat {
    size_t cur_total = 0;
    size_t hist_total = 0;
    size_t peak = 0;
    size_t mem_each_ = 0;
  };

  /**
   * A singleton class for TiledArray object allocation statistics and recycle
   */
  class TAManager {
    friend std::ostream& operator<<(std::ostream& out, const TAManager& manager);

  private:
    size_t peak_mem_ = 0;
    size_t cur_mem_ = 0;
    std::map<char, TA::TiledRange1> rangeTypes_;
    std::map<char, std::pair<size_t, size_t>> blockRangeTypes_;
    std::map<std::string, TAStat> rTAstat_;
    std::map<std::string, TAStat> cTAstat_;
    std::map<std::string, std::vector<TA::TArray<double>>> rTAs_;
    std::map<std::string, std::vector<TA::TArray<dcomplex>>> cTAs_;

    TAManager() = default;

  public:
    TAManager(TAManager const &) = delete;
    TAManager &operator=(TAManager const &) = delete;
    static TAManager& get() {
      static TAManager manager;
      return manager;
    }

    void addRangeType(char key, const TA::TiledRange1 &tr) {
      if (rangeTypes_.count(key))
        CErr(std::string("Key ") + key + std::string(" already exist in TAManager!"));
      rangeTypes_[key] = tr;
    }

    void addBlockRangeType(char key, size_t begin, size_t end) {
      if (blockRangeTypes_.count(key))
        CErr(std::string("Key ") + key + std::string(" already exist in TAManager!"));
      blockRangeTypes_[key] = std::pair(begin, end);
    }

    const TA::TiledRange1& getRange(char key) const {
      if (not rangeTypes_.count(key))
        CErr(std::string("Key ") + key + std::string(" does not exist in TAManager!"));
      return rangeTypes_.at(key);
    }

    std::vector<std::pair<size_t, size_t>> toBlockRange(std::string ranges) const;

    TA::TiledRange toRange(std::string ranges) const;

    size_t elem_per_TA(std::string ranges) const;

    template<typename MatsT>
    TA::TArray<MatsT> malloc_fresh(std::string ranges);

    template<typename MatsT>
    TA::TArray<MatsT> malloc(std::string ranges);

    template<typename MatsT>
    void free(std::string ranges, TA::TArray<MatsT> &&ta, bool discard = false);

    void discard_cache();

    void reset(bool reset_range_types = true);

  }; // class TAManager

}; // namespace ChronusQ
#endif
