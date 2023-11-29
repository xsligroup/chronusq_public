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

#include <cxxapi/options.hpp>

namespace ChronusQ {

  /**
   * \brief CQIntermediates class serves as a base class for all the intermediates
   */
  class CQIntermediates {
  public:
    virtual void saveToBin(const std::string &filename, const std::string &prefix) = 0;
    virtual void readFromBin(const std::string &filename, const std::string &prefix) = 0;
  };

  struct CQJob {
    JobType jobType; // Job type
    std::string storagePrefix; // Prefix for storage
    std::map<std::string, std::shared_ptr<CQIntermediates>> prerequisites;// Prerequisites
    std::map<std::string, std::shared_ptr<CQIntermediates>> computes;// Computes
    std::map<std::string, std::shared_ptr<CQIntermediates>> saves;// Saves

    CQJob() = delete;
    /**
     * \brief Constructor for CQJob
     */
    CQJob(JobType jobType) : jobType(jobType) {}

  };

}; // namespace ChronusQ