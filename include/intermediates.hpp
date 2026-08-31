/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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

#include <string>
#include <unordered_map>
#include <any>
#include <memory>


namespace ChronusQ {
  class CQIntermediates {
  private:
    // Private constructor to prevent direct instantiation
    CQIntermediates() {}

    // Data storage using an unordered_map
    std::unordered_map<std::string, std::any> dataMap;

  public:
    // Deleted copy constructor and assignment operator to prevent copies
    CQIntermediates(const CQIntermediates &) = delete;
    CQIntermediates(CQIntermediates &&) = delete;
    CQIntermediates &operator=(const CQIntermediates &) = delete;
    CQIntermediates &operator=(CQIntermediates &&) = delete;

    // Static method to access the singleton instance
    static CQIntermediates &getInstance() {
      static CQIntermediates instance;  // Guaranteed to be destroyed and instantiated correctly
      return instance;
    }

    void clear() {
      dataMap.clear();
    }

    bool hasData(const std::string &name) {
      return dataMap.find(name) != dataMap.end();
    }

    // Add data with a unique name
    template<typename T>
    void addData(const std::string &name, std::shared_ptr<T> data, bool raiseIfOverwrite = false) {
      if (raiseIfOverwrite and dataMap.find(name) != dataMap.end()) {
        CErr("Data with this name already exists.");
      }
      dataMap[name] = data;
    }

    // Get data by name
    template<typename T>
    std::shared_ptr<T> getData(const std::string &name) {
      auto it = dataMap.find(name);
      if (it != dataMap.end()) {
        try {
          return std::any_cast<std::shared_ptr<T>>(it->second);
        } catch (const std::bad_any_cast &) {
          CErr("Data of " + name + " does not match requested type.");
        }
      }
      CErr("Data with name " + name + " not found.");
      return nullptr;  // Return an empty shared_ptr on failure
    }

    // Erase data by name
    void eraseData(const std::string &name) {
      dataMap.erase(name);
    }
  };

}; // namespace ChronusQ