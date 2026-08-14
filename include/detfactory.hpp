/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2018 Li Research Group (University of Washington)
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

#define DETFACTORY_HPP

#include <chronusq_sys.hpp>
#include <detfactory/determinants.hpp>
#include <detfactory/excitationlist.hpp>

namespace ChronusQ {
  
/*
 *  Data Structure for categorical excitations
 * 
 *  Built as
 *  braCat <- exList1 <- ( auxCat <- exList2 ... <- ) ketCat
 *
 *  only bra and ket categorical ids are saved
 *  
 */ 
template <size_t NExL, size_t NEx>
struct DetsCatExcitation {
  
  const std::string term = ""; 
  std::array<std::shared_ptr<const NewExcitationList>, NExL> exLists = {};
  std::pair<size_t, size_t> categoricalIndices = {};
  std::array<size_t, NEx * 2> exSpaces = {};
  double symmetryFactor = 1.0;

}; // struct DetsCatExcitation


using DetsCat1eExcitation = DetsCatExcitation<1, 1>;
using DetsCat2eRIExcitation = DetsCatExcitation<2, 2>;

class DeterminantFactory {

 protected:  
  MPI_Comm comm_;

#ifdef CQ_ENABLE_MPI
  // DetFactoryTaskScheduler taskScheduler_;
#endif  

  // actual working active spaces
  std::vector<ActiveSpaceParameters> activeSpaces_;
  std::shared_ptr<CategoricalSpace> braCategoricalSpace_;
  std::shared_ptr<CategoricalSpace> ketCategoricalSpace_;
  
  // keep a copy for fast reference
  size_t nTotalCorrE_; 
  std::vector<size_t> nOrbitals_; 
  std::vector<size_t> orbOffs_; 
  
  // Excitation List Storage, only one copy is stored
  // Key for full list: nEx-ExSpaceInvolvedInKetDetCat
  // such as 1e-(ne,no) or 1e-(ne1,no1)-(ne2, no2) optional(-R)
  std::unordered_map<std::string, std::shared_ptr<NewExcitationList>> exLists_;
  
  // Contraction Tasks
  // Key will be the h1e(t,u) or g2e(t,u,w,v), here indices are spaces
  std::unordered_set<std::string> oneEExTerms_;
  std::unordered_set<std::string> twoEExTerms_;
  std::vector<DetsCat1eExcitation> oneEExcitations_; 
  std::vector<DetsCat2eRIExcitation> twoEExcitations_; 
  
  // helper functions for generating Exlists
  std::shared_ptr<NewExcitationList> constructFullCD1eExList(
    const DeterminantGroup &tGroup, size_t tSpace, const DeterminantGroup &uGroup, size_t uSpace);
  
  void constructOneEExcitation(const std::shared_ptr<const DeterminantCategory> &, size_t,
                               size_t, size_t, size_t, std::string);
  
  std::string constructTwoEExcitationWithExRI(const std::shared_ptr<const DeterminantCategory> &, size_t,
    size_t, size_t, size_t, size_t, const std::shared_ptr<const DeterminantCategory>&, size_t);

public:
  
  // Disable default constructor, 
  // and use default copy and move constructors
  DeterminantFactory() = delete;
  DeterminantFactory(const DeterminantFactory& other) = default;
  DeterminantFactory(DeterminantFactory&& other) = default;
  DeterminantFactory(MPI_Comm comm,
                     size_t nTotalCorrE,
                     const std::vector<ActiveSpaceParameters>& actS):
      comm_(comm), activeSpaces_(actS), nTotalCorrE_(nTotalCorrE) {
    for(const auto& s : activeSpaces_) {
      nOrbitals_.push_back(s.nOrbitals);
    }
    orbOffs_ = {0ul};
    for(auto i = 0ul; i < nOrbitals_.size() - 1; ++i) {
      orbOffs_.push_back(orbOffs_.back() + nOrbitals_[i]);   
    }
  } // constructor

  ~DeterminantFactory() { };
  
  // setters
  void setBraCategoricalSpace(std::shared_ptr<CategoricalSpace> categoricalSpace) { 
    braCategoricalSpace_ = categoricalSpace; 
  }
  void setKetCategoricalSpace(std::shared_ptr<CategoricalSpace> categoricalSpace) { 
    ketCategoricalSpace_ = categoricalSpace; 
  }
  void useSameBraDetsAsKetDets() { setBraCategoricalSpace(ketCategoricalSpace_); }
  
  // getters
  MPI_Comm MPIComm() const { return comm_; }
  std::shared_ptr<const CategoricalSpace> braCategoricalSpace() const { return braCategoricalSpace_; }
  std::shared_ptr<const CategoricalSpace> ketCategoricalSpace() const { return ketCategoricalSpace_; }
  const std::vector<ActiveSpaceParameters>& activeSpaces() const { return activeSpaces_; }
  const std::unordered_set<std::string>& oneEExTerms() const { return oneEExTerms_; }
  const std::unordered_set<std::string>& twoEExTerms() const { return twoEExTerms_; }
  const std::vector<DetsCat1eExcitation>& oneEExcitations() const { return oneEExcitations_; }
  const std::vector<DetsCat2eRIExcitation>& twoEExcitations() const { return twoEExcitations_; }
  const std::vector<size_t>& nOrbitalsInEachSpace() const { return nOrbitals_; }
  const std::vector<size_t>& orbOffsInEachSpace() const { return orbOffs_; }
  size_t nTotalCorrElectrons() const { return nTotalCorrE_; }  
  
  // ensure that it's using the same active space
  std::shared_ptr<CategoricalSpace> buildEmptyDetsSpace() const {
    return std::make_shared<CategoricalSpace>(activeSpaces_);
  }

  // establish connections between bra and ket categories
  // enforce no symmetry for 1e part and permutational symmetry for 2e part
  void generateComputingGraph(bool usingExRI = true);

  // estimate memory requirement
  void estimateMemoryRequirement() {
    double totalMem = 0.0;
    for (auto & l: exLists_) {
      totalMem += l.second->storageSize();
    }
    std::cout << "Total Memory Requirement for Storing Excitation Lists: " << std::fixed<<totalMem/1e9 << " GB" << std::endl;
    if (totalMem > CQMemManager::get().max_avail_allocatable<char>(1,totalMem)) CErr();
  }

  // virtual function to do precomputations
  void computeExcitationList() { 
    for (auto & l: exLists_) {
      // std::cout << "Compute Excitation List " << l.first << " (" << l.second->storageSize()/1e9 << "GB )"<< std::endl;
      l.second->computeExcitationList(); 
    }
  }

  void output(std::ostream & out, const std::string & s = "") const; 

}; // class DetFactory

} // namespace ChronusQ

// add implementations
#include <detfactory/impl.hpp>
