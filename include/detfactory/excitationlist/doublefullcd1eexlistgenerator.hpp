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

#ifndef DETFACTORY_EXCITATIONLIST
#error This file may only be included from detfactory/excitationlist.hpp
#endif

#include <tensor/tensorlooper.hpp>

namespace ChronusQ {
  
/*
 * \brief the DoubleFullCD1eExListGenerator class for
 * 
 *   <Kpqrs|E_pq|Jpqrs> <Jpqrs|E_rs|Lpqrs>
 *
 *   provide loop for Jpqrs, so will use <Jpq|E_qp|Kpq> list instead
 */ 
class DoubleFullCD1eExListGenerator {

 protected:
  
  const FullCD1eExList& qpExList_;
  const FullCD1eExList& rsExList_;
   
  // this two have same dimensions
  std::vector<size_t> JDims_;
  std::unordered_set<size_t> exSpaces_;
  std::array<size_t, 4> pqrsSpaces_;

  size_t pPos_;
  size_t qPos_;
  size_t rPos_;
  size_t sPos_;
  
  // used for internal loop
  size_t JpExOff_; 
  size_t JqExOff_; 
  size_t JrExOff_; 
  size_t JsExOff_;

 public:
  //Constructors
  DoubleFullCD1eExListGenerator() = delete;
  DoubleFullCD1eExListGenerator(
    const FullCD1eExList& qpExList,
    const FullCD1eExList& rsExList,
    const std::array<size_t, 4>& pqrsSpaces);
    
  DoubleFullCD1eExListGenerator(const DoubleFullCD1eExListGenerator & other) = default; 
  DoubleFullCD1eExListGenerator(DoubleFullCD1eExListGenerator && other)      = default;
  ~DoubleFullCD1eExListGenerator() { };  
  
  // getters
  size_t qpNNZ() const { return qpExList_.nNonZeroExcitations(); }
  size_t rsNNZ() const { return rsExList_.nNonZeroExcitations(); }
  size_t totalDimension() const {
    size_t JDim = 1ul;
    for (const auto& n : JDims_) JDim *= n;
    return JDim;
  }    
  const std::vector<size_t>& dimensions() const {  return JDims_; }
  const std::unordered_set<size_t>& excitationSpaces() const {  return exSpaces_; }
  const std::array<size_t, 4>& pqrsSpaces() const {  return pqrsSpaces_; }
  
  template <class Visitor>
  void visitAllExcitations(const std::vector<size_t>& LExOff, 
      const std::vector<size_t>& KExOff, Visitor visitor) const {
    visitExcitations(0ul, totalDimension(), LExOff, KExOff, std::forward<Visitor>(visitor));
  };

#ifdef CQ_ENABLE_SPARSE
  template <class Visitor>
  void visitAllSparseExcitations(const std::vector<size_t>& LExOff,
      const std::vector<size_t>& KExOff, Visitor visitor) const {
    visitSparseExcitations(0ul, totalDimension(), LExOff, KExOff, std::forward<Visitor>(visitor));
  };
#endif

  // Major interface for building contraction loops
  template <class Visitor>
  void visitExcitations(size_t JBegin, size_t JEnd,
                        const std::vector<size_t>& ketExOff,
                        const std::vector<size_t>& braExOff, Visitor visitor) const;

#ifdef CQ_ENABLE_SPARSE
  template <class Visitor>
  void visitSparseExcitations(size_t JBegin, size_t JEnd,
                        const std::vector<size_t>& ketExOff,
                        const std::vector<size_t>& braExOff, Visitor visitor) const;
#endif

}; // class DoubleFullCD1eExListGenerator 

inline DoubleFullCD1eExListGenerator::DoubleFullCD1eExListGenerator(
    const FullCD1eExList& qpExList,
    const FullCD1eExList& rsExList,
    const std::array<size_t, 4>& pqrsSpaces):
  qpExList_(qpExList), rsExList_(rsExList), pqrsSpaces_(pqrsSpaces) {

  // create an ordered map
  std::map<size_t, size_t> exSpaces_JDims_map;
  const auto& JqpGroups = qpExList.braCategory()->detGroups();
  const auto& JrsGroups = rsExList.braCategory()->detGroups();
  if (JqpGroups.size() > 1 ) {
    exSpaces_JDims_map.try_emplace(pqrsSpaces[0], JqpGroups[1].nDeterminants());
  }
  exSpaces_JDims_map.try_emplace(pqrsSpaces[1], JqpGroups[0].nDeterminants());
  exSpaces_JDims_map.try_emplace(pqrsSpaces[2], JrsGroups[0].nDeterminants());
  if (JrsGroups.size() > 1 ) {
    exSpaces_JDims_map.try_emplace(pqrsSpaces[3], JrsGroups[1].nDeterminants());
  }

  // copy over to exSpaces and JDims
  size_t cur = 0ul;

#ifdef DEBUG_DoubleOneEExList
  std::vector<size_t> tmp;
#endif

  for (const auto& it : exSpaces_JDims_map) {
    exSpaces_.emplace(it.first);
    JDims_.push_back(it.second);
    
#ifdef DEBUG_DoubleOneEExList
    tmp.push_back(it.first);
#endif

    // create map from total J to qpExList and rsExList
    if (it.first == pqrsSpaces_[0]) pPos_ = cur;    
    if (it.first == pqrsSpaces_[1]) qPos_ = cur;    
    if (it.first == pqrsSpaces_[2]) rPos_ = cur;    
    if (it.first == pqrsSpaces_[3]) sPos_ = cur;
    ++cur;
  }
  
  JpExOff_ = (pPos_ == qPos_) ? 0ul: JDims_[qPos_];
  JqExOff_ = 1ul; 
  JrExOff_ = 1ul; 
  JsExOff_ = (rPos_ == sPos_) ? 0ul: JDims_[rPos_];
  
#ifdef DEBUG_DoubleOneEExList
  std::cout << " exSpaces_ = "; 
  for (auto & i: exSpaces_) std::cout << i << " ";
  std::cout << std::endl;
  std::cout << " JDims_ = "; 
  for (auto & i: JDims_) std::cout << i << " ";
  std::cout << std::endl;
  std::cout << " p q r s ExSpaces= " << tmp[pPos_] 
            << " " << tmp[qPos_] 
            << " " << tmp[rPos_] 
            << " " << tmp[sPos_] << std::endl;
#endif      
  
  return;
} // DoubleOneEExListNZOperator constructor
  
template <class Visitor>
void DoubleFullCD1eExListGenerator::visitExcitations(
    size_t JBegin, size_t JEnd,
    const std::vector<size_t>& ketExOff,
    const std::vector<size_t>& braExOff,
    Visitor visitor) const {
  
  // sanity check on sizes
  // std::cout << "ketExOff.size() = " << ketExOff.size() << std::endl;
  // std::cout << "braExOff.size() = " << braExOff.size() << std::endl;
  // std::cout << "JDims.size() = "  << JDims_.size() << std::endl;
  // std::cout << "braExOff = ";
  // for (const auto& off : braExOff) std::cout << off << " ";
  // std::cout << std::endl;
  // std::cout << "ketExOff = ";
  // for (const auto& off : ketExOff) std::cout << off << " ";
  // std::cout << std::endl;
  
  assert(ketExOff.size() == JDims_.size());
  assert(braExOff.size() == JDims_.size());

  // TensorLooper returns the relative address of p,q,r,s in their own space
  auto TL = constructTensorLooper(JDims_);
  const auto& TLIndices = *TL;
  const auto& Jp = TLIndices[pPos_];
  const auto& Jq = TLIndices[qPos_];
  const auto& Jr = TLIndices[rPos_];
  const auto& Js = TLIndices[sPos_];
  
  // K, L offsets regarding to K, L addresses in a category
  const size_t LrExOff = ketExOff[rPos_];
  const size_t LsExOff = (rPos_ == sPos_) ? 0ul: ketExOff[sPos_];
  const size_t LpExOff = (pPos_ == rPos_ or pPos_ == sPos_) ? 0ul: ketExOff[pPos_];
  const size_t LqExOff = (qPos_ == pPos_ or qPos_ == rPos_ or qPos_ == sPos_) ? 0ul: ketExOff[qPos_];

  const size_t KpExOff = (pPos_ == qPos_) ? 0ul: braExOff[pPos_];
  const size_t KqExOff = braExOff[qPos_];
  const size_t KrExOff = (rPos_ == pPos_ or rPos_ == qPos_) ? 0ul: braExOff[rPos_];
  const size_t KsExOff = (sPos_ == rPos_ or sPos_ == pPos_ or sPos_ == qPos_) ? 0ul: braExOff[sPos_];
  
  // std::cout << "braExOff in p, q, r, s space: " << KpExOff << ", " << KqExOff << ", " << KrExOff << ", " << KsExOff << std::endl;
  // std::cout << "ketExOff in p, q, r, s space: " << LpExOff << ", " << LqExOff << ", " << LrExOff << ", " << LsExOff << std::endl;
  
  auto qpExGen = qpExList_.generator({KqExOff, KpExOff}, {0ul, 0ul});
  auto rsExGen = rsExList_.generator({LrExOff, LsExOff}, {0ul, 0ul});
  
  //const auto catJqp = qpExList_.braCategory();
  //const auto catJrs = rsExList_.braCategory();
  
  // main outer loop,  handles the excitation part
  // TODO: OPENMP parallelization here
  const auto& JIndex = TL->index(); 
  const auto& rsExs = rsExGen->excitations();
  const auto& qpExs = qpExGen->excitations();
  for (TL->setIndex(JBegin); JIndex < JEnd; TL->increment()) {   
    
    // computing the working ExLists
    const auto Jqp = Jp * JpExOff_ + Jq * JqExOff_;
    const auto Jrs = Jr * JrExOff_ + Js * JsExOff_;
    const auto Lqp = Jp * LpExOff + Jq * LqExOff;  
    const auto Krs = Jr * KrExOff + Js * KsExOff;  
    
    // std::cout << " JEx = " << TL->index() << std::endl;
    // std::cout << " Jqp = " << Jqp << ", " 
    //           << catJqp->addressToString(Jqp) 
    //           << std::endl;
    // std::cout << " Jrs = " << Jrs << ", " 
    //           << catJrs->addressToString(Jrs)
    //           << std::endl;


    // std::cout << "Build Contractions for Jrs = " << Jrs << ", Jqp = " << Jqp 
    //           << ", Lqp = " << Lqp << ", Krs = " << Krs << std::endl;
    rsExGen->updateExcitations(Jrs, Lqp); 
    qpExGen->updateExcitations(Jqp, Krs);
    
    //std::cout << "Contraction Build finished" << std::endl;

    visitor(JIndex, qpExs,  rsExs);
  } //  main loop

} // DoubleFullCD1eExListGenerator::visitExcitations

#ifdef CQ_ENABLE_SPARSE
template <class Visitor>
void DoubleFullCD1eExListGenerator::visitSparseExcitations(
    size_t JBegin, size_t JEnd,
    const std::vector<size_t>& ketExOff,
    const std::vector<size_t>& braExOff,
    Visitor visitor) const {

  assert(ketExOff.size() == JDims_.size());
  assert(braExOff.size() == JDims_.size());

  // TensorLooper returns the relative address of p,q,r,s in their own space
  auto TL = constructTensorLooper(JDims_);
  const auto& TLIndices = *TL;
  const auto& Jp = TLIndices[pPos_];
  const auto& Jq = TLIndices[qPos_];
  const auto& Jr = TLIndices[rPos_];
  const auto& Js = TLIndices[sPos_];

  // K, L offsets regarding to K, L addresses in a category
  const size_t LrExOff = ketExOff[rPos_];
  const size_t LsExOff = (rPos_ == sPos_) ? 0ul: ketExOff[sPos_];
  const size_t LpExOff = (pPos_ == rPos_ or pPos_ == sPos_) ? 0ul: ketExOff[pPos_];
  const size_t LqExOff = (qPos_ == pPos_ or qPos_ == rPos_ or qPos_ == sPos_) ? 0ul: ketExOff[qPos_];

  const size_t KpExOff = (pPos_ == qPos_) ? 0ul: braExOff[pPos_];
  const size_t KqExOff = braExOff[qPos_];
  const size_t KrExOff = (rPos_ == pPos_ or rPos_ == qPos_) ? 0ul: braExOff[rPos_];
  const size_t KsExOff = (sPos_ == rPos_ or sPos_ == pPos_ or sPos_ == qPos_) ? 0ul: braExOff[sPos_];

  auto qpExGen = qpExList_.generator({KqExOff, KpExOff}, {0ul, 0ul});
  auto rsExGen = rsExList_.generator({LrExOff, LsExOff}, {0ul, 0ul});


  // main outer loop, handles the excitation part with dynamic load balancing
  const auto& JIndex = TL->index();
  const auto& rsExs = rsExGen->excitations();
  const auto& qpExs = qpExGen->excitations();
    
  for (TL->setIndex(JBegin); JIndex < JEnd; TL->increment()) {
    // computing the working ExLists
    const auto Jqp = Jp * JpExOff_ + Jq * JqExOff_;
    const auto Jrs = Jr * JrExOff_ + Js * JsExOff_;
    const auto Lqp = Jp * LpExOff + Jq * LqExOff;
    const auto Krs = Jr * KrExOff + Js * KsExOff;

    // std::cout << " JEx = " << TL->index() << std::endl;
    // std::cout << " Jqp = " << Jqp << ", "
    //           << catJqp->addressToString(Jqp)
    //           << std::endl;
    // std::cout << " Jrs = " << Jrs << ", "
    //           << catJrs->addressToString(Jrs)
    //           << std::endl;


    // std::cout << "Build Contractions for Jrs = " << Jrs << ", Jqp = " << Jqp
    //           << ", Lqp = " << Lqp << ", Krs = " << Krs << std::endl;
    rsExGen->updateExcitations(Jrs, Lqp);
    qpExGen->updateExcitations(Jqp, Krs);

    //std::cout << "Contraction Build finished" << std::endl;

    visitor(JIndex, qpExs,  rsExs);
  } //  main loop

} // DoubleFullCD1eExListGenerator::visitSparseExcitations
#endif

} // namespace ChronusQ
