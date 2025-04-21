/*
*  This file is part of the Chronus Quantum (ChronusQ) software package
*
*  Copyright (C) 2014-2020 Li Research Group (University of Washington)
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
#include <coupledcluster/MBExpansion.hpp>

namespace ChronusQ{

    
  /// copy the values of another tensor 
  template <typename MatsT>
  MBTensor<MatsT>& MBTensor<MatsT>:: operator=(const MBTensor<MatsT>& other) {
    if (this != &other) {
      vir_index_range_ = other.vir_index_range_;
      occ_index_range_ = other.occ_index_range_;
      shape_ = other.shape_;
      ta_labels_ = other.ta_labels_;
      V_ = other.V_.clone();
      name_ = other.name_;
      rank_ = other.rank_;
      coefficient_ = other.coefficient_;
      index_of_vo_groups_ = other.index_of_vo_groups_;
      dims_by_type_ = other.dims_by_type_;
      dim_ = other.dim_;
      outOfBound_ = other.outOfBound_;
    }
    return *this;
  }
    
  /// copy constructor
  template <typename MatsT>
  MBTensor<MatsT>::MBTensor(const MBTensor<MatsT>& other) {
    vir_index_range_ = other.vir_index_range_;
    occ_index_range_ = other.occ_index_range_;
    shape_ = other.shape_;
    ta_labels_ = other.ta_labels_;
    name_ = other.name_;
    rank_ = other.rank_;
    coefficient_ = other.coefficient_;
    V_ = TAManager::get().malloc<MatsT>(shape_);
    V_ = other.V_.clone();
    index_of_vo_groups_ = other.index_of_vo_groups_;
    dims_by_type_ = other.dims_by_type_;
    dim_ = other.dim_;
    outOfBound_ = other.outOfBound_;
  }
  template MBTensor<dcomplex>::MBTensor(const MBTensor<dcomplex>& other);

  /// move constructor
  template <typename MatsT>
  MBTensor<MatsT>::MBTensor(MBTensor<MatsT>&& other) {
    std::swap(V_, other.V_);
    std::swap(vir_index_range_, other.vir_index_range_);
    std::swap(occ_index_range_, other.occ_index_range_);
    std::swap(shape_, other.shape_);
    std::swap(ta_labels_, other.ta_labels_);
    std::swap(name_, other.name_);
    std::swap(rank_, other.rank_);
    std::swap(coefficient_, other.coefficient_);
    std::swap(index_of_vo_groups_, other.index_of_vo_groups_);
    std::swap(dims_by_type_, other.dims_by_type_);
    std::swap(dim_, other.dim_);
    std::swap(outOfBound_, other.outOfBound_);
  }

  /// get space of the TA object
  template <typename MatsT>
  size_t MBTensor<MatsT>::size() const{
    TAManager &TAmanager = TAManager::get();
    if (shape_.length() == 0) return 0;
    std::unordered_map<char, int> letterCount;
    size_t size_vir = 1;
    for (char letter : vir_index_range_) {
        if (letterCount.find(letter) != letterCount.end()) {
            size_vir *= TAmanager.getRange(letter).extent() - letterCount[letter];
            letterCount[letter]++;
            size_vir /= letterCount[letter];
        } else {
            letterCount[letter] = 1;
            size_vir *= TAmanager.getRange(letter).extent();
        }
    }

    size_t size_occ = 1.0;
    letterCount.clear();
    for (char letter : occ_index_range_) {
        if (letterCount.find(letter) != letterCount.end()) {
            size_occ *= TAmanager.getRange(letter).extent() - letterCount[letter];
            letterCount[letter]++;
            size_occ /= letterCount[letter];
        } else {
            letterCount[letter] = 1;
            size_occ *= TAmanager.getRange(letter).extent();
        }
    }

    return size_occ * size_vir;

  }

  /// take conjugate of the tensor
  template <typename MatsT>
  void MBTensor<MatsT>::conjugate(){
    if (std::is_same<MatsT, double>::value) return;
    //if (vir_index_range_.length() != occ_index_range_.length()) 
    //  CErr("MBTensor:: tensor has wrong rank, cannot be conjugated.");
    V_(ta_labels_) = conj(V_(ta_labels_));
  }

  /// calculate the dot product of the tensor with other, if conjA is true, do inner product
  template <typename MatsT>
  MatsT MBTensor<MatsT>::dot(const MBTensor<MatsT> &other, bool conjA) const {
    if (rank_ != other.rank_) 
       CErr("MBTensor::dot product error, rank mismatch.");
    if (conjA){
      return V_(ta_labels_).inner_product(other.V_(ta_labels_)).get();
      TA::get_default_world().gop.fence();    
    } else {
      return V_(ta_labels_).dot(other.V_(ta_labels_)).get();
      TA::get_default_world().gop.fence();    
    }
    TA::get_default_world().gop.fence();
    return MatsT(0);
  }

  /// add alpha times X into itself 
  template <typename MatsT>
  void MBTensor<MatsT>::axpy(MatsT alpha, const MBTensor<MatsT> &X) {
    if (rank_ != X.rank_) 
       CErr("MBTensor::axpy error, rank mismatch.");
    V_(ta_labels_) += alpha * X.V_(ta_labels_);
  }

  template <typename MatsT>
  void MBTensor<MatsT>::scale(MatsT factor) {
    V_(ta_labels_) = factor * V_(ta_labels_);
  }

  /// average symmetrical elements in the tensor object with 2 fold symmetry, eg, A_ijxx=-A_jixx
  template <typename MatsT>
  void MBTensor<MatsT>::symmetrize_two_indices(int ii, int jj){
    int i = 2 * ii;
    int j = 2 * jj;
    std::string replace_to = ta_labels_;
    std::swap(replace_to[i], replace_to[j]);
    V_(ta_labels_) -= V_(replace_to);
    V_(ta_labels_) = 0.5 * V_(ta_labels_);
  }
    
  /// average symmetrical elements in the tensor object with 3 fold symmetry, where 6 elements are symmetrical
  ///eg. A_ijkxxx = -A_jikxx = A_jkixxx
  template <typename MatsT>
  void MBTensor<MatsT>::symmetrize_three_indices(int ii, int jj, int kk){
    int i = 2 * ii;
    int j = 2 * jj;
    int k = 2 * kk;
    std::string replace_to_2_1 = ta_labels_;
    std::string replace_to_3_1 = ta_labels_;
    std::string replace_to_3_2 = ta_labels_;
    std::swap(replace_to_2_1[i], replace_to_2_1[j]); // ijk -> jik
    std::swap(replace_to_3_1[i], replace_to_3_1[k]); // ijk -> kji
    std::swap(replace_to_3_2[j], replace_to_3_2[k]); // ijk -> ikj
    V_(ta_labels_) -= V_(replace_to_3_1) + V_(replace_to_3_2);
    V_(ta_labels_) -= V_(replace_to_2_1);
    V_(ta_labels_) = 1.0/6.0 * V_(ta_labels_);
  }
    
  /// average symmetrical elements in the tensor object with 4 fold symmetry, where 24 elements are symmetrical
  /// eg. A_ijklxxx=-A_jiklxx etc
  template <typename MatsT>
  void MBTensor<MatsT>::symmetrize_four_indices(int ii, int jj, int kk, int ll){
    int i = 2 * ii;
    int j = 2 * jj;
    int k = 2 * kk;
    int l = 2 * ll;
    std::string replace_to_2_1 = ta_labels_;
    std::string replace_to_3_1 = ta_labels_;
    std::string replace_to_3_2 = ta_labels_;
    std::string replace_to_4_1 = ta_labels_;
    std::string replace_to_4_2 = ta_labels_;
    std::string replace_to_4_3 = ta_labels_;
    std::swap(replace_to_2_1[i], replace_to_2_1[j]); // ijkl -> jikl
    std::swap(replace_to_3_1[i], replace_to_3_1[k]); // ijkl -> kjil
    std::swap(replace_to_3_2[j], replace_to_3_2[k]); // ijkl -> ikjl
    std::swap(replace_to_4_1[i], replace_to_4_1[l]); // ijkl -> lijk
    std::swap(replace_to_4_2[j], replace_to_4_2[l]); // ijkl -> ilkj
    std::swap(replace_to_4_3[k], replace_to_4_3[l]); // ijkl -> ijlk
    V_(ta_labels_) -= V_(replace_to_4_1) + V_(replace_to_4_2) + V_(replace_to_4_3);
    V_(ta_labels_) -= V_(replace_to_3_1) + V_(replace_to_3_2);
    V_(ta_labels_) -= V_(replace_to_2_1);
    V_(ta_labels_) = 1.0/24.0 * V_(ta_labels_);
  }
  /// average symmetrical elements in the tensor object 
  template <typename MatsT>
  void MBTensor<MatsT>::enforceSymmetry() {
    // symmetrize virtual space
    std::map<char, std::vector<int>> locations;
    for (int i = 0; i < vir_index_range_.length(); i++) {
      locations[vir_index_range_[i]].push_back(i);
    }
    for (auto space : locations) {
      switch (space.second.size()) {
        case 4: 
          symmetrize_four_indices(space.second[0], space.second[1], space.second[2], space.second[3]); 
          break;
        case 3: 
          symmetrize_three_indices(space.second[0], space.second[1], space.second[2]); 
          break;
        case 2: 
          symmetrize_two_indices(space.second[0], space.second[1]); 
          break;
        case 1: 
          break;
        default: 
          CErr("MBTensor:: enforceSymmetry with this rank is not yet implemented.");
      }
    }

    // symmetrize occupied space
    locations.clear();
    for (int i = 0; i < occ_index_range_.length(); i++) {
      locations[occ_index_range_[i]].push_back(i + vir_index_range_.length());
    }
    for (auto space : locations) {
      switch (space.second.size()) {
        case 4: 
          symmetrize_four_indices(space.second[0], space.second[1], space.second[2], space.second[3]); 
          break;
        case 3: 
          symmetrize_three_indices(space.second[0], space.second[1], space.second[2]); 
          break;
        case 2: 
          symmetrize_two_indices(space.second[0], space.second[1]); 
          break;
        case 1: 
          break;
        default: 
          CErr("MBTensor:: enforceSymmetry with this rank is not yet implemented.");
      }
    }
  }

  template <typename MatsT>
  std::vector<size_t> MBTensor<MatsT>::findTensorIndices(size_t idx) const{
    TAManager &TAmanager = TAManager::get();
    std::vector<size_t> pqrs(rank_, SIZE_MAX);
    std::vector<size_t> idxByType(dims_by_type_.size());
    size_t elementIndex = idx;
    for (int i = 0; i < dims_by_type_.size(); i++) {
        idxByType[i] = elementIndex % dims_by_type_[i];
        elementIndex /= dims_by_type_[i];
    }
    
    for (int i = 0; i < dims_by_type_.size() ; i++) {
        switch (index_of_vo_groups_[i].size()) {
          case 1: {
            pqrs[index_of_vo_groups_[i][0]] = idxByType[i];
            break;
          } 
          case 2: {
            size_t q = static_cast<size_t>(sqrt(2*idxByType[i] + 0.25) + 0.5);
            size_t p = idxByType[i] - q * (q-1) / 2;
            pqrs[index_of_vo_groups_[i][0]] = p;
            pqrs[index_of_vo_groups_[i][1]] = q;
            break;
          } 
          case 3: {
            size_t r = static_cast<size_t>(std::cbrt(6.0*idxByType[i]));
            while (idxByType[i] * 6 >= r * r * r - r) {
              r++;
            }
            size_t pq = idxByType[i] - r * (r-1) * (r-2) / 6;
            size_t q = static_cast<size_t>(sqrt(2*pq + 0.25) + 0.5);
            size_t p = pq - q * (q-1) / 2;
            pqrs[index_of_vo_groups_[i][0]] = p;
            pqrs[index_of_vo_groups_[i][1]] = q;
            pqrs[index_of_vo_groups_[i][2]] = r;
            
            break;
          } 
          case 4: {
            size_t s = static_cast<size_t>(std::sqrt(std::sqrt(24.0*idxByType[i])));
            while (idxByType[i] * 24 >= s * (s+1) * (s-1) * (s-2)) {
              s++;
            }
            size_t pqr = idxByType[i] - s * (s-1) * (s-2) * (s-3) / 24;

            size_t r = static_cast<size_t>(std::cbrt(6.0*pqr));
            while (pqr * 6 >= r * r * r - r) {
              r++;
            }
            size_t pq = idxByType[i] - r * (r-1) * (r-2) / 6;
            size_t q = static_cast<size_t>(sqrt(2*pq + 0.25) + 0.5);
            size_t p = pq - q * (q-1) / 2;
            pqrs[index_of_vo_groups_[i][0]] = p;
            pqrs[index_of_vo_groups_[i][1]] = q;
            pqrs[index_of_vo_groups_[i][2]] = r;
            pqrs[index_of_vo_groups_[i][3]] = s;
            break;
          } 
        }

    }
    for (auto p : pqrs) {
      if ( p == SIZE_MAX){
        CErr("MBTensor::findTensorIndices Something is wrong in calculation the tensor indicies.\n"
             "Does Your tensor have more than four (4) upper or lower indicies with the same space?\n" 
             "(eg. R^ooo_vvvvv has 5 v indices, and the function cannot handle it)  \n");
      }
    }
    return pqrs;
  }


  // Function to count the number of swaps needed to transform vector2 into vector1
  template <typename MatsT>
  int MBTensor<MatsT>::countSwaps(const std::vector<std::pair<size_t, char>>& vector1, const std::vector<std::pair<size_t, char>>& vector2) {
    int n = vector1.size();
    int swapCount = 0;

    // Create a copy of vector2
    std::vector<std::pair<size_t, char>> tempVector = vector2;

    // Iterate through each element in vector1
    for (int i = 0; i < n; ++i) {
        // Find the index of the current element in vector2
        auto it = std::find(tempVector.begin()+i, tempVector.end(), vector1[i]);

        // If the element is not found in vector2, return -1 indicating invalid transformation
        if (it == tempVector.end()) {
            return -1;
        }

        // Calculate the number of swaps needed to bring the current element to its correct position
        int index = std::distance(tempVector.begin(), it);
        if (i != index) {
            swapCount++;
            // Swap the current element with the element at index i in vector2
            std::swap(tempVector[i], tempVector[index]);
        }
    }

    return swapCount;
  }

  template <typename MatsT>
  void MBTensor<MatsT>::setSymmetricalElem(std::vector<size_t> pqrs, MatsT elem) {
    if (pqrs.size() != shape_.size())
       CErr("MBTensor::setSymmetrticalElem Tensor indices must correspond to the correct rank.");
   
	// attach space info to pqrs label 
    std::vector<std::pair<size_t, char>> indices;
    indices.reserve(pqrs.size());  
    for (size_t i = 0; i < pqrs.size(); ++i) {
        indices.emplace_back(pqrs[i], shape_[i]);  
    }

    std::vector<std::pair<size_t, char>> indices_save(indices);
    std::sort(indices.begin(), indices.end());
	//set all symmetrical elements with the correct signs
    do {
        bool match = true;
		// do not permute indices that belong to different shapes
		for (size_t i = 0; i < indices.size(); i++) {
		    if (indices[i].second != shape_[i] ) {
			  match = false;
			  break;
			}
		}
		if (match) {            
		    // find the correct sign
            MatsT val = (countSwaps(indices, indices_save) % 2) ? -elem : elem; 
			// extract the index into a vector(idx) from the vector of pairs(indices)
    		std::vector<size_t> idx(indices.size());
			std::transform(indices.begin(), indices.end(), idx.begin(),
                   [](const std::pair<size_t, char>& p) { return p.first; });
			// assign one element to the TA object
    		TA::foreach_inplace(V_, [&idx, val](TA::Tensor<MatsT> &tile) {
    		    const auto& lobound = tile.range().lobound();
    		    const auto& upbound = tile.range().upbound();

    		    bool in_bounds = true;
    		    for (size_t i = 0; i < idx.size(); ++i) {
    		        if (!(lobound[i] <= idx[i] && idx[i] < upbound[i])) {
    		            in_bounds = false;
    		            break;
    		        }
    		    }

    		    if (in_bounds) {
    		        tile[idx] = val;
    		    }
    		});
        }
    } while ( std::next_permutation(indices.begin(), indices.end()) );

    TA::get_default_world().gop.fence();
  }


  template <typename MatsT>
  MBExpansion<MatsT>::MBExpansion(std::vector<std::string> tensor_info1, bool memset):
    tensor_info(tensor_info1){
    if ( tensor_info.size() % 3 != 0)
      CErr("MBExpansion:: must construct MBExpansion with vir_space_range, occ_space_range, space_name.");
    
    //type_ = CC_TENSOR_TYPE::EE;
    //// takes advantage that DIP one body virtual space is 0. 
    //// Note we need to modify this part for DEA, but not for CVS-EE or active DIP
    //if (tensor_info[0].size() == 0) {
    //  type_ = CC_TENSOR_TYPE::DIP;
    //}

    //bool contain_active_space = false;
    //for (int i = 0; i < tensor_info.size(); i += 3) {
    //  // if all virtual ranges are all the same, like vvoo, where vv are the same
    //  if (tensor_info[i].find_first_not_of(vLabel_) != std::string::npos)
    //    contain_active_space = true;
    //  // if all occ ranges are all the same, like vvoo, unlike vvcV, where cV are different.
    //  if (tensor_info[i+1].find_first_not_of(oLabel_) != std::string::npos)
    //    contain_active_space = true;
    //}
    //if (contain_active_space) {
    //  CErr("MBExpansion:: The existing MBExpansion does not work for vectors that contain active space. Write setElem, toRaw, fromRaw, vLabel, oLabel to make it work!");
    //}

    initialize(tensor_info, memset);
    size_t offset = 0;
    for (size_t i = 0; i < V_.size(); i++) {
      tensor_offsets.push_back(offset);
      offset += V_[i].size();
    }

  }

  template <typename MatsT>
  MBExpansion<MatsT>::MBExpansion(const MBExpansion<MatsT> &other) {
    for (int i = 0; i < other.V_.size(); i++) {
      V_.emplace_back(other.V_[i].vir_index_range(), other.V_[i].occ_index_range(), other.V_[i].name());
    }
    operator=(other);
  }

  template <typename MatsT>
  MBExpansion<MatsT>::MBExpansion(MBExpansion<MatsT> &&other) {
    swap(other);
  }

  template <typename MatsT>
  MBExpansion<MatsT>::~MBExpansion() {
  }

  template <typename MatsT>
  MBExpansion<MatsT>& MBExpansion<MatsT>::operator=(const MBExpansion<MatsT>& other) {
    if (this != &other) {
      V0_ = other.V0_;
      //V_.reserve(other.V_.size());
      for ( int i = 0; i < other.V_.size(); i++) {
        V_[i] = other.V_[i];
        //V_.emplace_back(MBTensor<MatsT>(Vi));
      }
      tensor_info = other.tensor_info;
      tensor_offsets = other.tensor_offsets;
    }
    return *this;
  }

  template <typename MatsT>
  void MBExpansion<MatsT>::swap(MBExpansion<MatsT> &other) {
    std::swap(V0_, other.V0_);
    std::swap(V_, other.V_);
    std::swap(tensor_info, other.tensor_info);
    std::swap(tensor_offsets, other.tensor_offsets);
  }


  template <typename MatsT>
  void MBExpansion<MatsT>::initialize(std::vector<std::string> tensor_info, bool memset){

    V0_ = 0.0;

    for (int i = 0; i < tensor_info.size(); i += 3) {
      V_.emplace_back(tensor_info[i], tensor_info[i+1], tensor_info[i+2], memset);
    }
  }
  template <typename MatsT>
  MatsT MBExpansion<MatsT>::dot(const MBExpansion<MatsT> &other, bool conjA) const {

    MatsT dotProduct = 0.0;
    if (V_.size() != other.V_.size())
      CErr("MBExpansion::dot product error, number of tensors mismatch.");
    if (conjA) {
      dotProduct += SmartConj(V0_) * other.V0_;
    } else {
      dotProduct += V0_ * other.V0_;
    }
    for (int i = 0; i < V_.size(); i++){
      dotProduct += V_[i].coefficient() * V_[i].dot(other.V_[i], conjA);
      TA::get_default_world().gop.fence();    
    }
    TA::get_default_world().gop.fence();

    return dotProduct;
  }

  template <typename MatsT>
  void MBExpansion<MatsT>::axpy(MatsT alpha, const MBExpansion<MatsT> &X) {
    if (V_.size() != X.V_.size())
      CErr("MBExpansion::axpy product error, number of tensors mismatch.");
    V0_ += alpha * X.V0_;
    for (int i = 0; i < V_.size(); i++){
        V_[i].axpy(alpha, X.V_[i]);
    }
  }

  template <typename MatsT>
  double MBExpansion<MatsT>::norm() const {
    return std::sqrt(std::real(dot(*this)));
  }

  template <typename MatsT>
  double MBExpansion<MatsT>::absmax() const {
    std::vector<double> max;
    max.push_back(std::abs(V0_));
    for (auto V: V_) 
      max.push_back(TA::abs_max(V.data()).get());
    return *std::max_element(max.begin(), max.end());
  }

  template <typename MatsT>
  void MBExpansion<MatsT>::scale(MatsT factor) {
    V0_ *= factor;
    for (int i = 0; i < V_.size(); i++){
      V_[i].scale(factor);
    }
  }

  template <typename MatsT>
  void MBExpansion<MatsT>::conjugate() {
    if (std::is_same<MatsT, double>::value) return;
    V0_ = std::conj(V0_);
    for (int i = 0; i < V_.size(); i++){
      V_[i].conjugate();
    }
  }

  template <typename MatsT>
  void MBExpansion<MatsT>::normalize() {
    scale(1.0/norm());
  }

  template <typename MatsT>
  void MBExpansion<MatsT>::enforceSymmetry() {
    for (int i = 0; i < V_.size(); i++){
      V_[i].enforceSymmetry(); 
    }
  }


  template <typename MatsT>
  void MBExpansion<MatsT>::projectOut(const MBExpansion<MatsT> &other, bool normalized) {
    MatsT coef = dot(other);
    if (not normalized)
      coef /= other.dot(other);
      TA::get_default_world().gop.fence();    
    axpy(-coef, other);
  }

  template <typename MatsT>
  void MBExpansion<MatsT>::setElem(size_t idx, MatsT elem) {


    // find the tensor and adjust the idx to the index within this tensor
    size_t id_tensor;
    for (auto it = tensor_offsets.begin(); it != tensor_offsets.end() ; it++) {
        if (*it <= idx) {
            auto next_it = std::next(it);
            if ((next_it == tensor_offsets.end()) || (next_it != tensor_offsets.end() && *next_it > idx)) {
                idx -= *it; 
                id_tensor = it - tensor_offsets.begin();
                break;
            }
        }
    }

    std::vector<size_t> pqrs = V_[id_tensor].findTensorIndices(idx);
    V_[id_tensor].setSymmetricalElem(pqrs, elem);

  }
  template <typename MatsT>
  void MBExpansion<MatsT>::setZeroBodyElem(MatsT elem) {
    V0_ = elem;
  }

  template <typename MatsT>
  void MBExpansion<MatsT>::toRaw(MatsT *raw, bool includeZeroBody) const {

    std::fill_n(raw, length(includeZeroBody), 0.0);
    MatsT *raw1 = raw + (includeZeroBody ? 1 : 0);

    for (int i = 0; i < V_.size(); i++) {
        switch (V_[i].rank()) {
            case 1:
                toRankOneElements(V_[i], raw1 + tensor_offsets[i]);
                break;
            case 2:
                toRankTwoElements(V_[i], raw1 + tensor_offsets[i]);
                break;
            case 3:
                toRankThreeElements(V_[i], raw1 + tensor_offsets[i]);
                break;
            case 4:
                toRankFourElements(V_[i], raw1 + tensor_offsets[i]);
                break;
            case 5:
                toRankFiveElements(V_[i], raw1 + tensor_offsets[i]);
                break;
            case 6:
                toRankSixElements(V_[i], raw1 + tensor_offsets[i]);
                break;
            default:
                CErr("MBExpansion::toRaw Rank not implemented");
        }
    }
    TA::get_default_world().gop.fence();
    MatsT *raw1_copy = CQMemManager::get().malloc<MatsT>(length(false));
    std::copy_n(raw1, length(false), raw1_copy);
    std::fill_n(raw1, length(false), MatsT(0.0));
    MPIAllReduce(raw1_copy, length(false), raw1, MPI_COMM_WORLD);
    CQMemManager::get().free(raw1_copy);
    
    TA::get_default_world().gop.fence();

    if (includeZeroBody)
      raw[0] = V0_;

  }

  template <typename MatsT>
  void MBExpansion<MatsT>::toRankOneElements(const MBTensor<MatsT> &V, MatsT * raw1) const{

    TA::foreach_inplace(const_cast<TArray &>(V.data()), [&V, raw1](TA::Tensor<MatsT> &tile){
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
        size_t idx = V.toCompoundIdx(x[0], true);
        if (V.isInBound(idx))
          raw1[idx] = tile[x];
        }
    });
  }

  template <typename MatsT>
  void MBExpansion<MatsT>::toRankTwoElements(const MBTensor<MatsT> &V, MatsT * raw1) const{

    TA::foreach_inplace(const_cast<TArray &>(V.data()), [&V, raw1](TA::Tensor<MatsT> &tile){
      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
          size_t idx = V.toCompoundIdx(std::vector<size_t>{x[0], x[1]}, true);
          if (V.isInBound(idx)) {
            double signAB = V.sign(std::vector<size_t>{x[0], x[1]});
            raw1[idx] = tile[x] * signAB;
          }
        }
    });
  }

  template <typename MatsT>
  void MBExpansion<MatsT>::toRankThreeElements(const MBTensor<MatsT> &V, MatsT * raw2) const{

    TA::foreach_inplace(const_cast<TArray &>(V.data()), [&V, raw2](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2]) {
            size_t idx = V.toCompoundIdx(std::vector<size_t>{x[0], x[1], x[2]}, true);
            if (V.isInBound(idx)) {
              double signABI = V.sign(std::vector<size_t>{x[0], x[1], x[2]});
              raw2[idx] = tile[x] * signABI;
            }
          }
    });
  }

  template <typename MatsT>
  void MBExpansion<MatsT>::toRankFourElements(const MBTensor<MatsT> &V, MatsT * raw2) const{

    TA::foreach_inplace(const_cast<TArray &>(V.data()), [&V, raw2](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
              size_t idx = V.toCompoundIdx(std::vector<size_t>{x[0], x[1], x[2], x[3]}, true);
              if (V.isInBound(idx)) {
                double signABIJ = V.sign(std::vector<size_t>{x[0], x[1], x[2], x[3]});
                raw2[idx] = tile[x] * signABIJ;
              }
            }
    });
  }
  template <typename MatsT>
  void MBExpansion<MatsT>::toRankFiveElements(const MBTensor<MatsT> &V, MatsT * raw2) const{

    TA::foreach_inplace(const_cast<TArray &>(V.data()), [&V, raw2](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) 
              for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4]) {
                size_t idx = V.toCompoundIdx(std::vector<size_t>{x[0], x[1], x[2], x[3], x[4]}, true);
                if (V.isInBound(idx)) {
                  double signABIJ = V.sign(std::vector<size_t>{x[0], x[1], x[2], x[3], x[4]});
                  raw2[idx] = tile[x] * signABIJ;
                }
              }
    });
  }
  template <typename MatsT>
  void MBExpansion<MatsT>::toRankSixElements(const MBTensor<MatsT> &V, MatsT * raw2) const{

    TA::foreach_inplace(const_cast<TArray &>(V.data()), [&V, raw2](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0,0,0,0,0,0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) 
              for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4]) 
                for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]) {
                  size_t idx = V.toCompoundIdx(std::vector<size_t>{x[0], x[1], x[2], x[3], x[4], x[5]}, true);
                  if (V.isInBound(idx)) {
                    double signABIJ = V.sign(std::vector<size_t>{x[0], x[1], x[2], x[3], x[4], x[5]});
                    raw2[idx] = tile[x] * signABIJ;
                  }
                }
    });
  }

  template <typename MatsT>
  void MBExpansion<MatsT>::fromRaw(const MatsT *raw, bool hasZeroBody) {
      partlyFromRaw(raw, V_.size(), hasZeroBody);
  }

  template <typename MatsT>
  void MBExpansion<MatsT>::partlyFromRaw(const MatsT *raw, size_t nTensor, bool hasZeroBody) {
    const MatsT *raw1 = raw;
    if (hasZeroBody) {
      raw1++;
      V0_ = raw[0];
    }
    MPIBCast(const_cast<MatsT*>(raw), length(hasZeroBody), 0, MPI_COMM_WORLD);
    TA::get_default_world().gop.fence();

    for (int i = 0; i < nTensor; i++) {
        switch (V_[i].rank()) {
            case 1:
                fromRankOneElements(V_[i], raw1 + tensor_offsets[i]);
                break;
            case 2:
                fromRankTwoElements(V_[i], raw1 + tensor_offsets[i]);
                break;
            case 3:
                fromRankThreeElements(V_[i], raw1 + tensor_offsets[i]);
                break;
            case 4:
                fromRankFourElements(V_[i], raw1 + tensor_offsets[i]);
                break;
            case 5:
                fromRankFiveElements(V_[i], raw1 + tensor_offsets[i]);
                break;
            case 6:
                fromRankSixElements(V_[i], raw1 + tensor_offsets[i]);
                break;
            default:
                CErr("MBExpansion::fromRaw Rank not implemented");
        }
    }
    TA::get_default_world().gop.fence();
  }
    
  template <typename MatsT>
  void MBExpansion<MatsT>::fromRankOneElements(MBTensor<MatsT> &V, const MatsT * raw1){

    TA::foreach_inplace(V.data(), [&V, raw1](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
          size_t idx = V.toCompoundIdx(x[0], false);
          tile[x] = V.isInBound(idx) ? raw1[idx] : 0.0;
        }
    });
  }
  template <typename MatsT>
  void MBExpansion<MatsT>::fromRankTwoElements(MBTensor<MatsT> &V, const MatsT * raw1){

    TA::foreach_inplace(V.data(), [&V, raw1](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      std::size_t x[] = {0, 0};
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
          double signAB = V.sign(std::vector<size_t>{x[0],x[1]});
          if (signAB == 0.0) {
            tile[x] = 0.0;
          }
          else {
            size_t idx = V.toCompoundIdx(std::vector<size_t>{x[0], x[1]}, false);
            tile[x] = V.isInBound(idx) ? signAB * raw1[idx] : 0.0;
          }
        }
    });
  }

  template <typename MatsT>
  void MBExpansion<MatsT>::fromRankThreeElements(MBTensor<MatsT> &V, const MatsT * raw2){
    TA::foreach_inplace(V.data(), [&V, raw2](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      size_t x[] = {0,0,0};
      size_t a,b,i;
      double signABI;
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2]) {
              a = x[0];
              b = x[1];
              i = x[2];
              signABI = V.sign(std::vector<size_t>{a,b,i});
              if (signABI == 0.0)
                tile[x] = 0.0;
              else {
                size_t idx = V.toCompoundIdx(std::vector<size_t>{a, b, i}, false);
                tile[x] = V.isInBound(idx) ? signABI * raw2[idx] : 0.0;
              }
            }
    });
  }
  template <typename MatsT>
  void MBExpansion<MatsT>::fromRankFourElements(MBTensor<MatsT> &V, const MatsT * raw2){
    TA::foreach_inplace(V.data(), [&V, raw2](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      size_t x[] = {0,0,0,0};
      size_t a,b,i,j;
      double signABIJ;
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
              a = x[0];
              b = x[1];
              i = x[2];
              j = x[3];
              signABIJ = V.sign(std::vector<size_t>{a,b,i,j});
              if (signABIJ == 0.0)
                tile[x] = 0.0;
              else {
                size_t idx = V.toCompoundIdx(std::vector<size_t>{a, b, i, j}, false);
                tile[x] = V.isInBound(idx) ? signABIJ * raw2[idx] : 0.0;
              }
            }
    });
  }
  template <typename MatsT>
  void MBExpansion<MatsT>::fromRankFiveElements(MBTensor<MatsT> &V, const MatsT * raw2){
    TA::foreach_inplace(V.data(), [&V, raw2](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      size_t x[] = {0,0,0,0,0};
      size_t a,b,i,j,k;
      double signABIJ;
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
              for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4]) {
                a = x[0];
                b = x[1];
                i = x[2];
                j = x[3];
                k = x[4];
                signABIJ = V.sign(std::vector<size_t>{a,b,i,j,k});
                if (signABIJ == 0.0)
                  tile[x] = 0.0;
                else {
                  size_t idx = V.toCompoundIdx(std::vector<size_t>{a, b, i, j, k}, false);
                  tile[x] = V.isInBound(idx) ? signABIJ * raw2[idx] : 0.0;
                }
              }
    });
  }
  template <typename MatsT>
  void MBExpansion<MatsT>::fromRankSixElements(MBTensor<MatsT> &V, const MatsT * raw2){

      TA::foreach_inplace(V.data(), [&V, raw2](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      const auto& upbound = tile.range().upbound();

      size_t x[] = {0,0,0,0,0,0};
      size_t a,b,c,i,j,k;
      double signABIJ;
      for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
        for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1])
          for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
            for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3])
              for(x[4] = lobound[4]; x[4] < upbound[4]; ++x[4]) 
                for(x[5] = lobound[5]; x[5] < upbound[5]; ++x[5]) {
                  a = x[0];
                  b = x[1];
                  c = x[2];
                  i = x[3];
                  j = x[4];
                  k = x[5];
                  signABIJ = V.sign(std::vector<size_t>{a,b,c,i,j,k});
                  if (signABIJ == 0.0)
                    tile[x] = 0.0;
                  else {
                    size_t idx = V.toCompoundIdx(std::vector<size_t>{a, b, c, i, j, k}, false);
                    tile[x] = V.isInBound(idx) ? signABIJ * raw2[idx] : 0.0;
                  }
                }
    });
  }



  template <typename MatsT>
  void MBExpansionSet<MatsT>::initialize(std::vector<std::string> &tensor_info, size_t nVec) {
    vecs_.reserve(nVec);
    //std::vector<std::string> tmp;
    //tmp.push_back(std::string({vLabel_}));
    //tmp.push_back(std::string({oLabel_}));
    //tmp.push_back(std::string("OneBody"));
    //tmp.push_back(std::string({vLabel_,vLabel_}));
    //tmp.push_back(std::string({oLabel_,oLabel_}));
    //tmp.push_back(std::string("TwoBody"));
    for (size_t i = 0; i < nVec; i++)
      vecs_.emplace_back(tensor_info);
  }

  template <typename MatsT>
  void MBExpansionSet<MatsT>::multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                              MatsT alpha, MatsT const *B, int64_t ldb,
                                              MatsT beta, SolverVectors<MatsT> &C, size_t shiftC) const {

    if (transB != blas::Op::NoTrans)
      CErr("Transpose of B matrix NYI in MBExpansionSet::multiply_matrix");

    C.scale(beta, shiftC, n);
    
    tryDowncastReferenceTo<MBExpansionSet<MatsT>>(C,
                                                  [&] (auto& CRef, size_t extraShiftC) {
          shiftC += extraShiftC;
          for (size_t j = 0; j < n; j++) {
            MBExpansion<MatsT> &C_vec = CRef.get(j + shiftC);
            for (size_t i = 0; i < k; i++) {
              C_vec.axpy(alpha * B[i + j * ldb], get(i + shiftA));
            }
          }
        }
    );
  }

  template <typename MatsT>
  void MBExpansionSet<MatsT>::dot_product(size_t shiftA, const SolverVectors<MatsT> &B, size_t shiftB,
                                          int64_t m, int64_t n, MatsT *C, int64_t ldc, bool conjA) const {

    tryDowncastReferenceTo<MBExpansionSet<MatsT>>(B,
                                                  [&] (auto& BRef, size_t extraShiftB) {
          shiftB += extraShiftB;
          for (size_t i = 0; i < m; i++) {
            const MBExpansion<MatsT> &A_vec = get(i + shiftA);
            for (size_t j = 0; j < n; j++) {
              const MBExpansion<MatsT> &B_vec = BRef.get(j + shiftB);
              C[i + j * ldc] = A_vec.dot(B_vec, conjA);
            }
          }
        }
    );
  }

  template <typename MatsT>
  void MBExpansionSet<MatsT>::set_data(size_t shiftA, size_t nVec, const SolverVectors<MatsT> &B, size_t shiftB, bool moveable) {
    if (moveable) {
      swap_data(shiftA, nVec, const_cast<SolverVectors<MatsT>&>(B), shiftB);
      return;
    }

    tryDowncastReferenceTo<MBExpansionSet<MatsT>>(B,
                                                  [&] (auto& BRef, size_t extraShiftB) {
          shiftB += extraShiftB;
          for (size_t i = 0; i < nVec; i++) {
            get(i + shiftA) = BRef.get(i + shiftB);
          }
        }
    );
  }

  template <typename MatsT>
  void MBExpansionSet<MatsT>::swap_data(size_t shiftA, size_t nVec, SolverVectors<MatsT> &B, size_t shiftB) {

    tryDowncastReferenceTo<MBExpansionSet<MatsT>>(B,
                                                  [&] (auto& BRef, size_t extraShiftB) {
          shiftB += extraShiftB;
          for (size_t i = 0; i < nVec; i++) {
            get(i + shiftA).swap(BRef.get(i + shiftB));
          }
        }
    );
  }

  template <typename MatsT>
  void MBExpansionSet<MatsT>::scale(MatsT scalar, size_t shiftA, size_t nVec) {
    if (nVec == 0) return;

    this->sizeCheck(shiftA + nVec, "MBExpansionSet<MatsT>::scale");

    for (size_t i = 0; i < nVec; i++) {
      get(i + shiftA).scale(scalar);
    }
  }

  template <typename MatsT>
  void MBExpansionSet<MatsT>::conjugate(size_t shiftA, size_t nVec) {
    if (std::is_same<MatsT, double>::value) return;
    if (nVec == 0) return;

    this->sizeCheck(shiftA + nVec, "MBExpansionSet<MatsT>::conjugate");

    for (size_t i = 0; i < nVec; i++) {
      get(i + shiftA).conjugate();
    }
  }

  template <typename MatsT>
  void MBExpansionSet<MatsT>::axpy(size_t shiftY, size_t nVec, MatsT alpha, const SolverVectors<MatsT> &X, size_t shiftX) {

    tryDowncastReferenceTo<MBExpansionSet<MatsT>>(X,
                                                  [&] (auto& XRef, size_t extraShiftX) {
          shiftX += extraShiftX;
          for (size_t i = 0; i < nVec; i++) {
            get(i + shiftY).axpy(alpha, XRef.get(i + shiftX));
          }
        }
    );
  }

  template <typename MatsT>
  void MBExpansionSet<MatsT>::trsm(size_t shift, int64_t n, MatsT alpha, MatsT const *A, int64_t lda) {
    CErr("MBExpansionSet::trsm NYI.");
    abort();
  }

  template <typename MatsT>
  int MBExpansionSet<MatsT>::QR(size_t shift, size_t nVec, MatsT *R, int LDR) {
    CErr("MBExpansionSet::QR NYI.");
    abort();
  }

  template <typename MatsT>
  double MBExpansionSet<MatsT>::norm2F(size_t shift, size_t nVec) const {
    if (nVec == 0) return 0.0;

    this->sizeCheck(shift + nVec, "MBExpansionSet<MatsT>::norm2F");

    double norm = 0.0;

    for (size_t i = 0; i < nVec; i++) {
      double norm_i = get(i + shift).norm();
      norm += norm_i * norm_i;
    }

    return std::sqrt(norm);
  }

  template <typename MatsT>
  double MBExpansionSet<MatsT>::maxNormElement(size_t shift, size_t nVec) const {
    if (nVec == 0) return 0.0;

    this->sizeCheck(shift + nVec, "MBExpansionSet<MatsT>::maxNormElement");

    double absmax = 0.0;

    for (size_t i = 0; i < nVec; i++) {
      double absmax_i = get(i + shift).absmax();
      absmax = std::max(absmax, absmax_i);
    }

    return absmax;
  }

  template <typename MatsT>
  RawVectors<MatsT> MBExpansionSet<MatsT>::toRaw(
      MPI_Comm c, bool includeZeroBody, size_t shift, size_t nVec) const {
    if (nVec == 0) return RawVectors<MatsT>(c, length(includeZeroBody), nVec);

    if (nVec == std::numeric_limits<size_t>::max()) {
      this->sizeCheck(shift, "MBExpansionSet<MatsT>::toRaw");
      nVec = size() - shift;
    } else
      this->sizeCheck(shift + nVec, "MBExpansionSet<MatsT>::toRaw");


    RawVectors<MatsT> raw(c, length(includeZeroBody), nVec);
    MatsT* rawPtr = raw.getPtr();
    if (MPIRank(c) != 0)
      rawPtr = CQMemManager::get().malloc<MatsT>(nVec * length(includeZeroBody));

    for (size_t i = 0; i < nVec; i++)
      get(shift + i).toRaw(rawPtr + i * length(includeZeroBody), includeZeroBody);


    if (MPIRank(c) != 0)
      CQMemManager::get().free(rawPtr);

    return raw;

  }

  template <typename MatsT>
  RawVectors<MatsT> MBExpansionSet<MatsT>::toRaw(
      MPI_Comm c, const EOMCCBase<MatsT> &eom,
      bool includeZeroBody, size_t shift, size_t nVec) const {
    if (nVec == 0) return RawVectors<MatsT>(c, length(eom, includeZeroBody), nVec);

    if (nVec == std::numeric_limits<size_t>::max()) {
      this->sizeCheck(shift, "MBExpansionSet<MatsT>::toRaw");
      nVec = size() - shift;
    } else
      this->sizeCheck(shift + nVec, "MBExpansionSet<MatsT>::toRaw");


    RawVectors<MatsT> raw(c, length(eom, includeZeroBody), nVec);
    MatsT* rawPtr = raw.getPtr();
    if (MPIRank(c) != 0)
      rawPtr = CQMemManager::get().malloc<MatsT>(nVec * length(eom, includeZeroBody));

    for (size_t i = 0; i < nVec; i++)
      get(shift + i).toRaw(rawPtr + i * length(eom, includeZeroBody), includeZeroBody);


    if (MPIRank(c) != 0)
      CQMemManager::get().free(rawPtr);

    return raw;

  }

  template <typename MatsT>
  void MBExpansionSet<MatsT>::fromRaw(MPI_Comm c,
                                      const RawVectors<MatsT> &raw, const EOMCCBase<MatsT> &eom,
                                      bool hasZeroBody, size_t shiftThis, size_t shiftRaw, size_t nVec) {
    if (nVec == 0) return;

    if (nVec == std::numeric_limits<size_t>::max()) {
      this->sizeCheck(shiftThis, "MBExpansionSet<MatsT>::fromRaw");
      raw.sizeCheck(shiftRaw, "MBExpansionSet<MatsT>::fromRaw");
      nVec = std::min(size() - shiftThis, raw.size() - shiftRaw);
    } else {
      this->sizeCheck(shiftThis + nVec, "MBExpansionSet<MatsT>::fromRaw");
      raw.sizeCheck(shiftRaw + nVec, "MBExpansionSet<MatsT>::fromRaw");
    }


    MatsT* rawPtr = const_cast<MatsT*>(raw.getPtr(shiftRaw));
    if (MPIRank(c) != 0)
      rawPtr = CQMemManager::get().malloc<MatsT>(nVec * length(eom, hasZeroBody));

    TA::get_default_world().gop.broadcast(rawPtr, nVec * length(eom, hasZeroBody), 0);

    for (size_t i = 0; i < nVec; i++)
      get(shiftThis + i).fromRaw(rawPtr + i * length(eom, hasZeroBody), hasZeroBody);


    if (MPIRank(c) != 0)
      CQMemManager::get().free(rawPtr);

  }

  template <typename MatsT>
  void MBExpansionSet<MatsT>::fromRaw(MPI_Comm c, const RawVectors<MatsT> &raw,
                                      bool hasZeroBody, size_t shiftThis, size_t shiftRaw, size_t nVec) {
    if (nVec == 0) return;

    if (nVec == std::numeric_limits<size_t>::max()) {
      this->sizeCheck(shiftThis, "MBExpansionSet<MatsT>::fromRaw");
      raw.sizeCheck(shiftRaw, "MBExpansionSet<MatsT>::fromRaw");
      nVec = std::min(size() - shiftThis, raw.size() - shiftRaw);
    } else {
      this->sizeCheck(shiftThis + nVec, "MBExpansionSet<MatsT>::fromRaw");
      raw.sizeCheck(shiftRaw + nVec, "MBExpansionSet<MatsT>::fromRaw");
    }


    MatsT* rawPtr = const_cast<MatsT*>(raw.getPtr(shiftRaw));
    if (MPIRank(c) != 0)
      rawPtr = CQMemManager::get().malloc<MatsT>(nVec * length(hasZeroBody));

    TA::get_default_world().gop.broadcast(rawPtr, nVec * length(hasZeroBody), 0);

    for (size_t i = 0; i < nVec; i++)
      get(shiftThis + i).fromRaw(rawPtr + i * length(hasZeroBody), hasZeroBody);


    if (MPIRank(c) != 0)
      CQMemManager::get().free(rawPtr);

  }

  template <typename MatsT>
  void MBExpansionSet<MatsT>::writeToBinaryFile(const std::string& saveEntryName){
        dcomplex * r_vector;
        size_t Hbar_dim = vecs_[0].length();
        r_vector = CQMemManager::get().malloc<dcomplex>(Hbar_dim);
        for(size_t i = 0; i < vecs_.size(); i++) {
          get(i).toRaw(r_vector, false);
          TA::get_default_world().gop.fence();
          std::string numberedEntryName = saveEntryName+std::to_string(i);
          if(this->savFile_.exists()) savFile_.safeWriteData(numberedEntryName, r_vector, {Hbar_dim});
          TA::get_default_world().gop.fence();
        }
        CQMemManager::get().free(r_vector);

  }
  template <typename MatsT>
  double MBExpansionSetDebug<MatsT>::compareDebug(size_t shift, size_t nVec) {
    if (nVec == 0) return 0.0;

    if (nVec == std::numeric_limits<size_t>::max()) {
      eomccSet_.sizeCheck(shift, "MBExpansionSetDebug<MatsT>::compareDebug");
      rawSet_.sizeCheck(shift, "MBExpansionSetDebug<MatsT>::compareDebug");
      nVec = size() - shift;
    } else {
      eomccSet_.sizeCheck(shift + nVec, "MBExpansionSetDebug<MatsT>::compareDebug");
      rawSet_.sizeCheck(shift + nVec, "MBExpansionSetDebug<MatsT>::compareDebug");
    }

    std::cout<<"TA result "<<eomccSet_.norm2F(0,nVec)<< " vs Raw result " << rawSet_.norm2F(0, nVec)<<" ";
    return rawSet_.norm2F(0, nVec) - eomccSet_.norm2F(0, nVec);

  }

  template <typename MatsT>
  void MBExpansionSetDebug<MatsT>::multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                                   MatsT alpha, MatsT const *B, int64_t ldb,
                                                   MatsT beta, SolverVectors<MatsT> &C, size_t shiftC) const {

    tryDowncastReferenceTo<MBExpansionSetDebug<MatsT>>(C,
                                                       [&] (auto& C_debug, size_t extraShiftC) {
          shiftC += extraShiftC;
          eomccSet_.multiply_matrix(shiftA, transB, n, k, alpha, B, ldb, beta, C_debug.getEOMCCSet(), shiftC);

          TA::get_default_world().gop.fence();
          rawSet_.multiply_matrix(shiftA, transB, n, k, alpha, B, ldb, beta, C_debug.getRawSet(), shiftC);

          std::cout << "MBExpansionSetDebug::multiply_matrix error = "
                    << C_debug.compareDebug(shiftC, n) << std::endl;
        }
    );
  }

  template <typename MatsT>
  void MBExpansionSetDebug<MatsT>::dot_product(size_t shiftA, const SolverVectors<MatsT> &B, size_t shiftB,
                                               int64_t m, int64_t n, MatsT *C, int64_t ldc, bool conjA) const {
    if (m * n == 0) return;

    tryDowncastReferenceTo<MBExpansionSetDebug<MatsT>>(B,
                                                       [&] (auto& B_debug, size_t extraShiftB) {
          shiftB += extraShiftB;
          eomccSet_.dot_product(shiftA, B_debug.eomccSet_, shiftB, m, n, C, ldc, conjA);

          MatsT *C_ref = CQMemManager::get().malloc<MatsT>(m * n);

          TA::get_default_world().gop.fence();
          rawSet_.dot_product(shiftA, B_debug.rawSet_, shiftB, m, n, C_ref, m, conjA);

          if (ldc == m)
            blas::axpy(m * n, -1.0, C, 1, C_ref, 1);
          else
            for (size_t i = 0; i < n; i++)
              blas::axpy(m, -1.0, C + i * ldc, 1, C_ref + i * m, 1);

          std::cout << "MBExpansionSetDebug::dot_product error = "
                    << blas::nrm2(m * n, C_ref, 1) << std::endl;
        }
    );
  }

  template <typename MatsT>
  void MBExpansionSetDebug<MatsT>::set_data(size_t shiftA, size_t nVec, const SolverVectors<MatsT> &B, size_t shiftB, bool moveable) {
    if (moveable) {
      swap_data(shiftA, nVec, const_cast<SolverVectors<MatsT>&>(B), shiftB);
      return;
    }

    tryDowncastReferenceTo<MBExpansionSetDebug<MatsT>>(B,
                                                       [&] (auto& B_debug, size_t extraShiftB) {
          shiftB += extraShiftB;
          eomccSet_.set_data(shiftA, nVec, B_debug.eomccSet_, shiftB);

          TA::get_default_world().gop.fence();
          rawSet_.set_data(shiftA, nVec, B_debug.rawSet_, shiftB);

          std::cout << "MBExpansionSetDebug::set_data error = "
                    << compareDebug(shiftA, nVec) << std::endl;
        }
    );
  }

  template <typename MatsT>
  void MBExpansionSetDebug<MatsT>::swap_data(size_t shiftA, size_t nVec, SolverVectors<MatsT> &B, size_t shiftB) {

    tryDowncastReferenceTo<MBExpansionSetDebug<MatsT>>(B,
                                                       [&] (auto& B_debug, size_t extraShiftB) {
          shiftB += extraShiftB;
          
          eomccSet_.swap_data(shiftA, nVec, B_debug.eomccSet_, shiftB);

          TA::get_default_world().gop.fence();
          rawSet_.swap_data(shiftA, nVec, B_debug.rawSet_, shiftB);

          std::cout << "MBExpansionSetDebug::swap error = "
          << compareDebug(shiftA, nVec) << std::endl;
        }
    );

  }

  template <typename MatsT>
  void MBExpansionSetDebug<MatsT>::scale(MatsT scalar, size_t shiftA, size_t nVec) {

    eomccSet_.scale(scalar, shiftA, nVec);

    TA::get_default_world().gop.fence();
    rawSet_.scale(scalar, shiftA, nVec);

    std::cout << "MBExpansionSetDebug::scale error = "
              << compareDebug(shiftA, nVec) << std::endl;
  }

  template <typename MatsT>
  void MBExpansionSetDebug<MatsT>::conjugate(size_t shiftA, size_t nVec) {

    eomccSet_.conjugate(shiftA, nVec);

    TA::get_default_world().gop.fence();
    rawSet_.conjugate(shiftA, nVec);

    std::cout << "MBExpansionSetDebug::conjugate error = "
              << compareDebug(shiftA, nVec) << std::endl;
  }

  template <typename MatsT>
  void MBExpansionSetDebug<MatsT>::axpy(size_t shiftY, size_t nVec, MatsT alpha, const SolverVectors<MatsT> &X, size_t shiftX) {

    tryDowncastReferenceTo<MBExpansionSetDebug<MatsT>>(X,
                                                       [&] (auto& X_debug, size_t extraShiftX) {
          shiftX += extraShiftX;
          eomccSet_.axpy(shiftY, nVec, alpha, X_debug.eomccSet_, shiftX);

          TA::get_default_world().gop.fence();
          rawSet_.axpy(shiftY, nVec, alpha, X_debug.rawSet_, shiftX);

          std::cout << "MBExpansionSetDebug::axpy error = "
                    << compareDebug(shiftY, nVec) << std::endl;
        }
    );
  }

  template <typename MatsT>
  size_t MBExpansionSetDebug<MatsT>::GramSchmidt(size_t shift, size_t Mold, size_t Mnew,
                                                 size_t NRe, double eps) {
    double vecDiffBefore = compareDebug(shift, Mold + Mnew);

    std::cout << "MBExpansionSetDebug::GramSchmidt before error = "
              << vecDiffBefore << std::endl;

    size_t iOrtho = SolverVectors<MatsT>::GramSchmidt(shift, Mold, Mnew, NRe, eps);

    double vecDiffAfter = compareDebug(shift, Mold + Mnew);

    std::cout << "MBExpansionSetDebug::GramSchmidt after error = "
              << vecDiffAfter << std::endl;

    return iOrtho;

  }

  template <typename MatsT>
  void MBExpansionSetDebug<MatsT>::trsm(size_t shift, int64_t n, MatsT alpha, MatsT const *A, int64_t lda) {
    eomccSet_.trsm(shift, n, alpha, A, lda);

    TA::get_default_world().gop.fence();
    rawSet_.trsm(shift, n, alpha, A, lda);

    std::cout << "MBExpansionSetDebug::trsm error = "
              << compareDebug(shift, n) << std::endl;

  }

  template <typename MatsT>
  int MBExpansionSetDebug<MatsT>::QR(size_t shift, size_t nVec, MatsT *R, int LDR) {
    int iOrtho = eomccSet_.QR(shift, nVec, R, LDR);

    TA::get_default_world().gop.fence();
    int iOrthoRaw = rawSet_.QR(shift, nVec, R, LDR);

    std::cout << "MBExpansionSetDebug::QR error = "
              << compareDebug(shift, nVec) << std::endl;

    if (iOrtho != iOrthoRaw)
      std::cout << "MBExpansionSetDebug::QR iOrtho differs: eomccSet = "
                << iOrtho << ", rawSet = " << iOrthoRaw << std::endl;

    return iOrtho;
  }

  template <typename MatsT>
  double MBExpansionSetDebug<MatsT>::norm2F(size_t shift, size_t nVec) const {

    double norm = eomccSet_.norm2F(shift, nVec);

    TA::get_default_world().gop.fence();
    double normRaw = rawSet_.norm2F(shift, nVec);

    std::cout << "MBExpansionSetDebug::norm2F error = "
              << std::abs(norm - normRaw) << std::endl;

    return norm;
  }

  template <typename MatsT>
  double MBExpansionSetDebug<MatsT>::maxNormElement(size_t shift, size_t nVec) const {

    double absmax = eomccSet_.maxNormElement(shift, nVec);

    TA::get_default_world().gop.fence();
    double absmaxRaw = rawSet_.maxNormElement(shift, nVec);

    std::cout << "MBExpansionSetDebug::maxNormElement error = "
              << std::abs(absmax - absmaxRaw) << std::endl;

    return absmax;
  }

}
