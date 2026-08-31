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
#ifdef CQ_HAS_TA
#include <itersolver/solvervectors.hpp>
#include <coupledcluster/TAManager.hpp>
#include <iostream>

namespace ChronusQ {

  template <typename MatsT>
  class EOMCCBase;

  enum class MBTensorSymmetry { ANTI_SYMMETRIC, RCCSD };

  /// One single tensor of coupled cluster, such as T1, T2, L1, L2, R1, R2.
  template <typename MatsT>
  class MBTensor {
    using TArray = TA::TArray<MatsT>;
    TArray V_;
    std::string vir_index_range_; // eg: v, vv
    std::string occ_index_range_; // eg: c, cVV, co
    std::string shape_; // eg: vvoo
    std::string ta_labels_ = ""; // eg: a,b,i,j
    std::string name_; // eg: OneBody, TwoBody etc
    int rank_;
    MBTensorSymmetry symmetry_; // This controls the symmetry of the tensor
    double coefficient_; // coefficient for many body contribution, eg, 0.25 for vvoo type tensor
    std::vector<std::vector<int>> index_of_vo_groups_; //space index grouped by space. Eg. vvoo results in {{0,1},{2,3}}, vooo results in {{0},{1,2,3}}, lhcc results in {{0},{1},{2,3}}
    std::vector<size_t> dim_; // number of elements of each rank, lhcc results in {nL, nH, nC, nC}
    std::vector<size_t> dims_by_type_; // number of elements of each type. Eg vvoo results in {nV*(nV-1)/2, nO*(nO-1)/2}, lhcc results in {nL, nH, nC*(nC-1)/2}
    size_t outOfBound_ = 1;
    public:
    MBTensor(std::string vir_index_range, std::string occ_index_range, std::string name,
             bool memset = false, MBTensorSymmetry symmetry = MBTensorSymmetry::ANTI_SYMMETRIC)
          : name_(name), symmetry_(symmetry) {
        vir_index_range_ = vir_index_range.substr(0, vir_index_range.length()); 
        occ_index_range_ = occ_index_range.substr(0, occ_index_range.length());
        shape_ = vir_index_range_ + occ_index_range_;
        rank_ = shape_.length();
        
        //TiledArray labels
        char vir_start_label  = 'a';
        char occ_start_label  = 'i';
		std::string ta_labels = "";
        for (int i = 0; i < vir_index_range_.length(); i++) {
          ta_labels.append(1, vir_start_label + i);
        }
        for (int i = 0; i < occ_index_range_.length(); i++) {
          ta_labels.append(1, occ_start_label + i);
        }
	    for (size_t i = 0; i < ta_labels.length(); ++i) {
	        ta_labels_ += ta_labels[i];
	        if (i != ta_labels.length() - 1) {
	            ta_labels_ += ',';
	        }
	    }

        // TA memory allocation
        if (not V_.is_initialized())
          V_ = TAManager::get().malloc<MatsT>(shape_);

        // TA memory set to 0
        if (memset) {
          V_(ta_labels_) = 0.0 * V_(ta_labels_);
        }

        //coefficients
        coefficient_ = 1.0;
        std::map<char, int> locations;
        for (int i = 0; i < vir_index_range_.length(); i++) {
          locations[vir_index_range_[i]]++;
        }
        for (auto space : locations) {
          coefficient_ *= 1.0 / (double)factorial(space.second);
        }
		locations.clear();
        for (int i = 0; i < occ_index_range_.length(); i++) {
          locations[occ_index_range_[i]]++;
        }
        for (auto space : locations) {
          coefficient_ *= 1.0 / (double)factorial(space.second);
        }
	  	for (int i = 0; i < shape_.length(); ) {
	  	  std::vector<int> group;
	  	  for (int j = i; j < shape_.size(); j++) {
	  	    if (shape_[j] == shape_[i]){
	  	      group.push_back(j);
	  	    } else {
	  	      break;
	  	    }
	  	  }
	  	  index_of_vo_groups_.push_back(group);
	  	  i += group.size();
    	}
        for (auto group: index_of_vo_groups_){
          size_t this_dim = 1;
          auto iter = group.begin();
          int count = 0;
          while(iter != group.end()){
            size_t dim_this_rank = TAManager::get().getRange(shape_[*iter]).extent();
            this_dim *= dim_this_rank - count;
            dim_.push_back(dim_this_rank);
            count++;
            this_dim /= count;
            iter++;
          }
          dims_by_type_.push_back(this_dim);
          outOfBound_ *= this_dim;
        }
        outOfBound_++;
      if (rank() == 4 and symmetry_ == MBTensorSymmetry::RCCSD) {
        if (shape_ != "vvoo") {
          CErr("MBTensor::setSymmetricalElem for RCCSD type tensor only implemented for vvoo shape.");
        }
        size_t nV = dim_[0], nO = dim_[2];
        size_t nOV = nV * nO;
        outOfBound_ = nOV * (nOV + 1) / 2;
        coefficient_ = 1.0;
      }
      if (rank() == 2 and symmetry_ == MBTensorSymmetry::RCCSD) {
        coefficient_ = 2.0;
      }
    }

    MBTensor() = delete;

    /// construct tensor with the shape and value of another
    MBTensor(const MBTensor<MatsT>& other);
    
    /// move constructor
    MBTensor(MBTensor<MatsT>&& other);

    /// copy the values of another tensor 
    MBTensor<MatsT>& operator=(const MBTensor<MatsT>& other);

    ~MBTensor(){
        if (V_) {
          TAManager::get().free(shape_, std::move(V_));
        }
    }

    /// get coefficient to calculate contributing many body effect. Eg: 0.25 for vvoo type tensor
    double coefficient() const {return coefficient_;}

    /// get rank of tensor
    size_t rank() const {return shape_.length();}

    /// get symmetry of the tensor
    MBTensorSymmetry symmetry() const {return symmetry_;}

    /// get shape of tensor, eg vvoo
    const std::string shape() const {return shape_;}

    /// get shape of virtual space indices, eg vv for tensor shape vvoo
    const std::string vir_index_range() const {return vir_index_range_;}

    /// get shape of occupied space indices, eg oo for tensor shape vvoo
    const std::string occ_index_range() const {return occ_index_range_;}

    /// get name of tensor, name is constructed when tensor is constructed
    const std::string name() const {return name_;}

    /// get address of the TA object
    TArray & data(){return V_;}
    const TArray & data() const {return V_;}

    /// get space of the TA object
    size_t size() const ;

    /// get dimension of the tensor, eg. for vvoo type tensor, return {nV, nV, nO, nO}
    const std::vector<size_t>& dim() const {return dim_;}

    /// take conjugate of the tensor
    void conjugate();

    /// calculate the dot product of the tensor with other, if conjA is true, do inner product
    MatsT dot(const MBTensor<MatsT> &other, bool conjA) const;
    MatsT dot_rank4_RCCSD(const MBTensor<MatsT> &other, bool conjA) const;

    /// add alpha times X into itself 
    void axpy(MatsT alpha, const MBTensor<MatsT> &X);

    /// scale all tensor values by a factor
    void scale(MatsT factor) ;

    /// average symmetrical elements in the tensor object 
    void enforceSymmetry();

    void print(std::ostream& out) const{
      out << name_ << std::endl << V_ << std::endl;
    }
    private:
    /// average symmetrical elements in the tensor object with RCCSD symmetry, eg, A_ijab = A_jiba
    void symmetrize_rank4_RCCSD();

    /// average symmetrical elements in the tensor object with 2 fold symmetry, eg, A_ijxx=-A_jixx
    void symmetrize_two_indices(int ii, int jj);
      
    /// average symmetrical elements in the tensor object with 3 fold symmetry, where 6 elements are symmetrical
    ///eg. A_ijkxxx = -A_jikxx = A_jkixxx
    void symmetrize_three_indices(int ii, int jj, int kk);
      
    /// average symmetrical elements in the tensor object with 4 fold symmetry, where 24 elements are symmetrical
    /// eg. A_ijklxxx=-A_jiklxx etc
    void symmetrize_four_indices(int ii, int jj, int kk, int ll);

    public:
    /// given the tensor indices and the value, set all symmetrical and antisymmetrical elements of the tensor
    void setSymmetricalElem(const std::vector<size_t> indices, MatsT elem);
    /// given composed indices, find out the tensor indices, assume ij is arranged in order
    /// in order of   01, 02, 12, 03, 13, 23, 04 ...
    /// as opposed to 01, 02, 03, 04, ..., 0N, 12, 13, ...
    std::vector<size_t> findTensorIndices(size_t idx) const;
    bool isInBound(size_t idx) const { return idx < outOfBound_; }
      private:
    /// given i and j (etc) that belong to the same space and i < j, calculate the composed index of ij. (i=ijkl[0], j=ijkl[1]) 
    size_t toCompoundIdxOfSameSpace(std::vector<size_t>& ijkl) const {
      switch(ijkl.size()) {
        case 1: 
          return ijkl[0];
        case 2: 
          return ijkl[1]*(ijkl[1]-1)/2 + ijkl[0];
        case 3: 
          return ijkl[2]*(ijkl[2]-1)*(ijkl[2]-2)/6 + ijkl[1]*(ijkl[1]-1)/2 + ijkl[0];
        case 4: 
          return ijkl[3]*(ijkl[3]-1)*(ijkl[3]-2)*(ijkl[3]-3)/24 + ijkl[2]*(ijkl[2]-1)*(ijkl[2]-2)/6 + ijkl[1]*(ijkl[1]-1)/2 + ijkl[0];
        default:
          CErr("MBTensor::toCompoundIdxOfSameSpace Tensors with more then 4 indices from the same space are not implemented.");
      }
      return 0;
    }
    public:
    double sign(std::vector<size_t> pqrs) const{
      if (rank() == 4 and symmetry_ == MBTensorSymmetry::RCCSD) {
        return 1.0;
      }
      if (pqrs.size() != shape_.size())
       CErr("MBTensor::sign Tensor indices must correspond to the correct rank.");

      size_t swapCount = 0;
      size_t n = pqrs.size();

      //bubble sort
      for (auto i = 0; i < n - 1; ++i) {
          for (auto j = i + 1; j < n; ++j) {
              if (shape_[i] == shape_[j] && pqrs[i] == pqrs[j]) return 0.0;
              if (shape_[i] == shape_[j] && pqrs[i] > pqrs[j]) {
                  std::swap(pqrs[i], pqrs[j]);
                  ++swapCount;
              }
          }
      }
      return (swapCount % 2) ? -1.0 : 1.0; 
    }
    size_t toCompoundIdx(size_t i, bool i_less_j_only) const {
                return i;
    }
    size_t toCompoundIdx(std::vector<size_t> pqrs, bool i_less_j_only) const {
      if (pqrs.size() != shape_.size())
       CErr("MBTensor::toCompoundIdx Tensor indices must correspond to the correct rank.");
      size_t n = pqrs.size();

      if (n == 4 and symmetry_ == MBTensorSymmetry::RCCSD) {
        if (shape_ != "vvoo") {
          CErr("MBTensor::toCompoundIdx for RCCSD type tensor only implemented for vvoo shape.");
        }
        if (pqrs[0] > pqrs[1] or (pqrs[0] == pqrs[1] and pqrs[2] > pqrs[3])) {
          if (i_less_j_only) return outOfBound_;
          std::swap(pqrs[0], pqrs[1]);
          std::swap(pqrs[2], pqrs[3]);
        }
        // size_t ik = pqrs[0] * dim_[2] + pqrs[2];
        size_t nO = dim_[2];
        size_t jl = pqrs[1] * nO + pqrs[3];
        return jl*(jl+1)/2 + pqrs[0] * nO + pqrs[2]; // jl*(jl+1)/2 + ik
      }

      //bubble sort
      for (auto i = 0; i < n - 1; ++i) {
          for (auto j = i + 1; j < n; ++j) {
              if (shape_[i] == shape_[j] && pqrs[i] == pqrs[j]) return outOfBound_;
              if (shape_[i] == shape_[j] && pqrs[i] > pqrs[j]) {
                  if (i_less_j_only) return outOfBound_;
                  std::swap(pqrs[i], pqrs[j]);
              }
          }
      }
      size_t totalCompoundIdx = 0;
      for (auto i = 0; i < index_of_vo_groups_.size(); ++i) {
          std::vector<size_t> pqrsSameSpace(pqrs.begin() + index_of_vo_groups_[i][0], pqrs.begin() + index_of_vo_groups_[i][0] + index_of_vo_groups_[i].size());

          size_t currentSpaceCompoundIdx = toCompoundIdxOfSameSpace(pqrsSameSpace);
          size_t accumulating_dim = 1;
          for (auto j = 0; j < i; ++j) {
            accumulating_dim *= dims_by_type_[j];
          }
          totalCompoundIdx += currentSpaceCompoundIdx * accumulating_dim;
      }
      return totalCompoundIdx;
    }

    private:
    int countSwaps(const std::vector<std::pair<size_t, char>>& vector1, const std::vector<std::pair<size_t, char>>& vector2);

  };

  enum class CC_TENSOR_TYPE { EE, DIP };

  /// CCVector, contains various ranks of tensors for this CC calculation. 
  /// A typical example of a CCVector is T, where it contains T1 and T2 etc. Also are ground state Lambda, one R vector, etc.
  template <typename MatsT>
  class MBExpansion {
  protected:
    using TArray = TA::TArray<MatsT>;

    //char vLabel_ = ' ';
    //char oLabel_ = ' ';
    MatsT V0_ = 0.0;
    std::vector<std::string> tensor_info;
    std::vector<size_t> tensor_offsets;
    //bool contain_active_space;

    //CC_TENSOR_TYPE type_;
  public:
    std::vector<MBTensor<MatsT>> V_;

    /// construct MBExpansion with instructions of each tensor. The variable tensor_info should be
    /// virtural_space_range_tensor_1, occupied_space_range_tensor_1, name_tensor_1
    /// virtural_space_range_tensor_2, occupied_space_range_tensor_2, name_tensor_2
    MBExpansion(std::vector<std::string> tensor_info, bool memset = false, MBTensorSymmetry symmetry = MBTensorSymmetry::ANTI_SYMMETRIC);
    /// copy constructor
    MBExpansion(const MBExpansion<MatsT>& other);
    MBExpansion(MBExpansion<MatsT>&&);
    ~MBExpansion();

    size_t length(bool includeZeroBody = false) const {
      size_t size = 0;
      for (int i = 0; i < V_.size(); i++) 
        size += V_[i].size();
      size += includeZeroBody ? 1 : 0;
      return size;
    }

    TArray & get_tensor(std::string search_name) {
      for (int i = 0; i < V_.size(); i++) {
          if (V_[i].name() == search_name) {
              return V_[i].data();
          }
      }
      CErr("MBExpansion:: cannot find tensor of the name " + search_name);
      throw std::runtime_error("MBExpansion:: cannot find tensor of the name " + search_name);
    }  
    const TArray & get_tensor(std::string search_name) const {
      for (int i = 0; i < V_.size(); i++) {
          if (V_[i].name() == search_name) {
              return V_[i].data();
          }
      }
      CErr("MBExpansion:: cannot find tensor of the name " + search_name);
      throw std::runtime_error("MBExpansion:: cannot find tensor of the name " + search_name);
    }  


    MBExpansion<MatsT>& operator=(const MBExpansion<MatsT>&);

    void swap(MBExpansion<MatsT> &other);

    void initialize(std::vector<std::string> tensor_info, bool memset, MBTensorSymmetry symmetry = MBTensorSymmetry::ANTI_SYMMETRIC);
 
    MatsT dot(const MBExpansion<MatsT> &other, bool conjA = true) const;

    void axpy(MatsT alpha, const MBExpansion<MatsT> &X);

    double norm() const;

    double absmax() const;

    void scale(MatsT factor);

    void conjugate();

    void normalize();

    //void enforceTwoBodySymmetry();
    void enforceSymmetry();

    void projectOut(const MBExpansion<MatsT> &other, bool normalized = false);

    MatsT& zeroBody() { return V0_; }
    const MatsT& zeroBody() const { return V0_; }

    void setElem(size_t idx, MatsT elem);
    //void setElemEE(size_t idx, MatsT elem);
    //void setElemDIP(size_t idx, MatsT elem);
    void setZeroBodyElem(MatsT elem);
    //void setOneBodyElemEE(MBTensor<MatsT> &V, size_t a, size_t i, MatsT elem);
    //void setOneBodyElemDIP(MBTensor<MatsT> &V, size_t a, size_t i, MatsT elem);
    //void setTwoBodyElemEE(MBTensor<MatsT> &V, size_t a, size_t b, size_t i, size_t j, MatsT elem);
    //void setTwoBodyElemDIP(MBTensor<MatsT> &V, size_t a, size_t b, size_t i, size_t j, MatsT elem);
    //void setThreeBodyElemEE(MBTensor<MatsT> &V, size_t a, size_t b, size_t c, size_t i, size_t j, size_t k, MatsT elem);

    void print(std::ostream& out, std::string str) const{

      out << str << "[zero body]: " << V0_ << std::endl;
      for (auto V:V_){
          out << str << std::endl;
          V.print(out);
      }
    }

    void toRaw(MatsT *raw, bool includeZeroBody = false) const;
  private:
    void toRankOneElements(const MBTensor<MatsT> &V, MatsT *raw) const;
    void toRankTwoElements(const MBTensor<MatsT> &V, MatsT *raw) const;
    void toRankThreeElements(const MBTensor<MatsT> &V, MatsT *raw) const;
    void toRankFourElements(const MBTensor<MatsT> &V, MatsT *raw) const;
    void toRankFiveElements(const MBTensor<MatsT> &V, MatsT *raw) const;
    void toRankSixElements(const MBTensor<MatsT> &V, MatsT *raw) const;

  public:
    void fromRaw(const MatsT *raw, bool hasZeroBody = false);
    void partlyFromRaw(const MatsT *raw, size_t nTensor, bool hasZeroBody = false);
    std::vector<size_t> setOneBodyElem(std::vector<size_t> pq, MatsT val);
    std::vector<size_t> setTwoBodyElem(std::vector<size_t> pqrs, MatsT val);
  private:
    void fromRankOneElements(MBTensor<MatsT> &V, const MatsT *raw);
    void fromRankTwoElements(MBTensor<MatsT> &V, const MatsT *raw);
    void fromRankThreeElements(MBTensor<MatsT> &V, const MatsT *raw);
    void fromRankFourElements(MBTensor<MatsT> &V, const MatsT *raw);
    void fromRankFiveElements(MBTensor<MatsT> &V, const MatsT *raw);
    void fromRankSixElements(MBTensor<MatsT> &V, const MatsT *raw);
  };

  template <typename MatsT>
  class MBExpansionSet : public SolverVectors<MatsT> {
  protected:
    using TArray = TA::TArray<MatsT>;

    //char vLabel_ = 'v';
    //char oLabel_ = 'o';
    
    std::vector<MBExpansion<MatsT>> vecs_;

    void initialize(std::vector<std::string> &tensor_info, size_t nVec, MBTensorSymmetry symmetry);
   
    SafeFile & savFile_; 
  public:
    MBExpansionSet(std::vector<std::string> tensor_info, size_t nVec, SafeFile & savFile,
                   MBTensorSymmetry symmetry = MBTensorSymmetry::ANTI_SYMMETRIC): savFile_(savFile){
      initialize(tensor_info, nVec, symmetry);
    }

    virtual size_t length() const override {
      if (vecs_.size() > 0)
        return vecs_[0].length();
      else
        return 0;
    }
    size_t length(bool includeZeroBody) const {
      return length() + (includeZeroBody ? 1 : 0);
    }
    size_t length(const EOMCCBase<MatsT> &eom, bool includeZeroBody) const {
      return eom.getHbarDim() + (includeZeroBody ? 1 : 0);
    }

    virtual size_t size() const override { return vecs_.size(); }

    // Get element
    virtual MatsT get(size_t i, size_t j) const override {
      CErr("Get element in MBExpansionSet object is invalid.");
      return 0.0;
    }
    // Set element
    virtual void set(size_t i, size_t j, MatsT value) override {
      vecs_[j].setElem(i, value);
    }

    // Get MBExpansion
    MBExpansion<MatsT>& get(size_t i) {
      this->sizeCheck(i, "in MBExpansionSet<MatsT>::get");
      return vecs_[i];
    }
    const MBExpansion<MatsT>& get(size_t i) const {
      this->sizeCheck(i, "in MBExpansionSet<MatsT>::get");
      return vecs_[i];
    }
    // Set MBExpansion
    virtual void set(size_t i, const MBExpansion<MatsT>& other) {
      vecs_[i] = other;
    }

    using SolverVectors<MatsT>::clear;
    void clear(size_t shift, size_t nVec) override {
      if (nVec == 0) return;
      this->sizeCheck(shift + nVec, "in MBExpansionSet<MatsT>::clear");
      for (size_t i = shift; i < vecs_.size(); i++)
        vecs_[i].scale(0.0);
    }

    using SolverVectors<MatsT>::print;
    void print(std::ostream& out, std::string str, size_t shift, size_t nVec) const override {
      if (nVec == 0) {
        out << std::endl << str + ": " << std::endl;
        return;
      }
      this->sizeCheck(shift + nVec, "in MBExpansionSet<MatsT>::print");
      for (size_t i = 0; i < nVec; i++) {
        get(i + shift).print(out, str + "(" + std::to_string(i + shift) + ")");
      }
    }

    RawVectors<MatsT> toRaw(MPI_Comm c,
                            bool includeZeroBody = false, size_t shift = 0,
                            size_t nVec = std::numeric_limits<size_t>::max()) const;

    RawVectors<MatsT> toRaw(MPI_Comm c, const EOMCCBase<MatsT> &eom,
                            bool includeZeroBody = false, size_t shift = 0,
                            size_t nVec = std::numeric_limits<size_t>::max()) const;

    void fromRaw(MPI_Comm c, const RawVectors<MatsT> &raw, bool hasZeroBody = false,
                 size_t shiftThis = 0, size_t shiftRaw = 0, size_t nVec = std::numeric_limits<size_t>::max());

    void fromRaw(MPI_Comm c, const RawVectors<MatsT> &raw,
                 const EOMCCBase<MatsT> &eom, bool hasZeroBody = false,
                 size_t shiftThis = 0, size_t shiftRaw = 0, size_t nVec = std::numeric_limits<size_t>::max());

    void somFromRaw(MPI_Comm c, const RawVectors<MatsT> &raw,
                 const EOMCCBase<MatsT> &eom, bool hasZeroBody = false,
                 size_t shiftThis = 0, size_t shiftRaw = 0, size_t nTensor = 0, size_t nVec = std::numeric_limits<size_t>::max());

    virtual void multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                 MatsT alpha, MatsT const *B, int64_t ldb,
                                 MatsT beta, SolverVectors<MatsT> &C, size_t shiftC) const override;

    virtual void dot_product(size_t shiftA, const SolverVectors<MatsT> &B, size_t shiftB,
                             int64_t m, int64_t n, MatsT *C, int64_t ldc, bool conjA = true) const override;

    virtual void set_data(size_t shiftA, size_t nVec, const SolverVectors<MatsT> &B, size_t shiftB, bool moveable = false) override;

    virtual void swap_data(size_t shiftA, size_t nVec, SolverVectors<MatsT> &B, size_t shiftB) override;

    using SolverVectors<MatsT>::scale;
    virtual void scale(MatsT scalar, size_t shift, size_t nVec) override;

    using SolverVectors<MatsT>::conjugate;
    virtual void conjugate(size_t shift, size_t nVec) override;

    virtual void axpy(size_t shiftY, size_t nVec, MatsT alpha, const SolverVectors<MatsT> &X, size_t shiftX) override;

//    virtual size_t GramSchmidt(size_t shift, size_t Mold, size_t Mnew,
//                               size_t NRe = 0, double eps = 1e-12) override;

    virtual void trsm(size_t shift, int64_t n, MatsT alpha, MatsT const *A, int64_t lda) override;

    virtual int QR(size_t shift, size_t nVec, MatsT *R = nullptr, int LDR = 0) override;

    using SolverVectors<MatsT>::norm2F;
    virtual double norm2F(size_t shift, size_t nVec) const override;

    using SolverVectors<MatsT>::maxNormElement;
    virtual double maxNormElement(size_t shift, size_t nVec) const override;

    virtual void writeToBinaryFile(const std::string& saveEntryName);

    virtual ~MBExpansionSet() {}
    
  };



  template <typename MatsT>
  class MBExpansionSetDebug : public SolverVectors<MatsT> {

  protected:
    using TArray = TA::TArray<MatsT>;
    MBExpansionSet<MatsT> eomccSet_;
    RawVectors<MatsT> rawSet_;

  public:
    MBExpansionSetDebug(std::vector<std::string> tensor_info, size_t nVec, SafeFile & savFile,
                        MPI_Comm c) :
                          eomccSet_(tensor_info, nVec, savFile),
                          rawSet_(c, eomccSet_.length(), nVec){}

    double compareDebug(size_t shift = 0, size_t nVec = std::numeric_limits<size_t>::max());

    MBExpansionSet<MatsT>& getEOMCCSet() {
      return eomccSet_;
    }

    RawVectors<MatsT>& getRawSet() {
      return rawSet_;
    }

    virtual size_t length() const override { return eomccSet_.length(); }
    virtual size_t size() const override { return eomccSet_.size(); }

    // Get element
    virtual MatsT get(size_t i, size_t j) const override {
      return eomccSet_.get(i, j);
    }
    // Set element
    virtual void set(size_t i, size_t j, MatsT value) override {
      eomccSet_.set(i, j, value);
      rawSet_.set(i, j, value);
    }

    using SolverVectors<MatsT>::clear;
    void clear(size_t shift, size_t nVec) override {
      eomccSet_.clear(shift, nVec);
      rawSet_.clear(shift, nVec);
    }

    using SolverVectors<MatsT>::print;
    void print(std::ostream& out, std::string str, size_t shift, size_t nVec) const override {
      eomccSet_.print(out, str, shift, nVec);
      rawSet_.print(out, str, shift, nVec);
    }

    virtual void multiply_matrix(size_t shiftA, blas::Op transB, int64_t n, int64_t k,
                                 MatsT alpha, MatsT const *B, int64_t ldb,
                                 MatsT beta, SolverVectors<MatsT> &C, size_t shiftC) const override;

    virtual void dot_product(size_t shiftA, const SolverVectors<MatsT> &B, size_t shiftB,
                             int64_t m, int64_t n, MatsT *C, int64_t ldc, bool conjA = true) const override;

    virtual void set_data(size_t shiftA, size_t nVec, const SolverVectors<MatsT> &B, size_t shiftB, bool moveable = false) override;

    virtual void swap_data(size_t shiftA, size_t nVec, SolverVectors<MatsT> &B, size_t shiftB) override;

    using SolverVectors<MatsT>::scale;
    virtual void scale(MatsT scalar, size_t shift, size_t nVec) override;

    using SolverVectors<MatsT>::conjugate;
    virtual void conjugate(size_t shift, size_t nVec) override;

    virtual void axpy(size_t shiftY, size_t nVec, MatsT alpha, const SolverVectors<MatsT> &X, size_t shiftX) override;

    virtual size_t GramSchmidt(size_t shift, size_t Mold, size_t Mnew,
                               size_t NRe = 0, double eps = 1e-12) override;

    virtual void trsm(size_t shift, int64_t n, MatsT alpha, MatsT const *A, int64_t lda) override;

    virtual int QR(size_t shift, size_t nVec, MatsT *R = nullptr, int LDR = 0) override;

    using SolverVectors<MatsT>::norm2F;
    virtual double norm2F(size_t shift, size_t nVec) const override;

    using SolverVectors<MatsT>::maxNormElement;
    virtual double maxNormElement(size_t shift, size_t nVec) const override;
    
    virtual void writeToBinaryFile(const std::string& saveEntryName) {
                    CErr("MBExpansionSetDebug: writeToBinaryFile NYI");
    }

    virtual ~MBExpansionSetDebug() {}
  };

}; // namespace ChronusQ
#endif
