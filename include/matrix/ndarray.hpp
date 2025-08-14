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
#include <memmanager.hpp>
#include <cerr.hpp>
#include <util/matout.hpp>
#include <cqlinalg.hpp>

namespace ChronusQ {

  template <typename MatsT>
  bool hasNaN(MatsT *ptr, size_t N) {
    for (size_t i = 0; i < N; i++) {
      if (std::isnan(std::real(ptr[i]))) return true;
      if (std::isnan(std::imag(ptr[i]))) return true;
    }
    return false;
  }

  namespace cqmatrix {

    /**
     * @brief NDArray class
     *          A class for a general N-dimensional array
     *          Data is stored in column-major
     */
    template <typename MatsT>
    class NDArray {

      template <typename MatsU>
      friend class NDArray;

    protected:
      std::vector<size_t> dims_;
      MatsT *ptr_ = nullptr;     ///< Raw ndarray storage

    private:
      std::vector<size_t> strides_;

      void computeStrides() {
        strides_.resize(dims_.size());
        strides_[0] = 1;
        for (size_t i = 1; i < dims_.size(); ++i) {
          strides_[i] = strides_[i - 1] * dims_[i - 1];
        }
      }

      // Helper function to recursively print the NDArray
      void outputRecursive(std::ostream &out, size_t dim, size_t offset, size_t depth) const {
        if (dim == dims_.size() - 1) {
          out << std::string(depth, ' ') << "[";
          for (size_t i = 0; i < dims_[dim]; ++i) {
            if (i > 0) {
              out << ", ";
            }
            if (i > 0 && i % 5 == 0) {
              out << "\n" << std::string(depth + 1, ' ');
            }
            out << std::setw(10) << ptr_[offset + i];
          }
          out << "]";
        } else {
          out << std::string(depth, ' ') << "[";
          for (size_t i = 0; i < dims_[dim]; ++i) {
            if (i > 0) {
              out << ",\n";
            }
            outputRecursive(out, dim + 1, offset + i * strides_[dim], depth + 1);
          }
          out << "]";
        }
      }

    public:
      size_t nElements() const {
        return std::accumulate(dims_.begin(), dims_.end(), 1, std::multiplies<size_t>());
      }

      const std::vector<size_t>& dimensions() const { return dims_; }

      template <typename MatU>
      bool isSameDimension(const NDArray<MatU>& other) const {
        return dims_ == other.dims_;
      }

      bool isMatrix() const {
        return dims_.size() == 2;
      }

      // Constructor
      NDArray() = delete;
      /**
       * @brief Construct a new NDArray object
       *
       * @param dims Dimensions of the NDArray
       */
      NDArray(const std::vector<size_t> &dims):
          dims_(dims) {
        malloc();
      }

      NDArray( const NDArray &other ):
          NDArray(other.dims_) {
        std::copy_n(other.ptr_, nElements(), ptr_);
      }
      template <typename MatsU>
      NDArray( const NDArray<MatsU> &other, int = 0 ):
          NDArray(other.dims_) {
        if (std::is_same<MatsU, dcomplex>::value
            and std::is_same<MatsT, double>::value)
          CErr("Cannot create a Real NDArray from a Complex one.");
        std::copy_n(other.ptr_, nElements(), ptr_);
      }
      NDArray( NDArray &&other ):
          dims_(std::move(other.dims_)), ptr_(other.ptr_), strides_(std::move(other.strides_)) {
        other.ptr_ = nullptr;
      }

      // Column-major element access operator
      MatsT& operator[](const std::vector<size_t>& indices) {
        if (indices.size() != dims_.size()) {
          CErr("NDArray: Number of indices must match the number of dimensions.");
        }

        size_t index = 0;
        size_t multiplier = 1;

        for (size_t i = 0; i < dims_.size(); ++i) {
          if (indices[i] >= dims_[i]) {
            CErr("NDArray: Index out of bounds.");
          }
          index += indices[i] * strides_[i];
          multiplier *= dims_[i];
        }

        return ptr_[index];
      }

      const MatsT& operator[](const std::vector<size_t>& indices) const {
        if (indices.size() != dims_.size()) {
          CErr("NDArray: Number of indices must match the number of dimensions.");
        }

        size_t index = 0;
        size_t multiplier = 1;

        for (size_t i = 0; i < dims_.size(); ++i) {
          if (indices[i] >= dims_[i]) {
            CErr("NDArray: Index out of bounds.");
          }
          index += indices[i] * strides_[i];
          multiplier *= dims_[i];
        }

        return ptr_[index];
      }

      // Column-major element access operator
      template <typename... Indices>
      MatsT& operator()(Indices... indices) {
        std::vector<size_t> indices_vec = {static_cast<size_t>(indices)...};
        return operator[](indices_vec);
      }

      template <typename... Indices>
      const MatsT& operator()(Indices... indices) const {
        std::vector<size_t> indices_vec = {static_cast<size_t>(indices)...};
        return operator[](indices_vec);
      }

      NDArray& operator=( const NDArray &other ) {
        if (this != &other) {
          if (dims_ != other.dims_) {
            dims_ = other.dims_;
            malloc();
          }
          std::copy_n(other.ptr_, nElements(), ptr_);
        }
        return *this;
      }
      NDArray& operator=( NDArray &&other ) {
        if (this != &other) {
          dims_ = std::move(other.dims_);
          std::swap(ptr_, other.ptr_);
          strides_ = std::move(other.strides_);
        }
        return *this;
      }

      void resize(const std::vector<size_t> &dims) {
        size_t newN = std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<size_t>());
        if (nElements() != newN) {
          dims_ = dims;
          malloc();
        } else if (dims_ != dims) {
          dims_ = dims;
          computeStrides();
        }
      }

      // NDArray direct access
      MatsT* pointer() { return ptr_; }
      const MatsT* pointer() const { return ptr_; }

      NDArray<double> real_part() {
        NDArray<double> realMat(dims_);
        GetMatRE('N', 1, nElements(), 1., pointer(), 1, realMat.pointer(), 1);
        return realMat;
      }

      void clear() {
        std::fill_n(ptr_, nElements(), MatsT(0.));
      }

      void output(std::ostream &out, const std::string &s = "",
                  bool printFull = false) const {
        std::string matStr;
        if (s == "")
          matStr = "NDArray";
        else
          matStr = "NDArray[" + s + "]";
        matStr += " (";
        for (size_t i = 0; i < dims_.size(); ++i) {
          matStr += std::to_string(dims_[i]);
          if (i < dims_.size() - 1) {
            matStr += ",";
          }
        }
        matStr += ")";
        matStr += (printFull ? ":" : "");
        out << matStr << std::endl;
        if (printFull)
          outputRecursive(out, 0, 0, 0);
      }

      void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {

#ifdef CQ_ENABLE_MPI
        // BCast matrix to all MPI processes
        if( MPISize(comm) > 1 ) {
          std::cerr  << "  *** Scattering an ndarray ***\n";
          size_t oldN = nElements();
          size_t nDims = dims_.size();
          MPIBCast(nDims,root,comm);
          bool strideChanged = false;
          if (nDims != dims_.size()) {
            dims_.resize(nDims);
            strideChanged = true;
          }

          std::vector<size_t> dims_bcast = dims_;
          MPIBCast(dims_bcast.data(), dims_bcast.size(), root, comm);

          if (dims_bcast != dims_) {
            dims_ = dims_bcast;
            strideChanged = true;
          }

          if (oldN != nElements())
            malloc();
          else if (strideChanged)
            computeStrides();

          MPIBCast(ptr_, nElements(), root, comm);
        }
#endif

      }

      bool hasNaN() const {
        return ChronusQ::hasNaN(ptr_, nElements());
      }

      void malloc() {
#pragma omp critical
        {
          if (ptr_) CQMemManager::get().unsafe_free(ptr_);
          size_t N = nElements();
          if (N != 0) {
            try { ptr_ = CQMemManager::get().unsafe_malloc<MatsT>(N); }
            catch (...) {
              std::cout << std::fixed;
              std::cout << "Insufficient memory for the full NDArray ("
                        << (N / 1e9) * sizeof(double) << " GB)" << std::endl;
              std::cout << std::endl << CQMemManager::get() << std::endl;
              throw std::bad_alloc();
            }
          }
        }
        computeStrides();
      }

      ~NDArray() {
        if(ptr_) CQMemManager::get().free(ptr_);
      }

    }; // class NDArray

    template <typename MatsT>
    std::ostream& operator<<(std::ostream&, const NDArray<MatsT>&);

  } // namespace cqmatrix
} // namespace ChronusQ
