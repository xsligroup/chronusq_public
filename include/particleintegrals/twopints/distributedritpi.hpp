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

#include <particleintegrals/twopints.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cxxapi/output.hpp>
#include <particleintegrals/twopints/eri3j.hpp>
namespace ChronusQ {


template <typename MatsT, typename IntsT>
class DistributedRITPIContraction : public RITPIContraction<MatsT, IntsT> {

  template <typename MatsU, typename IntsU>
  friend class DistributedRITPIContraction;

  public:
    // Constructors
    DistributedRITPIContraction() = delete;

    DistributedRITPIContraction(std::shared_ptr<TwoPInts<IntsT>> tpi):
      RITPIContraction<MatsT,IntsT>(tpi) {}

    template <typename MatsU>
    DistributedRITPIContraction(const DistributedRITPIContraction<MatsU, IntsT>& other, int dummy = 0):
      RITPIContraction<MatsT,IntsT>(other.ints_) {}
    template <typename MatsU>
    DistributedRITPIContraction(DistributedRITPIContraction<MatsU, IntsT>&& other, int dummy = 0):
      RITPIContraction<MatsT,IntsT>(std::move(other.ints_)) {}

    DistributedRITPIContraction(const DistributedRITPIContraction& other):
      RITPIContraction<MatsT,IntsT>(other, 0) {}
    DistributedRITPIContraction(DistributedRITPIContraction&& other):
      RITPIContraction<MatsT,IntsT>(std::move(other), 0 ) {}

    // Computation interfaces
    virtual void JContract(MPI_Comm comm, TwoBodyContraction<MatsT>& C) const override;
    void JContractSplitNB(MPI_Comm comm, TwoBodyContraction<MatsT>& C, const std::shared_ptr<DistributedERI3J<IntsT>>& eri3j) const;
    void JContractSplitNBRI(MPI_Comm comm, TwoBodyContraction<MatsT>& C, const std::shared_ptr<DistributedERI3J<IntsT>>& eri3j) const;
    
    virtual void KContract(MPI_Comm comm, TwoBodyContraction<MatsT>& C) const override;
    void KContractSplitNB(MPI_Comm comm, TwoBodyContraction<MatsT>& C, const std::shared_ptr<DistributedERI3J<IntsT>>& eri3j) const;
    void KContractSplitNBRI(MPI_Comm comm, TwoBodyContraction<MatsT>& C, const std::shared_ptr<DistributedERI3J<IntsT>>& eri3j) const;
    
    void KCoefContract(MPI_Comm, size_t nO, MatsT *X, MatsT *AX) const;

    virtual ~DistributedRITPIContraction() {}

  protected:
    // Helper functions for Distributed K contractions
    void printOwnedKBlocks(const std::map<std::pair<int, int>, MatsT*>& blocks,
                        const DistributedERI3J<IntsT>& eri3j, int NB, int size) const;

    void placeKBlocksNTT(const std::map<std::pair<int, int>, MatsT*>& blocks,
                        MatsT* AX_local_NTT, const DistributedERI3J<IntsT>& eri3j, int NB, int size) const;

    void insertBlockIntoFullK(MatsT* AX, int i, int j, const MatsT* blockData, 
                        const DistributedERI3J<IntsT>& eri3j, int NB, int size) const;
    
    void KContractSplitNBRI_real_impl(MPI_Comm comm, const double *Xr, double *AXr, double *Ktemp, const DistributedERI3J<double>& eri3j) const;
  };



  template <typename MatsT, typename IntsT>
  class DistributedAsymmRITPIContraction : public RITPIContraction<MatsT, IntsT> {

    template <typename MatsU, typename IntsU>
    friend class DistributedAsymmRITPIContraction;

    public:
      DistributedAsymmRITPIContraction() = delete;

      DistributedAsymmRITPIContraction(std::shared_ptr<TwoPInts<IntsT>> tpi):
        RITPIContraction<MatsT,IntsT>(tpi) {}

      template <typename MatsU>
      DistributedAsymmRITPIContraction(const DistributedAsymmRITPIContraction<MatsU, IntsT>& other, int dummy = 0):
        RITPIContraction<MatsT,IntsT>(other.ints_) {}
      template <typename MatsU>
      DistributedAsymmRITPIContraction(DistributedAsymmRITPIContraction<MatsU, IntsT>&& other, int dummy = 0):
        RITPIContraction<MatsT,IntsT>(std::move(other.ints_)) {}

      DistributedAsymmRITPIContraction(const DistributedAsymmRITPIContraction& other):
        RITPIContraction<MatsT,IntsT>(other, 0) {}
      DistributedAsymmRITPIContraction(DistributedAsymmRITPIContraction&& other):
        RITPIContraction<MatsT,IntsT>(std::move(other), 0) {}

      // Computation interfaces
      virtual void JContract(MPI_Comm comm, TwoBodyContraction<MatsT>& C) const override;
      void JContractSplitNB(MPI_Comm comm, TwoBodyContraction<MatsT>& C, const std::shared_ptr<DistributedERI3J<IntsT>>& eri3j1, const std::shared_ptr<DistributedERI3J<IntsT>>& eri3j2) const;
      void JContractSplitNBRI(MPI_Comm comm, TwoBodyContraction<MatsT>& C, const std::shared_ptr<DistributedERI3J<IntsT>>& eri3j1, const std::shared_ptr<DistributedERI3J<IntsT>>& eri3j2) const;
      virtual void KContract(MPI_Comm comm, TwoBodyContraction<MatsT>& C) const override{
        CErr("K Contraction for (ee|pp) integrals are not valid");}

      virtual ~DistributedAsymmRITPIContraction() {}

  };
}; // namespace ChronusQ
