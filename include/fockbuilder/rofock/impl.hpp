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

#include <fockbuilder/rofock.hpp>
#include <cqlinalg.hpp>

namespace ChronusQ {

  /**   
   *  \brief Forms the Roothaan's effective Fock matrix from general fock
   *
   *  Reference: J. Phys. Chem. A, Vol. 114, No. 33, 2010;
   *             Mol. Phys. 1974, 28, 819–828.
   *
   */
  template <typename MatsT, typename IntsT>
  void ROFock<MatsT,IntsT>::rohfFock(SingleSlater<MatsT,IntsT> &ss) {

    size_t NB = ss.basisSet().nBasis;
    size_t NB2 = NB*NB;

    ROOT_ONLY(ss.comm);

    //construct focka and fockb
    std::vector<cqmatrix::Matrix<MatsT>> SCR = ss.fockMatrix->template spinGatherToBlocks<MatsT>(false);

    //construct projectors for closed, open, virtual
    //pc = dmb * S
    //po = (dma - dmb) * S
    //pv = I - dma * S
    MatsT* pc  = CQMemManager::get().malloc<MatsT>(NB*NB);
    MatsT* po  = CQMemManager::get().malloc<MatsT>(NB*NB);
    MatsT* pv  = CQMemManager::get().malloc<MatsT>(NB*NB);
    cqmatrix::Matrix<MatsT> tmp(NB);
    MatsT* tmp2 = CQMemManager::get().malloc<MatsT>(NB*NB);

    //overlap matrix
    tmp = 0.5 * (ss.onePDM->S() - ss.onePDM->Z());
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),tmp.pointer(),NB,
         cqmatrix::Matrix<MatsT>(ss.aoints_->overlap->matrix()).pointer(),NB,
         MatsT(0.),pc,NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),ss.onePDM->Z().pointer(),NB,
         cqmatrix::Matrix<MatsT>(ss.aoints_->overlap->matrix()).pointer(),NB,
         MatsT(0.),po,NB);
    tmp = 0.5 * (ss.onePDM->S() + ss.onePDM->Z());
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(-1.),tmp.pointer(),NB,
         cqmatrix::Matrix<MatsT>(ss.aoints_->overlap->matrix()).pointer(),NB,
         MatsT(0.),pv,NB);
    for(auto j = 0; j < NB; j++) pv[j*NB+j] = MatsT(1.) + pv[j*NB+j];
    /*
     * construct Roothaan's effective fock
        ======== ======== ====== =========
        space     closed   open   virtual
        ======== ======== ====== =========
        closed      Fc      Fb     Fc
        open        Fb      Fc     Fa
        virtual     Fc      Fa     Fc
        ======== ======== ====== =========
        where Fc = (Fa + Fb) / 2
     */
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(0.5),ss.fockMatrix->S().pointer(),NB,
         pc,NB,MatsT(0.),tmp2,NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(0.5),pc,NB,tmp2,NB,MatsT(0.),tmp.pointer(),NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),pv,NB,tmp2,NB,MatsT(1.),tmp.pointer(),NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(0.5),ss.fockMatrix->S().pointer(),NB,
         po,NB,MatsT(0.),tmp2,NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(0.5),po,NB,tmp2,NB,MatsT(1.),tmp.pointer(),NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(0.5),ss.fockMatrix->S().pointer(),NB,
         pv,NB,MatsT(0.),tmp2,NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(0.5),pv,NB,tmp2,NB,MatsT(1.),tmp.pointer(),NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),SCR[1].pointer(),NB,pc,NB,MatsT(0.),tmp2,NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),po,NB,tmp2,NB,MatsT(1.),tmp.pointer(),NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),SCR[0].pointer(),NB,pv,NB,MatsT(0.),tmp2,NB);
    blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,NB,NB,NB,MatsT(1.),po,NB,tmp2,NB,MatsT(1.),tmp.pointer(),NB);
    MatAdd('C','N',NB,NB,MatsT(1.),tmp.pointer(),NB,MatsT(1.),
           tmp.pointer(),NB,ss.fockMatrix->S().pointer(),NB);

    CQMemManager::get().free(pc,po,pv,tmp2);
  }; // ROFock<MatsT, IntsT>::rohfFock

  /**
   *  \brief Forms the Roothaan's effective Fock matrix for a single slater 
   *  determinant using the 1PDM.
   *
   *  \param [in] increment Whether or not the Fock matrix is being
   *  incremented using a previous density
   *
   *  Populates / overwrites fock strorage in SingleSlater &ss
   */
  template <typename MatsT, typename IntsT>
  void ROFock<MatsT,IntsT>::formFock(SingleSlater<MatsT,IntsT> &ss,
    EMPerturbation &pert, bool increment, double xHFX) {

    // General fock build
    FockBuilder<MatsT,IntsT>::formFock(ss, pert, increment, xHFX);

    // ROHF fock build
    rohfFock(ss);

  } // ROFock<MatsT,IntsT>::formFock


}; // namespace ChronusQ
