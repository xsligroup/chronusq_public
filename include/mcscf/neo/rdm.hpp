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

#include <mcscf.hpp>
#include <util/matout.hpp>
#include <util/print.hpp>
#include <cxxapi/output.hpp>
#include <mointstransformer/moranges.hpp>

namespace ChronusQ {

    template <typename MatsT, typename IntsT>
    void NEOMCSCF<MatsT,IntsT>::computeOneRDM(size_t i)
    {
        this->ciBuilder->computeOneRDM(*this,this->CIVecs[i],this->oneRDM[i]);
    }

    template <typename MatsT, typename IntsT>
    void NEOMCSCF<MatsT,IntsT>::computeOneRDM()
    {
        for(size_t i = 0; i < this->NStates; i++)
        {
            computeOneRDM(i);
        }
    }

    template <typename MatsT, typename IntsT>
    void NEOMCSCF<MatsT,IntsT>::computePOneRDM(size_t i)
    {
//        std::shared_ptr<NEOCASCI<MatsT,IntsT>> NEOCI = std::make_shared<NEOCASCI<MatsT,Ints>>(dynamic_cast(this->ciBuilder));
        NEOCIBuilder->computePOneRDM(*this,this->CIVecs[i],this->PoneRDM[i]);
    }

    template <typename MatsT, typename IntsT>
    void NEOMCSCF<MatsT,IntsT>::computePOneRDM()
    {
        for(size_t i = 0; i < this->NStates; i++)
        {
            computePOneRDM(i);
        }
    }

    // TODO: The two functions below are specialized to doubles to get this working quickly
    // These functions should probably be extracted elsewhere, or made more general members
    // of the parent classes of (NEO)MCSCF

    template <typename MatsT, typename IntsT>
    std::vector<std::shared_ptr<cqmatrix::Matrix<double>>> NEOMCSCF<MatsT,IntsT>::getOnePDM()
    {
        // Get the orbital offsets for the state
        size_t nInact = this->ewfn_->MOPartition.nInact;
        size_t nCorrO = this->ewfn_->MOPartition.nCorrO;
        size_t nAO = this->ewfn_->ref_.mo[0].nRows();
        double * MO = (double*)this->ewfn_->ref_.mo[0].pointer() + nAO * nInact;

        std::vector<std::shared_ptr<cqmatrix::Matrix<double>>> PDMs;
        PDMs.reserve(this->NStates);

        cqmatrix::Matrix<double> SCR(nAO);
        for(size_t i = 0; i < this->NStates; i++)
        {
            cqmatrix::Matrix<double> PDM(nAO);
            double * rdm = (double*)this->oneRDM[i].pointer();
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,nAO,nCorrO,nCorrO,1.0,MO,nAO,rdm,nCorrO,0.0,SCR.pointer(),nAO);
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,nAO,nAO,nCorrO,1.0,SCR.pointer(),nAO,MO,nAO,0.0,PDM.pointer(),nAO);
            PDMs.emplace_back(std::make_shared<cqmatrix::Matrix<double>>(PDM));
        }
        return PDMs;
    }

    template <typename MatsT, typename IntsT>
    std::vector<std::shared_ptr<cqmatrix::Matrix<double>>> NEOMCSCF<MatsT,IntsT>::getPOnePDM()
    {
        // Get the orbital offsets for the state
        size_t nInact = this->pwfn_->MOPartition.nInact;
        size_t nCorrO = this->pwfn_->MOPartition.nCorrO;
        size_t nAO = this->pwfn_->ref_.mo[0].nRows();
        double * MO = (double*)this->pwfn_->ref_.mo[0].pointer() + nAO * nInact;

        std::vector<std::shared_ptr<cqmatrix::Matrix<double>>> PDMs;
        PDMs.reserve(this->NStates);

        cqmatrix::Matrix<double> SCR(nAO);
        for(size_t i = 0; i < this->NStates; i++)
        {
            cqmatrix::Matrix<double> PDM(nAO);
            double * rdm = (double*)this->PoneRDM[i].pointer();
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,nAO,nCorrO,nCorrO,double(1.0),MO,nAO,rdm,nCorrO,double(0.0),SCR.pointer(),nAO);
            blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::Trans,nAO,nAO,nCorrO,double(1.0),SCR.pointer(),nAO,MO,nAO,double(0.0),PDM.pointer(),nAO);
            PDMs.emplace_back(std::make_shared<cqmatrix::Matrix<double>>(PDM));
        }
        return PDMs;
    }

}; // namespace ChronusQ
