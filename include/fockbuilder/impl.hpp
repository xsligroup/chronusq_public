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

#include <fockbuilder.hpp>
#include <memory>
#include <util/timer.hpp>
#include <cqlinalg.hpp>
#include <matrix.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <particleintegrals/gradints/incore.hpp>
#include <particleintegrals/gradints/direct.hpp>
#include <fockbuilder/rofock/impl.hpp>
#include <quantum/properties.hpp>
#include <fockbuilder/neofock.hpp>
#include <fockbuilder/interparticlefock.hpp>
#include <fockbuilder/fourcompfock/impl.hpp>
#include <fockbuilder/fourcompfock/batchgd.hpp>
#include <fockbuilder/matrixfock.hpp>

#include <typeinfo>

//#define USE_ONEPDMGRAD

namespace ChronusQ {

  /**
   *  Constructs a FockBuilder object from another of a another (possibly the
   *  same) type by copy.
   *
   *  \param [in] other FockBuilder object to copy
   */
  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  FockBuilder<MatsT,IntsT>::FockBuilder(const FockBuilder<MatsU,IntsT> &other):
      hamiltonianOptions_(other.hamiltonianOptions_){}

  /**
   *  Constructs a FockBuilder object from another of a another (possibly the
   *  same) by move.
   *
   *  \warning Deallocates the passed FockBuilder object
   *
   *  \param [in] other FockBuilder object to move
   */
  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  FockBuilder<MatsT,IntsT>::FockBuilder(FockBuilder<MatsU,IntsT> &&other):
      hamiltonianOptions_(other.hamiltonianOptions_){}

  /**
   *  \brief Forms the Hartree-Fock perturbation tensor
   *
   *  Populates / overwrites GD storage (and JScalar and K storage)
   */
  template <typename MatsT, typename IntsT>
  void FockBuilder<MatsT,IntsT>::formGD(SingleSlater<MatsT,IntsT> &ss,
    EMPerturbation &pert, bool increment, double xHFX, bool HerDen) {
  
    // Decide list of onePDMs to use
    cqmatrix::PauliSpinorMatrices<MatsT> &contract1PDM
        = increment ? *ss.deltaOnePDM : *ss.onePDM;
    
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> 
      onePDMs, coulombMatrices, exchangeMatrices, twoeHs;
    
    // setup pointers 
    if (increment) onePDMs.push_back(ss.deltaOnePDM);
    else onePDMs.push_back(ss.onePDM);
     
    exchangeMatrices.push_back(ss.exchangeMatrix);
    twoeHs.push_back(ss.twoeH);
    
    coulombMatrices.push_back(
      std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(
      ss.coulombMatrix->nRows(), false, false)
    );

    formRawGDInBatches(ss, pert, increment, xHFX, HerDen, onePDMs, coulombMatrices, exchangeMatrices, twoeHs);
    
    * ss.coulombMatrix = coulombMatrices[0]->S(); 
    
  } // FockBuilder::formGD 
  
  /**
   *  \brief Forms the Hartree-Fock perturbation tensor
   *
   *  actual work 
   */
  template <typename MatsT, typename IntsT>
  void FockBuilder<MatsT,IntsT>::formRawGDInBatches(SingleSlater<MatsT,IntsT> &ss,
    EMPerturbation &pert, bool increment, double xHFX, bool HerDen, 
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & onePDMs, 
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & coulombMatrices, 
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & exchangeMatrices,
    std::vector<std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>>> & twoeHs) {

    size_t NB = ss.basisSet().nBasis;
    size_t nBatch = onePDMs.size();
    bool computeCoulomb  = coulombMatrices.size() > 0;
    bool computeExchange = (std::abs(xHFX) > 1e-12) and exchangeMatrices.size() > 0;
    bool computeTwoeHs   = twoeHs.size() > 0;
    
    if (not computeCoulomb and not computeExchange and not computeTwoeHs) {
     CErr("Nothing specified to compute in FockBuilder::formRawGDInBatches");
    }  

    if( ss.nC == 4 )
      CErr("4C formGD is implemented in class FourCompFock.");

    // Zero out J and K[i]
    if(not increment) {
      for (auto i = 0ul; i < nBatch; i++) {
        if (computeCoulomb)  coulombMatrices[i]->clear();
        if (computeExchange) exchangeMatrices[i]->clear();
        if (computeTwoeHs)   twoeHs[i]->clear();
      }
    }

    std::vector<TwoBodyContraction<MatsT>> contract;
    
    if (computeCoulomb or computeTwoeHs) {
      auto & coulombContainers = computeCoulomb ? coulombMatrices: twoeHs;
      for (auto i = 0ul; i < nBatch; i++) {
        contract.push_back(
          {onePDMs[i]->S().pointer(), coulombContainers[i]->S().pointer(), HerDen, COULOMB}
        );
      }
    }

    if(ss.TPI->printContractionTiming)
        std::cout << "      " << std::string(ss.particle.charge>0 ? "Protonic" : "Electronic") << " Subsystem Symm-Contraction Timing: " << std::endl;

    auto beginContract = tick();

    auto distributed_ritpi = std::dynamic_pointer_cast<DistributedRITPIContraction<MatsT, IntsT>>(ss.TPI);
    // Determine how many (if any) exchange terms to calculate
    if( std::abs(xHFX) > 1e-12 and not increment and ss.nC == 1 and
        (std::dynamic_pointer_cast<InCoreRITPIContraction<MatsT, IntsT>>(ss.TPI) or
         (distributed_ritpi && distributed_ritpi->canUseKCoef()))) {

      auto ritpi_incore = std::dynamic_pointer_cast<InCoreRITPIContraction<MatsT, IntsT>>(ss.TPI);
      auto ritpi_dist   = std::dynamic_pointer_cast<DistributedRITPIContraction<MatsT, IntsT>>(ss.TPI);
      auto riKCoeffBegin = tick();
      bool isRoot = (MPIRank(ss.comm) == 0);

      const bool useCholeskyMOs = !ss.denEqCoeff_;
      //const bool useCholeskyMOs = true;
      decltype(ss.getCholeskyMOs()) choleskyMOs;
      if (useCholeskyMOs && isRoot)
      choleskyMOs = ss.getCholeskyMOs();

      // Contract one spin block: handles MO selection, broadcast, dispatch
      auto contractSpin = [&](int spin) -> cqmatrix::Matrix<MatsT> {
      size_t nO = 0;
      MatsT* mo = nullptr;
      if (isRoot) {
        if (useCholeskyMOs) {
          nO = choleskyMOs[spin]->nColumns();
          mo = choleskyMOs[spin]->pointer();
        } else {
          nO = (spin == 0) ? ss.nOA : ss.nOB;
          mo = ss.mo[spin].pointer();
        }
      }

      cqmatrix::Matrix<MatsT> Kblock(NB);

      if (ritpi_incore) {
        if (isRoot && nO > 0)
          ritpi_incore->KCoefContract(ss.comm, nO, mo, Kblock.pointer());
        else
          Kblock.clear();
      } else {
      #ifdef CQ_ENABLE_MPI
        MPIBCast(&nO, 1, 0, ss.comm);
      #endif
        if (nO > 0) {
          auto* mo_buf = CQMemManager::get().malloc<MatsT>(NB * nO);
          if (isRoot) std::copy_n(mo, NB * nO, mo_buf);
      #ifdef CQ_ENABLE_MPI
          MPIBCast(mo_buf, NB * nO, 0, ss.comm);
      #endif
          ritpi_dist->KCoefContract(ss.comm, nO, mo_buf, Kblock.pointer());
          CQMemManager::get().free(mo_buf);
        } else {
          Kblock.clear();
        }
      }
      return Kblock;
      };

      auto AAblock = contractSpin(0);
      auto BBblock = ss.iCS ? cqmatrix::Matrix<MatsT>(NB) : contractSpin(1);

      // Only root assembles spin-block exchange matrices
      if (isRoot) {
      for (auto i = 0ul; i < nBatch; i++)
        *exchangeMatrices[i] = ss.iCS
          ? cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(AAblock)
          : cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(AAblock, BBblock);
      }

      if (ss.TPI->printContractionTiming) {
        int wRank = 0, wSize = 1;
#ifdef CQ_ENABLE_MPI
        if (ritpi_dist) { MPI_Comm_rank(ss.comm, &wRank); MPI_Comm_size(ss.comm, &wSize); }
#endif
        ss.TPI->printTiming("K-Coeff-Contraction duration (s): ", tock(riKCoeffBegin), ss.comm, wRank, wSize);
      }

    } else if(computeExchange) {

      // Use density to do K contraction
      for (auto i = 0ul; i < nBatch; i++) { 
        contract.push_back(
            {onePDMs[i]->S().pointer(), exchangeMatrices[i]->S().pointer(), HerDen, EXCHANGE}
        );
        // If this is a protonic ss, then beta block is 0, then 2S=AA+BB should be equal to 2Z=AA-BB
        // We can avoid doing K contraction for Z component
        if (exchangeMatrices[i]->hasZ() and ss.particle.charge<0)
          contract.push_back(
            {onePDMs[i]->Z().pointer(), exchangeMatrices[i]->Z().pointer(), HerDen, EXCHANGE}
          );
        // TangDD: do full K contraction for 2C proton
        if (exchangeMatrices[i]->hasXY() and ss.particle.charge>0)
          contract.push_back(
            {onePDMs[i]->Z().pointer(), exchangeMatrices[i]->Z().pointer(), HerDen, EXCHANGE}
          ); 
        if (exchangeMatrices[i]->hasXY()) {
          contract.push_back(
            {onePDMs[i]->Y().pointer(), exchangeMatrices[i]->Y().pointer(), HerDen, EXCHANGE}
          );
          contract.push_back(
            {onePDMs[i]->X().pointer(), exchangeMatrices[i]->X().pointer(), HerDen, EXCHANGE}
          );
        }
      }

    }

    ss.TPI->twoBodyContract(ss.comm, contract, pert);

    // Copy the S component of protonic K matrix to its Z component
    if(computeExchange and ss.particle.charge>0)
      // TangDD: do full K contraction for 2C proton 
      if (ss.nC == 1)
        for (auto i = 0ul; i < nBatch; i++)  std::copy_n(exchangeMatrices[i]->S().pointer(),NB*NB,exchangeMatrices[i]->Z().pointer());

    if(ss.TPI->printContractionTiming)
        std::cout << "        " << std::left << std::setw(38) << "Cholesky-Symm-Contraction duration = " << tock(beginContract) << " s " << std::endl;

    ROOT_ONLY(ss.comm); // Return if not root (J/K only valid on root process)

    if (computeTwoeHs) {
      for (auto i = 0ul; i < nBatch; i++) { 
        // G[D] += 2*J[D]
        if (computeCoulomb) {
          *twoeHs[i] += 2.0 * *coulombMatrices[i];
        } else {
          *twoeHs[i] *= 2.0;
        }

        // Form GD: G[D] = 2.0*J[D] - K[D]
        if (computeExchange) {
          *twoeHs[i] -= xHFX * *exchangeMatrices[i];
        } 
      }
    }
#if 0
  //printJ(std::cout);
    printK(std::cout);
  //printGD(std::cout);
#endif

  } // FockBuilder::formGD

  /*******************************************************************************/
  /* Compute memory requirement for build 4C GD in Batches                       */
  /* Returns:                                                                    */
  /*   size_t SCR size needed for one batch                                      */
  /*   IMPORTANT HERE: size are all in MatsT                                     */
  /*******************************************************************************/
  template <typename MatsT, typename IntsT>
  size_t FockBuilder<MatsT,IntsT>::formRawGDSCRSizePerBatch(SingleSlater<MatsT,IntsT> &ss,
    bool computeExchange, bool HerDen) const {
  
      size_t SCRSize  = 0ul;
  
      if( std::dynamic_pointer_cast<GTODirectTPIContraction<MatsT,IntsT>>(ss.TPI) ) {
        
        GTODirectTPIContraction<MatsT,IntsT> &ERICon =
            *std::dynamic_pointer_cast<GTODirectTPIContraction<MatsT,IntsT>>(ss.TPI);
        
        size_t nConPerBatch = computeExchange ? 5: 1;

        SCRSize += ERICon.directScaffoldNewSCRSize() * nConPerBatch;
      } 
      
      return SCRSize;
  } // FockBuilder::formRawGDSCRSizePerBatch

  /**
   *  \brief Forms the Fock matrix for a single slater determinant using
   *  the 1PDM.
   *
   *  \param [in] increment Whether or not the Fock matrix is being
   *  incremented using a previous density
   *
   *  Populates / overwrites fock strorage in SingleSlater &ss
   */
  template <typename MatsT, typename IntsT>
  void FockBuilder<MatsT,IntsT>::formFock(SingleSlater<MatsT,IntsT> &ss,
    EMPerturbation &pert, bool increment, double xHFX) {

    auto GDStart = tick(); // Start time for G[D]

    // Form G[D]
    // See: https://doi.org/10.1016/j.cplett.2005.01.115
    // regarding the option to NOT include the proton two body term
    // TODOAL: Fix the logic here for hard-coded ss.particle.charge != -1.0
    if(ss.particle.charge != -1.0 && (this->hamiltonianOptions_.ignoreProtonTwoBody || ss.nO == 1))
    {
      ss.twoeH->clear();
    }
    else
    {
      formGD(ss,pert,increment,xHFX);
    }

    ss.GDDur = tock(GDStart); // G[D] Duraction
    //std::cout <<"formGD time = "<<ss.GDDur <<std::endl;

    ROOT_ONLY(ss.comm);

    // Form Fock
    *ss.fockMatrix = *ss.coreH + *ss.twoeH;


    // Add in the electric field contributions
    // FIXME: the magnetic field contribution should go here as well to allow for RT
    // manipulation

    if( pert_has_type(pert,Electric) ) {

      auto dipAmp = pert.getDipoleAmp(Electric);

      if (ss.nC == 4) {
        std::vector<cqmatrix::PauliSpinorMatrices<dcomplex>>
          lenElectric4C = *(ss.aoints_->lenElectric->gather4CDipole());
        for (auto i = 0; i < 3; i++)
          ss.fockMatrix->S() -= 2.0 * dipAmp[i] * lenElectric4C[i].S();
      } else {
        for (auto i = 0; i < 3; i++)
          ss.fockMatrix->S() -= 2. * dipAmp[i] * (*ss.aoints_->lenElectric)[i]->matrix();
      }
    }

#if 0
    ss.printFock(std::cout);
#endif
  }



  template <typename MatsT, typename IntsT>
  void MatrixFock<MatsT,IntsT>::formFock(SingleSlater<MatsT,IntsT> &ss,
                                         EMPerturbation &pert, bool increment, double xHFX) {

    *ss.fockMatrix = fockMatrix;

    ROOT_ONLY(ss.comm);
    *ss.twoeH = fockMatrix - *ss.coreH;

  }

  /**
   *  \brief The pointer convertor. This static function converts
   *  the underlying polymorphism correctly to hold a different
   *  type of matrices. It is called when the corresponding
   *  SingleSlater object is being converted.
   */
  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  std::shared_ptr<FockBuilder<MatsU,IntsT>>
  FockBuilder<MatsT,IntsT>::convert(const std::shared_ptr<FockBuilder<MatsT,IntsT>>& fb) {

    if (not fb) return nullptr;

    const std::type_info &tID(typeid(*fb));

    if (tID == typeid(ROFock<MatsT,IntsT>)) {
      return std::make_shared<ROFock<MatsU,IntsT>>(
               *std::dynamic_pointer_cast<ROFock<MatsT,IntsT>>(fb));
    }
    else if (tID == typeid(NEOFockBuilder<MatsT,IntsT>)) {
      return std::make_shared<NEOFockBuilder<MatsU,IntsT>>(
               *std::dynamic_pointer_cast<NEOFockBuilder<MatsT,IntsT>>(fb));
    } 
    else if (tID == typeid(FourCompFock<MatsT,IntsT>)) {
      return std::make_shared<FourCompFock<MatsU,IntsT>>(
          *std::dynamic_pointer_cast<FourCompFock<MatsT,IntsT>>(fb));

    } else if (tID == typeid(MatrixFock<MatsT,IntsT>)) {
      return std::make_shared<MatrixFock<MatsU,IntsT>>(
          *std::dynamic_pointer_cast<MatrixFock<MatsT,IntsT>>(fb));

    } else {
      return std::make_shared<FockBuilder<MatsU,IntsT>>(
               *std::dynamic_pointer_cast<FockBuilder<MatsT,IntsT>>(fb));
    }

  } // FockBuilder<MatsT,IntsT>::convert

  
  template <typename MatsT, typename IntsT>
  std::vector<double> FockBuilder<MatsT,IntsT>::getGDGrad(
    SingleSlater<MatsT,IntsT>& ss, EMPerturbation& pert, double xHFX) {

    size_t NB = ss.basisSet().nBasis;
    size_t nGrad = 3*ss.molecule().nAtoms;

    bool hasXY = ss.exchangeMatrix->hasXY();
    bool hasZ = ss.exchangeMatrix->hasZ();

    if( not ss.aoints_->gradERI )
      CErr("Gradient ERI missing in FockBuilder::getGDGrad!");

    GradInts<TwoPInts,IntsT>& gradERI = *ss.aoints_->gradERI;

    // Form contraction
    // TODO: There's gotta be a better way to do this...
    std::unique_ptr<GradContractions<MatsT,IntsT>> contract = nullptr;
    if ( std::dynamic_pointer_cast<InCore4indexTPI<IntsT>>(gradERI[0]) ) {
      contract = std::make_unique<InCore4indexGradContraction<MatsT,IntsT>>(gradERI);
    }
    else if ( std::dynamic_pointer_cast<DirectTPI<IntsT>>(gradERI[0]) ) {
      contract = std::make_unique<DirectGradContraction<MatsT,IntsT>>(gradERI);
    }
    else
      CErr("Gradients of RI NYI!");

    // Create contraction list
    std::vector<std::vector<TwoBodyContraction<MatsT>>> cList;

    std::vector<cqmatrix::Matrix<MatsT>> JList;
    std::vector<cqmatrix::PauliSpinorMatrices<MatsT>> KList;

    JList.reserve(nGrad);
    KList.reserve(nGrad);

    for( auto iGrad = 0; iGrad < nGrad; iGrad++ ) {
      std::vector<TwoBodyContraction<MatsT>> tempCont;

      // Coulomb
      JList.emplace_back(NB);
      JList.back().clear();
      tempCont.push_back(
         {ss.onePDM->S().pointer(), JList.back().pointer(), true, COULOMB}
      );

      // Exchange
      if( std::abs(xHFX) > 1e-12 ) {

        KList.emplace_back(NB, hasXY, hasZ);
        KList.back().clear();

        tempCont.push_back(
          {ss.onePDM->S().pointer(), KList.back().S().pointer(), true, EXCHANGE}
        );

        if (hasZ) {
          tempCont.push_back(
            {ss.onePDM->Z().pointer(), KList.back().Z().pointer(), true, EXCHANGE}
          );
        }
        if (hasXY) {
          tempCont.push_back(
            {ss.onePDM->Y().pointer(), KList.back().Y().pointer(), true, EXCHANGE}
          );
          tempCont.push_back(
            {ss.onePDM->X().pointer(), KList.back().X().pointer(), true, EXCHANGE}
          );
        }
      }


      cList.push_back(tempCont);
    }

    // Contract to J/K
    contract->gradTwoBodyContract(MPI_COMM_WORLD, true, cList, pert);

    // Contract to gradient
    std::vector<double> gradient;
    cqmatrix::PauliSpinorMatrices<MatsT> twoEGrad(NB, hasXY, hasZ);

    for( auto iGrad = 0; iGrad < nGrad; iGrad++ ) {

      if( std::abs(xHFX) > 1e-12 ){
        // Scale K by alpha
        twoEGrad = -xHFX * KList[iGrad];
        // G[S] = 2 * J[S] + alpha * K[S]
        twoEGrad.S() += 2. * JList[iGrad];
      } else{
        twoEGrad.S() = 2. * JList[iGrad];
      }

      double gradVal = ss.template computeOBProperty<SCALAR>(
        twoEGrad.S().pointer()
      );
      if( hasZ )
        gradVal += ss.template computeOBProperty<MZ>(
          twoEGrad.Z().pointer()
        );
      if( hasXY ) {
        gradVal += ss.template computeOBProperty<MY>(
          twoEGrad.Y().pointer()
        );
        gradVal += ss.template computeOBProperty<MX>(
          twoEGrad.X().pointer()
        );
      }
      gradient.push_back(0.25*gradVal);

    }

    return gradient;

  } // FockBuilder::getGDGrad

  
  // Pulay contribution
  //
  // We have two options here: 
  // 1. Compute the energy-weighted density matrix W and compute W * dS/dR 
  //    (Int. J. Quantum Chem., Quant. Chem. Symp., S13 (1979) 225-41. DOI: 10.1002/qua.560160825)
  // 2. Compute dV/dR and use that to compute the gradient 
  //    (J. Chem. Phys. 22 August 2005; 123 (8): 084106. DOI: 10.1063/1.2008258)
  template <typename MatsT, typename IntsT>
  std::vector<double> FockBuilder<MatsT,IntsT>::getPulayGrad(
    SingleSlater<MatsT,IntsT>& ss, bool equil, bool useW) {

    std::vector<double> pulayGrad;

    size_t NB = ss.basisSet().nBasis;
    size_t nGrad = 3*ss.molecule().nAtoms;
    size_t nSp = ss.fockMatrix->nComponent();
    bool hasXY = ss.exchangeMatrix->hasXY();
    bool hasZ = ss.exchangeMatrix->hasZ();

    if( useW ) {

      ss.formEWDM(equil);
      for( size_t iGrad = 0; iGrad < nGrad; iGrad++ ) {
        double gradVal = std::real(blas::dot(NB*NB, ss.W->S().pointer(), 1, (*ss.aoints_->gradOverlap)[iGrad]->pointer(), 1));
        if (hasZ) {
          gradVal += std::real(blas::dot(NB*NB, ss.W->Z().pointer(), 1, (*ss.aoints_->gradOverlap)[iGrad]->pointer(), 1));
        }
        if (hasXY) {
          gradVal += std::real(blas::dot(NB*NB, ss.W->Y().pointer(), 1, (*ss.aoints_->gradOverlap)[iGrad]->pointer(), 1));
          gradVal += std::real(blas::dot(NB*NB, ss.W->X().pointer(), 1, (*ss.aoints_->gradOverlap)[iGrad]->pointer(), 1));
        }
        pulayGrad.push_back(-gradVal);
      }

    } else {

      if (hasXY)
        CErr("Computing Pulay gradient with dV/dR not yet implemented for XY components!");

      // S^{-1/2}
      auto orthoForward = ss.orthoAB->forwardPointer();
      //auto orthoForward = orthoSpinor->forwardPointer();

      // Allocate
      cqmatrix::Matrix<MatsT> vdv(NB);
      cqmatrix::Matrix<MatsT> dvv(NB);
      cqmatrix::PauliSpinorMatrices<MatsT> SCR(NB, hasXY, hasZ);

      // allocate one-PDM gradient matrices
      if (ss.onePDMGrad.size() == 0) {
        ss.onePDMGrad.reserve(nGrad);
        for( size_t iGrad = 0; iGrad < nGrad; iGrad++ ) 
        ss.onePDMGrad.emplace_back(std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, hasXY, hasZ));
      }

      // XXX: This requires copying the overlap gradients, but it is for
      //      copying to MatsT != IntsT
      std::vector<cqmatrix::Matrix<MatsT>> gradOverlap;
      gradOverlap.reserve(nGrad);
      for( size_t iGrad = 0; iGrad < nGrad; iGrad++ ) {
        gradOverlap.emplace_back((*ss.aoints_->gradOverlap)[iGrad]->matrix());
      }

      // Calculate dV
      std::vector<cqmatrix::Matrix<MatsT>> gradOrtho;
      gradOrtho.reserve(nGrad);
      for( size_t iGrad = 0; iGrad < nGrad; iGrad++ ) {
        gradOrtho.emplace_back(NB);
      }
      ss.orthoAB->getOrthogonalizationGradients(gradOrtho, gradOverlap);

      for( size_t iGrad = 0; iGrad < nGrad; iGrad++ ) {

        // Form VdV and dVV
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
          NB,NB,NB,MatsT(1.),orthoForward->pointer(),NB,
          gradOrtho[iGrad].pointer(),NB,MatsT(0.),vdv.pointer(),NB);
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
          NB,NB,NB,MatsT(1.),gradOrtho[iGrad].pointer(),NB,
          orthoForward->pointer(),NB,MatsT(0.),dvv.pointer(),NB);

        // Form FVdV and dVVF for the non-xc part of F
        for( auto iSp = 0; iSp < nSp; iSp++ ) {
          auto comp = static_cast<cqmatrix::PAULI_SPINOR_COMPS>(iSp);
          
          cqmatrix::Matrix<MatsT> nonXC_F(NB);
          //nonXC_F = (*coreH)[comp] + (*twoeH)[comp];
          nonXC_F = (*ss.fockMatrix)[comp];

          
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
            NB,NB,NB,MatsT(1.),nonXC_F.pointer(),NB,
            vdv.pointer(),NB,MatsT(0.),SCR[comp].pointer(),NB);
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
            NB,NB,NB,MatsT(1.),dvv.pointer(),NB,
            nonXC_F.pointer(),NB,MatsT(1.),SCR[comp].pointer(),NB);
        }
  
#ifdef USE_ONEPDMGRAD
        // Compute 1PDM Gradient and Save: dP/dR = -(VdV * P + P * dVV)
        // S part:
        // Compute VdV * P and negate the result
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
          NB,NB,NB,MatsT(-1.),vdv.pointer(),NB, // Note the change here from 1. to -1.
          this->onePDM->S().pointer(),NB,MatsT(0.),onePDMGrad[iGrad]->S().pointer(),NB);

        // Compute P * dVV, add to VdV * P with sign change, resulting in -(P * dVV + VdV * P)
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
          NB,NB,NB,MatsT(-1.),this->onePDM->S().pointer(),NB, // Note the change here from 1. to -1.
          dvv.pointer(),NB,MatsT(1.),onePDMGrad[iGrad]->S().pointer(),NB); 

        // Mz part:
        if( hasZ ){
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
            NB,NB,NB,MatsT(-1.),vdv.pointer(),NB, 
            this->onePDM->Z().pointer(),NB,MatsT(0.),onePDMGrad[iGrad]->Z().pointer(),NB);
          blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
            NB,NB,NB,MatsT(-1.),this->onePDM->Z().pointer(),NB,
            dvv.pointer(),NB,MatsT(1.),onePDMGrad[iGrad]->Z().pointer(),NB); 
        }
#endif

        // Trace
        double gradVal = ss.template computeOBProperty<SCALAR>(
          SCR.S().pointer()
        );

        if( hasZ )
          gradVal += ss.template computeOBProperty<MZ>(
            SCR.Z().pointer()
        );
        if( hasXY ) {
          gradVal += ss.template computeOBProperty<MY>(
            SCR.Y().pointer()
          );
          gradVal += ss.template computeOBProperty<MX>(
            SCR.X().pointer()
          );
        }

        pulayGrad.push_back(-0.5*gradVal);
      }
    }

    return pulayGrad;

  } // FockBuilder::getPulayGrad

}; // namespace ChronusQ
