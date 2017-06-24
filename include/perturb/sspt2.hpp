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

#include <perturb.hpp>

#include <util/timer.hpp>
#include <util/matout.hpp>
#include <util/threads.hpp>

//#define _DEBUG_PT2_impl

#define PT2_LOOP_INIT() \
  size_t nC = refMCwfn->referenceWaveFunction().nC; \
  auto ciBuilder = dynamic_cast<RASCI<MatsT,IntsT>&>(*this->ciBuilder); \
  auto  RASString = std::dynamic_pointer_cast<RASStringManager>(this->detStr); \
  std::vector<size_t> nActO = this->MOPartition.nActOs; \
  std::vector<int> fCat = this->MOPartition.fCat; \
  size_t nCat = RASString->nCategory(); \
  size_t mxHole = RASString->maxHole(); \
  size_t mxElec = RASString->maxElectron(); \
  std::vector<std::vector<size_t>> LCat = RASString->LCategory(); \
  auto & ptFock = *ptFock_;


namespace ChronusQ {

  typedef std::vector<std::vector<int>> int_matrix;

  /**
   * 
   *  \brief build PT fock matrix (dimension nCorrO X nCorrO);
   *            use state-specific 1RDM or state-averaged 1RDM;
   *            currently building fock for correlated space.
   *
   */
  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::buildFock(cqmatrix::Matrix<MatsT> & ptFock,
                cqmatrix::Matrix<MatsT> & oneRDM) {

    ProgramTimer::tick("Fock build");
    auto & mopart = this->MOPartition;
    auto & nActO = mopart.nActOs;
    size_t nCorrO = mopart.nCorrO;

    auto & hCore = *(this->moints->template getIntegral<OnePInts, MatsT>
                                                    ("hCore_Correlated_Space")); 
    auto & moERI = *(this->moints->template getIntegral<InCore4indexTPI, MatsT>
                                                    ("ERI_Correlated_Space"));

    size_t u, v, w, t;
    for (auto u = 0ul; u < nCorrO; u++)
    for (auto v = 0ul; v < nCorrO; v++) {
      ptFock(u, v) = hCore(u, v);
      for (auto w = 0ul; w < nActO[0]+nActO[1]; w++)
      for (auto t = 0ul; t < nActO[0]+nActO[1]; t++) {
        if (w >= nActO[0] and t >= nActO[0])
          ptFock(u, v) += oneRDM(w - nActO[0], t - nActO[0])*(moERI(u, v, w, t)
                            - moERI(u, t, w, v));
        else if (w == t)
          ptFock(u, v) += (moERI(u, v, w, t) - moERI(u, t, w, v));
      }
    }

    ProgramTimer::tock("Fock build");    
#ifdef _DEBUG_PT2_impl
    prettyPrintSmart(std::cout,"LL PT2 fock -- ", ptFock.pointer(), nCorrO, nCorrO, nCorrO);
#endif

  } // PERTURB::buildFock

  /**
   *  
   *   \brief Compute the zero-order Energy <0|F|0>
   *              make use of the oneRDM E0 = \sum_uv f_uv * D_uv;
   *              E0 does not include fock in frozen core orbital space.
   *
   */
  template <typename MatsT, typename IntsT>
  MatsT PERTURB<MatsT,IntsT>::computeZeroE(cqmatrix::Matrix<MatsT> & ptFock,
                  cqmatrix::Matrix<MatsT> & oneRDM) {

    auto & mopart = this->MOPartition;
    auto & nActO = mopart.nActOs;

    MatsT E = MatsT(0.);
    for (auto u = 0ul; u < nActO[0]+nActO[1]; u++)
    for (auto v = 0ul; v < nActO[0]+nActO[1]; v++) {
      if (u >= nActO[0] and v >= nActO[0])
        E += ptFock(u, v) * oneRDM(u - nActO[0], v - nActO[0]);
      else if (u == v)
        E += ptFock(u, v);
    }

    std::cout<<"zero-order energy: "<< std::setprecision(12)<<E<<std::endl;
    return E;

  } // PERTURB::computeZeroE


  /**
   * 
   *  \brief Full LHS: (F - E0 * I)
   *    F_IJ = <I|F|J> = sum_(uv) f_uv * <I|E_(uv)|J>
   *    E0 is the reference zero-order energy
   *
   */
  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::buildFull(cqmatrix::Matrix<MatsT> &LHS, MatsT E0) {

    ProgramTimer::tick("Full LHS build");
    PT2_LOOP_INIT(); // check top for details

    // Initialize LHS
    std::fill_n(LHS.pointer(),SDsize*SDsize, MatsT(0.));

    std::shared_ptr<const ExcitationList> LrsList;
    size_t iCatJ, iCatL;
    std::vector<size_t> LCatAlt(3, 0);
    std::vector<size_t> LCatAlt2(3, 0);
    std::vector<size_t> nDetAlt(3, 0);
    std::vector<size_t> nDetJ(3, 0);
    size_t r, s, l, l1, l2, j, j1, j2, j3, iNZ, j23off;
    double sign;

    for (auto rBlk = 0ul; rBlk < 3; rBlk++)
    for (auto sBlk = 0ul; sBlk < 3; sBlk++) {
      //std::cout<<"rBlk: "<<rBlk<<", sBlk: "<<sBlk<<std::endl;
      if (nActO[rBlk] == 0 or nActO[sBlk] == 0) continue;
      iCatJ = 0;
      for (auto iEJ = 0ul; iEJ <= mxElec; iEJ++)
      for (auto iHJ = 0ul; iHJ <= mxHole; iHJ++, iCatJ++) {
        //std::cout<<"LCat"<<iCatJ<<" : "<<LCat[iCatJ][3]<<std::endl;
        if (iCatJ == 0 or LCat[iCatJ][3] == 0) continue;
        for (iCatL = 1; iCatL < nCat; iCatL++) {
          LrsList = RASString->find_exlist(rBlk, sBlk, iCatJ, iCatL);
          if (LrsList) break;
        }
        if (fCat[iCatL] < 0 or LrsList == nullptr) continue;
        nDetJ[0] = LCat[iCatJ][1];
        nDetJ[1] = LCat[iCatJ][2]/LCat[iCatJ][1];
        nDetJ[2] = LCat[iCatJ][3]/LCat[iCatJ][2];


        if (rBlk == sBlk) {// diagonal
          int rsoff = 0;
          for (auto i = 0ul; i < rBlk; i++) rsoff += nActO[i];
          ciBuilder.reassembleCatDet1RAS(LCat[iCatJ], nDetJ, rBlk, LCatAlt, nDetAlt);
          #pragma omp parallel for schedule(static, 1) default(shared) \
                private(j1, j2, j3, iNZ, l1, l, j23off, j, sign, r, s)
          for (j3 = 0; j3 < nDetAlt[2]; j3++)
          for (j2 = 0; j2 < nDetAlt[1]; j2++)
          for (j1 = 0; j1 < nDetAlt[0]; j1++) {
            j23off = j3 * LCatAlt[2] + j2 * LCatAlt[1];
            j = j1 * LCatAlt[0] + j23off;
            const int * exList_j1 = LrsList->pointerAtDet(j1);
            for (iNZ = 0; iNZ < LrsList->nNonZero(); iNZ++, exList_j1+=4) {
              UNPACK_EXCITATIONLIST_4(exList_j1, r, s, l1, sign);
              l = l1 * LCatAlt[0] + j23off; // <L|rs|J>
              LHS(l+fCat[iCatL]-CIsize, j+fCat[iCatJ]-CIsize) -= ptFock(r+rsoff, s+rsoff) * sign;
            }
          } // End openMP


        } else { // off-diagonal
          int roff = 0;
          int soff = 0;
          for (auto i = 0ul; i < rBlk; i++) roff += nActO[i];
          for (auto i = 0ul; i < sBlk; i++) soff += nActO[i];
          ciBuilder.reassembleCatDet2RAS(LCat[iCatJ], nDetJ, rBlk, sBlk,
                                LCatAlt, nDetAlt, LCatAlt2, LCat[iCatL]);
          #pragma omp parallel for schedule(static, 1) default(shared) \
                private(j1, j2, j3, iNZ, l1, l, l2, j, sign, r, s)
          for (j3 = 0; j3 < nDetAlt[2]; j3++)
          for (j1 = 0; j1 < nDetAlt[0]; j1++)
          for (j2 = 0; j2 < nDetAlt[1]; j2++) {
            j = j1*LCatAlt[0] + j2*LCatAlt[1] + j3*LCatAlt[2];
            const int * exList_j = LrsList->pointerAtDet(j1, j2);
            for (iNZ = 0; iNZ < LrsList->nNonZero(); iNZ++, exList_j+=5) {
              UNPACK_EXCITATIONLIST_5(exList_j, r, s, l1, l2, sign);
              l = l1*LCatAlt2[0] + l2*LCatAlt2[1] + j3*LCatAlt2[2];
              LHS(l+fCat[iCatL]-CIsize, j+fCat[iCatJ]-CIsize) -= ptFock(r+roff, s+soff) * sign;
            }
          } // End openMP
        }
      }
    }
    for (auto I = 0ul; I < SDsize; I++) LHS(I, I) += E0 - computeShift(I,E0);
#ifdef _DEBUG_PT2_impl
    prettyPrintSmart(std::cout,"LL PT2 full LHS -- ", LHS.pointer(), SDsize, SDsize, SDsize);
#endif

    ProgramTimer::tock("Full LHS build");

  } // PERTRUB::buildFull


  /**
   *  
   *  \brief Compute the diagonal element of LHS as preconditioner
   *          here diagLHS does not include the zero-order energy part
   *
   */
  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::buildDiag() {

    ProgramTimer::tick("LHS diagonal build");
    PT2_LOOP_INIT(); //check top for details
    size_t nCorrE = RASString->nElectron();

    // initialize diagLHS
    std::fill_n(diagLHS, SDsize, MatsT(0.));
    std::vector<size_t> j(3, 0);
    std::vector<size_t> nEras(3, 0);
    std::vector<size_t> nDras(3, 0);
    std::vector<size_t> elecPos(nCorrE, 0);
    int_matrix empty = {{}};

    size_t iCatJ = 0;
    for (auto iEJ = 0ul; iEJ <= mxElec; iEJ++)
    for (auto iHJ = 0ul; iHJ <= mxHole; iHJ++, iCatJ++) {

      if (iCatJ == 0) continue;
      nEras[0] = nActO[0] - iHJ;
      nEras[2] = iEJ;
      if (nEras[0] + nEras[2] > nCorrE or nActO[1] + nEras[0] + nEras[2] < nCorrE)
        continue;

      nEras[1] = nCorrE - nEras[0] - nEras[2];
      nDras[0] = LCat[iCatJ][1];
      nDras[1] = LCat[iCatJ][2]/LCat[iCatJ][1];
      nDras[2] = LCat[iCatJ][3]/LCat[iCatJ][2];

      MatsT* diagLHS_J = diagLHS + fCat[iCatJ] - CIsize;

      std::vector<int_matrix> addr_array;
      std::vector<int_matrix> de_addr_array;
      for (auto i = 0ul; i < 3; i++){
        if (nActO[i] == 0) {
          addr_array.push_back(empty);
          de_addr_array.push_back(empty);
        } else {
        addr_array.push_back(RASString->buildAddressingArray(nEras[i], nActO[i]));
        de_addr_array.push_back(RASString->buildDeAddressingArray(addr_array[i]));
        }
      }

      for (j[2] = 0; j[2] < nDras[2]; j[2]++)
      for (j[1] = 0; j[1] < nDras[1]; j[1]++)
      for (j[0] = 0; j[0] < nDras[0]; j[0]++) {
        size_t j_addr = j[0] * LCat[iCatJ][0] + j[1] * LCat[iCatJ][1] + j[2] * LCat[iCatJ][2];
        size_t orb_offset = 0;
        size_t e_offset = 0;
        for (auto i = 0ul; i < 3; i++){
          if (nActO[i] == 0) continue;
          std::vector<size_t> eP = RASString->address2DetString(j[i], addr_array[i],
                                            de_addr_array[i]);
          for (auto pos = 0ul; pos < nEras[i]; pos++)
            elecPos[pos+e_offset] = eP[pos] + orb_offset;
          orb_offset += nActO[i];
          e_offset += nEras[i];
        }

        // Update ptFock part
        for (auto pos_i = 0ul; pos_i < nCorrE; pos_i++) {
          diagLHS_J[j_addr] -= ptFock(elecPos[pos_i], elecPos[pos_i]);
        }
      }
    }
#ifdef _DEBUG_PT2_impl
    if (nC != 1) {
      auto &LHS = *LHS_;
      for (auto i = 0ul; i < SDsize; i++) {
        if (std::abs(LHS(i, i)-diagLHS[i]) >= 1e-12)
          std::cout<<"diag wrong"<<std::endl; 
      }
    }
#endif

    ProgramTimer::tock("LHS diagonal build");

  } //PERTURB::buildDiag


  /**
   *  
   *  \brief RHS: H * P_0
   *    H is the second-quantized Hamiltonian
   *    C is the CI coefficients of the reference state spanned over the new space
   *    Will use RASCI::buildSigma function
   */
  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::buildRHS(MatsT * RHS, size_t i) {

    ProgramTimer::tick("RHS build");
    size_t NDet = this->NDet;
    std::fill_n(this->CIVecs[i], NDet, MatsT(0.));
    std::copy_n(refMCwfn->CIVecs[SoI_[i]], CIsize, this->CIVecs[i]);

    MatsT * scr = CQMemManager::get().malloc<MatsT>(NDet);

    this->ciBuilder->buildSigma(dynamic_cast<MCWaveFunction<MatsT,IntsT>&>(*this), 1,
                                    this->CIVecs[i], scr);

    for (auto i = CIsize; i < NDet; i++) RHS[i-CIsize] = scr[i];
#ifdef _DEBUG_PT2_impl
    prettyPrintSmart(std::cout,"LL PT2 RHS --", RHS, SDsize, 1, SDsize);
#endif

    CQMemManager::get().free(scr);

    ProgramTimer::tock("RHS build");

  } // PERTRUB::buildRHS

  /**
   * 
   *  \brief Compute HV to be used in CHV
   *    using RASCI sigma build:
   *    ptV: 1st perturbed vector for state N; S is external determinant.
   *
   */
  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::computeHV(MatsT * HV, MatsT * C) {

    ProgramTimer::tick("compute HV");
    size_t NDet = this->NDet;
    MatsT * P = CQMemManager::get().malloc<MatsT>(NDet);

    for (auto i = 0ul; i < CIsize; i++) P[i] = MatsT(0.);
    for (auto i = CIsize; i < NDet; i++) P[i] = C[i];
    this->ciBuilder->buildSigma(dynamic_cast<MCWaveFunction<MatsT,IntsT>&>(*this), 1, P, HV);

    CQMemManager::get().free(P);
    ProgramTimer::tock("compute HV");

  } // PERTRUB::computeHV


  /**
   * 
   *  \brief Compute the second-order perturbed energy
   *        CHV(M,N) = \sum_IS (P0_MI)'* H_IS * ptV_SN
   *        P0: CI vector for state M; I is cas determinant.
   *        ptV: 1st perturbed vector for state N; S is external determinant.
   *
   */
  template <typename MatsT, typename IntsT>
  MatsT PERTURB<MatsT,IntsT>::computeCHV(MatsT * C, MatsT * HV) {

    // CHV = (P0)'* HV
    MatsT CHV = blas::dot(CIsize,C,1,HV,1);
    //std::cout<<"E2 raw: "<< std::setprecision(12)<<CHV<<std::endl;

    return CHV;

  } // PERTURB::computeCHV

  /**
   *  
   *  \brief Compute the imaginary shift in the i-th LHS diagonal term
   *         shift_i = levelShift + (imaginaryShift)^2 / (diagLHS[i] - E_0[state])
   */
  template <typename MatsT, typename IntsT>
  MatsT PERTURB<MatsT,IntsT>::computeShift(size_t i, MatsT E0) {

    MatsT shift = MatsT(PTopts.levelShift);
    shift += MatsT(std::pow(PTopts.imaginaryShift, 2)) / (diagLHS[i] - E0);
    //std::cout<<"i: "<<i<<", imaginary shift: "<<shift<<std::endl;
    return shift;

  } // PERTURB::computeShift

  /**
   *  
   *  \brief Compute shift correction based on Hylleraas functional
   *         ShiftCorrection = <Psi_1|shift|Psi_1>
   *            levelShift => computeCC
   *            imaginaryShift => sum_ij C_i^H <D_i|shift|D_j>C_j
   */
  template <typename MatsT, typename IntsT>
  MatsT PERTURB<MatsT,IntsT>::computeShiftCorrection(size_t state) {

    MatsT *C_pt = this->CIVecs[state] + CIsize;
    MatsT correction = MatsT(0.);

    if (PTopts.imaginaryShift) {
      MatsT *scr = CQMemManager::get().malloc<MatsT>(this->SDsize);
      // scr_i = C_i * shift;
      for (auto i = 0ul; i < SDsize; i++) {
        scr[i] = C_pt[i] * computeShift(i, E0_[state]);
      }
      correction +=  blas::dot(SDsize,C_pt,1,scr,1);
      CQMemManager::get().free(scr);
    }
    else if (PTopts.levelShift)
      correction += MatsT(PTopts.levelShift) * blas::dot(SDsize,C_pt,1,C_pt,1);

    return correction;

  } // PERTURB::computeShiftCorrection



  /**
   *  
   *  \brief forming vector AV used in GMRES
   *      trial vector: V
   *      LHS: A
   *      AV[I] = sum_J <I|F-E_0|J> * V[J].
   *
   */
  template <typename MatsT, typename IntsT>
  void PERTURB<MatsT,IntsT>::buildAV(MatsT * V, MatsT * AV, MatsT E0) {

    //std::cout<<"LL PT2 iterative AV build."<<std::endl;
    ProgramTimer::tick("LHS sigma build");
    PT2_LOOP_INIT(); // see top for details

    // Initialize AV
    std::fill_n(AV, SDsize, MatsT(0.));
    std::shared_ptr<const ExcitationList> LrsList;
    size_t iCatJ, iCatL;
    std::vector<size_t> LCatAlt(3, 0);
    std::vector<size_t> LCatAlt2(3, 0);
    std::vector<size_t> nDetAlt(3, 0);
    std::vector<size_t> nDetJ(3, 0);
    size_t r, s, l, l1, l2, j, j1, j2, j3, iNZ, j23off;
    double sign;

    for (auto rBlk = 0ul; rBlk < 3; rBlk++)
    for (auto sBlk = 0ul; sBlk < 3; sBlk++) {
      if (nActO[rBlk] == 0 or nActO[sBlk] == 0) continue;
      iCatJ = 0;
      for (auto iEJ = 0ul; iEJ <= mxElec; iEJ++)
      for (auto iHJ = 0ul; iHJ <= mxHole; iHJ++, iCatJ++) {      
        if (iCatJ == 0 or LCat[iCatJ][3] == 0) continue;
        for (iCatL = 1; iCatL < nCat; iCatL++) {
          LrsList = RASString->find_exlist(rBlk, sBlk, iCatJ, iCatL);
          if (LrsList) break;
        }
        if (fCat[iCatL] < 0 or LrsList == nullptr) continue;
        nDetJ[0] = LCat[iCatJ][1];
        nDetJ[1] = LCat[iCatJ][2]/LCat[iCatJ][1];
        nDetJ[2] = LCat[iCatJ][3]/LCat[iCatJ][2];

        if (rBlk == sBlk) {// diagonal
          int rsoff = 0;
          for (auto i = 0ul; i < rBlk; i++) rsoff += nActO[i];
          ciBuilder.reassembleCatDet1RAS(LCat[iCatJ], nDetJ, rBlk, LCatAlt, nDetAlt);
          #pragma omp parallel for schedule(static, 1) default(shared) \
            private(j1, j2, j3, iNZ, l1, l, j23off, j, sign, r, s)
          for (j3 = 0; j3 < nDetAlt[2]; j3++)
          for (j2 = 0; j2 < nDetAlt[1]; j2++)
          for (j1 = 0; j1 < nDetAlt[0]; j1++) {
            j23off = j3 * LCatAlt[2] + j2 * LCatAlt[1];
            j = j1 * LCatAlt[0] + j23off;
            const int * exList_j1 = LrsList->pointerAtDet(j1);
            for (iNZ = 0; iNZ < LrsList->nNonZero(); iNZ++, exList_j1+=4) {
              UNPACK_EXCITATIONLIST_4(exList_j1, r, s, l1, sign);
              l = l1 * LCatAlt[0] + j23off; // <L|rs|J>
              AV[l+fCat[iCatL]-CIsize] -= ptFock(r+rsoff, s+rsoff) * sign * V[j+fCat[iCatJ]-CIsize];
            }
          } // End openMP
        } else { // off-diagonal
          int roff = 0;
          int soff = 0;
          for (auto i = 0ul; i < rBlk; i++) roff += nActO[i];
          for (auto i = 0ul; i < sBlk; i++) soff += nActO[i];
          ciBuilder.reassembleCatDet2RAS(LCat[iCatJ], nDetJ, rBlk, sBlk,
                                LCatAlt, nDetAlt, LCatAlt2, LCat[iCatL]);
          #pragma omp parallel for schedule(static, 1) default(shared) \
            private(j1, j2, j3, iNZ, l1, l, l2, j, sign, r, s)
          for (j3 = 0; j3 < nDetAlt[2]; j3++)
          for (j1 = 0; j1 < nDetAlt[0]; j1++)
          for (j2 = 0; j2 < nDetAlt[1]; j2++) {
            j = j1*LCatAlt[0] + j2*LCatAlt[1] + j3*LCatAlt[2];
            const int * exList_j = LrsList->pointerAtDet(j1, j2);
            for (iNZ = 0; iNZ < LrsList->nNonZero(); iNZ++, exList_j+=5) {
              UNPACK_EXCITATIONLIST_5(exList_j, r, s, l1, l2, sign);
              l = l1*LCatAlt2[0] + l2*LCatAlt2[1] + j3*LCatAlt2[2];
              AV[l+fCat[iCatL]-CIsize] -= ptFock(r+roff, s+soff) * sign * V[j+fCat[iCatJ]-CIsize];
            }
          } // End openMP
        }
      }
    }
    for (auto i = 0; i < SDsize; i++) AV[i] += (E0 - computeShift(i,E0)) * V[i];
    ProgramTimer::tock("LHS sigma build");

  } // PERTURB::buildAV


} // namespace ChronusQ

