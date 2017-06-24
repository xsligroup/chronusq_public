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

#include <cibuilder/rasci.hpp>
#include <detstringmanager.hpp>
#include <cqlinalg.hpp>
#include <cqlinalg/matfunc.hpp>

#include <util/timer.hpp>
#include <util/matout.hpp>
#include <util/threads.hpp>

//#define _DEBUG_CIBUILDER_RASCI_IMPL

namespace ChronusQ {

  /**
   * \brief Re-assemble LCatAlt and nDetAlt for 1 RAS space
   */
  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::reassembleCatDet1RAS(std::vector<size_t> & LCat,
        std::vector<size_t> & nDet, size_t e, std::vector<size_t> & LCatAlt,
        std::vector<size_t> & nDetAlt) {

    LCatAlt[0] = LCat[e];
    nDetAlt[0] = nDet[e];
    size_t ii = 0;
    for (auto i = 0ul; i < 3; i++)
      if (i != e) {
        ii += 1;
        LCatAlt[ii] = LCat[i];
        nDetAlt[ii] = nDet[i];
      }
  } // RASCI::reassembleCatDet1RAS

  /**
   * \brief Re-assemble LCatAlt and nDetAlt for 2 RAS spaces
   */
  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::reassembleCatDet2RAS(std::vector<size_t> & LCat,
        std::vector<size_t> & nDet, size_t e1, size_t e2, std::vector<size_t> & LCatAlt,
        std::vector<size_t> & nDetAlt, std::vector<size_t> & Jrs,
        std::vector<size_t> & LCat2) {

    if (LCat2.size() > 0) {
      Jrs[0] = LCat2[e1];
      Jrs[1] = LCat2[e2];
    } else {
      Jrs[0] = e1;
      Jrs[1] = e2;
    }
    LCatAlt[0] = LCat[e1];
    LCatAlt[1] = LCat[e2];
    nDetAlt[0] = nDet[e1];
    nDetAlt[1] = nDet[e2];
    for (auto i = 0ul; i < 3; i++)
      if (i != e1 and i != e2) {
        if (LCat2.size() > 0) Jrs[2] = LCat2[i];
        else Jrs[2] = i;
        LCatAlt[2] = LCat[i];
        nDetAlt[2] = nDet[i];
      }
  } // RASCI::reassembleCatDet2RAS

  /**
   * \brief Build intermediates D blockwise for GHF RAS Sigma -- Diagonal
   *        using blockwise excitation list as:
   *        D[rs][J] = sum_L <J|E_rs|L> C_L
   *
   */
  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::SigmaD_diag(MatsT * C, MatsT * DBlk, size_t nNZof,
        size_t nNZB, std::shared_ptr<const ExcitationList> exList, size_t rsB,
        std::vector<size_t> & LCat,
        std::vector<size_t> & nDet) {

    size_t l1, l, j23off, j, p, q, j1, j2, j3, iNZ;
    double sign;

    // re-assemble LCat and nDet into LCatAlt and nDetAlt
    std::vector<size_t> LCatAlt(3, 0);
    std::vector<size_t> nDetAlt(3, 0);
    reassembleCatDet1RAS(LCat, nDet, rsB, LCatAlt, nDetAlt);

    #pragma omp parallel for schedule(static, 1) default(shared) private(j1, j2, j3, iNZ, l1, l, j23off, j, sign, p, q)
    for (j3 = 0; j3 < nDetAlt[2]; j3++)
    for (j2 = 0; j2 < nDetAlt[1]; j2++)
    for (j1 = 0; j1 < nDetAlt[0]; j1++) {
      j23off = j3 * LCatAlt[2] + j2 * LCatAlt[1];
      j = j1 * LCatAlt[0] + j23off;
      const int * exList_j1 = exList->pointerAtDet(j1) + nNZof * 4;
      for (iNZ = 0; iNZ < nNZB; iNZ++, exList_j1+=4) {
        UNPACK_EXCITATIONLIST_4(exList_j1, p, q, l1, sign);
        l = l1 * LCatAlt[0] + j23off;
        DBlk[iNZ+j*nNZB] += sign * C[l];
      }
    }
    // End openMP

#ifdef _DEBUG_CIBUILDER_RASCI_IMPL
    prettyPrintSmart(std::cout,"LL RASCI Sigma D diag --", DBlk, nNZB, LCat[3], nNZB);
#endif

  } //RASCI::SigmaD_diag

  /**
   * \brief Build intermediates D blockwise for GHF RAS Sigma -- Off-diagonal
   *        using blockwise excitation list as:
   *        D[rs][J] = sum_L <J|E_rs|L> C_L
   *
   */
  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::SigmaD_offdiag(MatsT * C, MatsT * DBlk, size_t nNZof,
        size_t nNZB, std::shared_ptr<const ExcitationList> exList, size_t rB, size_t sB,
        std::vector<size_t> & LCat1, std::vector<size_t> & LCat2,
        std::vector<size_t> & nDet) {

    // re-assemble LCatAlt and NDAlt
    std::vector<size_t> LCatAlt1(3, 0);
    std::vector<size_t> LCatAlt2(3, 0);
    std::vector<size_t> nDetAlt(3, 0);
    reassembleCatDet2RAS(LCat1, nDet, sB, rB, LCatAlt1, nDetAlt, LCatAlt2, LCat2);

    size_t j, l, l1, l2, p, q, j1, j2, j3, iNZ;
    double sign;

    #pragma omp parallel for schedule(static, 1) default(shared) private(j1, j2, j3, iNZ, l1, l, l2, j, sign, p, q)
    for (j3 = 0; j3 < nDetAlt[2]; j3++)
    for (j1 = 0; j1 < nDetAlt[0]; j1++)
    for (j2 = 0; j2 < nDetAlt[1]; j2++) {
      j = j1*LCatAlt1[0] + j2*LCatAlt1[1] + j3*LCatAlt1[2];
      const int * exList_j = exList->pointerAtDet(j1, j2) + nNZof * 5;
      for (iNZ = 0; iNZ < nNZB; iNZ++, exList_j+=5) {
        UNPACK_EXCITATIONLIST_5(exList_j, p, q, l1, l2, sign);
        l = l1*LCatAlt2[0] + l2*LCatAlt2[1] + j3*LCatAlt2[2];
        DBlk[iNZ+j*nNZB] += sign * C[l];
      }
    }
    // End openMP

#ifdef _DEBUG_CIBUILDER_RASCI_IMPL
    prettyPrintSmart(std::cout,"LL RASCI Sigma D off-diag --", DBlk, nNZB, LCat1[3], nNZB);
#endif

  } // RASCI::SigmaD_offdiag

  /**
   * \brief Update Sigma 1e for GHF RAS -- diagonal
   *        Using block excitationlist as <J|E_sr|L>
   *        Sigma_1e[K] = \sum_pq hCorP_pq  * D[pq, K];
   *
   */
  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::Sigma1e_diag(MatsT * Sigma, MatsT * Sscr, MatsT * DBlk,
        MatsT * Dscr, OnePInts<MatsT> & hCoreP, MatsT * hCPscr, size_t nNZof,
        size_t nNZB, std::shared_ptr<const ExcitationList> exList, size_t rsB,
        std::vector<size_t> & LCat,
        std::vector<size_t> & nDet, int rsoff) {

    // re-assemble LCatAlt and NDAlt
    std::vector<size_t> LCatAlt(3, 0);
    std::vector<size_t> nDetAlt(3, 0);
    reassembleCatDet1RAS(LCat, nDet, rsB, LCatAlt, nDetAlt);

    size_t n23 = LCat[3] / nDetAlt[0];
    size_t j1off, j23, j13off, j, r, s, l1;
    double sign;

    #pragma omp for schedule(static, 1)
    for (auto j1 = 0ul; j1 < nDetAlt[0]; j1++) {
      j1off = j1 * LCatAlt[0];
      // Gathering D, implicitly J_rs = j1;
      j23 = 0;
      for (auto j3 = 0ul; j3 < nDetAlt[2]; j3++) {
        j13off = j1off + j3 * LCatAlt[2];
        for (auto j2 = 0ul; j2 < nDetAlt[1]; j2++, j23++) {
          j = j13off + j2 * LCatAlt[1];
          std::copy_n(DBlk+j*nNZB, nNZB, Dscr+j23*nNZB);
        }
      }
      // Gathering HCorePrime
      const int * exList_j1 = exList->pointerAtDet(j1) + nNZof * 4;
      for (auto iNZ = 0ul; iNZ < nNZB; iNZ++, exList_j1+=4) {
        UNPACK_EXCITATIONLIST_4(exList_j1, s, r, l1, sign);
        s += rsoff;
        r += rsoff;
        hCPscr[iNZ] = hCoreP(r, s);
      }
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
        1,n23,nNZB,MatsT(1.),hCPscr,1,Dscr,nNZB,
        MatsT(0.),Sscr,1);
      // Scattering back to Sigma
      j23 = 0;
      for (auto j3 = 0ul; j3 < nDetAlt[2]; j3++) {
        j13off = j1off + j3 * LCatAlt[2];
        for (auto j2 = 0ul; j2 < nDetAlt[1]; j2++, j23++) {
          j = j13off + j2 * LCatAlt[1];
          Sigma[j] += Sscr[j23];
        }
      }
    }
    // End omp for

  } //RASCI::Sigma1e_diag

  /**
   * \brief Update Sigma 1e for GHF RAS -- off-diagonal
   *        Using block excitationlist as <J|E_sr|L>
   *        Sigma_1e[K] = \sum_pq hCorP_pq  * D[pq, K];
   *
   */
  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::Sigma1e_offdiag(MatsT * Sigma, MatsT * Sscr, MatsT * DBlk,
        MatsT * Dscr, OnePInts<MatsT> & hCoreP, MatsT * hCPscr, size_t nNZof,
        size_t nNZB, std::shared_ptr<const ExcitationList> exList, size_t rB, size_t sB,
        std::vector<size_t> & LCat, std::vector<size_t> & nDet,
        int roff, int soff) {

    // locate the RAS space involved in excitation
    size_t exBlk1 = std::min(rB, sB);
    size_t exBlk2 = std::max(rB, sB);

    // re-assemble LCatAlt and nDetAlt
    std::vector<size_t> Jrs(3, 0);
    std::vector<size_t> LCatAlt(3, 0);
    std::vector<size_t> nDetAlt(3, 0);
    std::vector<size_t> LCat2 = {};
    reassembleCatDet2RAS(LCat, nDet, exBlk1, exBlk2, LCatAlt, nDetAlt, Jrs, LCat2);

    std::vector<size_t> Jras(3, 0);
    size_t j, j12off, s, r, l1, l2;
    double sign;

    #pragma omp for schedule(static, 1)
    for (auto j1 = 0ul; j1 < nDetAlt[0]; j1++)
    for (auto j2 = 0ul; j2 < nDetAlt[1]; j2++) {
      j12off = j2 * LCatAlt[1] + j1 * LCatAlt[0];
      Jras[Jrs[0]] = j1;
      Jras[Jrs[1]] = j2;
      // Gathering D
      for (auto j3 = 0ul; j3 < nDetAlt[2]; j3++) {
        j = j12off + j3 * LCatAlt[2];
        std::copy_n(DBlk+j*nNZB, nNZB, Dscr+j3*nNZB);
      }
      // Gathering HCorePrime
      const int * exList_j = exList->pointerAtDet(Jras[sB], Jras[rB]) + nNZof * 5;
      for (auto iNZ = 0ul; iNZ < nNZB; iNZ++, exList_j+=5) {
        UNPACK_EXCITATIONLIST_5(exList_j, s, r, l1, l2, sign);
        s += soff;
        r += roff;
        hCPscr[iNZ] = hCoreP(r, s);
      }
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
        1,nDetAlt[2],nNZB,MatsT(1.),hCPscr,1,
        Dscr,nNZB,MatsT(0.),Sscr,1);
      // Scattering back to Sigma
      for (auto j3 = 0ul; j3 < nDetAlt[2]; j3++) {
        j = j12off + j3 * LCatAlt[2];
        Sigma[j] += Sscr[j3];
      }
    }

  } // RASCI::Sigma1e_offdiag

  /**
   * \brief Build intermediates G blockwise for GHF RAS Sigma
   *        diagonal rs -- diagonal pq
   *        Using block excitationlist as <J|E_sr|L>
   *        G(pq, J) = (pq|rs) * D(rs, J)
   *        
   */
  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::SigmaG_dd(MatsT * GBlk, MatsT * Gscr, MatsT * DBlk,
        MatsT * Dscr, MatsT * ERIscr, InCore4indexTPI<MatsT> & moERI,
        size_t nNZrsof, size_t nNZrsB, size_t nNZpqof, size_t nNZpqB,
        std::shared_ptr<const ExcitationList> pqList, std::shared_ptr<const ExcitationList> rsList,
        size_t rsB, size_t pqB,
        std::vector<size_t> & LCat, std::vector<size_t> & nDet, int rsoff, int pqoff) {

    // determine excitation cases:
    //   1. one RAS space is involved in excitation
    //   2. two RAS spaces are involved in excitation

    size_t exBlk1, exBlk2, j, p, q, r, s;
    double sign;
    std::vector<size_t> Jrs(3, 0);
    std::vector<size_t> LCatAlt(3, 0);
    std::vector<size_t> nDetAlt(3, 0);

    // Case 1: one RAS space is involved in excitation
    if (rsB == pqB) {
      // locate the RAS space involved in excitation
      exBlk1 = rsB;
      // re-assemble LCatAlt and nDetAlt
      reassembleCatDet1RAS(LCat, nDet, exBlk1, LCatAlt, nDetAlt);

      size_t n23 = LCat[3]/nDetAlt[0];
      size_t j1off, j23, j13off,l1;

      #pragma omp for schedule(static, 1)
      for (auto j1 = 0ul; j1 < nDetAlt[0]; j1++) {
        // Gathering D, implicitly jrs = jpq = j1
        j1off = j1 * LCatAlt[0];
        j23 = 0;
        for (auto j3 = 0ul; j3 < nDetAlt[2]; j3++) {
          j13off = j1off + j3 * LCatAlt[2];
          for (auto j2 = 0ul; j2 < nDetAlt[1]; j2++, j23++) {
            j = j13off + j2 * LCatAlt[1];
            std::copy_n(DBlk+j*nNZrsB, nNZrsB, Dscr+j23*nNZrsB);
          }
        }
        // Gathering moERI
        const int * exList_pq = pqList->pointerAtDet(j1) + nNZpqof * 4;
        for (auto iNZpq = 0ul; iNZpq < nNZpqB; iNZpq++, exList_pq+=4) {
          UNPACK_EXCITATIONLIST_4(exList_pq, p, q, l1, sign);
          p += pqoff;
          q += pqoff;
          const int * exList_rs = rsList->pointerAtDet(j1) + nNZrsof * 4;
          for (auto iNZrs = 0ul; iNZrs < nNZrsB; iNZrs++, exList_rs+=4) {
            UNPACK_EXCITATIONLIST_4(exList_rs, s, r, l1, sign);
            s += rsoff;
            r += rsoff;
            ERIscr[iNZrs+iNZpq*nNZrsB] = moERI(r, s, p, q);
          }
        }
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
            nNZpqB,n23,nNZrsB,MatsT(1.),ERIscr,nNZrsB,
            Dscr,nNZrsB,MatsT(0.),Gscr,nNZpqB);
        // Scattering back to GBlk
        j23 = 0;
        for (auto j3 = 0ul; j3 < nDetAlt[2]; j3++) {
          j13off = j1off + j3 * LCatAlt[2];
          for (auto j2 = 0ul; j2 < nDetAlt[1]; j2++, j23++) {
            j = j13off + j2 * LCatAlt[1];
            MatAdd('N','N',nNZpqB,1,MatsT(1.),Gscr+j23*nNZpqB,nNZpqB,
                MatsT(1.),GBlk+j*nNZpqB,nNZpqB,GBlk+j*nNZpqB,nNZpqB);
          }
        }
      }
    // Case 2: two RAS spaces are involved in excitation
    //         implicitly pBlk < rBlk
    } else {
      // locate the RAS space involved in excitation
      exBlk1 = std::min(rsB, pqB);
      exBlk2 = std::max(rsB, pqB);

      // re-assemble LCatAlt and nDetAlt
      std::vector<size_t> LCat2 = {};
      reassembleCatDet2RAS(LCat, nDet, exBlk1, exBlk2, LCatAlt, nDetAlt, Jrs, LCat2);

      std::vector<size_t> Jras(3, 0);
      size_t j12off, l1;

      #pragma omp for schedule(static, 1)
      for (auto j1 = 0ul; j1 < nDetAlt[0]; j1++)
      for (auto j2 = 0ul; j2 < nDetAlt[1]; j2++) {
        j12off = j2 * LCatAlt[1] + j1 * LCatAlt[0];
        Jras[Jrs[0]] = j1;
        Jras[Jrs[1]] = j2;
        // Gathering D
        for (auto j3 = 0ul; j3 < nDetAlt[2]; j3++) {
          j = j12off + j3 * LCatAlt[2];
          std::copy_n(DBlk+j*nNZrsB, nNZrsB, Dscr+j3*nNZrsB);
        }
        // Gathering moERI, implicitly pBlk < rBlk
        const int * exList_pq = pqList->pointerAtDet(Jras[pqB]) + nNZpqof * 4;
        for (auto iNZpq = 0ul; iNZpq < nNZpqB; iNZpq++, exList_pq+=4) {
          UNPACK_EXCITATIONLIST_4(exList_pq, p, q, l1, sign);
          p += pqoff;
          q += pqoff;
          const int * exList_rs = rsList->pointerAtDet(Jras[rsB]) + nNZrsof * 4;
          for (auto iNZrs = 0ul; iNZrs < nNZrsB; iNZrs++, exList_rs+=4) {
            UNPACK_EXCITATIONLIST_4(exList_rs, s, r, l1, sign);
            s += rsoff;
            r += rsoff;
            ERIscr[iNZrs+iNZpq*nNZrsB] = moERI(r, s, p, q) - moERI(p, s, r, q);
          }
        }

        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
            nNZpqB,nDetAlt[2],nNZrsB,MatsT(1.),ERIscr,nNZrsB,
            Dscr,nNZrsB,MatsT(0.),Gscr,nNZpqB);
        // Scattering back to GBlk
        for (auto j3 = 0ul; j3 < nDetAlt[2]; j3++) {
          j = j12off + j3 * LCatAlt[2];
          MatAdd('N','N',nNZpqB,1,MatsT(1.),Gscr+j3*nNZpqB,nNZpqB,
                MatsT(1.),GBlk+j*nNZpqB,nNZpqB,GBlk+j*nNZpqB,nNZpqB);
        }
      }
    }
#ifdef _DEBUG_CIBUILDER_RASCI_IMPL
    prettyPrintSmart(std::cout,"LL RASCI Sigma G diag-diag --", GBlk, nNZpqB, LCat[3], nNZpqB);
#endif
  } // RASCI::SigmaG_dd

  /**
   * \brief Build intermediates G blockwise for GHF RAS Sigma
   *        diagonal rs -- off-diagonal pq
   *        Using block excitationlist as <J|E_sr|L>
   *        G(pq, J) = (pq|rs) * D(rs, J)
   *
   */
  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::SigmaG_do(MatsT * GBlk, MatsT * Gscr, MatsT * DBlk,
        MatsT * Dscr, MatsT * ERIscr, InCore4indexTPI<MatsT> & moERI, size_t nNZrsof,
        size_t nNZrsB, size_t nNZpqof, size_t nNZpqB,
        std::shared_ptr<const ExcitationList> pqList,
        std::shared_ptr<const ExcitationList> rsList, size_t rsB, size_t pB, size_t qB,
        std::vector<size_t> & LCat, std::vector<size_t> & nDet,
        int rsoff, int poff, int qoff) {

    std::vector<size_t> exRAS(3, 0);
    std::vector<size_t> Jrs(3, 0);
    std::vector<size_t> LCatAlt(3, 0);
    std::vector<size_t> nDetAlt(3, 0);
    std::vector<size_t> Jras(3, 0);
    size_t exBlk1, exBlk2, j, p, q, r, s, l1, l2;
    double sign;

    // determine excitation cases:
    //  2. two RAS spaces are involved in excitation
    //  3. three RAS spaces are involved in excitation
    exRAS[rsB] = 1;
    exRAS[pB] = 1;
    exRAS[qB] = 1;
    int eCase = std::accumulate(exRAS.begin(), exRAS.end(), 0);
    // Case 2: two RAS spaces are involved in excitation
    if (eCase == 2) {
      // locate the RAS space involved in excitation
      exBlk1 = std::min({rsB, pB, qB});
      exBlk2 = std::max({rsB, pB, qB});
      std::vector<size_t> LCat2 = {};
      reassembleCatDet2RAS(LCat, nDet, exBlk1, exBlk2, LCatAlt, nDetAlt, Jrs, LCat2);

      #pragma omp for schedule(static, 1)
      for (auto j1 = 0ul; j1 < nDetAlt[0]; j1++)
      for (auto j2 = 0ul; j2 < nDetAlt[1]; j2++) {
        size_t j12off = j2 * LCatAlt[1] + j1 * LCatAlt[0];
        Jras[Jrs[0]] = j1;
        Jras[Jrs[1]] = j2;
        // Gathering D
        for (auto j3 = 0ul; j3 < nDetAlt[2]; j3++) {
          j = j12off + j3 * LCatAlt[2];
          std::copy_n(DBlk+j*nNZrsB, nNZrsB, Dscr+j3*nNZrsB);
        }
        // Scattering moERI
        const int * exList_pq = pqList->pointerAtDet(Jras[pB], Jras[qB]) + nNZpqof * 5;
        for (auto iNZpq = 0ul; iNZpq < nNZpqB; iNZpq++, exList_pq+=5) {
          UNPACK_EXCITATIONLIST_5(exList_pq, p, q, l1, l2, sign);
          p += poff;
          q += qoff;
          const int * exList_rs = rsList->pointerAtDet(Jras[rsB]) + nNZrsof * 4;
          for (auto iNZrs = 0ul; iNZrs < nNZrsB; iNZrs++, exList_rs+=4) {
            UNPACK_EXCITATIONLIST_4(exList_rs, s, r, l1, sign);
            s += rsoff;
            r += rsoff;
            ERIscr[iNZrs+iNZpq*nNZrsB] = moERI(r, s, p, q);
          }
        }
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
            nNZpqB,nDetAlt[2],nNZrsB,MatsT(1.),ERIscr,nNZrsB,
            Dscr,nNZrsB,MatsT(0.),Gscr,nNZpqB);
        // Scattering back to GBlk
        for (auto j3 = 0ul; j3 < nDetAlt[2]; j3++) {
          j = j12off + j3 * LCatAlt[2];
          MatAdd('N','N',nNZpqB,1,MatsT(1.),Gscr+j3*nNZpqB,nNZpqB,
                MatsT(1.),GBlk+j*nNZpqB,nNZpqB,GBlk+j*nNZpqB,nNZpqB);
        }
      }
    // Case 3: three RAS spaces are involved in excitation
    //         implicitly that pBlk != rBlk, qBlk != sBlk
    } else if (eCase == 3) {
      // LCatAlt == LCat, nDetAlt == nDet
      #pragma omp for schedule(static, 1)
      for (auto j3 = 0ul; j3 < nDet[2]; j3++)
      for (auto j2 = 0ul; j2 < nDet[1]; j2++)
      for (auto j1 = 0ul; j1 < nDet[0]; j1++) {
        Jras[0] = j1;
        Jras[1] = j2;
        Jras[2] = j3;
        j = j1 * LCat[0] + j2 * LCat[1] + j3 * LCat[2];
        // assemble moERI
        const int * exList_pq = pqList->pointerAtDet(Jras[pB], Jras[qB]) + nNZpqof * 5;
        for (auto iNZpq = 0ul; iNZpq < nNZpqB; iNZpq++, exList_pq+=5) {
          UNPACK_EXCITATIONLIST_5(exList_pq, p, q, l1, l2, sign);
          p += poff;
          q += qoff;
          const int * exList_rs = rsList->pointerAtDet(Jras[rsB]) + nNZrsof * 4;
          for (auto iNZrs = 0ul; iNZrs < nNZrsB; iNZrs++, exList_rs+=4) {
            UNPACK_EXCITATIONLIST_4(exList_rs, s, r, l1, sign);
            s += rsoff;
            r += rsoff;
            ERIscr[iNZrs+iNZpq*nNZrsB] = moERI(r, s, p, q) - moERI(p, s, r, q);
          }
        }
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
            nNZpqB,1,nNZrsB,MatsT(1.),ERIscr,nNZrsB,
            DBlk+j*nNZrsB,nNZrsB,MatsT(1.),GBlk+j*nNZpqB,nNZpqB);
      }
    }
#ifdef _DEBUG_CIBUILDER_RASCI_IMPL
    prettyPrintSmart(std::cout,"LL RASCI Sigma G diag-off --", GBlk, nNZpqB, LCat[3], nNZpqB);
#endif
  } // RASCI::SigmaG_do

  /**
   * \brief Build intermediates G blockwise for GHF RAS Sigma
   *        off-diagonal rs -- diagonal pq
   *        Using block excitationlist as <J|E_sr|L>
   *        G(pq, J) = (pq|rs) * D(rs, J)
   *
   */
  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::SigmaG_od(MatsT * GBlk, MatsT * Gscr, MatsT * DBlk,
        MatsT * Dscr, MatsT * ERIscr, InCore4indexTPI<MatsT> & moERI, size_t nNZrsof,
        size_t nNZrsB, size_t nNZpqof, size_t nNZpqB,
        std::shared_ptr<const ExcitationList> pqList,
        std::shared_ptr<const ExcitationList> rsList, size_t rB, size_t sB, size_t pqB,
        std::vector<size_t> & LCat, std::vector<size_t> & nDet,
        int roff, int soff, int pqoff) {

    std::vector<size_t> exRAS(3, 0);
    std::vector<size_t> Jrs(3, 0);
    std::vector<size_t> LCatAlt(3, 0);
    std::vector<size_t> nDetAlt(3, 0);
    std::vector<size_t> Jras(3, 0);
    size_t exBlk1, exBlk2, j, p, q, r, s, l1, l2;
    double sign;

    // determine excitation cases:
    //  2. two RAS spaces are involved in excitation
    //  3. three RAS spaces are involved in excitation
    exRAS[rB] = 1;
    exRAS[sB] = 1;
    exRAS[pqB] = 1;
    int eCase = std::accumulate(exRAS.begin(), exRAS.end(), 0);
    // Case 2: two RAS spaces are involved in excitation
    if (eCase == 2) {
      // locate the RAS space involved in excitation
      exBlk1 = std::min({rB, sB, pqB});
      exBlk2 = std::max({rB, sB, pqB});
      std::vector<size_t> LCat2 = {};
      reassembleCatDet2RAS(LCat, nDet, exBlk1, exBlk2, LCatAlt, nDetAlt, Jrs, LCat2);

      #pragma omp for schedule(static, 1)
      for (auto j1 = 0ul; j1 < nDetAlt[0]; j1++)
      for (auto j2 = 0ul; j2 < nDetAlt[1]; j2++) {
        size_t j12off = j2 * LCatAlt[1] + j1 * LCatAlt[0];
        Jras[Jrs[0]] = j1;
        Jras[Jrs[1]] = j2;
        // Gathering D
        for (auto j3 = 0ul; j3 < nDetAlt[2]; j3++) {
          j = j12off + j3 * LCatAlt[2];
          std::copy_n(DBlk+j*nNZrsB, nNZrsB, Dscr+j3*nNZrsB);
        }
        // Scattering moERI
        const int * exList_pq = pqList->pointerAtDet(Jras[pqB]) + nNZpqof * 4;
        for (auto iNZpq = 0ul; iNZpq < nNZpqB; iNZpq++, exList_pq+=4) {
          UNPACK_EXCITATIONLIST_4(exList_pq, p, q, l1, sign);
          p += pqoff;
          q += pqoff;
          const int * exList_rs = rsList->pointerAtDet(Jras[sB], Jras[rB]) + nNZrsof * 5;
          for (auto iNZrs = 0ul; iNZrs < nNZrsB; iNZrs++, exList_rs+=5) {
            UNPACK_EXCITATIONLIST_5(exList_rs, s, r, l1, l2, sign);
            s += soff;
            r += roff;
            ERIscr[iNZrs+iNZpq*nNZrsB] = moERI(r, s, p, q);
          }
        }
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
            nNZpqB,nDetAlt[2],nNZrsB,MatsT(1.),ERIscr,nNZrsB,
            Dscr,nNZrsB,MatsT(0.),Gscr,nNZpqB);
        // Scattering back to GBlk
        for (auto j3 = 0ul; j3 < nDetAlt[2]; j3++) {
          j = j12off + j3 * LCatAlt[2];
          MatAdd('N','N',nNZpqB,1,MatsT(1.),Gscr+j3*nNZpqB,nNZpqB,
                MatsT(1.),GBlk+j*nNZpqB,nNZpqB,GBlk+j*nNZpqB,nNZpqB);
        }
      }
    // Case 3: three RAS spaces are involved in excitation
    //         implicitly that pBlk != rBlk, qBlk != sBlk
    } else if (eCase == 3) {
      // LCatAlt == LCat, nDetAlt == nDet
      #pragma omp for schedule(static, 1)
      for (auto j3 = 0ul; j3 < nDet[2]; j3++)
      for (auto j2 = 0ul; j2 < nDet[1]; j2++)
      for (auto j1 = 0ul; j1 < nDet[0]; j1++) {
        Jras[0] = j1;
        Jras[1] = j2;
        Jras[2] = j3;
        j = j1 * LCat[0] + j2 * LCat[1] + j3 * LCat[2];
        // assemble moERI
        const int * exList_pq = pqList->pointerAtDet(Jras[pqB]) + nNZpqof * 4;
        for (auto iNZpq = 0ul; iNZpq < nNZpqB; iNZpq++, exList_pq+=4) {
          UNPACK_EXCITATIONLIST_4(exList_pq, p, q, l1, sign);
          p += pqoff;
          q += pqoff;
          const int * exList_rs = rsList->pointerAtDet(Jras[sB], Jras[rB]) + nNZrsof * 5;
          for (auto iNZrs = 0ul; iNZrs < nNZrsB; iNZrs++, exList_rs+=5) {
            UNPACK_EXCITATIONLIST_5(exList_rs, s, r, l1, l2, sign);
            s += soff;
            r += roff;
            ERIscr[iNZrs+iNZpq*nNZrsB] = moERI(r, s, p, q) - moERI(p, s, r, q);
          }
        }
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
            nNZpqB,1,nNZrsB,MatsT(1.),ERIscr,nNZrsB,
            DBlk+j*nNZrsB,nNZrsB,MatsT(1.),GBlk+j*nNZpqB,nNZpqB);
      }
    }
#ifdef _DEBUG_CIBUILDER_RASCI_IMPL
    prettyPrintSmart(std::cout,"LL RASCI Sigma G off-diag --", GBlk, nNZpqB, LCat[3], nNZpqB);
#endif
  } // RASCI::SigmaG_od

  /**
   * \brief Build intermediates G blockwise for GHF RAS Sigma
   *        off-diagonal rs -- off-diagonal pq
   *        Using block excitationlist as <J|E_sr|L>
   *        G(pq, J) = (pq|rs) * D(rs, J)
   *
   */
  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::SigmaG_oo(MatsT * GBlk, MatsT * Gscr, MatsT * DBlk,
        MatsT * Dscr, MatsT * ERIscr, InCore4indexTPI<MatsT> & moERI, size_t nNZrsof,
        size_t nNZrsB, size_t nNZpqof, size_t nNZpqB,
        std::shared_ptr<const ExcitationList> pqList,
        std::shared_ptr<const ExcitationList> rsList, size_t rB, size_t sB, size_t pB, size_t qB,
        std::vector<size_t> & LCat, std::vector<size_t> & nDet,
        int roff, int soff, int poff, int qoff) {

    std::vector<size_t> exRAS(3, 0);
    std::vector<size_t> Jrs(3, 0);
    std::vector<size_t> LCatAlt(3, 0);
    std::vector<size_t> nDetAlt(3, 0);
    std::vector<size_t> Jras(3, 0);
    size_t exBlk1, exBlk2, j, p, q, r, s, l1, l2;
    double sign;

    // determine excitation cases:
    //  2. two RAS spaces are involved in excitation
    //  3. three RAS spaces are involved in excitation
    exRAS[rB] = 1;
    exRAS[sB] = 1;
    exRAS[pB] = 1;
    exRAS[qB] = 1;
    int eCase = std::accumulate(exRAS.begin(), exRAS.end(), 0);
    // Case 2: two RAS spaces are involved in excitation
    if (eCase == 2) {
      // locate the RAS space involved in excitation
      exBlk1 = std::min({rB, sB, pB, qB});
      exBlk2 = std::max({rB, sB, pB, qB});
      std::vector<size_t> LCat2 = {};
      reassembleCatDet2RAS(LCat, nDet, exBlk1, exBlk2, LCatAlt, nDetAlt, Jrs, LCat2);

      #pragma omp for schedule(static, 1)
      for (auto j1 = 0ul; j1 < nDetAlt[0]; j1++)
      for (auto j2 = 0ul; j2 < nDetAlt[1]; j2++) {
        size_t j12off = j2 * LCatAlt[1] + j1 * LCatAlt[0];
        Jras[Jrs[0]] = j1;
        Jras[Jrs[1]] = j2;
        // Gathering D
        for (auto j3 = 0ul; j3 < nDetAlt[2]; j3++) {
          j = j12off + j3 * LCatAlt[2];
          std::copy_n(DBlk+j*nNZrsB, nNZrsB, Dscr+j3*nNZrsB);
        }
        // Scattering moERI
        const int * exList_pq = pqList->pointerAtDet(Jras[pB], Jras[qB]) + nNZpqof * 5;
        for (auto iNZpq = 0ul; iNZpq < nNZpqB; iNZpq++, exList_pq+=5) {
          UNPACK_EXCITATIONLIST_5(exList_pq, p, q, l1, l2, sign);
          p += poff;
          q += qoff;
          const int * exList_rs = rsList->pointerAtDet(Jras[sB], Jras[rB]) + nNZrsof * 5;
          for (auto iNZrs = 0ul; iNZrs < nNZrsB; iNZrs++, exList_rs+=5) {
            UNPACK_EXCITATIONLIST_5(exList_rs, s, r, l1, l2, sign);
            s += soff;
            r += roff;
            if (not (pB != rB and qB != sB))
              ERIscr[iNZrs+iNZpq*nNZrsB] = moERI(r, s, p, q);
            else ERIscr[iNZrs+iNZpq*nNZrsB] = moERI(r, s, p, q) - moERI(p, s, r, q);
          }
        }
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
            nNZpqB,nDetAlt[2],nNZrsB,MatsT(1.),ERIscr,nNZrsB,
            Dscr,nNZrsB,MatsT(0.),Gscr,nNZpqB);
        // Scattering back to GBlk
        for (auto j3 = 0ul; j3 < nDetAlt[2]; j3++) {
          j = j12off + j3 * LCatAlt[2];
          MatAdd('N','N',nNZpqB,1,MatsT(1.),Gscr+j3*nNZpqB,nNZpqB,
                MatsT(1.),GBlk+j*nNZpqB,nNZpqB,GBlk+j*nNZpqB,nNZpqB);
        }
      }
    // Case 3: three RAS spaces are involved in excitation
    //         implicitly that pBlk != rBlk, qBlk != sBlk
    } else if (eCase == 3) {
      // LCatAlt == LCat, nDetAlt == nDet
      #pragma omp for schedule(static, 1)
      for (auto j3 = 0ul; j3 < nDet[2]; j3++)
      for (auto j2 = 0ul; j2 < nDet[1]; j2++)
      for (auto j1 = 0ul; j1 < nDet[0]; j1++) {
        Jras[0] = j1;
        Jras[1] = j2;
        Jras[2] = j3;
        j = j1 * LCat[0] + j2 * LCat[1] + j3 * LCat[2];
        // assemble moERI
        const int * exList_pq = pqList->pointerAtDet(Jras[pB], Jras[qB]) + nNZpqof * 5;
        for (auto iNZpq = 0ul; iNZpq < nNZpqB; iNZpq++, exList_pq+=5) {
          UNPACK_EXCITATIONLIST_5(exList_pq, p, q, l1, l2, sign);
          p += poff;
          q += qoff;
          const int * exList_rs = rsList->pointerAtDet(Jras[sB], Jras[rB]) + nNZrsof * 5;
          for (auto iNZrs = 0ul; iNZrs < nNZrsB; iNZrs++, exList_rs+=5) {
            UNPACK_EXCITATIONLIST_5(exList_rs, s, r, l1, l2, sign);
            s += soff;
            r += roff;
            if (not (pB != rB and qB != sB))
              ERIscr[iNZrs+iNZpq*nNZrsB] = moERI(r, s, p, q);
            else ERIscr[iNZrs+iNZpq*nNZrsB] = moERI(r, s, p, q) - moERI(p, s, r, q);
          }
        }
        blas::gemm(blas::Layout::ColMajor, blas::Op::Trans, blas::Op::NoTrans,
            nNZpqB,1,nNZrsB,MatsT(1.),ERIscr,nNZrsB,
            DBlk+j*nNZrsB,nNZrsB,MatsT(1.),GBlk+j*nNZpqB,nNZpqB);
      }
    }
#ifdef _DEBUG_CIBUILDER_RASCI_IMPL
    prettyPrintSmart(std::cout,"LL RASCI Sigma G off-off --", GBlk, nNZpqB, LCat[3], nNZpqB);
#endif
  } // RASCI::SigmaG_oo

  /**
   * \brief Update Sigma 2e for GHF RAS -- diagonal
   *        Using block excitation list as <J|E_sr|L> 
   *        Sigma_2e[K] = \sum_pq 0.5 * <K|E_pq|J> * G(pq,J)
   *
   */
  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::Sigma2e_diag(MatsT * Sigma, MatsT * GBlk, size_t nNZof,
        size_t nNZB, std::shared_ptr<const ExcitationList> exList, size_t pqB,
        std::vector<size_t> & LCat,
        std::vector<size_t> & nDet, double symmFc) {

    // re-assemble LCatAlt and NDAlt
    std::vector<size_t> LCatAlt(3, 0);
    std::vector<size_t> nDetAlt(3, 0);
    reassembleCatDet1RAS(LCat, nDet, pqB, LCatAlt, nDetAlt);

    size_t j23off, k1, k, j, p, q, j1, j2, j3, iNZ;
    double sign;

    #pragma omp parallel for schedule(static, 1) default(shared) private(j1, j2, j3, iNZ, k1, k, j23off, j, sign, p, q)
    for (j3 = 0; j3 < nDetAlt[2]; j3++)
    for (j2 = 0; j2 < nDetAlt[1]; j2++)
    for (j1 = 0; j1 < nDetAlt[0]; j1++) {
      j23off = j3 * LCatAlt[2] + j2 * LCatAlt[1];
      j = j1 * LCatAlt[0] + j23off;
      const int * exList_j1 = exList->pointerAtDet(j1) + nNZof * 4;
      for (iNZ = 0; iNZ < nNZB; iNZ++, exList_j1+=4) {
        UNPACK_EXCITATIONLIST_4(exList_j1, p, q, k1, sign);
        sign = sign * symmFc;
        k = k1 * LCatAlt[0] + j23off;
        Sigma[k] += sign * GBlk[iNZ + j * nNZB];
      }
    }
    // End openMP
  } // RASCI::Sigma2e_diag

  /**
   * \brief Update Sigma 2e for GHF RAS -- off-diagonal
   *        Using block excitation list as <J|E_sr|L> 
   *        Sigma_2e[K] = \sum_pq 0.5 * <K|E_pq|J> * G(pq,J)
   *
   */
  template <typename MatsT, typename IntsT>
  void RASCI<MatsT,IntsT>::Sigma2e_offdiag(MatsT * Sigma, MatsT * GBlk, size_t nNZof,
        size_t nNZB, std::shared_ptr<const ExcitationList> exList, size_t pB, size_t qB,
        std::vector<size_t> & LCat1, std::vector<size_t> & LCat2,
        std::vector<size_t> & nDet, double symmFc) {

    // re-assemble LCatAlt and NDAlt
    std::vector<size_t> LCatAlt1(3, 0);
    std::vector<size_t> LCatAlt2(3, 0);
    std::vector<size_t> nDetAlt(3, 0);
    reassembleCatDet2RAS(LCat1, nDet, pB, qB, LCatAlt1, nDetAlt, LCatAlt2, LCat2);

    size_t j, k, k1, k2, p, q, j1, j2, j3, iNZ;
    double sign;

    #pragma omp parallel for schedule(static, 1) default(shared) private(j1, j2, j3, iNZ, k1, k, k2, j, sign, p, q)
    for (j3 = 0; j3 < nDetAlt[2]; j3++)
    for (j1 = 0; j1 < nDetAlt[0]; j1++)
    for (j2 = 0; j2 < nDetAlt[1]; j2++) {
      j = j1*LCatAlt1[0] + j2*LCatAlt1[1] + j3*LCatAlt1[2];
      const int * exList_j = exList->pointerAtDet(j1, j2) + nNZof * 5;
      for (iNZ = 0; iNZ < nNZB; iNZ++, exList_j+=5) {
        UNPACK_EXCITATIONLIST_5(exList_j, p, q, k1, k2, sign);
        k = k1*LCatAlt2[0] + k2*LCatAlt2[1] + j3*LCatAlt2[2];
        sign = sign * symmFc;
        Sigma[k] += sign * GBlk[iNZ + j * nNZB];
      }
    }
    // End openMP
  } // RASCI::Sigma2e_offdiag

}; // namespace ChronusQ







