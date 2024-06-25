/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *
 *  This program is free software; you ca redistribute it and/or modify
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

#include <mcwavefunction.hpp>
#include <cibuilder.hpp>

namespace ChronusQ {

  /**
   *  \brief The RASCI Class. 
   */
  template <typename MatsT, typename IntsT>
  class RASCI: public CIBuilder<MatsT,IntsT> {
      
  private:
    // RASCI helper function
    void SigmaD_diag(MatsT *, MatsT *, size_t, size_t, std::shared_ptr<const ExcitationList>,
            size_t, std::vector<size_t> &, std::vector<size_t> &);
    void SigmaD_offdiag(MatsT *, MatsT *, size_t, size_t, std::shared_ptr<const ExcitationList>,
            size_t, size_t, std::vector<size_t> &, std::vector<size_t> &,
            std::vector<size_t> &);
    void Sigma1e_diag(MatsT * Sigma, MatsT * Sscr, MatsT * DBlk,
            MatsT * Dscr, OnePInts<MatsT> &, MatsT *, size_t, size_t,
            std::shared_ptr<const ExcitationList>,
            size_t, std::vector<size_t> &, std::vector<size_t> &, int);
    void Sigma1e_offdiag(MatsT *, MatsT *, MatsT *, MatsT *, OnePInts<MatsT> &,
            MatsT *, size_t, size_t, std::shared_ptr<const ExcitationList>, size_t, size_t,
            std::vector<size_t> &, std::vector<size_t> &, int, int);
    void SigmaG_dd(MatsT *, MatsT *, MatsT *,
        MatsT *, MatsT *, InCore4indexTPI<MatsT> &,
        size_t, size_t, size_t, size_t, std::shared_ptr<const ExcitationList>,
        std::shared_ptr<const ExcitationList>, size_t, size_t, std::vector<size_t> &,
        std::vector<size_t> &, int, int);
    void SigmaG_do(MatsT *, MatsT *, MatsT *,
        MatsT *, MatsT *, InCore4indexTPI<MatsT> &, size_t,
        size_t, size_t, size_t, std::shared_ptr<const ExcitationList>,
        std::shared_ptr<const ExcitationList>,
        size_t, size_t, size_t, std::vector<size_t> &,
        std::vector<size_t> &, int, int, int);
    void SigmaG_od(MatsT *, MatsT *, MatsT *,
        MatsT *, MatsT *, InCore4indexTPI<MatsT> &, size_t,
        size_t, size_t, size_t, std::shared_ptr<const ExcitationList>,
        std::shared_ptr<const ExcitationList>,
        size_t, size_t, size_t, std::vector<size_t> &,
        std::vector<size_t> &, int, int, int);
    void SigmaG_oo(MatsT *, MatsT *, MatsT *,
        MatsT *, MatsT *, InCore4indexTPI<MatsT> &, size_t,
        size_t, size_t, size_t, std::shared_ptr<const ExcitationList>,
        std::shared_ptr<const ExcitationList>,
        size_t, size_t, size_t, size_t, std::vector<size_t> &,
        std::vector<size_t> &, int, int, int, int);
    void Sigma2e_diag(MatsT *, MatsT *, size_t,
        size_t, std::shared_ptr<const ExcitationList>, size_t, std::vector<size_t> &,
        std::vector<size_t> &, double);
    void Sigma2e_offdiag(MatsT *, MatsT *, size_t,
        size_t, std::shared_ptr<const ExcitationList>, size_t, size_t, std::vector<size_t> &,
        std::vector<size_t> &, std::vector<size_t> &, double);
  
  public:
    // Constructors

    // Disable default constructor
    RASCI() = default; 
    
    // Same or Different type
    template <typename MatsU>
    RASCI(const RASCI<MatsU,IntsT> & other):
    CIBuilder<MatsT, IntsT>(other) {};

    template <typename MatsU>
    RASCI(RASCI<MatsU,IntsT> && other):
    CIBuilder<MatsT, IntsT>(other) {};
  
    // destructor
    ~RASCI() {};

    // RASCI helper function
    void reassembleCatDet1RAS(std::vector<size_t> &, std::vector<size_t> &,
            size_t, std::vector<size_t> &, std::vector<size_t> &);
    void reassembleCatDet2RAS(std::vector<size_t> &, std::vector<size_t> &,
            size_t, size_t, std::vector<size_t> &, std::vector<size_t> &,
            std::vector<size_t> &, std::vector<size_t> &);
    
    // Solving RASCI Functions
    void buildFullH(MCWaveFunction<MatsT, IntsT> &, MatsT *);
    
    void buildDiagH(MCWaveFunction<MatsT, IntsT> &, MatsT *);
	
    void buildSigma(MCWaveFunction<MatsT, IntsT> &, size_t, MatsT *, MatsT *);
    void buildMu(MCWaveFunction<MatsT, IntsT> &, size_t, MatsT *, MatsT *, EMPerturbation & pert);

    void computeOneRDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, cqmatrix::Matrix<MatsT> &);
    void computeTwoRDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, InCore4indexTPI<MatsT> &);
    void computeTDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, MatsT *, cqmatrix::Matrix<MatsT> &);
  }; // class RASCI

}; // namespace ChronusQ
