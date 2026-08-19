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
   * \brief A general CI builder
   */
  template <typename MatsT, typename IntsT>
  class CASHelper {

  public:
    // Sections where the actual CI building is done
    // Note for all of these the argument passed as MatsT* is simply added to
    // i.e., no zero-ing is done in these functions, it is the role of the
    // derived class to zero out blocks if needed
    void buildFullHOneParticle(MCWaveFunction<MatsT,IntsT>&,
                               MatsT*,
                               std::shared_ptr<const ExcitationList>,
                               const std::string,
                               const std::string);
    void buildFullHTwoParticle(MCWaveFunction<MatsT,IntsT>&,
                               MatsT*,
                               std::shared_ptr<const ExcitationList>,
                               std::shared_ptr<const ExcitationList>,
                               const std::string);
    void buildDiagHOneParticle(MCWaveFunction<MatsT,IntsT>&,
                               MatsT*,
                               std::shared_ptr<const ExcitationList>,
                               const std::string,
                               const std::string);
    void buildDiagHTwoParticle(MCWaveFunction<MatsT,IntsT>&,
                               MatsT*,
                               size_t,
                               size_t,
                               std::shared_ptr<const ExcitationList>,
                               std::shared_ptr<const ExcitationList>,
                               const std::string);
    void buildSigmaOneParticle(MCWaveFunction<MatsT,IntsT>&,
                               MatsT * C,
                               MatsT * sigma,
                               size_t,
                               size_t,
                               std::shared_ptr<const ExcitationList>,
                               const std::string,
                               const std::string);
    void buildSigmaTwoParticle(MCWaveFunction<MatsT,IntsT>&,
                               MatsT * C,
                               MatsT * sigma,
                               size_t,
                               size_t,
                               std::shared_ptr<const ExcitationList>,
                               std::shared_ptr<const ExcitationList>,
                               const std::string,
                               const double chargeproduct = 1.0);

    template<typename ... MatsArgs>
    void transposeVectors(size_t nVec,
                          size_t,
                          size_t,
                          MatsArgs...);

    void addBlockToMatrixBlockDiagonal(size_t,
                                       size_t,
                                       MatsT,
                                       MatsT *,
                                       MatsT *);
    void addVecToBlockedVector(size_t,
                               size_t,
                               MatsT,
                               MatsT*,
                               MatsT*);
    //void oneRDM(){};
    //void TwoRDM(){};
    void computeTDM(MCWaveFunction<MatsT,IntsT>&,
                    MatsT * Cm,
                    MatsT * Cn,
                    size_t,
                    std::shared_ptr<const ExcitationList>,
                    cqmatrix::Matrix<MatsT> &);

}; // class CASHelper


/**
 *  \brief The CASCI Class. 
 */
template <typename MatsT, typename IntsT>
class CASCI: public CASHelper<MatsT,IntsT>,
              public CIBuilder<MatsT,IntsT> {
    
public:
  // Constructors

  // Disable default constructor
  CASCI() = default;
  
  // Same or Different type
  template <typename MatsU>
  CASCI(const CASCI<MatsU,IntsT> & other):
  CIBuilder<MatsT, IntsT>(other) {};

  template <typename MatsU>
  CASCI(CASCI<MatsU,IntsT> && other):
  CIBuilder<MatsT, IntsT>(other) {};

  // destructor
  ~CASCI() {};

  // Solving CASCI Functions
  void buildFullH(MCWaveFunction<MatsT, IntsT> &, MatsT *);
  void buildDiagH(MCWaveFunction<MatsT, IntsT> &, MatsT *);
  void buildSigma(MCWaveFunction<MatsT, IntsT> &, size_t, MatsT *, MatsT *);
  void buildMu(MCWaveFunction<MatsT, IntsT> &, size_t, MatsT *, MatsT *, EMPerturbation & pert);
  
  void computeOneRDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, cqmatrix::Matrix<MatsT> &);
  void computeTwoRDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, InCore4indexTPI<MatsT> &);
  void computeTDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, MatsT *, cqmatrix::Matrix<MatsT> &);
}; // class CASCI

}; // namespace ChronusQ
