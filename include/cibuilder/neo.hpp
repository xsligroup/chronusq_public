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
#include <cibuilder/casci.hpp>
#include <cibuilder/rasci.hpp>
#include <cibuilder.hpp>

namespace ChronusQ {

  /**
   *  \brief The CASCI Class. 
   */
  template <typename MatsT, typename IntsT>
  class NEOCASCI: public CASCI<MatsT,IntsT> {
      
  public:
    // Constructors

    // Disable default constructor
    NEOCASCI() = default;
    
    // Same or Different type
    template <typename MatsU>
    NEOCASCI(const CASCI<MatsU,IntsT> & other):
    CASCI<MatsT, IntsT>(other) {};

    template <typename MatsU>
    NEOCASCI(CASCI<MatsU,IntsT> && other):
    CASCI<MatsT, IntsT>(other) {};
  
    // destructor
    ~NEOCASCI() {};

    // Solving CASCI Functions
    void buildFullH(MCWaveFunction<MatsT, IntsT> &, MatsT *);
    void buildDiagH(MCWaveFunction<MatsT, IntsT> &, MatsT *);
    void buildSigma(MCWaveFunction<MatsT, IntsT> &, size_t, MatsT *, MatsT *);
    
    void computeOneRDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, cqmatrix::Matrix<MatsT> &);
    void computePOneRDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, cqmatrix::Matrix<MatsT> &);
    //void computeTwoRDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, InCore4indexTPI<MatsT> &);
    void computeTDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, MatsT *, cqmatrix::Matrix<MatsT> &);
    void computePTDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, MatsT *, cqmatrix::Matrix<MatsT> &);

    // SMG 07/08/24
    // Two RDM aren't included since without doing NEO-CASSCF they aren't useful
    // They will be added back with the NEO-CASSCF update

    //// For NEO
    //void computePOneRDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, cqmatrix::Matrix<MatsT> &);
    //void computePTDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, MatsT *, cqmatrix::Matrix<MatsT> &);

  }; // class CASCI

}; // namespace ChronusQ
