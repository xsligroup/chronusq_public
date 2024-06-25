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

#include <chronusq_sys.hpp>
#include <mcwavefunction.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>

// CI Headers

//#define DEBUG_CI

namespace ChronusQ {
  
  /* 
   * Brief Definition of CIBuilder Class
   */
  
  template <typename MatsT, typename IntsT>
  class CIBuilder { 
    
  public:
    
    // default Constructor
    CIBuilder() = default;
  
    // Different type
	template <typename MatsU>
    CIBuilder(const CIBuilder<MatsU,IntsT> &);
	
	template <typename MatsU>
    CIBuilder(CIBuilder<MatsU,IntsT> &&);
	
	~CIBuilder() = default;

    // Pointer convertor
    template <typename MatsU>
    static std::shared_ptr<CIBuilder<MatsU,IntsT>>
    convert(const std::shared_ptr<CIBuilder<MatsT,IntsT>>&);

	// Virtual Solving CI Functions
    virtual void buildFullH(MCWaveFunction<MatsT, IntsT> &, MatsT *) = 0;
    virtual void buildDiagH(MCWaveFunction<MatsT, IntsT> &, MatsT *) = 0;
    virtual void buildSigma(MCWaveFunction<MatsT, IntsT> &, size_t, MatsT *, MatsT *) = 0;
    virtual void buildMu(MCWaveFunction<MatsT, IntsT> &, size_t, MatsT *, MatsT *, EMPerturbation&) = 0;
    
    virtual void computeOneRDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, cqmatrix::Matrix<MatsT> &) = 0;
    virtual void computeTwoRDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, InCore4indexTPI<MatsT> &) = 0;
    virtual void computeTDM(MCWaveFunction<MatsT, IntsT> &, MatsT *, MatsT *, cqmatrix::Matrix<MatsT> &) = 0;
  }; // class CIBuilder
  

}; // namespace ChronusQ

// Include declaration for specialization of CIBuilder
#include <cibuilder/casci.hpp>
#include <cibuilder/rasci.hpp>


