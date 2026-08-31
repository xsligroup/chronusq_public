/*
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *
 *  Copyright (C) 2014-2026 Li Research Group (University of Washington)
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

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  class InterParticleFockBase {

    protected:

    FockBuilder<MatsT, IntsT>* upstream = nullptr;
    cqmatrix::Matrix<MatsT>* outMat = nullptr;
    SingleSlater<MatsT,IntsT>* aux_ss = nullptr;

    public:

    virtual FockBuilder<MatsT,IntsT>* getIntraParticleUpstream() {
      if( auto p = dynamic_cast<InterParticleFockBase<MatsT,IntsT>*>(upstream) ) {
        return p->getIntraParticleUpstream();
      }
      else {
        return upstream;
      }
    }

    void setAux(SingleSlater<MatsT,IntsT>* ss) {
      aux_ss = ss;
    }

    void setOutput(cqmatrix::Matrix<MatsT>* out) {
      outMat = out;
    }

    void setUpstream(FockBuilder<MatsT,IntsT>* up) {
      upstream = up;
    }

    cqmatrix::Matrix<MatsT> * getOutMat(){return outMat;}

    // Getters
    FockBuilder<MatsT,IntsT>* getUpstream(){ return upstream; }
  };

  template <typename MatsT, typename IntsT>
  class InterParticleFockBuilder:
    public FockBuilder<MatsT,IntsT>,
    public InterParticleFockBase<MatsT,IntsT>
  {

    template<typename MatsU, typename IntsU>
    friend class InterParticleFockBuilder;

    protected:

    std::shared_ptr<TPIContractions<MatsT,IntsT>> contraction = nullptr;
    GradInts<TwoPInts,IntsT>* gradTPI = nullptr;

    public:

    //
    // Constructors
    // XXX: Constructors do NOT populate the protected members - not even the
    //      copy/move (since they can't convert types)
    //
    InterParticleFockBuilder() = delete;
    InterParticleFockBuilder(HamiltonianOptions hamiltonianOptions) :
      FockBuilder<MatsT,IntsT>(hamiltonianOptions) { }

    // Other type constructors
    template <typename MatsU>
    InterParticleFockBuilder(const InterParticleFockBuilder<MatsU,IntsT> &other ) : 
      FockBuilder<MatsT,IntsT>( dynamic_cast<const FockBuilder<MatsU,IntsT>&>(other) )
      { }
    template <typename MatsU>
    InterParticleFockBuilder(FockBuilder<MatsU,IntsT> &&other ) :
      FockBuilder<MatsT,IntsT>( dynamic_cast<FockBuilder<MatsU,IntsT>&&>(other) )
      { }

    // Setters
    void setContraction(std::shared_ptr<TPIContractions<MatsT,IntsT>> cont) {
      contraction = cont;
    }

    void setPrintContractionTiming(bool arg) {
      contraction->printContractionTiming = arg;
    }

    void setGradientIntegrals(GradInts<TwoPInts,IntsT>* tpi) {
      gradTPI = tpi;
    }

    // Inter-particle Coulomb interaction
    void formInterParticleCoulomb(SingleSlater<MatsT,IntsT>&, bool increment = false);

    std::vector<double> formInterParticleCoulombGrad(SingleSlater<MatsT,IntsT>&,
      EMPerturbation&, double xHFX);

    // Interface method
    virtual void formFock(SingleSlater<MatsT,IntsT>&, EMPerturbation&,
      bool increment = false, double xHFX = 1.);

    virtual std::vector<double> getGDGrad(SingleSlater<MatsT,IntsT>&,
      EMPerturbation&, double xHFX = 1.);

  };

}

#include <fockbuilder/interparticlefock/impl.hpp>
