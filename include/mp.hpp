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
#include <cerr.hpp>
#include <util/files.hpp>
#include <itersolver.hpp>
#include <posthartreefock.hpp>
#include <mcscf.hpp>
#include <mointstransformer.hpp>

namespace ChronusQ {


  class MP2Base
  {
    protected:

      // Total MP2 Energy
      double MP2_CorrE = 0.0;

    public:

      // Whether or not to construct MP2NOs
      bool makeMP2NOs=false;

      // Main MP2 interface
      virtual void runMP2(EMPerturbation &) = 0;
  };

  template <typename MatsT, typename IntsT>
  class MP2 : public MP2Base, public PostHartreeFock<MatsT,IntsT>
  {
    protected:
      std::shared_ptr<InCore4indexTPI<MatsT>> T2;
      std::shared_ptr<InCore4indexTPI<MatsT>> moERI;

      // For Direct implementation, we can save on storage by only calculating
      // sub-blocks of the full moERI tensor, but require a second moERI object
      // so that we can store the exchange type terms
      std::shared_ptr<InCore4indexTPI<MatsT>> exmoERI;

      //std::shared_ptr<cqmatrix::Matrix<MatsT>> MP2_1RDM;
      std::vector<double> MP2_NO_occs;

    public:

      // Opposite and Same Spin energy contributions
      double eos = 0.0;
      double ess = 0.0;

      // Data memebers for convenience
      size_t NB;
      size_t nocc;
      size_t nvir;

      // Constructors
      template <typename MatsU>
      MP2(std::shared_ptr<SingleSlater<MatsU,IntsT>> ref):
        PostHartreeFock<MatsT,IntsT>(ref,0) 
        {
          // Guard against using MP2 for untested applications
          if(this->ref_->nC != 1)
            CErr("MP2 NYI for non-1C references!");
          if(this->ref_->nOA != this->ref_->nOB  && this->ref_->particle.charge < 0)
            CErr("MP2 NYI for non-closed shell references!");
          if(!std::is_same<IntsT,double>::value)
            CErr("MP2 for complex IntsT untested!");

          // For saving to bin
          this->savFile = this->ref_->savFile;

          NB = this->ref_->nAlphaOrbital();
          nocc = this->ref_->nOA;
          nvir = NB - nocc;
        };

      // Functions that are purely virtual in PostHartreeFock
      void run(EMPerturbation & pert)
      {
        runMP2(pert);
      };
      void computeTDM(size_t, size_t, std::shared_ptr<cqmatrix::Matrix<MatsT>>)
      {
        CErr("Doesn't make sense to call PostHartreeFock virtual compute TDM from MP2!");
      };
      
      // MP2 functionality
      void runMP2(EMPerturbation &);
      void saveState();

      void MP2EnergyEval(EMPerturbation &);
      // Helper functions for the energy evaluation
      void MP2DirectEnergyEval(EMPerturbation &);
      void MP2InCoreEnergyEval(EMPerturbation &);
      void MP2EnergyLooper(MatsT *, MatsT*,
                           size_t,size_t,size_t,size_t,
                           size_t,size_t,size_t,size_t);
      void make1RDM();
      void generateMP2NOs();

      // Dummy implementation for PostHartreeFock Base
      void compute2TDM(size_t, size_t, std::shared_ptr<InCore4indexTPI<MatsT>>){}; 
      void compute2RDM(size_t, size_t, std::shared_ptr<InCore4indexTPI<MatsT>>){};

      // Memory concerns
      // Frees the memory associated with the moERI and T2 intermediates
      // back to CQMemManager
      void dealloc()
      {
        T2 = nullptr;
        moERI = nullptr;
      }

  };

  template <typename MatsT, typename IntsT>
  class NEOMP2: public MP2Base, public PostHartreeFock<MatsT,IntsT>
  {
    protected:
      std::shared_ptr<NEOSS<MatsT,IntsT>> neoref_;
      std::vector<std::shared_ptr<MP2Base>> subsystems;

      std::shared_ptr<InCore4indexTPI<MatsT>> epERI;
      std::shared_ptr<InCore4indexTPI<MatsT>> MixedT2;

      // Quantities for convenience (and potential future generalization)
      std::vector<std::shared_ptr<MP2<MatsT,IntsT>>> refs;
      std::vector<size_t> NB;
      std::vector<size_t> nocc;
      std::vector<size_t> nvir;

    public:

      // NEO electron-proton and proton-proton
      double ep = 0.0;

      template <typename MatsU>
      NEOMP2(std::shared_ptr<NEOSS<MatsU,IntsT>> ref):
        neoref_(ref),
        PostHartreeFock<MatsT,IntsT>(std::dynamic_pointer_cast<SingleSlater<MatsT,IntsT>>(ref),0) 
        {
          // For saving to bin
          this->savFile = this->ref_->savFile;
          // Make subsystem MP2 objects
          subsystems.push_back(std::make_shared<MP2<MatsT,IntsT>>(neoref_->template getSubsystem<SingleSlater>(std::string("Electronic"))));
          subsystems.push_back(std::make_shared<MP2<MatsT,IntsT>>(neoref_->template getSubsystem<SingleSlater>(std::string("Protonic"))));
        };
      // Functions that are purely virtual in PostHartreeFock
      void run(EMPerturbation &)
      {
        CErr("Doesn't make sense to call PostHartreeFock virtual run from NEOMP2!  The MP2 interface is runMP2()");
      };
      void computeTDM(size_t, size_t, std::shared_ptr<cqmatrix::Matrix<MatsT>>)
      {
        CErr("Doesn't make sense to call PostHartreeFock virtual compute TDM from NEOMP2!");
      };

      // MP2 functionality 
      void runMP2(EMPerturbation &);
      void saveState();
      void addepto1RDM();

      // Dummy implementation for PostHartreeFock Base
      void compute2TDM(size_t, size_t, std::shared_ptr<InCore4indexTPI<MatsT>>){}; 
      void compute2RDM(size_t, size_t, std::shared_ptr<InCore4indexTPI<MatsT>>){};

  };



}; // namespace ChronusQ
