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
#include <memmanager.hpp>
#include <singleslater.hpp>
#include <singleslater/neoss.hpp>


// RT Headers
#include <realtime/enums.hpp>
#include <realtime/fields.hpp>



namespace ChronusQ {

  /**
   *  \brief A struct to store information pertinent to the time
   *  propagation procedure.
   */ 
  struct IntegrationScheme {

    IntegrationAlgorithm intAlg  = MMUT;         ///< Integration Algorithm
    PropagationStep      rstStep = ExplicitMagnus2; ///< Restart Step
    PropagatorAlgorithm  prpAlg  = Diagonalization; ///< exp(-iF) Algorithm

    double tMax    = 0.1;  ///< Max simulation time in AU
    double deltaT  = 0.01; ///< Time-step in AU

    size_t iRstrt  = 50; ///< Restart MMUT every N steps

    size_t iSave    = 50; ///< Save progress every N steps
    size_t restoreStep = 0;  ///< Restore propagation from this step

    size_t nSteps = 0; ///< Electronic steps to update tMax

    bool   includeSCFField = true;  ///< Whether to include the SCF field

  }; // struct IntegrationScheme

  /**
   *  \brief A struct to store information pertinent to the current
   *  state of the time propagation
   */ 
  struct IntegrationProgress {

    double  xTime = 0.; ///< Current time point
    size_t  iStep = 0;  ///< Step index of current time point
    double  stepSize;   ///< Current step size

    PropagationStep curStep;  ///< Current integration step

  };

  /**
   *  \brief A struct to store the property data obtained throughout the
   *  RealTime simulation
   */ 
  struct IntegrationData {

    std::vector<double> Time;
    std::vector<double> Energy;
    std::vector<std::array<double,3>> ElecDipole;

    // Field
    std::vector<std::array<double,3>> ElecDipoleField;
  };




  struct RealTimeBase {

    SafeFile savFile; ///< Data File

    IntegrationScheme intScheme;   ///< Integration scheme (MMUT, etc)
    TDEMPerturbation  pert;        ///< TD field perturbation
    EMPerturbation    scfPert;     ///< SCF Perturbation

    IntegrationProgress curState;  ///< Current state of the time propagation
    IntegrationData     data;      ///< Data collection

    int printLevel = 1; ///< Amount of printing in RT calc
    bool printDen = false; ///< Print density to out file
    bool printContractionTiming =false; ///< Print contraction timing during RT propagation
    size_t orbitalPopFreq = 0; ///< Amount of printing in RT calc
    
    bool restart   = false; ///< Restarting calc from bin file

    RealTimeBase()                     = delete;
    RealTimeBase(const RealTimeBase &) = delete;
    RealTimeBase(RealTimeBase &&)      = delete;


    RealTimeBase( CQMemManager &memManager): memManager_(memManager){ }


    // RealTimeBase procedural functions
    virtual void doPropagation()         = 0;
    virtual double totalEnergy()         = 0;
    virtual std::vector<double> getGrad(EMPerturbation&) = 0;
    virtual void formCoreH(EMPerturbation&)              = 0;
    virtual void updateAOProperties(double) = 0;
    virtual void createRTDataSets(size_t maxPoints) = 0;

    // Progress functions
    void printRTHeader();
    void printRTStep();
    void appendStepRecord();


    /**
     *  \brief Adds a field to the time-dependent electromagnetic
     *  perturbation.
     *
     *  Calls TDEMPerturbation::addField. See include/realtime/fields.hpp
     *  for proper documentation.
     */ 
    template <typename... Args>
    inline void addField(Args... args){ pert.addField(args...); }


    inline void setSCFPerturbation( EMPerturbation& scfp ) {
      scfPert = scfp;
    }

  protected:

    CQMemManager     &memManager_; ///< Memory manager

  };


  template <template <typename, typename> class _SSTyp, typename IntsT>
  class RealTime : public RealTimeBase {
   
    typedef dcomplex*                 oper_t;
    typedef std::vector<oper_t>       oper_t_coll;

    SingleSlaterBase         *reference_ = nullptr;  ///< Initial conditions
    _SSTyp<dcomplex,IntsT>    propagator_; ///< Total system with complex matrices 
    std::vector<SingleSlater<dcomplex, IntsT>*> systems_; ///< Objects for time propagation

    std::vector<std::shared_ptr<PauliSpinorSquareMatrices<dcomplex>>> DOSav;
    std::vector<std::shared_ptr<PauliSpinorSquareMatrices<dcomplex>>> UH;
    
  public:


    // Constructors

    // Disable default, copy and move constructors
    RealTime()                 = delete;
    RealTime(const RealTime &) = delete;
    RealTime(RealTime &&)      = delete;


    /**
     *  \brief RealTime Constructor.
     *
     *  Stores references to a "reference" SingleSlater object and
     *  CQMemManager and makes a copy of the reference into a complex
     *  SingleSlater object for the propagation.
     */ 
    template <typename RefMatsT>
    RealTime(_SSTyp<RefMatsT,IntsT> &reference) : 
      RealTimeBase(reference.memManager),
      reference_(&reference), propagator_(reference) { 

      alloc<RefMatsT>(); 

    }; // RealTime constructor
  
    inline double totalEnergy(){
      //propagator_.computeEnergy();
      return propagator_.totalEnergy;
    }

    inline void formCoreH(EMPerturbation &emPert) {
      return propagator_.formCoreH(emPert, false);
    }

    inline std::vector<double> getGrad(EMPerturbation &emPert) {
      return propagator_.getGrad(emPert,false,false);
    }

    // RealTime procedural functions
    // RealTime procedural functions
    void doPropagation(); // From RealTimeBase
    void propagateStep();
    void formPropagator(size_t);
    void formFock(bool,double,size_t);
    void updateAOProperties(double t);
    void propagateWFN(size_t);
    void saveState(EMPerturbation&);
    void restoreState(); 
    void createRTDataSets(size_t maxPoints);
    void orbitalPop();

    // Progress functions
    void printRTHeader();
    void printRTStep();
    void printStepSummary();
    void printStepDetail();
    void appendStepRecord();


    // Memory functions
    template <typename MatsT>
    void alloc();

  }; // class RealTime
  

}; // namespace ChronusQ

