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

#include <cerr.hpp>
#include <chronusq_sys.hpp>
#include <manybodywavefunction.hpp>
#include <cubegen.hpp>
#include <mcwavefunction.hpp>
#include <singleslater.hpp>
#include <singleslater/neoss.hpp>

// RT Headers
#include <orbitalmodifieroptions.hpp>
#include <realtime/realtimesingleslater/fields.hpp>

// RTMR Headers
#include <realtime/realtimemultislater/vectormanager.hpp>

#include <cubegenoptions.hpp>

namespace ChronusQ {
  /**
   *  \brief A struct to store information pertinent to the time
   *  propagation procedure.
   */ 
  struct RTPostHFIntegrationScheme {

    RealTimeAlgorithm    integrationAlgorithm = RealTimeAlgorithm::RTRungeKuttaOrderFour;         ///< Integration Algorithm
    bool nonhermitian_propagation = false; ///< Nonhermitian (TDCC) or Hermitian(TDCI) propagation

    double tMax    = 0.1;  ///< Max simulation time in AU
    double deltaT  = 0.01; ///< Time-step in AU

    size_t iRstrt  = 50; ///< Restart MMUT every N steps

    size_t iSave    = 50; ///< Save progress every N steps
    size_t iCube    = 0; ///< Save progress every N steps
    size_t restoreStep = 0;  ///< Restore propagation from this step

    CubeGenOptions cubeOptsRTMS; ///< Cube settings
    std::vector<std::shared_ptr<CubeGen>> rtcubes; ///< Cube files

    size_t StatePopFreq = 0; ///< frequency to calculate (and optionally save/print) MB populations
    size_t StatePopNStates = 0;///< Number of state to calculate MB populations
  
    size_t nSteps = 0; ///< Electronic steps to update tMax
  
    bool includeSCFField = true; ///< Whether to include the SCF field
  
    size_t RealTimeCorrelationFunctionFreq =
        0; ///< frequency to calculate(and optionally save/print) the real time
           ///< correlation function
    double RealTimeCorrelationFunctionStart =
        0.0; ///< the time at which to start the time correlation function
             ///< collection
  
    double transform_threshold = 1e-10; // manually set tight/ user option?
    bool time_independent_ham =
        false; // if true -> H is time independent (we apply mu explicitly
               // separately for the field) if false -> H(t) has the field folded
               // into it (as a one electron operator aka in HCore)
    bool time_independent_ham_transformed =
        false; // only used if time_independent_ham == true. If true then we have
               // done the transform once and we can skip the rest of the MO int
               // transforms
  
  }; // struct IntegrationScheme

  /**
   *  \brief A struct to store information pertinent to the current
   *  state of the time propagation
   */
  struct IntegrationProgress {
  
    double xTime = 0.; ///< Current time point
    size_t iStep = 0;  ///< Step index of current time point
    double stepSize;   ///< Current step size
  
    RealTimeAlgorithm curStep =
        RealTimeAlgorithm::RTRungeKuttaOrderFour; ///< Integration Algorithm
  };
  
  /**
   *  \brief A struct to store the property data obtained throughout the
   *  RealTime simulation
   */
  struct IntegrationData {
  
    std::vector<double> Time;
    std::vector<double> Energy;
    std::vector<std::array<double, 3>> ElecDipole;
  
    // Field
    std::vector<std::array<double, 3>> ElecDipoleField;
  
    std::vector<dcomplex> RealTimeCorrFunc;
  };
  
  class RealTimeBase {
  public:
    SafeFile savFile;                    ///< Data File
    RTPostHFIntegrationScheme intScheme; ///< Integration scheme (SSO, RK4, etc)
    TDEMPerturbation pert;               ///< TD field perturbation
    EMPerturbation scfPert;              ///< SCF Perturbation
  
    IntegrationProgress curState; ///< Current state of the time propagation
    IntegrationData data;         ///< Data collection
  
    int printLevel = 1;   ///< Amount of printing in RT calc
    bool restart = false; ///< Restarting calc from bin file
  
    RealTimeBase() = default;
    RealTimeBase(const RealTimeBase &) = delete;
    RealTimeBase(RealTimeBase &&) = delete;
  
    ~RealTimeBase(){};
  
    // RealTimeBase procedural functions
    virtual void doPropagation() = 0;
    virtual double totalEnergy() = 0;
    virtual void createRTDataSets(size_t maxPoints) = 0;
    virtual void run(bool firstStep, EMPerturbation &emPert) = 0;
  
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
    template <typename... Args> inline void addField(Args... args) {
      pert.addField(args...);
    }
  
    inline void setTDPerturbation(TDEMPerturbation &inpert) { pert = inpert; }
  
    inline void setSCFPerturbation(EMPerturbation &scfp) { scfPert = scfp; }
  };
  
  template <typename MatsT, typename IntsT>
  class RealTimeMultiSlaterBase : public RealTimeBase {
  public:
    // Disable default, copy and move constructors
    RealTimeMultiSlaterBase() = delete;
    RealTimeMultiSlaterBase(const RealTimeMultiSlaterBase &) = delete;
    RealTimeMultiSlaterBase(RealTimeMultiSlaterBase &&) = delete;
    ~RealTimeMultiSlaterBase() { dealloc(); };
    /**
     *  \brief RealTimeMultiSlaterBase Constructor.
     *
     */
    RealTimeMultiSlaterBase(
        std::shared_ptr<RealTimeMultiSlaterVectorManagerBase> vecManager_,
        RealTimeAlgorithm MRRTAlg)
        : vecManager(std::move(vecManager_)) {
      intScheme.integrationAlgorithm = MRRTAlg;
      this->alloc();
    }; // RealTimeMultiSlaterBase constructor
  
    std::shared_ptr<RealTimeMultiSlaterVectorManagerBase> vecManager;
  
    // Do we need to transform the field
    std::valarray<double> old_amp{std::numeric_limits<double>::max(),
                                  std::numeric_limits<double>::max(),
                                  std::numeric_limits<double>::max()};
    double total_energy = 0.0;
    std::array<double, 3> Dipole;
  
    // RealTimeMultiSlater procedural functions
    inline double totalEnergy() override { return this->total_energy; };
  
    void createRTDataSets(size_t maxPoints) override;
  
    void doPropagation() override;
    void propagateWFN(bool Start, bool Finish);
    void propagateWFN_SSO(bool Start, bool Finish);
    void propagateWFN_RK4(bool Start, bool Finish);
  
    void RealTimeCorrelationFunction();
    void saveState(EMPerturbation &);
    void restoreState();
    void run(bool firstStep, EMPerturbation &emPert) override {
      // Get correct time length
      if (!firstStep) {
        intScheme.restoreStep = curState.iStep;
        intScheme.tMax = intScheme.tMax + intScheme.nSteps * intScheme.deltaT;
      }
  
      doPropagation();
    }
  
    // Progress functions
    void printRTHeader();
    void printRTStep();
    void printStepSummary();
    void printStepDetail();
    void appendStepRecord();
  
    // Memory functions
    void alloc();
    void dealloc();
  
    virtual void statePop() = 0;
    virtual void buildSigma(std::shared_ptr<SolverVectors<MatsT>> cin,
                            std::shared_ptr<SolverVectors<MatsT>> sigma_out,
                            double t, bool do_left = false) = 0;
    virtual void buildMu(std::shared_ptr<SolverVectors<MatsT>> cin,
                         std::shared_ptr<SolverVectors<MatsT>> sigma_out,
                         double t, bool do_left = false) = 0;
    virtual void constShiftTotalE() = 0;
    virtual void genInitialState() = 0;
    virtual void calculateDipole() = 0;
    virtual void formHamiltonian(double) = 0;
  
    // Generate cube files
    virtual void genCubes() = 0;
  };
  
  template <typename MatsT, typename IntsT>
  class RealTimeCI : public RealTimeMultiSlaterBase<MatsT, IntsT> {
  
  public:
    std::shared_ptr<MCWaveFunction<MatsT, IntsT>> reference_;
  
    // Constructors
  
    // Disable default, copy and move constructors
    RealTimeCI() = delete;
    RealTimeCI(const RealTimeCI &) = delete;
    RealTimeCI(RealTimeCI &&) = delete;
  
    /**
     *  \brief RealTimeCI Constructor.
     *
     *  Stores references to a "reference" MultiSlater object and
     *  CQMemManager and makes a copy of the reference into a propagated
     * wavefunction object
     */
    RealTimeCI(std::shared_ptr<MCWaveFunction<MatsT, IntsT>> reference,
               std::shared_ptr<RealTimeMultiSlaterVectorManagerBase> vecManager_,
               RealTimeAlgorithm MRRTAlg)
        : reference_(reference),
          RealTimeMultiSlaterBase<MatsT, IntsT>(
              vecManager_, MRRTAlg){}; // RealTimeMultiSlater constructor
  
    void buildSigma(std::shared_ptr<SolverVectors<MatsT>> cin,
                    std::shared_ptr<SolverVectors<MatsT>> sigma_out,
                    double t, bool do_left = false) override;
    void buildMu(std::shared_ptr<SolverVectors<MatsT>> cin,
                 std::shared_ptr<SolverVectors<MatsT>> sigma_out,
                 double t, bool do_left = false) override;
    void genInitialState() override;
    void calculateDipole() override;
    void formHamiltonian(double) override;
    void statePop() override;
    void constShiftTotalE() override {
      this->total_energy = reference_->reference().molecule().nucRepEnergy +
                           reference_->InactEnergy;
    };
  
    // Generate cube files
    void genCubes();
  }; // class RealTimeMultiSlater

}; // namespace ChronusQ
