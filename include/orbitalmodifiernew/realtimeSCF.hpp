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

/*
 *    This is a proposed object for carrying out a composite SCF calculation
 *    (where you can use a combination of several methods). However, this
 *    object has not been implemented or tested in anyway. It is included
 *    here for someone in the future to use it if they want to. 
 *
 */
#include <orbitalmodifiernew.hpp>
#include <realtime/enums.hpp>
#include <realtime/fields.hpp>


namespace ChronusQ {


  /**
   *  \brief A struct to store information pertinent to the current
   *  state of the time propagation
   */
  struct IntegrationProgressNew {

    double  currentTime = 0.; ///< Current time point
    size_t  currentStep = 0;  ///< Step index of current time point
    double  currentDeltaT;   ///< Current step size
    size_t  maxSavePoints;   ///< Maximum number of save points
    size_t  lastSavePoint = 0;   ///< Index of last save point

    std::vector<double> time;
    std::vector<double> energy;
    std::vector<std::array<double,3>> electricDipole;

    // Field
    std::vector<std::array<double,3>> electricDipoleField;
  };

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
class RealTimeSCF : public OrbitalModifierNew<singleSlaterT, MatsT, IntsT> {

  TDSCFOptions     tdSCFOptions;   ///< Integration scheme (MMUT, etc)
  TDEMPerturbation  &tdEMPerturbation;        ///< TD field perturbation
  EMPerturbation    staticEMPerturbation;     ///< SCF Perturbation

  IntegrationProgressNew integrationProgress;  ///< Current state of the time propagation

  int printLevel = 1; ///< Amount of printing in RT calc
  bool printDen = false; ///< Print density to out file
  bool printContractionTiming =false; ///< Print contraction timing during RT propagation

  typedef MatsT*                 oper_t;
  typedef std::vector<oper_t>       oper_t_coll;

  SafeFile savFile; ///< Data File

  //singleSlaterT<dcomplex,IntsT>    singleSlaterSystem; ///< Total system with complex matrices
  //std::vector<SingleSlater<MatsT, IntsT>*> systems_; ///< Objects for time propagation

  //std::vector<std::vector<SquareMatrix<MatsT>>> onePDMSquareAOSave;
  std::vector<SquareMatrix<MatsT>> previousOnePDMSquareOrtho;
  std::vector<SquareMatrix<MatsT>> previousFockSquareOrtho;
  std::vector<SquareMatrix<MatsT>> unitarySquareOrtho;

public:

  // Constructors
  // Disable default, copy and move constructors
  RealTimeSCF()                 = delete;
  RealTimeSCF(const RealTimeSCF &) = delete;
  RealTimeSCF(RealTimeSCF &&)      = delete;

  /**
   *  \brief RealTime Constructor.
   *
   *  Stores references to a "reference" SingleSlater object and
   *  CQMemManager and makes a copy of the reference into a complex
   *  SingleSlater object for the propagation.
   */
  RealTimeSCF(TDSCFOptions sC, TDEMPerturbation &tdPert, singleSlaterT<MatsT,IntsT> &referenceSS, MPI_Comm comm, CQMemManager& mem):
          tdSCFOptions(sC), tdEMPerturbation(tdPert), OrbitalModifierNew<singleSlaterT,MatsT,IntsT>(referenceSS, comm, mem) {

    if(!std::is_same<MatsT, std::complex<double>>::value) {
      throw std::runtime_error("RealTimeSCF: MatsT must be dcomplex");
    }

    this->savFile = referenceSS.savFile;
  }
  inline std::vector<double> getGrad(EMPerturbation &emPert) {
    return this->singleSlaterSystem.getGrad(emPert, false, false);
  }

  // RealTime procedural functions
  // RealTime procedural functions
  void run(EMPerturbation&) override; // From RealTimeBase
  void getNewOrbitals(EMPerturbation&) {};
  void printRunHeader(EMPerturbation&) override;
  void printIteration(bool printDiff = false) override;

  void formPropagator(std::vector<std::shared_ptr<SquareMatrix<MatsT>>> fockSquareAO = {});
  void formFock(bool,double);
  void doPropagation();
  void saveState(EMPerturbation&);
  void restoreState();
  void createRTDataSets(size_t maxPoint = 0);

  // Progress functions
  void printStepSummary();
  void printStepDetail();
  void appendStepRecord();


  // Memory functions
  virtual void initialize(size_t maxPoints = 0) override;

  /**
 *  \brief Adds a field to the time-dependent electromagnetic
 *  perturbation.
 *
 *  Calls TDEMPerturbation::addField. See include/realtime/fields.hpp
 *  for proper documentation.
 */
  template <typename... Args>
  inline void addField(Args... args){ tdEMPerturbation.addField(args...); }


  inline void setSCFPerturbation( EMPerturbation& scfp ) {
    staticEMPerturbation = scfp;
  }


  void RTFormattedLineNew(std::ostream &, std::string);
  void RTFormattedLineNew(std::ostream &, std::string, double);
  void RTFormattedLineNew(std::ostream &, std::string, size_t);
  void RTFormattedLineNew(std::ostream &, std::string, std::string);
  void RTFormattedLineNew(std::ostream &, std::string, double, std::string);

};

template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
void RealTimeSCF<singleSlaterT,MatsT,IntsT>::initialize(size_t maxPoints) {

  // Create data set (on root process)
  this->createRTDataSets(maxPoints);
  std::cout<<"xsli test RealTimeSCF initialize 2"<<std::endl;

  std::vector<std::shared_ptr<SquareMatrix<MatsT>>> onePDM = this->singleSlaterSystem.getOnePDM();
  for (auto &d: onePDM) previousOnePDMSquareOrtho.emplace_back(d->memManager(), d->dimension());
  std::vector<std::shared_ptr<SquareMatrix<MatsT>>> fock = this->singleSlaterSystem.getFock();
  for (auto &f: fock) {
    previousFockSquareOrtho.emplace_back(f->memManager(), f->dimension());
    unitarySquareOrtho.emplace_back(f->memManager(), f->dimension());
  }

  // XSLI: number of max steps should be calculated when the tdSCFOptions is set up
  tdSCFOptions.maxSteps = (size_t) ((tdSCFOptions.tMax + tdSCFOptions.deltaT / 4) / tdSCFOptions.deltaT);

}

};   // NameSpace ChronusQ
