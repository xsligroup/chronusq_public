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

#include <singleslater.hpp>
#include <singleslater/multiparticless.hpp>

/*
 *     Brief: This header defines the interface between OrbitalModifier Objects and
 *            the objects that call them. The objects inherit the SCFInterface
 *            base class
 *
 */

namespace ChronusQ {

  // Assign input types to alias for ease of use
//template<typename MatsT>
//using vecMORef = std::vector<std::reference_wrapper<cqmatrix::Matrix<MatsT>>>;
//using vecEPtr  = std::vector<double*>;
template<typename MatsT>
using vecShrdPtrMat = std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>>;
template<typename MatsT>
using vecShrdPtrOrtho = std::vector<std::shared_ptr<Orthogonalization<MatsT>>>;

struct OrbitalModifierNewBase {

  OrbitalModifierNewBase(MPI_Comm mpiComm) : mpiComm(mpiComm) {};
  ~OrbitalModifierNewBase() = default;
  OrbitalModifierNewBase(const OrbitalModifierNewBase&) = delete;
  OrbitalModifierNewBase& operator=(const OrbitalModifierNewBase&) = delete;
  OrbitalModifierNewBase(OrbitalModifierNewBase&&) = delete;
  OrbitalModifierNewBase& operator=(OrbitalModifierNewBase&&) = delete;

  // the whole optimization/simulation
  virtual void run(EMPerturbation&) = 0;
  virtual void initialize(size_t maxPoints = 0) = 0;

protected:
    MPI_Comm mpiComm;         ///< MPI Communication

};
/*
 *   Brief: Abstract Base class for the modifyOrbitals object. This allows
 *          us to abstract the algorithms that modify the orbitals into one
 *          interface. This is the interface that the object that owns OrbitalModifier
 *          runs the OrbitalModifier algorithms.
 */
template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
class OrbitalModifierNew: public OrbitalModifierNewBase {

  public:
    singleSlaterT<MatsT,IntsT> &singleSlaterSystem;
    // Orthogonal Fock/Density Matrices
    std::vector<cqmatrix::Matrix<MatsT>> fockSquareOrtho;
    std::vector<cqmatrix::Matrix<MatsT>> onePDMSquareOrtho;
    std::vector<cqmatrix::Matrix<MatsT>> onePDMSquareAO;

    OrbitalModifierNew() = delete;
    OrbitalModifierNew(singleSlaterT<MatsT,IntsT> &ss, MPI_Comm mpiComm):
    singleSlaterSystem(ss),
    OrbitalModifierNewBase(mpiComm) {
      // Allocate ortho Fock and Den
      vecShrdPtrMat<MatsT> fock = this->singleSlaterSystem.getFock();
      for( auto& f : fock ) fockSquareOrtho.emplace_back(f->nRows());

      vecShrdPtrMat<MatsT> onePDM = this->singleSlaterSystem.getOnePDM();
      for( auto& d : onePDM ) {
        onePDMSquareOrtho.emplace_back(d->nRows());
        onePDMSquareAO.emplace_back(d->nRows());
      }
    }

    // Reset Function
    virtual void reInitialize(){};

    // getNewOrbitals performs only a single step in the optimization/simulation
    virtual void getNewOrbitals(EMPerturbation& pert) = 0;

    // Printing functions
    virtual void printRunHeader(EMPerturbation&) const = 0;
    virtual void printIteration(bool printDiff = true)  = 0;

    void ao2orthoFock(std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> fockSquareAO = {});
    void ao2orthoDen(std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> onePDMSquareAO = {});
    void ortho2aoDen(std::vector<cqmatrix::Matrix<MatsT>> onePDMSquareOrtho = {});
    void diagOrthoFock();
    void ortho2aoMOs();
    void ao2orthoMOs();
};

};   // Namespace ChronusQ

#include <orbitalmodifiernew/orbitaloptimizernew.hpp>
#include <orbitalmodifiernew/conventionalSCFnew.hpp>
#include <orbitalmodifiernew/impl.hpp>
//#include <orbitalmodifiernew/newtonRaphsonSCF.hpp>
#include <orbitalmodifiernew/realtimeSCF.hpp>

namespace ChronusQ {

/**
 *  \brief Construct a ConventionalSCFNew driver for the concrete
 *  single-slater type of \p ss (MultiParticleSS, HartreeFock, or KohnSham).
 *
 *  \param [in] sC  SCF controls governing the run
 *  \param [in] ss  Reference single-slater object
 *  \returns A shared_ptr to the SCF driver, or nullptr if \p ss is not a
 *           supported concrete type.
 */
template <typename MatsT, typename IntsT>
std::shared_ptr<OrbitalModifierNewBase> buildConventionalSCF(
    SCFControls sC, SingleSlater<MatsT, IntsT> &ss) {

  sC.printLevel    = ss.printLevel;
  sC.refLongName_  = ss.refLongName_;
  sC.refShortName_ = ss.refShortName_;

  if (auto *p = dynamic_cast<MultiParticleSS<MatsT, IntsT> *>(&ss))
    return std::make_shared<ConventionalSCFNew<MultiParticleSS, MatsT, IntsT>>(
        sC, *p, ss.comm);
  if (auto *p = dynamic_cast<HartreeFock<MatsT, IntsT> *>(&ss))
    return std::make_shared<ConventionalSCFNew<HartreeFock, MatsT, IntsT>>(
        sC, *p, ss.comm);
  if (auto *p = dynamic_cast<KohnSham<MatsT, IntsT> *>(&ss))
    return std::make_shared<ConventionalSCFNew<KohnSham, MatsT, IntsT>>(
        sC, *p, ss.comm);

  return nullptr;
}

/**
 *  \brief Type-erased overload of buildConventionalSCF that dispatches over
 *  the supported MatsT/IntsT instantiations.
 *
 *  \returns A shared_ptr to the SCF driver, or nullptr if \p ss is not a
 *           supported concrete type.
 */
inline std::shared_ptr<OrbitalModifierNewBase> buildConventionalSCF(
    SCFControls sC, SingleSlaterBase &ss) {

  if (auto *p = dynamic_cast<SingleSlater<double, double> *>(&ss))
    return buildConventionalSCF(sC, *p);
  if (auto *p = dynamic_cast<SingleSlater<dcomplex, double> *>(&ss))
    return buildConventionalSCF(sC, *p);
  if (auto *p = dynamic_cast<SingleSlater<dcomplex, dcomplex> *>(&ss))
    return buildConventionalSCF(sC, *p);

  return nullptr;
}

/**
 *  \brief Construct a RealTimeSCF driver for the concrete single-slater
 *  type of \p ss (MultiParticleSS, HartreeFock, or KohnSham).
 *
 *  Real-time propagation evolves a complex density, so only the dcomplex
 *  MatsT instantiations are supported.
 *
 *  \param [in] tdSCOptions  TD-SCF controls governing the run
 *  \param [in] tdPert       TD field perturbation
 *  \param [in] ss           Reference single-slater object
 *  \returns A shared_ptr to the RT driver, or nullptr if \p ss is not a
 *           supported concrete type.
 */
inline std::shared_ptr<OrbitalModifierNewBase> buildRealTimeSCF(
    TDSCFOptions &tdSCOptions, TDEMPerturbation &tdPert, SingleSlaterBase &ss) {

  if (auto *p = dynamic_cast<MultiParticleSS<dcomplex, double> *>(&ss))
    return std::make_shared<RealTimeSCF<MultiParticleSS, dcomplex, double>>(
        tdSCOptions, tdPert, *p, ss.comm);
  if (auto *p = dynamic_cast<MultiParticleSS<dcomplex, dcomplex> *>(&ss))
    return std::make_shared<RealTimeSCF<MultiParticleSS, dcomplex, dcomplex>>(
        tdSCOptions, tdPert, *p, ss.comm);
  if (auto *p = dynamic_cast<HartreeFock<dcomplex, double> *>(&ss))
    return std::make_shared<RealTimeSCF<HartreeFock, dcomplex, double>>(
        tdSCOptions, tdPert, *p, ss.comm);
  if (auto *p = dynamic_cast<HartreeFock<dcomplex, dcomplex> *>(&ss))
    return std::make_shared<RealTimeSCF<HartreeFock, dcomplex, dcomplex>>(
        tdSCOptions, tdPert, *p, ss.comm);
  if (auto *p = dynamic_cast<KohnSham<dcomplex, double> *>(&ss))
    return std::make_shared<RealTimeSCF<KohnSham, dcomplex, double>>(
        tdSCOptions, tdPert, *p, ss.comm);
  if (auto *p = dynamic_cast<KohnSham<dcomplex, dcomplex> *>(&ss))
    return std::make_shared<RealTimeSCF<KohnSham, dcomplex, dcomplex>>(
        tdSCOptions, tdPert, *p, ss.comm);

  return nullptr;
}

};   // Namespace ChronusQ
