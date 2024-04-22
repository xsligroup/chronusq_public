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

#include <singleslater.hpp>

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
      for( auto& f : fock ) fockSquareOrtho.emplace_back(f->dimension());

      vecShrdPtrMat<MatsT> onePDM = this->singleSlaterSystem.getOnePDM();
      for( auto& d : onePDM ) {
        onePDMSquareOrtho.emplace_back(d->dimension());
        onePDMSquareAO.emplace_back(d->dimension());
      }
    }


    // getNewOrbitals performs only a single step in the optimization/simulation
    virtual void getNewOrbitals(EMPerturbation& pert) = 0;

    // Printing functions
    virtual void printRunHeader(EMPerturbation&)   = 0;
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
