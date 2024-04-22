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
#include <basisset.hpp>
#include <molecule.hpp>
#include <fields.hpp>
#include <hamiltonianoptions.hpp>

#include <type_traits>

namespace ChronusQ {

  /**
   *  The operator types to evaluate integrals.
   */
  enum OPERATOR {
    ELECTRON_REPULSION,
    EP_ATTRACTION,
    OVERLAP,
    KINETIC,
    NUCLEAR_POTENTIAL,
    LEN_ELECTRIC_MULTIPOLE,
    VEL_ELECTRIC_MULTIPOLE,
    MAGNETIC_MULTIPOLE,
    MAGNETIC_4COMP_rVr,
    MAGNETIC_4COMP_PVrprVP,
    MAGNETIC_4COMP_PVrmrVP
  };

  /**
   *  \brief Templated class to handle the evaluation and storage of
   *  one particle integral matrix in a finite basis set.
   */
  class ParticleIntegrals {

  protected:
    size_t NB;

  public:

    // Constructor
    ParticleIntegrals() = delete;
    ParticleIntegrals( const ParticleIntegrals & ) = default;
    ParticleIntegrals( ParticleIntegrals && ) = default;
    ParticleIntegrals(size_t nb): NB(nb) {}

    size_t nBasis() const{ return NB; }

    // Computation interfaces
    /// Evaluate AO Integrals according to a basis set.
    virtual void computeAOInts(BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) = 0;

    /// Evaluate AO Integrals according to two basis sets
    virtual void computeAOInts(BasisSet&, BasisSet&, Molecule&, EMPerturbation&,
        OPERATOR, const HamiltonianOptions&) = 0;

    virtual void clear() = 0;

    virtual void output(std::ostream&, const std::string& = "",
                        bool printFull = false) const = 0;

    virtual void broadcast(MPI_Comm comm = MPI_COMM_WORLD, int root = 0) {

#ifdef CQ_ENABLE_MPI
      // BCast ParticleIntegrals to all MPI processes
      if( MPISize(comm) > 1 ) {
        std::cerr  << "  *** Scattering a ParticleIntegrals object ***\n";
        MPIBCast(NB,root,comm);
      }
#endif

    }

    template <typename TransT>
    static std::shared_ptr<ParticleIntegrals> transform(
        const ParticleIntegrals&, char TRANS, const TransT* T, int NT, int LDT);

    virtual ~ParticleIntegrals() {}

  }; // class ElectronIntegrals

  std::ostream& operator<<(std::ostream &out,
                           const ParticleIntegrals &ints);

}; // namespace ChronusQ
