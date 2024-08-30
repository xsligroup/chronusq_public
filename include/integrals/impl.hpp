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

#include <integrals.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/incoreritpi.hpp>
#include <particleintegrals/onepints/relativisticints.hpp>
#include <matrix.hpp>

namespace ChronusQ {

  /**
   *  \brief Allocate, compute  and store the 1-e integrals +
   *  orthonormalization matricies over the given CGTO basis.
   *
   *  Computes:
   *    Overlap + length gauge Electric Multipoles
   *    Kinetic energy matrix
   *    Nuclear potential energy matrix
   *    Core Hamiltonian (T + V)
   *    Orthonormalization matricies (Lowdin / Cholesky)
   *
   */
  template <typename IntsT>
  void Integrals<IntsT>::computeAOOneP(
      Molecule &mol, BasisSet &basis, EMPerturbation &emPert,
      const std::vector<std::pair<OPERATOR,size_t>> &ops,
      const HamiltonianOptions &options) {

    const std::array<std::string,3> dipoleList =
      { "X","Y","Z" };
    const std::array<std::string,6> quadrupoleList =
      { "XX","XY","XZ","YY","YZ","ZZ" };
    const std::array<std::string,10> octupoleList =
      { "XXX","XXY","XXZ","XYY","XYZ","XZZ","YYY",
        "YYZ","YZZ","ZZZ" };
    const std::array<std::string,9> quadrupoleListAsymm =
      { "XX","XY","XZ","YX","YY","YZ","ZX","ZY","ZZ" };

    size_t NB = basis.nBasis;

    bool electron = (options.particle.charge < 0.);

    std::string prefix = electron ? "INTS/" : "PINTS/";

#if 0
      std::cout << std::endl<<"Molecular Geometry in ComputeOneP: (Bohr)"<<std::endl;
      size_t i = 0;
      for( Atom& atom : mol.atoms ) {

        std::cout << std::right <<"AtomicNumber = " << std::setw(4) << atom.atomicNumber 
                  << std::right <<"  X= "<< std::setw(16) << atom.coord[0]
                  << std::right <<"  Y= "<< std::setw(16) << atom.coord[1]
                  << std::right <<"  Z= "<< std::setw(16) << atom.coord[2] <<std::endl;
        i += 3;
 
      }

      for(auto& shell: basis.shells)
        for(i = 0; i < 3; i++)
          std::cout << shell.O[i] << std::endl;
#endif

    for (const std::pair<OPERATOR,size_t> &op : ops)
      switch (op.first) {
      case OVERLAP:
        overlap = std::make_shared<OnePInts<IntsT>>(NB);
        overlap->computeAOInts(basis, mol, emPert, OVERLAP, options);
        if( savFile.exists() )
          savFile.safeWriteData(prefix + "OVERLAP", overlap->pointer(), {NB,NB});
        break;

      case KINETIC:
        kinetic = std::make_shared<OnePInts<IntsT>>(NB);
        kinetic->computeAOInts(basis, mol, emPert, KINETIC, options);
        if( savFile.exists() )
          savFile.safeWriteData(prefix + "KINETIC", kinetic->pointer(), {NB,NB});
        break;

      case NUCLEAR_POTENTIAL:
        // Use Libint for point nuclei, in-house for Gaussian nuclei
        if (options.OneEScalarRelativity)
          potential = std::make_shared<OnePRelInts<IntsT>>(
              NB, options.OneESpinOrbit);
        else
          potential = std::make_shared<OnePInts<IntsT>>(NB);
        potential->computeAOInts(basis, mol, emPert, NUCLEAR_POTENTIAL, options);
        if( savFile.exists() ) {
          std::string potentialTag = options.finiteWidthNuc ? "_FINITE_WIDTH" : "";
          savFile.safeWriteData(prefix + "POTENTIAL" + potentialTag,potential->pointer(),{NB,NB});
        }
        break;

      case LEN_ELECTRIC_MULTIPOLE:
        lenElectric = std::make_shared<MultipoleInts<IntsT>>(NB, op.second, true);
        // If 4C, we convert OnePInts stored in MultipoleInts object to be OnePRelInts (only handles dipole)
        if(op.second == 1 and options.OneEScalarRelativity) lenElectric->convert2OnePRelInts(NB, options.OneESpinOrbit);
        lenElectric->computeAOInts(basis, mol, emPert, LEN_ELECTRIC_MULTIPOLE, options);

        // If 4C, we gather all components of dipole integralsa (only handles dipole)
        if(op.second == 1 and options.OneEScalarRelativity) lenElectric4C = lenElectric->gather4CDipole();

        if( savFile.exists() ) {
          // Length Gauge electric dipole
          for(auto i = 0; i < 3; i++)
            savFile.safeWriteData(prefix + "ELEC_DIPOLE_LEN_" +
              dipoleList[i], (*lenElectric)[dipoleList[i]]->pointer(), {NB,NB} );

          // Length Gauge electric quadrupole
          if(op.second >= 2)
            for(auto i = 0; i < 6; i++)
              savFile.safeWriteData(prefix + "ELEC_QUADRUPOLE_LEN_" +
                quadrupoleList[i], (*lenElectric)[quadrupoleList[i]]->pointer(), {NB,NB} );

          // Length Gauge electric octupole
          if(op.second >= 3)
            for(auto i = 0; i < 10; i++)
              savFile.safeWriteData(prefix + "/ELEC_OCTUPOLE_LEN_" +
                octupoleList[i], (*lenElectric)[octupoleList[i]]->pointer(), {NB,NB} );
        }
        break;

      case VEL_ELECTRIC_MULTIPOLE:
        velElectric = std::make_shared<MultipoleInts<IntsT>>(NB, op.second, true);
        velElectric->computeAOInts(basis, mol, emPert, VEL_ELECTRIC_MULTIPOLE, options);
        if( savFile.exists() ) {
          // Velocity Gauge electric dipole
          for(auto i = 0; i < 3; i++)
            savFile.safeWriteData(prefix + "ELEC_DIPOLE_VEL_" +
              dipoleList[i], (*velElectric)[dipoleList[i]]->pointer(), {NB,NB} );

          // Velocity Gauge electric quadrupole
          if(op.second >= 2)
            for(auto i = 0; i < 6; i++)
              savFile.safeWriteData(prefix + "ELEC_QUADRUPOLE_VEL_" +
                quadrupoleList[i], (*velElectric)[quadrupoleList[i]]->pointer(), {NB,NB} );

          // Velocity Gauge electric octupole
          if(op.second >= 3)
            for(auto i = 0; i < 10; i++)
              savFile.safeWriteData(prefix + "ELEC_OCTUPOLE_VEL_" +
                octupoleList[i], (*velElectric)[octupoleList[i]]->pointer(), {NB,NB} );
        }
        break;

      case MAGNETIC_MULTIPOLE:
        magnetic = std::make_shared<MultipoleInts<IntsT>>(NB, op.second, false);
        magnetic->computeAOInts(basis, mol, emPert, MAGNETIC_MULTIPOLE, options);
        if( savFile.exists() ) {
          // Magnetic Dipole
          for(auto i = 0; i < 3; i++)
            savFile.safeWriteData(prefix + "MAG_DIPOLE_" +
              dipoleList[i], (*magnetic)[dipoleList[i]]->pointer(), {NB,NB} );

          // Magnetic Quadrupole
          if(op.second >= 2)
            for(auto i = 0; i < 6; i++)
              savFile.safeWriteData(prefix + "MAG_QUADRUPOLE_" +
                quadrupoleList[i], (*magnetic)[i+3]->pointer(), {NB,NB} );
        }
        break;

      // Calculate additional integrals, if using GIAO + X2C
      case MAGNETIC_4COMP_rVr:
        if (options.x2cType==X2C_TYPE::ONEE and options.basisType==ChronusQ::COMPLEX_GIAO) {
          rVr = std::make_shared<VectorInts<IntsT>>(NB, 2, true);
          rVr->computeAOInts(basis, mol, emPert, MAGNETIC_4COMP_rVr, options);
        } else {
          CErr("Case MAGNETIC_4COMP only for GIAO + X2CCH!",std::cout);
        }
        break;
      case MAGNETIC_4COMP_PVrprVP:
        if (options.x2cType==X2C_TYPE::ONEE and options.basisType==ChronusQ::COMPLEX_GIAO) {
          PVrprVP = std::make_shared<VectorInts<IntsT>>(NB, 2, false);
          PVrprVP->computeAOInts(basis, mol, emPert, MAGNETIC_4COMP_PVrprVP, options); 
        } else {
          CErr("Case MAGNETIC_4COMP only for GIAO + X2CCH!",std::cout);
        }
        break;
      case MAGNETIC_4COMP_PVrmrVP:
        if (options.x2cType==X2C_TYPE::ONEE and options.basisType==ChronusQ::COMPLEX_GIAO) {
          PVrmrVP = std::make_shared<VectorInts<IntsT>>(NB, 2, false);
          PVrmrVP->computeAOInts(basis, mol, emPert, MAGNETIC_4COMP_PVrmrVP, options);  
        } else {
          CErr("Case MAGNETIC_4COMP only for GIAO + X2CCH!",std::cout);
        }
        break;

      case ELECTRON_REPULSION:
        CErr("Electron repulsion integrals are not implemented in computeAOOneP,"
             " they are implemented in TwoPInts",std::cout);
        break;

        default:
          break;
      }

  }; // AOIntegrals<IntsT>::computeAOOneP



  // Computes the integrals necessary for the gradients
  template <typename IntsT>
  void Integrals<IntsT>::computeGradInts(
    Molecule &mol, BasisSet &basis, EMPerturbation &emPert,
    const std::vector<std::pair<OPERATOR,size_t>> &ops,
    const HamiltonianOptions &options) {

    size_t NB = basis.nBasis;
    size_t NAt = mol.nAtoms;

    auto computeOneE = [&](std::shared_ptr<GradInts<OnePInts,IntsT>>& p, OPERATOR o) {
      if (p == nullptr)
        p = std::make_shared<GradInts<OnePInts,IntsT>>(NB, NAt);
      else
        p->clear();

      p->computeAOInts(basis, mol, emPert, o, options);
    };

    for ( auto& op: ops ) {


      switch (op.first) {
        case OVERLAP:
          computeOneE(gradOverlap, op.first);
          break;

        case KINETIC:
          computeOneE(gradKinetic, op.first);
          break;

        case NUCLEAR_POTENTIAL:
          if ( options.OneEScalarRelativity )
            CErr("Relativistic gradients not yet implemented!");
          computeOneE(gradPotential, op.first);
          break;

        case MAGNETIC_MULTIPOLE:
        case LEN_ELECTRIC_MULTIPOLE:
        case VEL_ELECTRIC_MULTIPOLE:
          CErr("Gradients of multipoles are not yet implemented");
          break;

        case ELECTRON_REPULSION:
          if ( gradERI == nullptr )
            CErr("ERI gradients must be allocated outside of computeGradInts!");

          gradERI->computeAOInts(basis, mol, emPert, op.first, options);
          break;
        default:
          break;
      }

      // TODO: Write code to save gradient integrals to binary file. Should be
      //   toggleable with default off to keep binary file to reasonable size.

    }

    if(options.includeTau){
      // Note: we allocate NAtom*3 number of S0a matrices to use existing infrastructure, 
      //       but we only use the three matrices that correspond the derivatives on ket:
      //       <ϕ_P | d/dx ϕ_Q>, 
      //       <ϕ_P | d/dy ϕ_Q>, 
      //       <ϕ_P | d/dz ϕ_Q>
      //       we places these three matrices as the first three elements of S0a
      if (S0a == nullptr)
        S0a = std::make_shared<GradInts<OnePInts,IntsT>>(NB, NAt);
      else
        S0a->clear();
      
      std::vector<IntsT*> S0aPtrs(3*NAt, nullptr);
      for (auto i = 0; i < 3*NAt; i++) S0aPtrs[i] = (*S0a)[i]->pointer();

      OnePInts<IntsT>::OnePDriverLibint(
        libint2::Operator::overlap, mol, basis, S0aPtrs, options.particle, 1, 1
      );
      //for (size_t i = 0; i < S0aPtrs.size(); i++)
      //  prettyPrintSmart(std::cout,"ovlp0a " +  std::to_string(i), S0aPtrs[i],NB,NB,NB);
    }

  }; // AOIntegrals<IntsT>::computeGradInts


}; // namespace ChronusQ
