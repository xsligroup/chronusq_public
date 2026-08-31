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

#include <integrals/impl.hpp>
#include <particleintegrals/inhouseaointegral.hpp>
#include <cqlinalg.hpp>
#include <cqlinalg/svd.hpp>
#include <cqlinalg/blasutil.hpp>
#include <physcon.hpp>
#include <util/matout.hpp>
#include <util/threads.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Core>


// Debug directives
//#define _DEBUGORTHO
//#define _DEBUGERI
//#define _DEBUGGIAOERI //SS
//#define _DEBUGGIAOONEE //SS 


namespace ChronusQ {

  typedef std::vector<libint2::Shell> shell_set;

  template <>
  template <size_t NOPER, bool SYMM, typename F>
  void OnePInts<dcomplex>::OnePDriverLocal(
      const F &obFunc, const Molecule &mol, BasisSet &basis, std::vector<dcomplex*> mats,
      OPERATOR op, const HamiltonianOptions &options, size_t deriv) {

    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();

    shell_set& shells = basis.shells;
    // Determine the number of basis functions for the passed shell set
    size_t NB = std::accumulate(shells.begin(),shells.end(),0,
      [](size_t init, libint2::Shell &sh) -> size_t {
        return init + sh.size();
      }
    );

    size_t NBSQ = NB*NB;

    // Determine the maximum angular momentum of the passed shell set
    int maxL = std::max_element(shells.begin(), shells.end(),
      [](libint2::Shell &sh1, libint2::Shell &sh2){
        return sh1.contr[0].l < sh2.contr[0].l;
      }
    )->contr[0].l;

    // Determine the maximum contraction depth of the passed shell set
    int maxPrim = std::max_element(shells.begin(), shells.end(),
      [](libint2::Shell &sh1, libint2::Shell &sh2){
        return sh1.alpha.size() < sh2.alpha.size();
      }
    )->alpha.size();

    std::vector<
      Eigen::Map<
        Eigen::Matrix<dcomplex,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> 
      > 
    > matMaps;

    for( auto i = 0; i < mats.size(); i++ ) {
      std::fill_n(mats[i],NBSQ,0.);  
      matMaps.emplace_back(mats[i],NB,NB);
    }

//    if(basisType == REAL_GTO)
      // pre compute all the shellpair data
//      auto pair_to_use = genShellPairs(shells,std::log(std::numeric_limits<double>::lowest()));

    #pragma omp parallel
    {
      int thread_id = GetThreadID();

    size_t n1,n2,atom1,atom2;
    // Loop over unique shell pairs
    for(size_t s1(0), bf1_s(0), s12(0); s1 < shells.size(); bf1_s+=n1, s1++){ 
      n1 = shells[s1].size(); // Size of Shell 1
      atom1 = basis.mapSh2Cen[s1]; // Index of atom for Shell 1
    for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++, s12++) {
      n2 = shells[s2].size(); // Size of Shell 2
      atom2 = basis.mapSh2Cen[s2]; // Index of atom for Shell 2

      // Round Robbin work distribution
      #ifdef _OPENMP
      if( s12 % nthreads != thread_id ) continue;
      #endif

      libint2::ShellPair pair_to_use;
      pair_to_use.init(shells[s1],shells[s2],-1000);

      auto buff = obFunc(pair_to_use, shells[s1],shells[s2]);

      // Number of matrice must match
      if (NOPER>0)
        assert(buff.size() == NOPER);
      else if (NOPER==0) {
        // If NOPER is 0, check by catagory
        if (op == NUCLEAR_POTENTIAL)
          assert(buff.size() == 6+3*mol.atomsC.size());
        else
          CErr("Use pre-defined NOPER!");
      } else
        CErr("Number of Components Undefined!");

      // Place integral blocks into their respective matricies
      auto add_shellset_to_mat = [&](size_t iOp, size_t iMat) {

        // If the integrals were screened, do nothing (currently no screening implemented for GIAO)
        //if(buff[iOp] == nullptr) return;

        Eigen::Map<
          const Eigen::Matrix<
            dcomplex,
            Eigen::Dynamic,Eigen::Dynamic,  
            Eigen::RowMajor>>
          bufMat(&buff[iOp][0],n1,n2);

        matMaps[iMat].block(bf1_s,bf2_s,n1,n2) += bufMat.template cast<dcomplex>();
      };

      switch (deriv) {

        case 0: {
          for(auto iMat = 0; iMat < buff.size(); iMat++){
            add_shellset_to_mat(iMat, iMat);
          }
        }
        break; // case deriv == 0

        case 1: {
          size_t result_idx = 0;

          // Map For gradients
          // S and T: Ax, Ay, Az, Bx, By, Bz 
            
          // First the bra and ket
          for (auto xyz = 0; xyz < 3; xyz++, result_idx++)
            add_shellset_to_mat(result_idx, 3*atom1 + xyz);

          for (auto xyz = 0; xyz < 3; xyz++, result_idx++)
            add_shellset_to_mat(result_idx, 3*atom2 + xyz);

          // Gradient of operator
          // V: Ax, Ay, Az, Bx, By, Bz, Cx, Cy, Cz 
          if (op == NUCLEAR_POTENTIAL) {
            auto nAtoms = mol.atomsC.size();
            for (auto iAt = 0; iAt < nAtoms; iAt++) {
              for ( auto xyz = 0; xyz < 3; xyz++, result_idx++) {
                add_shellset_to_mat(result_idx, 3*mol.atomsC[iAt]+ xyz);
              }
            }
          }
        }
        break; // case deriv == 1

        // S0a
        case 10: {

          assert(op == TAUS0a);
          size_t result_idx = 0;

          // Map For gradients
          // S and T: Ax, Ay, Az 
          for (auto xyz = 0; xyz < 3; xyz++, result_idx++)
            add_shellset_to_mat(result_idx, 3*atom1 + xyz);
        }
        break; // case deriv == 10

        // Len1
        case 11: {
          size_t result_idx = 0;

          // X, Y, Z
          for (auto icomp=0; icomp<3; icomp++) {
            result_idx = 0;
            for (auto xyz = 0; xyz < 3; xyz++, result_idx++)
              add_shellset_to_mat(3*result_idx+icomp, 3*(3*atom1 + xyz)+icomp);
            for (auto xyz = 0; xyz < 3; xyz++, result_idx++)
              add_shellset_to_mat(3*result_idx+icomp, 3*(3*atom2 + xyz)+icomp);
          }
        }
        break;

        // Len2
        case 12: {
          size_t result_idx = 0;

          // XX, XY, XZ, YY, YZ, ZZ
          for (auto icomp=0; icomp<6; icomp++) {
            result_idx = 0;
            for (auto xyz = 0; xyz < 3; xyz++, result_idx++)
              add_shellset_to_mat(6*result_idx+icomp, 6*(3*atom1 + xyz)+icomp);
            for (auto xyz = 0; xyz < 3; xyz++, result_idx++)
              add_shellset_to_mat(6*result_idx+icomp, 6*(3*atom2 + xyz)+icomp);
          }
        }
        break;
        
        default: {
          CErr("Required Gradient NYI in GIAO!",std::cout);
        }
        break;

      } // switch deriv

    } // Loop over s2 <= s1
    } // Loop over s1

    } // end of omp

    // Symmetrize the matricies 
    // XXX: USES EIGEN
    // FIXME: not SYMM -> creates a temporary
    for(auto nMat = 0; nMat < matMaps.size(); nMat++) {
      if(SYMM) {
        for(auto i = 0  ; i < NB; ++i)
        for(auto j = i+1; j < NB; ++j)
          matMaps[nMat](i,j) = std::conj(matMaps[nMat](j,i));
      } else {
        for(auto i = 0  ; i < NB; ++i)
        for(auto j = i+1; j < NB; ++j)
          matMaps[nMat](i,j) = - std::conj(matMaps[nMat](j,i));
      }
    }

    // Future modify matrices if particle is proton
    // If engine is K, scale it by 1/m
    if (op == KINETIC) {
      for(auto nMat = 0; nMat < matMaps.size(); nMat++)
        for(auto i = 0  ; i < NB; ++i)
        for(auto j = 0  ; j < NB; ++j)
          matMaps[nMat](i,j) *= 1.0 / options.particle.mass;
    }
    // If engine is V or q<r...r>, scale it by charge
    // X2C + NEO currently not considered 
    if(op == NUCLEAR_POTENTIAL or
       op == LEN_ELECTRIC_MULTIPOLE or
       op == VEL_ELECTRIC_MULTIPOLE) {
      for(auto nMat = 0; nMat < matMaps.size(); nMat++)
        for(auto i = 0  ; i < NB; ++i)
        for(auto j = 0  ; j < NB; ++j)
          matMaps[nMat](i,j) *= -1.0 * options.particle.charge;
    }

  }; // OnePInts::OnePDriverLocal


  template <>
  void OnePInts<dcomplex>::computeAOInts(BasisSet &basis, Molecule &mol,
      EMPerturbation &emPert, OPERATOR op, const HamiltonianOptions &options) {
    if (options.basisType == REAL_GTO)
      CErr("Real GTOs are not allowed in OneEInts<dcomplex>",std::cout);
    if (options.basisType == COMPLEX_GTO)
      CErr("Complex GTOs NYI in OneEInts<dcomplex>",std::cout);
    if (op == NUCLEAR_POTENTIAL and (options.OneEScalarRelativity or options.OneESpinOrbit))
      CErr("Relativistic integrals are implemented in OnePRelInts",std::cout);

    auto magAmp = emPert.getDipoleAmp(Magnetic);
    std::vector<dcomplex*> tmp(1, pointer());

    switch (op) {
    case OVERLAP:
      OnePInts<dcomplex>::OnePDriverLocal<1,true>(
          std::bind(&ComplexGIAOIntEngine::computeGIAOOverlapS,
                    std::placeholders::_1, std::placeholders::_2,
                    std::placeholders::_3, &magAmp[0],options.particle.charge),
          mol, basis, tmp, op, options, 0);
      break;
    case KINETIC:
      OnePInts<dcomplex>::OnePDriverLocal<1,true>(
          std::bind(&ComplexGIAOIntEngine::computeGIAOKineticT,
                    std::placeholders::_1, std::placeholders::_2,
                    std::placeholders::_3, &magAmp[0],options.particle.charge),
          mol, basis, tmp, op, options, 0);
      break;
    case NUCLEAR_POTENTIAL:
      options.finiteWidthNuc ?
      OnePInts<dcomplex>::OnePDriverLocal<1,true>(
          [&](libint2::ShellPair& pair, libint2::Shell& sh1, 
              libint2::Shell& sh2) -> std::vector<std::vector<dcomplex>> { 
            return ComplexGIAOIntEngine::computeGIAOPotentialV(
                mol.chargeDist,pair,sh1,sh2,&magAmp[0],mol.retainCNuc(),options.particle.charge);
            }, mol, basis, tmp, op, options, 0) :
      OnePInts<dcomplex>::OnePDriverLocal<1,true>(
          [&](libint2::ShellPair& pair, libint2::Shell& sh1,
              libint2::Shell& sh2) -> std::vector<std::vector<dcomplex>> {
            return ComplexGIAOIntEngine::computeGIAOPotentialV(
                pair,sh1,sh2,&magAmp[0],mol.retainCNuc(),options.particle.charge);
            }, mol, basis, tmp, op, options, 0);
      break;
    case ELECTRON_REPULSION:
      CErr("Electron repulsion integrals are not implemented in OnePInts,"
           " they are implemented in TwoEInts",std::cout);
      break;
    case LEN_ELECTRIC_MULTIPOLE:
    case VEL_ELECTRIC_MULTIPOLE:
    case MAGNETIC_MULTIPOLE:
      CErr("Requested operator is not implemented in OnePInts,"
           " it is implemented in MultipoleInts",std::cout);
      break;
    default:
      CErr("Requested operator is not implemented in OneEInts.");
      break;
    }

  };

  template <>
  void VectorInts<dcomplex>::computeAOInts(BasisSet &basis, Molecule &mol,
      EMPerturbation &emPert, OPERATOR op, const HamiltonianOptions &options) {
    if (options.basisType == REAL_GTO)
      CErr("Real GTOs are not allowed in VectorInts<dcomplex>",std::cout);
    if (options.basisType == COMPLEX_GTO)
      CErr("Complex GTOs NYI in VectorInts<dcomplex>",std::cout);
    // TangDD: Magnetic 4component integrals are placed here. lifting
    //if (options.OneEScalarRelativity or options.OneESpinOrbit)
    //  CErr("Relativistic integrals are implemented in OnePRelInts",std::cout);

    auto magAmp = emPert.getDipoleAmp(Magnetic);

    switch (op) {
    case OVERLAP:
    case KINETIC:
    case NUCLEAR_POTENTIAL:
      CErr("Requested operator is not implemented in VectorInts,"
           " it is implemented in OneEInts",std::cout);
      break;
    case ELECTRON_REPULSION:
      CErr("Electron repulsion integrals are not implemented in VectorInts,"
           " they are implemented in TwoEInts",std::cout);
      break;
    case LEN_ELECTRIC_MULTIPOLE:
      switch (order()) {
      case 1:
        OnePInts<dcomplex>::OnePDriverLocal<3,true>(
            std::bind(&ComplexGIAOIntEngine::computeGIAOEDipoleE1_len,
                      std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3, &magAmp[0],options.particle.charge),
            mol, basis, pointers(), op, options, 0);
        break;
      case 2:
        OnePInts<dcomplex>::OnePDriverLocal<6,true>(
            std::bind(&ComplexGIAOIntEngine::computeGIAOEQuadrupoleE2_len,
                      std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3, &magAmp[0],options.particle.charge),
            mol, basis, pointers(), op, options, 0);
        break;
      case 3:
        OnePInts<dcomplex>::OnePDriverLocal<10,true>(
            std::bind(&ComplexGIAOIntEngine::computeGIAOEOctupoleE3_len,
                      std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3, &magAmp[0],options.particle.charge),
            mol, basis, pointers(), op, options, 0);
        break;
      default:
        CErr("Requested operator is NYI in VectorInts.",std::cout);
        break;
      }
      break;
    case VEL_ELECTRIC_MULTIPOLE:
      CErr("Requested operator is NYI in VectorInts.",std::cout);
        break;
    case MAGNETIC_MULTIPOLE:
      switch (order()) {
      case 1:
        OnePInts<dcomplex>::OnePDriverLocal<3,false>(
            std::bind(&ComplexGIAOIntEngine::computeGIAOAngularL,
                      std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3, &magAmp[0],options.particle.charge),
            mol, basis, pointers(), op, options, 0);
        break;
      default:
        CErr("Requested operator is NYI in VectorInts.",std::cout);
        break;
      }
      break;

    // GIAO + X2C
    // Calculate rVr and pVAAVp integrals
    case MAGNETIC_4COMP_rVr:
      switch (order()) {
      case 2:
        if (options.finiteWidthNuc)
          OnePInts<dcomplex>::OnePDriverLocal<6,true>(
            [&](libint2::ShellPair& pair, libint2::Shell& sh1,
                libint2::Shell& sh2) -> std::vector<std::vector<dcomplex>> {  
              return ComplexGIAOIntEngine::computeGIAOrVr(
                mol.chargeDist, pair,sh1,sh2,&magAmp[0],mol.retainCNuc(),options.particle.charge); 
              }, mol, basis, pointers(), op, options, 0);
        else
          OnePInts<dcomplex>::OnePDriverLocal<6,true>(
            [&](libint2::ShellPair& pair, libint2::Shell& sh1,
                libint2::Shell& sh2) -> std::vector<std::vector<dcomplex>> {  
              return ComplexGIAOIntEngine::computeGIAOrVr(
                pair,sh1,sh2,&magAmp[0],mol.retainCNuc(),options.particle.charge); 
              }, mol, basis, pointers(), op, options, 0);
        break;
      default:
        CErr("Requested operator is not implemented in VectorInts.");
        break;
      }
      break;
    case MAGNETIC_4COMP_PVrprVP:
      switch (order()) {
      case 2:
        if (options.finiteWidthNuc) 
          OnePInts<dcomplex>::OnePDriverLocal<9,false>(
            [&](libint2::ShellPair& pair, libint2::Shell& sh1,
                libint2::Shell& sh2) -> std::vector<std::vector<dcomplex>> {  
              return ComplexGIAOIntEngine::computeGIAOpVrprVp(
                mol.chargeDist, pair,sh1,sh2,&magAmp[0],mol.retainCNuc(),options.particle.charge); 
              }, mol, basis, pointers(), op, options, 0);
        else
          OnePInts<dcomplex>::OnePDriverLocal<9,false>(
            [&](libint2::ShellPair& pair, libint2::Shell& sh1,
                libint2::Shell& sh2) -> std::vector<std::vector<dcomplex>> {  
              return ComplexGIAOIntEngine::computeGIAOpVrprVp(
                pair,sh1,sh2,&magAmp[0],mol.retainCNuc(),options.particle.charge); 
              }, mol, basis, pointers(), op, options, 0);
        break;
      default:
        CErr("Requested operator is not implemented in VectorInts.");
        break;
      }
      break;
    case MAGNETIC_4COMP_PVrmrVP:
      switch (order()) {
      case 2:
        if (options.finiteWidthNuc)  
          OnePInts<dcomplex>::OnePDriverLocal<9,true>(
            [&](libint2::ShellPair& pair, libint2::Shell& sh1,
                libint2::Shell& sh2) -> std::vector<std::vector<dcomplex>> {  
              return ComplexGIAOIntEngine::computeGIAOpVrmrVp(
                mol.chargeDist, pair,sh1,sh2,&magAmp[0],mol.retainCNuc(),options.particle.charge); 
              }, mol, basis, pointers(), op, options, 0);
        else
          OnePInts<dcomplex>::OnePDriverLocal<9,true>(
            [&](libint2::ShellPair& pair, libint2::Shell& sh1,
                libint2::Shell& sh2) -> std::vector<std::vector<dcomplex>> {  
              return ComplexGIAOIntEngine::computeGIAOpVrmrVp(
                pair,sh1,sh2,&magAmp[0],mol.retainCNuc(),options.particle.charge); 
              }, mol, basis, pointers(), op, options, 0);
        break;
      default:
        CErr("Requested operator is not implemented in VectorInts.");
        break;
      }
      break;

      default:
        break;
    }

  };

  template <>
  void MultipoleInts<dcomplex>::computeAOInts(BasisSet &basis, Molecule &mol,
      EMPerturbation &emPert, OPERATOR op, const HamiltonianOptions &options) {
    if (options.basisType == REAL_GTO)
      CErr("Real GTOs are not allowed in MultipoleInts<dcomplex>",std::cout);
    if (options.basisType == COMPLEX_GTO)
      CErr("Complex GTOs NYI in MultipoleInts<dcomplex>",std::cout);
    //if (options.OneEScalarRelativity or options.OneESpinOrbit)
    //  CErr("Relativistic integrals are implemented in OnePRelInts",std::cout);

    auto magAmp = emPert.getDipoleAmp(Magnetic);

    switch (op) {
    case OVERLAP:
    case KINETIC:
    case NUCLEAR_POTENTIAL:
      CErr("Requested operator is not implemented in MultipoleInts,"
           " it is implemented in OneEInts",std::cout);
      break;
    case ELECTRON_REPULSION:
      CErr("Electron repulsion integrals are not implemented in MultipoleInts,"
           " they are implemented in TwoEInts",std::cout);
      break;
    case LEN_ELECTRIC_MULTIPOLE:
    case VEL_ELECTRIC_MULTIPOLE:
    case MAGNETIC_MULTIPOLE:
      for (VectorInts<dcomplex> &vInts: components_) {
        vInts.computeAOInts(basis, mol, emPert, op, options);
      }
      break;
    default:
      CErr("Requested operator is not implemented in MultipoleInts.");
      break;
    }

  };

  template <>
  void OnePRelInts<dcomplex>::computeAOInts(BasisSet &basis, Molecule &mol,
      EMPerturbation &emPert, OPERATOR op, const HamiltonianOptions &options) {
    if (options.basisType == REAL_GTO)
      CErr("Real GTOs are not allowed in OnePRelInts<dcomplex>",std::cout);
    if (options.basisType == COMPLEX_GTO)
      CErr("Complex GTOs NYI in OnePRelInts<dcomplex>",std::cout);

    auto magAmp = emPert.getDipoleAmp(Magnetic);
    std::vector<dcomplex*> tmp(1, pointer());

    // All Magnetic integrals are placed under potential pointer!
    // Finite Nuc width
    std::vector<dcomplex*> _potential(1, pointer());
    if (options.finiteWidthNuc)
      OnePInts<dcomplex>::OnePDriverLocal<1,true>(
          [&](libint2::ShellPair& pair, libint2::Shell& sh1, 
              libint2::Shell& sh2) -> std::vector<std::vector<dcomplex>> { 
            return ComplexGIAOIntEngine::computeGIAOPotentialV(
                mol.chargeDist,pair,sh1,sh2,&magAmp[0],mol.retainCNuc(),options.particle.charge);
            }, mol, basis, _potential, op, options, 0);
    else
      OnePInts<dcomplex>::OnePDriverLocal<1,true>(
          [&](libint2::ShellPair& pair, libint2::Shell& sh1, 
              libint2::Shell& sh2) -> std::vector<std::vector<dcomplex>> { 
            return ComplexGIAOIntEngine::computeGIAOPotentialV(
                pair,sh1,sh2,&magAmp[0],mol.retainCNuc(),options.particle.charge);
            }, mol, basis, _potential, op, options, 0);

    // Point nuclei is used when chargeDist is empty
    const std::vector<libint2::Shell> &chargeDist = options.finiteWidthNuc ?
        mol.chargeDist : std::vector<libint2::Shell>();

    // pVp Part
    if (options.OneESpinOrbit) {
      if (not hasSpinOrbit())
        CErr("computeAOInts: Requested spin-orbit integrals, "
             "but the OnePRelInts object does not contain spin-orbit components");

      OnePInts<dcomplex>::OnePDriverLocal<3,false>(
            [&](libint2::ShellPair& pair, libint2::Shell& sh1,
                libint2::Shell& sh2) -> std::vector<std::vector<dcomplex>> {
              return ComplexGIAOIntEngine::computeGIAOSL(chargeDist,
                  pair,sh1,sh2,&magAmp[0],mol.retainCNuc(),options.particle.charge);
              }, mol, basis, SOXYZPointers(), op, options, 0);       
    }

    std::vector<dcomplex*> _PVdP(1, scalar().pointer());
    OnePInts<dcomplex>::OnePDriverLocal<1,true>(
          [&](libint2::ShellPair& pair, libint2::Shell& sh1,
              libint2::Shell& sh2) -> std::vector<std::vector<dcomplex>> {
            return ComplexGIAOIntEngine::computeGIAOpVdotp(chargeDist,
                pair,sh1,sh2,&magAmp[0],mol.retainCNuc(),options.particle.charge);
            }, mol, basis, _PVdP, op, options, 0);

  };

  // --------------------- Gradient integrals
  template<>
  void GradInts<OnePInts,dcomplex>::computeAOInts(BasisSet& basis,
    Molecule& mol, EMPerturbation &emPert, OPERATOR op, const HamiltonianOptions& options)
  {

    if (options.basisType == REAL_GTO)
      CErr("Real GTOs are not allowed in OnePInts<dcomplex>",std::cout);
    if (options.basisType == COMPLEX_GTO)
      CErr("Complex GTOs NYI in OnePInts<dcomplex>",std::cout);

    auto magAmp = emPert.getDipoleAmp(Magnetic);

    std::vector<dcomplex*> gradPtrs(3*nAtoms_, nullptr);

    for (auto i = 0; i < 3*nAtoms_; i++) {
      gradPtrs[i] = components_[i]->pointer();
    }

    switch (op) {
    case OVERLAP:
      OnePInts<dcomplex>::OnePDriverLocal<6,true>(
          std::bind(&ComplexGIAOIntEngine::computeGIAOOverlapGradS,
                    std::placeholders::_1, std::placeholders::_2,
                    std::placeholders::_3, &magAmp[0],options.particle.charge),
          mol, basis, gradPtrs, op, options, 1);
      break;
    case KINETIC:
      OnePInts<dcomplex>::OnePDriverLocal<6,true>(
          std::bind(&ComplexGIAOIntEngine::computeGIAOKineticGradT,
                    std::placeholders::_1, std::placeholders::_2,
                    std::placeholders::_3, &magAmp[0],options.particle.charge),
          mol, basis, gradPtrs, op, options, 1);
      break;
    case NUCLEAR_POTENTIAL:
      //OnePInts<dcomplex>::OnePDriverLocal<9,true>(
      //    [&](libint2::ShellPair& pair, libint2::Shell& sh1, 
      //        libint2::Shell& sh2) -> std::vector<std::vector<dcomplex>> { 
      //      return ComplexGIAOIntEngine::computeGIAOPotentialGradV(
      //          mol.chargeDist,pair,sh1,sh2,&magAmp[0],mol.retainCNuc(),options.particle.charge);
      //      }, mol, basis, tmp, op, options, 0) :
      if (options.finiteWidthNuc) {
      std::cerr<<"no finite nuclei grdient yet"<<std::endl;
      } else {
        OnePInts<dcomplex>::OnePDriverLocal<0,true>(
            [&](libint2::ShellPair& pair, libint2::Shell& sh1,
                libint2::Shell& sh2) -> std::vector<std::vector<dcomplex>> {
              return ComplexGIAOIntEngine::computeGIAOPotentialGradV(
                  pair,sh1,sh2,&magAmp[0],mol.retainCNuc(),options.particle.charge);
              }, mol, basis, gradPtrs, op, options, 1);
      }
      break;
    case TAUS0a:
      OnePInts<dcomplex>::OnePDriverLocal<3,true>(
          std::bind(&ComplexGIAOIntEngine::computeGIAOOverlapGradS0a,
                    std::placeholders::_1, std::placeholders::_2,
                    std::placeholders::_3, &magAmp[0],options.particle.charge),
          mol, basis, gradPtrs, op, options, 10);
      break;
    case ELECTRON_REPULSION:
    case LEN_ELECTRIC_MULTIPOLE:
    case MAGNETIC_MULTIPOLE:
      CErr("Requested operator is not implemented... yet",std::cout);
      break;
    case VEL_ELECTRIC_MULTIPOLE:
      CErr("Requested operator is not planned in GradInts-OnePInts.",std::cout);
      break;
    default:
      CErr("Requested operator is not implemented in GradInts-OnePInts.",std::cout); 
      break;
    }
  };

  template<>
  void GradInts<VectorInts, dcomplex>::computeAOInts(BasisSet& basis,
    Molecule& mol, EMPerturbation &emPert, OPERATOR op, const HamiltonianOptions& options) {

    if (options.basisType == REAL_GTO)
      CErr("Real GTOs are not allowed in VectorInts<dcomplex>",std::cout);
    if (options.basisType == COMPLEX_GTO)
      CErr("Complex GTOs NYI in VectorInts<dcomplex>",std::cout);

    auto magAmp = emPert.getDipoleAmp(Magnetic);

    switch (op) {
    case OVERLAP:
    case KINETIC:
    case NUCLEAR_POTENTIAL:
    case ELECTRON_REPULSION:
    case VEL_ELECTRIC_MULTIPOLE:
      CErr("Requested operator is not implemented in GradInts-VectorInts.",std::cout); 
    case MAGNETIC_MULTIPOLE:
      switch (components_[0]->order()) {
      case 1: {
        // Gradient Placer [iGrad, Component]
        std::vector<dcomplex*> gradPtrs(3*nAtoms_*3, nullptr);
        for (auto i = 0; i < 3*nAtoms_; i++) {
          for (auto j = 0; j < 3; j++) { 
            gradPtrs[i*3+j] = components_[i]->pointers()[j];
          }
        }
        OnePInts<dcomplex>::OnePDriverLocal<18,false>(
            std::bind(&ComplexGIAOIntEngine::computeGIAOAngularGradL,
                      std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3, &magAmp[0],options.particle.charge),
            mol, basis, gradPtrs, op, options, 11);
        break;
      }
      default:
        CErr("Requested operator is NYI in VectorInts.",std::cout);
        break;
      }
      break;
    case LEN_ELECTRIC_MULTIPOLE:
      switch (components_[0]->order()) {
      case 2: {
        // Gradient Placer [iGrad, Component]
        std::vector<dcomplex*> gradPtrs(3*nAtoms_*6, nullptr);
        for (auto i = 0; i < 3*nAtoms_; i++) {
          for (auto j = 0; j < 6; j++) { 
            gradPtrs[i*6+j] = components_[i]->pointers()[j];
          }
        }
        OnePInts<dcomplex>::OnePDriverLocal<36,true>(
            std::bind(&ComplexGIAOIntEngine::computeGIAOEQuadrupoleGradE2_len,
                      std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3, &magAmp[0],options.particle.charge),
            mol, basis, gradPtrs, op, options, 12);
        break;
      }
      default:
        CErr("Requested operator is NYI in VectorInts.",std::cout);
        break;
      }
      break;
    default:
      CErr("Requested operator is not implemented in GradInts-VectorInts.",std::cout); 
      break;
    }

  };

  template<>
  void GradInts<MultipoleInts, dcomplex>::computeAOInts(BasisSet& basis,
    Molecule& mol, EMPerturbation &emPert, OPERATOR op, const HamiltonianOptions& options) {
    CErr("Gradients cannot be called directly through MultipoleInts. Use VectorInts instead.");
  };
  template <>
  void GradInts<OnePInts, dcomplex>::computeAOInts(BasisSet&, BasisSet&,
    Molecule&, EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("Two basis gradients not implemented for OnePInts");
  };
  template <>
  void GradInts<VectorInts, dcomplex>::computeAOInts(BasisSet&, BasisSet&,
    Molecule&, EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("Two basis gradients not implemented for VectorInts");
  };
  template <>
  void GradInts<MultipoleInts, dcomplex>::computeAOInts(BasisSet&, BasisSet&,
    Molecule&, EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("Two basis gradients not implemented for MultipoleInts");
  };

  template void Integrals<dcomplex>::computeAOOneP(
      Molecule&, BasisSet&, EMPerturbation&,
      const std::vector<std::pair<OPERATOR,size_t>>&,
      const HamiltonianOptions&);


  template void 
  Integrals<dcomplex>::computeGradInts(
      Molecule&, BasisSet&, EMPerturbation&,
      const std::vector<std::pair<OPERATOR,size_t>>&,
      const HamiltonianOptions&);
      //CErr("Requested operator is not implemented in MultipoleInts.");


}; // namespace ChronusQ

