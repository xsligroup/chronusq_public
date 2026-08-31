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

#include <integrals.hpp>
#include <particleintegrals/inhouseaointegral.hpp>
#include <cqlinalg.hpp>
#include <cqlinalg/svd.hpp>
#include <cqlinalg/blasutil.hpp>
#include <physcon.hpp>
#include <util/matout.hpp>
#include <util/timer.hpp>
#include <util/threads.hpp>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <iomanip>

#include <integrals/impl.hpp>
#include <libcint.hpp>

// Debug directives
//#define _DEBUGORTHO
//#define _DEBUGERI
//#define LIBCINT_FOR_FINITE_NUC

namespace ChronusQ {

  typedef std::vector<libint2::Shell> shell_set;

  /**
   *  \brief A general wrapper for 1-e (2 index) integral evaluation.
   *
   *  Currently computes 1-e integrals using Libint2. Shells sets are
   *  passed in order to be possibly general to the uncontracted basis.
   *  Handles all internal memory allocation including the evaluated matricies
   *  themselves
   *
   *  \param [in] op     Operator for which to calculate the 1-e integrals
   *  \param [in] shells Shell set for the integral evaluation
   *
   *  \returns    A vector of properly allocated pointers which store the
   *              1-e evaluations.
   *
   *  This function returns a vector of pointers as it sometimes makes sense
   *  to evaluate several matricies together if they are inimately related,
   *  namely the length gauge electric multipoles and the overlap.
   *
   *  z.B. op == libint2::Operator::emultipole3
   *
   *  The function will return a vector of 20 pointers in the following order
   *  { overlap, 
   *    dipole_x, dipole_y, dipole_z, 
   *    quadrupole_xx, quadrupole_xy, quadrupole_xz, quadrupole_yy,
   *      quadrupole_yz, quadrupole_zz,
   *    octupole_xxx, octupole_xxy, octupole_xxz, octupole_xyy,
   *      octupole_xyz, octupole_xzz, octupole_yyy, octupole_yyz,
   *      octupole_yzz, octupole_zzz
   *  }
   *
   *  z.B. op == libint2::Operator::kinetic
   *
   *  The function will return a vector of 1 pointer
   *
   *  { kinetic }
   */ 
  template <>
  void OnePInts<dcomplex>::OnePDriverLibint(libint2::Operator op,
      const Molecule &mol, const BasisSet& basis, std::vector<dcomplex*> mats,
      Particle p, size_t deriv, size_t S0a) {
    CErr("Only real GTOs are allowed",std::cout);
  };

  template <>
  void OnePInts<double>::OnePDriverLibint(libint2::Operator op,
      const Molecule &mol, const BasisSet& basis, std::vector<double*> mats, 
      Particle p, size_t deriv, size_t S0a) {

    const shell_set& shells = basis.shells;

    // Determine the number of basis functions for the passed shell set
    size_t NB = std::accumulate(shells.cbegin(),shells.cend(),0,
      [](size_t init, const libint2::Shell &sh) -> size_t {
        return init + sh.size();
      }
    );

    size_t NBSQ = NB*NB;


    // Determine the maximum angular momentum of the passed shell set
    int maxL = std::max_element(shells.cbegin(), shells.cend(),
      [](const libint2::Shell &sh1, const libint2::Shell &sh2){
        return sh1.contr[0].l < sh2.contr[0].l;
      }
    )->contr[0].l;

    // Determine the maximum contraction depth of the passed shell set
    int maxPrim = std::max_element(shells.cbegin(), shells.cend(),
      [](const libint2::Shell &sh1, const libint2::Shell &sh2){
        return sh1.alpha.size() < sh2.alpha.size();
      }
    )->alpha.size();

    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();

    // Create a vector of libint2::Engines for possible threading
    std::vector<libint2::Engine> engines(nthreads);

    // Initialize the first engine for the integral evaluation
    engines[0] = libint2::Engine(op,maxPrim,maxL,deriv);
    engines[0].set_precision(0.0);

    // If engine is K, prescale it by 1/m
    if (op == libint2::Operator::kinetic)
      engines[0].prescale_by(1.0 / p.mass);  

    // If engine is V, define nuclear charges (pseudo molecule is used for NEO)
    if(op == libint2::Operator::nuclear){
      std::vector<std::pair<double,std::array<double,3>>> q;
      for (auto ind : mol.atomsC) // loop over classical atoms
        q.push_back( { -1.0 * p.charge * mol.atoms[ind].nucCharge, mol.atoms[ind].coord } );

      engines[0].set_params(q);
      
    }

    // for multipoles, prescale it by charge
    if (op == libint2::Operator::emultipole1 or op == libint2::Operator::emultipole2 or op == libint2::Operator::emultipole3)
      engines[0].prescale_by(-1.0 * p.charge);

    // Copy over the engines to other threads if need be
    for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];


    std::vector<
      Eigen::Map<
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>
      > 
    > matMaps;
    for( auto i = 0; i < mats.size(); i++ ) {
      std::fill_n(mats[i],NBSQ,0.);
      matMaps.emplace_back(mats[i],NB,NB);
    }


    #pragma omp parallel
    {
      int thread_id = GetThreadID();

      const auto& buf_vec = engines[thread_id].results();
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

        // Compute the integrals       
        engines[thread_id].compute(shells[s1],shells[s2]);

        // adds the iOp result of the engine to the iMat matrix 
        //   For non-gradients, iOp and iMat should be the same
        //   For gradients, they can differ
        auto add_shellset_to_mat = [&](size_t iOp, size_t iMat) {


          // If the integrals were screened, do nothing
          if(buf_vec[iOp] == nullptr) return;

          // std::cout << "iOp: " << iOp << " iMat: " << iMat << std::endl;
          Eigen::Map<
            const Eigen::Matrix<
              double,
              Eigen::Dynamic,
              Eigen::Dynamic,
              Eigen::RowMajor
            >
          > bufMat(buf_vec[iOp],n1,n2);

          size_t _idx = 0;
          for ( auto r_idx = 0; r_idx < n1; r_idx++) {
            for (auto c_idx = 0; c_idx < n2; c_idx++, _idx++) {
              // std::cout << buf_vec[iOp][_idx] << " ";
            }
            // std::cout << std::endl;
          }

          matMaps[iMat].block(bf1_s, bf2_s, n1, n2) += bufMat;

        };

        // Compute S0a matrix
        if (S0a) {
          for (size_t xyz = 0; xyz < 3; xyz ++) {
            // "map" buffer to a const Eigen matrix, and copy it to the
            // corresponding blocks of the result
            Eigen::Map<
              const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
                Eigen::RowMajor>>
              bufMat(buf_vec[3+xyz],n1,n2);
            Eigen::Map<
              const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
                Eigen::RowMajor>>
              bufMat2(buf_vec[xyz],n1,n2);
              //bufMat(buf_vec[xyz],n1,n2);
            if (p.charge > 0) {
              //std::cout << "(" << s1 << "," << s2 << ")" << std::endl;
              matMaps[xyz].block(bf1_s,bf2_s,n1,n2) = bufMat;
              if (s1 != s2)
                matMaps[xyz].block(bf2_s,bf1_s,n2,n1) = bufMat2.transpose();
            }
          }
          //std::cout << "computed S0a" << std::endl;
          continue; 
        }


        // Place integral blocks into their respective matricies
        switch (deriv) {

          case 0:
            for(auto iMat = 0; iMat < buf_vec.size(); iMat++){
              add_shellset_to_mat(iMat, iMat);
            }
            break; // case deriv == 0

          case 1:
            // For gradients, libint returns first the gradients of the
            //   bra/ket, and then the gradients of the operator. We handle
            //   these separately.
            // e.g.
            //   For the (O1s|V|H1s) nuclear attraction gradients in H2O with
            //   atom indices: O:0, H:1, H:2, libint will return 15 derivative
            //   integrals.
            //   (3 cartesian indices * (2 shell centers + 3 nuclear centers))
            //   There are only 9 gradient integrals
            //   (3 cartesian indices * 3 nuclear centers)
            //
            //   The results will be mapped to their respective gradient
            //   integrals by:
            //
            //   | ======================================================== |
            //   |   Engine result    | Gradient integral |  iOps   | iMats |
            //   | ------------------ + ----------------- + ------- + ----- |
            //   | (d/dR0 O1s|V| H1s) | d/dR0 (O1s|V|H1s) | [0,2]   | [0,2] |
            //   | (O1s|V| d/dR1 H1s) | d/dR1 (O1s|V|H1s) | [3,5]   | [3,5] |
            //   | (O1s|d/dR0 V| H1s) | d/dR0 (O1s|V|H1s) | [6,8]   | [0,2] |
            //   | (O1s|d/dR1 V| H1s) | d/dR1 (O1s|V|H1s) | [9,11]  | [3,5] |
            //   | (O1s|d/dR2 V| H1s) | d/dR2 (O1s|V|H1s) | [12,14] | [6,8] |
            //   | ======================================================== |
            //
            // For geometry independent operators, libint will only return 6
            //   derivative integrals. (bra then ket)
            // std::cout << "(" << s1 << "," << s2 << ")" << std::endl;
            // for (auto& x: buf_vec) {
            //   std::cout << "**************************************" << std::endl;
            //   for ( auto i = 0 ; i < n1*n2 ; i++ ) {
            //     std::cout << i << ": " << x[i] << std::endl;
            //   }
            // }

            size_t result_idx = 0;
            
            // First the bra and ket
            for (auto xyz = 0; xyz < 3; xyz++, result_idx++)
              add_shellset_to_mat(result_idx, 3*atom1 + xyz);

            for (auto xyz = 0; xyz < 3; xyz++, result_idx++)
              add_shellset_to_mat(result_idx, 3*atom2 + xyz);

            // Gradient of operator
            if (op == libint2::Operator::nuclear) {
              auto nAtoms = mol.atomsC.size();
              for (auto iAt = 0; iAt < nAtoms; iAt++) {
                for ( auto xyz = 0; xyz < 3; xyz++, result_idx++) {
                  add_shellset_to_mat(result_idx, 3*mol.atomsC[iAt]+ xyz);
                }
              }
            }
            break; // case deriv == 1
        } // switch deriv

      } // Loop over s2 <= s1
      } // Loop over s1

    } // end OpenMP context


    // Symmetrize the matricies 
    if (S0a == 0) {
      for(auto nMat = 0; nMat < matMaps.size(); nMat++) 
        matMaps[nMat] = matMaps[nMat].template selfadjointView<Eigen::Lower>();
    }

  }; // OnePInts::OnePDriver


  /** 
   *  /brief Integral driver with two basis support
   *  Only op=overlap is actually tested, be careful when using other operators!
   *  Returns matrices with dimension (basisB.size(), basisA.size())
   **/
  template <>
  void OnePInts<dcomplex>::OnePDriverLibint(libint2::Operator op,
      const Molecule &mol, const BasisSet& basis_B, const BasisSet& basis_A, std::vector<dcomplex*> mats, 
      Particle p, size_t deriv) {
    CErr("OnePDriverLibint: Only real GTOs are allowed");
  };
  template <>
  void OnePInts<double>::OnePDriverLibint(libint2::Operator op,
      const Molecule &mol, const BasisSet& basis_B, const BasisSet& basis_A, std::vector<double*> mats, 
      Particle p, size_t deriv) {

    const shell_set& shells_A = basis_A.shells;
    const shell_set& shells_B = basis_B.shells;


    // Determine the number of basis functions for the passed shell set
    auto shell_counter = [](size_t init, const libint2::Shell &sh) -> size_t {
        return init + sh.size();
    };

    const size_t NB_A = std::accumulate(shells_A.cbegin(),shells_A.cend(), 0, shell_counter);
    const size_t NB_B = std::accumulate(shells_B.cbegin(),shells_B.cend(), 0, shell_counter);


    // Determine the maximum angular momentum of the passed shell set
    auto max_l_compare = [](const libint2::Shell &sh1, const libint2::Shell &sh2){
        return sh1.contr[0].l < sh2.contr[0].l;
    };

    const int maxL_A = std::max_element(shells_A.cbegin(), shells_A.cend(), max_l_compare)->contr[0].l;
    const int maxL_B = std::max_element(shells_B.cbegin(), shells_B.cend(), max_l_compare)->contr[0].l;
    const int maxL = std::max(maxL_A, maxL_B);

    // Determine the maximum contraction depth of the passed shell set
    auto max_prim_compare = [](const libint2::Shell &sh1, const libint2::Shell &sh2){
        return sh1.alpha.size() < sh2.alpha.size();
    };

    const int maxPrim_A = std::max_element(shells_A.cbegin(), shells_A.cend(), max_prim_compare)->alpha.size();
    const int maxPrim_B = std::max_element(shells_B.cbegin(), shells_B.cend(), max_prim_compare)->alpha.size();
    const int maxPrim = std::max(maxPrim_A, maxPrim_B);

    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();

    // Create a vector of libint2::Engines for possible threading
    std::vector<libint2::Engine> engines(nthreads);

    // Initialize the first engine for the integral evaluation
    engines[0] = libint2::Engine(op,maxPrim,maxL,deriv);
    engines[0].set_precision(0.0);

    // If engine is K, prescale it by 1/m
    if (op == libint2::Operator::kinetic)
      engines[0].prescale_by(1.0 / p.mass);  

    // If engine is V, define nuclear charges (pseudo molecule is used for NEO)
    if(op == libint2::Operator::nuclear){
      std::vector<std::pair<double,std::array<double,3>>> q;
      for (auto ind : mol.atomsC) // loop over classical atoms
        q.push_back( { -1.0 * p.charge * mol.atoms[ind].nucCharge, mol.atoms[ind].coord } );

      engines[0].set_params(q);
      
    }

    // for multipoles, prescale it by charge
    if (op == libint2::Operator::emultipole1 or op == libint2::Operator::emultipole2 or op == libint2::Operator::emultipole3)
      engines[0].prescale_by(-1.0 * p.charge);

    // Copy over the engines to other threads if need be
    for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];


    std::vector<
      Eigen::Map<
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>
      > 
    > matMaps;
    for( auto i = 0; i < mats.size(); i++ ) {
      std::fill_n(mats[i],NB_A*NB_B,0.);
      matMaps.emplace_back(mats[i],NB_B,NB_A);
    }


    #pragma omp parallel
    {
      int thread_id = GetThreadID();

      const auto& buf_vec = engines[thread_id].results();
      size_t n1,n2,atom1,atom2;

      // Loop over unique shell pairs
      // TODO: exploit symmetry
      for(size_t s1(0), bf1_s(0), s12(0); s1 < shells_B.size(); bf1_s+=n1, s1++){ 
        n1 = shells_B[s1].size(); // Size of Shell 1
        atom1 = basis_B.mapSh2Cen[s1]; // Index of atom for Shell 1
      for(size_t s2(0), bf2_s(0); s2 < shells_A.size(); bf2_s+=n2, s2++, s12++) {
        n2 = shells_A[s2].size(); // Size of Shell 2
        atom2 = basis_A.mapSh2Cen[s2]; // Index of atom for Shell 2

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s12 % nthreads != thread_id ) continue;
        #endif

        // Compute the integrals       
        engines[thread_id].compute(shells_B[s1],shells_A[s2]);

        // adds the iOp result of the engine to the iMat matrix 
        //   For non-gradients, iOp and iMat should be the same
        //   For gradients, they can differ
        auto add_shellset_to_mat = [&](size_t iOp, size_t iMat) {


          // If the integrals were screened, do nothing
          if(buf_vec[iOp] == nullptr) return;

          // std::cout << "iOp: " << iOp << " iMat: " << iMat << std::endl;
          Eigen::Map<
            const Eigen::Matrix<
              double,
              Eigen::Dynamic,
              Eigen::Dynamic,
              Eigen::RowMajor
            >
          > bufMat(buf_vec[iOp],n1,n2);

          matMaps[iMat].block(bf1_s, bf2_s, n1, n2) += bufMat;

        };

        // Place integral blocks into their respective matricies
        switch (deriv) {

          case 0:
            for(auto iMat = 0; iMat < buf_vec.size(); iMat++){
              add_shellset_to_mat(iMat, iMat);
            }
            break; // case deriv == 0

          case 1:
            // For gradients, libint returns first the gradients of the
            //   bra/ket, and then the gradients of the operator. We handle
            //   these separately.
            // e.g.
            //   For the (O1s|V|H1s) nuclear attraction gradients in H2O with
            //   atom indices: O:0, H:1, H:2, libint will return 15 derivative
            //   integrals.
            //   (3 cartesian indices * (2 shell centers + 3 nuclear centers))
            //   There are only 9 gradient integrals
            //   (3 cartesian indices * 3 nuclear centers)
            //
            //   The results will be mapped to their respective gradient
            //   integrals by:
            //
            //   | ======================================================== |
            //   |   Engine result    | Gradient integral |  iOps   | iMats |
            //   | ------------------ + ----------------- + ------- + ----- |
            //   | (d/dR0 O1s|V| H1s) | d/dR0 (O1s|V|H1s) | [0,2]   | [0,2] |
            //   | (O1s|V| d/dR1 H1s) | d/dR1 (O1s|V|H1s) | [3,5]   | [3,5] |
            //   | (O1s|d/dR0 V| H1s) | d/dR0 (O1s|V|H1s) | [6,8]   | [0,2] |
            //   | (O1s|d/dR1 V| H1s) | d/dR1 (O1s|V|H1s) | [9,11]  | [3,5] |
            //   | (O1s|d/dR2 V| H1s) | d/dR2 (O1s|V|H1s) | [12,14] | [6,8] |
            //   | ======================================================== |
            //
            // For geometry independent operators, libint will only return 6
            //   derivative integrals. (bra then ket)
            // std::cout << "(" << s1 << "," << s2 << ")" << std::endl;
            // for (auto& x: buf_vec) {
            //   std::cout << "**************************************" << std::endl;
            //   for ( auto i = 0 ; i < n1*n2 ; i++ ) {
            //     std::cout << i << ": " << x[i] << std::endl;
            //   }
            // }

            size_t result_idx = 0;
            
            // First the bra and ket
            for (auto xyz = 0; xyz < 3; xyz++, result_idx++)
              add_shellset_to_mat(result_idx, 3*atom1 + xyz);

            for (auto xyz = 0; xyz < 3; xyz++, result_idx++)
              add_shellset_to_mat(result_idx, 3*atom2 + xyz);

            // Gradient of operator
            if (op == libint2::Operator::nuclear) {
              auto nAtoms = mol.atomsC.size();
              for (auto iAt = 0; iAt < nAtoms; iAt++) {
                for ( auto xyz = 0; xyz < 3; xyz++, result_idx++) {
                  add_shellset_to_mat(result_idx, 3*mol.atomsC[iAt]+ xyz);
                }
              }
            }
            break; // case deriv == 1
        } // switch deriv

      } // Loop over s2 <= s1
      } // Loop over s1

    } // end OpenMP context

  }; // OnePInts::OnePDriver


  /**
   *  \brief A general wrapper for 1-e (2 index) integral evaluation by libcint
   *         Currently support overlap, kinetic, and bare nuclear potential
   *
   *
   *  \param [in] op     Operator for which to calculate the 1-e integrals
   *
   */
  template <>
  void OnePInts<dcomplex>::OnePDriverLibcint(OPERATOR, const Molecule&,
      const BasisSet&, const HamiltonianOptions&) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void OnePInts<double>::OnePDriverLibcint(OPERATOR op,
      const Molecule &molecule_, const BasisSet &originalBasisSet,
      const HamiltonianOptions &options) {

    if (originalBasisSet.forceCart)
      CErr("Libcint + cartesian GTO NYI.");

    BasisSet basisSet_ = originalBasisSet.groupGeneralContractionBasis();

    size_t buffSize = std::max_element(basisSet_.shells.begin(),
                                       basisSet_.shells.end(),
                                       [](libint2::Shell &a, libint2::Shell &b) {
                                         return a.size() < b.size();
                                       })->size();
    buffSize *= buffSize;

    int nAtoms = molecule_.nAtoms;
    int nShells = basisSet_.nShell;

    // ATM_SLOTS = 6; BAS_SLOTS = 8;
    int *atm = CQMemManager::get().template malloc<int>(nAtoms * ATM_SLOTS);
    int *bas = CQMemManager::get().template malloc<int>(nShells * BAS_SLOTS);
    double *env = CQMemManager::get().template malloc<double>(basisSet_.getLibcintEnvLength(molecule_));


    basisSet_.setLibcintEnv(molecule_, atm, bas, env, options.finiteWidthNuc);

    auto intFunc = &int1e_kin_sph;
    switch (op) {
    case OVERLAP:
      intFunc = &int1e_ovlp_sph;
      break;
    case KINETIC:
      intFunc = &int1e_kin_sph;
      break;
    case NUCLEAR_POTENTIAL:
      intFunc = &int1e_nuc_sph;
      break;
    default:
      CErr("Requested OPERATOR type not implemented in OnePDriverLibcint");
    }


    size_t cache_size = 0;
    for (int i = 0; i < nShells; i++) {
      size_t n;
      int shls[2]{i,i};
      n = intFunc(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
      cache_size = std::max(cache_size, n);
    }

    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();

    double *buffAll = CQMemManager::get().template malloc<double>(buffSize*nthreads);
    double *cacheAll = CQMemManager::get().template malloc<double>(cache_size*nthreads);

    clear();
    Eigen::Map<
      Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>
    > matMap(pointer(), NB, NB);

    #pragma omp parallel
    {
      int thread_id = GetThreadID();
      size_t n1,n2;
      int shls[2];
      double *buff = buffAll + buffSize * thread_id;
      double *cache = cacheAll + cache_size * thread_id;

      // Loop over unique shell pairs
      for(size_t s1(0), bf1_s(0), s12(0); s1 < basisSet_.nShell; bf1_s+=n1, s1++){
        n1 = basisSet_.shells[s1].size(); // Size of Shell 1
      for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++, s12++) {
        n2 = basisSet_.shells[s2].size(); // Size of Shell 2

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s12 % nthreads != thread_id ) continue;
        #endif

        shls[0] = int(s2);
        shls[1] = int(s1);

        // Compute the integrals
        if(intFunc(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;

        // Place integral blocks into their respective matricies
        Eigen::Map<
          const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
            Eigen::RowMajor>
        > bufMat(buff, n1, n2);

        matMap.block(bf1_s,bf2_s,n1,n2) = bufMat.template cast<double>();

      } // Loop over s2 <= s1
      } // Loop over s1

    } // end OpenMP context

    CQMemManager::get().free(cacheAll, buffAll, env, bas, atm);


    // Symmetrize the matricies
    matMap = matMap.template selfadjointView<Eigen::Lower>();


    // If engine is K, scale it by 1/m
    if (op == KINETIC)
      matrix() *= 1.0 / options.particle.mass;

    // If engine is V, scale it by charge
    if(op == NUCLEAR_POTENTIAL)
      matrix() *= -1.0 * options.particle.charge;

  }; // OnePInts::OnePDriverLibcint

  // Evaluate Cartesian components of (r x p) using Libcint.
  // Returned components follow the same real-valued convention used by computeAngularL.
  static void OnePDriverLibcintAngularL(const Molecule &molecule_,
      const BasisSet &originalBasisSet, const HamiltonianOptions &options,
      std::vector<double*> mats, size_t NB) {

    if (mats.size() != 3)
      CErr("OnePDriverLibcintAngularL requires 3 output matrices.");

    if (originalBasisSet.forceCart)
      CErr("Libcint + cartesian GTO NYI.");

    BasisSet basisSet_ = originalBasisSet.groupGeneralContractionBasis();

    size_t buffSize = std::max_element(basisSet_.shells.begin(),
                                       basisSet_.shells.end(),
                                       [](libint2::Shell &a, libint2::Shell &b) {
                                         return a.size() < b.size();
                                       })->size();
    buffSize *= buffSize * 3;

    int nAtoms = molecule_.nAtoms;
    int nShells = basisSet_.nShell;

    int *atm = CQMemManager::get().template malloc<int>(nAtoms * ATM_SLOTS);
    int *bas = CQMemManager::get().template malloc<int>(nShells * BAS_SLOTS);
    double *env = CQMemManager::get().template malloc<double>(basisSet_.getLibcintEnvLength(molecule_));

    basisSet_.setLibcintEnv(molecule_, atm, bas, env, options.finiteWidthNuc);

    int nthreads = GetNumThreads();
    double *buffAll = CQMemManager::get().template malloc<double>(buffSize*nthreads);

    std::vector<
      Eigen::Map<
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>
      >
    > matMaps;
    for (auto i = 0; i < mats.size(); i++) {
      std::fill_n(mats[i], NB*NB, 0.0);
      matMaps.emplace_back(mats[i], NB, NB);
    }

    #pragma omp parallel
    {
      int thread_id = GetThreadID();
      size_t n1,n2;
      int shls[2];
      double *buff = buffAll + buffSize * thread_id;

      for (size_t s1(0), bf1_s(0), s12(0); s1 < basisSet_.nShell; bf1_s += n1, s1++) {
        n1 = basisSet_.shells[s1].size();
        for (size_t s2(0), bf2_s(0); s2 <= s1; bf2_s += n2, s2++, s12++) {
          n2 = basisSet_.shells[s2].size();

          #ifdef _OPENMP
          if (s12 % nthreads != thread_id) continue;
          #endif

          shls[0] = int(s2);
          shls[1] = int(s1);

          if (cint1e_cg_irxp_sph(buff, shls, atm, nAtoms, bas, nShells, env) == 0)
            continue;

          for (size_t iXYZ = 0; iXYZ < 3; iXYZ++) {
            Eigen::Map<
              const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor>
            > bufMat(buff + iXYZ*n1*n2, n1, n2);
            matMaps[iXYZ].block(bf1_s,bf2_s,n1,n2) = bufMat;
          }
        }
      }
    }

    CQMemManager::get().free(buffAll, env, bas, atm);

    for (auto iXYZ = 0; iXYZ < 3; iXYZ++) {
      for (auto i = 0ul; i < NB; i++) {
        matMaps[iXYZ](i,i) = 0.0;
        for (auto j = i + 1; j < NB; j++)
          matMaps[iXYZ](i,j) = -matMaps[iXYZ](j,i);
      }
    }
  }


  /**
   *  \brief Computes relativistic nuclear potential integrals,
   *         including V, pVp, and pxVp by libcint
   *
   */
  template <>
  void OnePRelInts<dcomplex>::OnePRelDriverLibcint(const Molecule&,
      const BasisSet&, const HamiltonianOptions&) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void OnePRelInts<double>::OnePRelDriverLibcint(const Molecule &molecule_,
      const BasisSet &originalBasisSet, const HamiltonianOptions &options) {

    if (originalBasisSet.forceCart)
      CErr("Libcint + cartesian GTO NYI.");

    BasisSet basisSet_ = originalBasisSet.groupGeneralContractionBasis();

    size_t buffSize = std::max_element(basisSet_.shells.begin(),
                                       basisSet_.shells.end(),
                                       [](libint2::Shell &a, libint2::Shell &b) {
                                         return a.size() < b.size();
                                       })->size();
    buffSize *= buffSize;

    int nAtoms = molecule_.nAtoms;
    int nShells = basisSet_.nShell;

    // ATM_SLOTS = 6; BAS_SLOTS = 8;
    int *atm = CQMemManager::get().template malloc<int>(nAtoms * ATM_SLOTS);
    int *bas = CQMemManager::get().template malloc<int>(nShells * BAS_SLOTS);
    double *env = CQMemManager::get().template malloc<double>(basisSet_.getLibcintEnvLength(molecule_));


    basisSet_.setLibcintEnv(molecule_, atm, bas, env, options.finiteWidthNuc);

    clear();

    Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > VMap(pointer(), NB, NB);
    Eigen::Map< Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> > pVpMap(scalar().pointer(), NB, NB);

    std::vector< Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> > pxVpMaps;
    if (options.OneESpinOrbit) {

      if (not hasSpinOrbit())
        CErr("OnePRelDriverLibcint: Requested spin-orbit integrals, "
             "but the OnePRelInts object does not contain spin-orbit components");

      buffSize *= 3;
      pxVpMaps.reserve(3);
      for (double *ptr : SOXYZPointers())
        pxVpMaps.emplace_back(ptr, NB, NB);
    }

    size_t cache_size = 0;
    for (int i = 0; i < nShells; i++) {
      size_t n;
      int shls[2]{i,i};
      n = int1e_nuc_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
      cache_size = std::max(cache_size, n);
      n = int1e_pnucp_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
      cache_size = std::max(cache_size, n);
      if (options.OneESpinOrbit) {
        n = int1e_pnucxp_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
      }
    }

    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();

    double *buffAll = CQMemManager::get().template malloc<double>(buffSize*nthreads);
    double *cacheAll = CQMemManager::get().template malloc<double>(cache_size*nthreads);


    #pragma omp parallel
    {
      int thread_id = GetThreadID();
      size_t n1,n2;
      int shls[2];
      double *buff = buffAll + buffSize * thread_id;
      double *cache = cacheAll + cache_size * thread_id;

      // Loop over unique shell pairs
      for(size_t s1(0), bf1_s(0), s12(0); s1 < basisSet_.nShell; bf1_s+=n1, s1++){
        n1 = basisSet_.shells[s1].size(); // Size of Shell 1
      for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++, s12++) {
        n2 = basisSet_.shells[s2].size(); // Size of Shell 2

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s12 % nthreads != thread_id ) continue;
        #endif

        // Assign shells, note row-major in libcint
        shls[0] = int(s2);
        shls[1] = int(s1);

        // Place integral blocks into their respective matricies
        Eigen::Map<
          const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
            Eigen::RowMajor>
        > bufMat(buff, n1, n2);

        // Compute the bare potential integrals
        if(int1e_nuc_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)) {
          VMap.block(bf1_s,bf2_s,n1,n2) = bufMat.template cast<double>();
        }

        // Compute the pVp integrals
        if(int1e_pnucp_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)) {
          pVpMap.block(bf1_s,bf2_s,n1,n2) = bufMat.template cast<double>();
        }

        // Compute the pxVp integrals
        if(options.OneESpinOrbit and
           int1e_pnucxp_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)) {
          size_t n1n2 = n1*n2;
          // Place integral blocks into their respective matricies
          for(auto iMat = 0; iMat < 3; iMat++){
            Eigen::Map<
              const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
                Eigen::RowMajor>>
              bufMat(buff + iMat * n1n2,n1,n2);

            // Negetive sign reflects row-major to column-major switch
            pxVpMaps[iMat].block(bf1_s,bf2_s,n1,n2) = -bufMat.template cast<double>();
            pxVpMaps[iMat].block(bf2_s,bf1_s,n2,n1) = bufMat.transpose().template cast<double>();
          }
        }

      } // Loop over s2 <= s1
      } // Loop over s1

    } // end OpenMP context

    CQMemManager::get().free(cacheAll, buffAll, env, bas, atm);


    // Symmetrize the matricies
    VMap = VMap.template selfadjointView<Eigen::Lower>();
    pVpMap = pVpMap.template selfadjointView<Eigen::Lower>();

    // scale it by charge
    matrix() *= -1.0 * options.particle.charge;
    smallComponent_ *= -1.0 * options.particle.charge;

  }; // OnePRelInts::OnePRelDriverLibcint


  template <>
  template <size_t NOPER, bool SYMM, typename F>
  void OnePInts<double>::OnePDriverLocal(
      const F &obFunc, shell_set& shells, std::vector<double*> mats) {

    // Determine the number of basis functions for the passed shell set
    size_t NB = std::accumulate(shells.begin(),shells.end(),0,
      [](size_t init, libint2::Shell &sh) -> size_t {
        return init + sh.size();
      }
    );

    size_t NBSQ = NB*NB;

    std::vector<
      Eigen::Map<
        Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor> 
      > 
    > matMaps;

    for( auto i = 0; i < mats.size(); i++ ) {
      std::fill_n(mats[i],NBSQ,0.);  
      matMaps.emplace_back(mats[i],NB,NB);
    }

//    if(basisType == REAL_GTO)
      // pre compute all the shellpair data
//      auto pair_to_use = genShellPairs(shells,std::log(std::numeric_limits<double>::lowest()));
    
    // Loop over unique shell pairs



    auto start  = tick();

    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();
    #pragma omp parallel  
    {
    int thread_id = GetThreadID();

    size_t n1,n2;
    for(size_t s1(0), bf1_s(0), s12(0); s1 < shells.size(); bf1_s+=n1, s1++){ 
      n1 = shells[s1].size(); // Size of Shell 1
    for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++, s12++) {
      n2 = shells[s2].size(); // Size of Shell 2


        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s12 % nthreads != thread_id ) continue;
        #endif


      libint2::ShellPair pair_to_use;
      pair_to_use.init(shells[s1],shells[s2],-1000);

      auto buff = obFunc(pair_to_use, shells[s1],shells[s2]);

/*
#pragma omp critical
{
      std::cout<<"s1= "<<s1<<" s2 = "<<s2<<std::endl;
      for ( int elements = 0 ; elements < buff[0].size() ; elements++ ) {
        std::cout<<"buff["<<elements<<"]= "<<buff[0][elements]<<std::endl;
      } 
}  // critical 
*/
      assert(buff.size() == NOPER);

 /*     
      // Place integral blocks into their respective matricies
      for ( int iidx = 0 ; iidx < n1 ; iidx++ ) {
        for ( int jidx = 0 ; jidx < n2 ; jidx++ ) {
          for ( int icomp = 0 ; icomp < NOPER ; icomp++ ) {
            mats[icomp][(iidx+bf1_s)*NB+bf2_s+jidx] = buff[icomp][iidx*n2+jidx];
            std::cout<<"iidx+bf1_s= "<<iidx+bf1_s<<"  bf2_s+jidx= "<<bf2_s+jidx<<" elements = "<<iidx*n2+jidx<<" value "<<buff[icomp][iidx*n2+jidx]<<mats[icomp][(iidx+bf1_s)*NB+bf2_s+jidx]<<std::endl;
          }
        }
      }
*/

      for(auto iMat = 0; iMat < buff.size(); iMat++){
        Eigen::Map<
          const Eigen::Matrix<
            double,
            Eigen::Dynamic,Eigen::Dynamic,  
            Eigen::RowMajor>>
          bufMat(&buff[iMat][0],n1,n2);

        matMaps[iMat].block(bf1_s,bf2_s,n1,n2) = bufMat.template cast<double>();
      }


    } // Loop over s2 <= s1
    } // Loop over s1
    }   // omp
 
    double end = tock(start);
    //std::cout<<"onee driver time= "<<end<<std::endl;

    // Symmetrize the matricies 
    // XXX: USES EIGEN
    // FIXME: not SYMM -> creates a temporary
    for(auto nMat = 0; nMat < matMaps.size(); nMat++) {
      if(SYMM) matMaps[nMat] = matMaps[nMat].template selfadjointView<Eigen::Lower>();
      else {
        for(auto i = 0  ; i < NB; ++i)
        for(auto j = i+1; j < NB; ++j)
          matMaps[nMat](i,j) = - matMaps[nMat](j,i);
      }
    }

  }; // OnePInts::OnePDriverLocal

  template <>
  void OnePInts<double>::computeAOInts(BasisSet &basis, Molecule &mol,
      EMPerturbation&, OPERATOR op, const HamiltonianOptions &options) {

    if (options.basisType != REAL_GTO)
      CErr("Only Real GTOs are allowed in OnePInts<double>",std::cout);
    if (op == NUCLEAR_POTENTIAL and
        (options.OneEScalarRelativity or options.OneESpinOrbit))
      CErr("Relativistic nuclear potential is not implemented in OnePInts,"
           " they are implemented in OnePRelInts",std::cout);


    if (options.Libcint) {
      bool useLocalAngularImpl = false;
      switch (op) {
      case OVERLAP:
      case KINETIC:
      case NUCLEAR_POTENTIAL:
        OnePDriverLibcint(op, mol, basis, options);
        return;
      case SPIN_DOT_ANGULAR:
      case SPIN_VECTOR:
        // These operators are currently evaluated by local/libint-based routines.
        useLocalAngularImpl = true;
        break;
      case ANGULAR_MOMENTUM_SQUARED:
      case TOTAL_ANGULAR_MOMENTUM_SQUARED:
        CErr("Squared angular one-electron operators are deprecated; use non-squared operators and expectation-value contractions.", std::cout);
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
      if (!useLocalAngularImpl) return;
    }

    std::vector<double*> tmp(1, pointer());

    switch (op) {
    case OVERLAP:
      OnePDriverLibint(libint2::Operator::overlap,mol,basis,tmp,options.particle);
      break;
    case KINETIC:
      OnePDriverLibint(libint2::Operator::kinetic,mol,basis,tmp,options.particle);
      //output(std::cout,"",true);
      break;
    case NUCLEAR_POTENTIAL:
      if (options.finiteWidthNuc) {
#ifdef LIBCINT_FOR_FINITE_NUC
        std::cout << "Using Libcint for finite nuclear integrals." << std::endl;
        OnePDriverLibcint(op, mol, basis, options);
#else
        OnePDriverLocal<1,true>(
            [&](libint2::ShellPair& pair, libint2::Shell& sh1,
                libint2::Shell& sh2) -> std::vector<std::vector<double>> {
              return RealGTOIntEngine::computePotentialV(mol.chargeDist,
                  pair,sh1,sh2,mol);
              }, basis.shells,tmp);
#endif
      }
      else
        OnePDriverLibint(libint2::Operator::nuclear,mol,basis,tmp,options.particle);
      break;
    case ELECTRON_REPULSION:
      CErr("Electron repulsion integrals are not implemented in OnePInts,"
           " they are implemented in TwoEInts",std::cout);
      break;
    case ANGULAR_MOMENTUM_SQUARED:
      CErr("ANGULAR_MOMENTUM_SQUARED AO integral build disabled: use non-squared operators and expectation-value contractions.", std::cout);
      break;
    case SPIN_DOT_ANGULAR: {
      // Build S·L properly: SL = Sx*Lx + Sy*Ly + Sz*Lz
      // Spatial S components = overlap * 1/2
      size_t NB = this->nBasis();
      double *Lx = nullptr; double *Ly = nullptr; double *Lz = nullptr;
      double *Sx = nullptr; double *Sy = nullptr; double *Sz = nullptr;
      try {
        Lx = CQMemManager::get().malloc<double>(NB*NB);
        Ly = CQMemManager::get().malloc<double>(NB*NB);
        Lz = CQMemManager::get().malloc<double>(NB*NB);
        Sx = CQMemManager::get().malloc<double>(NB*NB);
        Sy = CQMemManager::get().malloc<double>(NB*NB);
        Sz = CQMemManager::get().malloc<double>(NB*NB);
      } catch(...) { CErr("Insufficient memory for SPIN_DOT_ANGULAR temporary", std::cout); }

      std::vector<double*> tmpL = {Lx,Ly,Lz};
      if (options.Libcint) {
        OnePDriverLibcintAngularL(mol, basis, options, tmpL, NB);
        for (auto *component : tmpL)
          for (size_t i = 0; i < NB * NB; ++i)
            component[i] *= -1.0;
      } else
        OnePInts<double>::OnePDriverLocal<3,false>(
            std::bind(&RealGTOIntEngine::computeAngularL,
                      std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3), basis.shells, tmpL);

      std::vector<double*> tmpS = {Sx,Sy,Sz};
      // Use the libint-based computeSL implementation to build S·L components
      // Pass an empty nuclear shell vector (non-finite-width case)
      std::vector<libint2::Shell> empty_chargeDist;
      OnePInts<double>::OnePDriverLocal<3,false>(
          [&](libint2::ShellPair& pair, libint2::Shell& sh1, libint2::Shell& sh2){
            return RealGTOIntEngine::computeSL(empty_chargeDist, pair, sh1, sh2, mol);
          }, basis.shells, tmpS);

      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> LxM(Lx,NB,NB);
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> LyM(Ly,NB,NB);
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> LzM(Lz,NB,NB);
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> SxM(Sx,NB,NB);
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> SyM(Sy,NB,NB);
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> SzM(Sz,NB,NB);
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> outM(pointer(),NB,NB);

      // Sum libint-computed SL components into the scalar S·L AO matrix
      outM = SxM + SyM + SzM;

      CQMemManager::get().free(Lx,Ly,Lz,Sx,Sy,Sz);
    }
      break;
    case SPIN_VECTOR:
      // Spin operator acts only in spin space; spatial part is overlap scaled by 1/2
      OnePDriverLocal<3,true>(
          [&](libint2::ShellPair& pair, libint2::Shell& sh1, libint2::Shell& sh2){
            auto ov = RealGTOIntEngine::computeOverlapS(pair, sh1, sh2);
            std::vector<std::vector<double>> out(3);
            out[0].resize(ov[0].size()); out[1].resize(ov[0].size()); out[2].resize(ov[0].size());
            for(size_t i=0;i<ov[0].size();++i){ out[0][i]=ov[0][i]*0.5; out[1][i]=ov[0][i]*0.5; out[2][i]=ov[0][i]*0.5; }
            return out;
          }, basis.shells, tmp);
      break;
    case TOTAL_ANGULAR_MOMENTUM_SQUARED:
      CErr("TOTAL_ANGULAR_MOMENTUM_SQUARED AO integral build disabled: use non-squared operators and expectation-value contractions.", std::cout);
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
  void VectorInts<double>::computeAOInts(BasisSet &basis, Molecule &mol,
      EMPerturbation&, OPERATOR op, const HamiltonianOptions &options) {
    if (options.basisType != REAL_GTO)
      CErr("Only Real GTOs are allowed in VectorInts<double>",std::cout);
    if ((options.OneEScalarRelativity or options.OneESpinOrbit) and
        op != SPIN_VECTOR and
        op != ANGULAR_MOMENTUM_VECTOR and
        op != TOTAL_ANGULAR_MOMENTUM_VECTOR)
      CErr("Relativistic multipole integrals are implemented in OnePRelInts",std::cout);

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
      CErr("Len Electric multipole integrals are not implemented in VectorInts,"
           " they are implemented in MultipoleInts",std::cout);
      break;
    case VEL_ELECTRIC_MULTIPOLE:
      switch (order()) {
      case 1:
        OnePInts<double>::OnePDriverLocal<3,false>(
            std::bind(&RealGTOIntEngine::computeEDipoleE1_vel,
                      std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3),
            basis.shells, pointers());
        break;
      case 2:
        OnePInts<double>::OnePDriverLocal<6,false>(
            std::bind(&RealGTOIntEngine::computeEQuadrupoleE2_vel,
                      std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3),
            basis.shells, pointers());
        break;
      case 3:
        OnePInts<double>::OnePDriverLocal<10,false>(
            std::bind(&RealGTOIntEngine::computeEOctupoleE3_vel,
                      std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3),
            basis.shells, pointers());
        break;
      default:
        CErr("Requested operator is NYI in VectorInts.",std::cout);
        break;
      }
      break;
    case MAGNETIC_MULTIPOLE:
      switch (order()) {
      case 1:
        OnePInts<double>::OnePDriverLocal<3,false>(
            std::bind(&RealGTOIntEngine::computeAngularL,
                      std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3),
            basis.shells, pointers());
        break;
      case 2:
        OnePInts<double>::OnePDriverLocal<9,false>(
            std::bind(&RealGTOIntEngine::computeMQuadrupoleM2_vel,
                      std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3),
            basis.shells, pointers());
        break;
      default:
        CErr("Requested operator is NYI in VectorInts.",std::cout);
        break;
      }
      break;
    case SPIN_VECTOR:
      OnePInts<double>::OnePDriverLocal<3,false>(
          [&](libint2::ShellPair& pair, libint2::Shell& sh1, libint2::Shell& sh2) {
            auto ov = RealGTOIntEngine::computeOverlapS(pair, sh1, sh2);
            std::vector<std::vector<double>> out(3, ov[0]);
            for (auto &comp : out)
              for (auto &val : comp)
                val *= 0.5;
            return out;
          }, basis.shells, pointers());
      break;
    case ANGULAR_MOMENTUM_VECTOR:
      if (options.Libcint){
        OnePDriverLibcintAngularL(mol, basis, options, pointers(), NB);
        for (auto *component : pointers())
          for (size_t i = 0; i < NB * NB; ++i)
            component[i] *= -1.0;
      } else {
        OnePInts<double>::OnePDriverLocal<3,false>(
            std::bind(&RealGTOIntEngine::computeAngularL,
                      std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3),
            basis.shells, pointers());
      }
      break;
    case TOTAL_ANGULAR_MOMENTUM_VECTOR: {
      // J = L + S where S (spatial part) = overlap * 1/2
      // First compute L into the output pointers
      if (options.Libcint)  {
        OnePDriverLibcintAngularL(mol, basis, options, pointers(), NB);
        for (auto *component : pointers())
          for (size_t i = 0; i < NB * NB; ++i)
            component[i] *= -1.0;
      } else {
        OnePInts<double>::OnePDriverLocal<3,false>(
            std::bind(&RealGTOIntEngine::computeAngularL,
                      std::placeholders::_1, std::placeholders::_2,
                      std::placeholders::_3),
            basis.shells, pointers());
      }

      // Now compute spatial S components (overlap scaled by 1/2) into temps and add
      size_t NB = this->nBasis();
      double *Sx = nullptr; double *Sy = nullptr; double *Sz = nullptr;
      try {
        Sx = CQMemManager::get().malloc<double>(NB*NB);
        Sy = CQMemManager::get().malloc<double>(NB*NB);
        Sz = CQMemManager::get().malloc<double>(NB*NB);
      } catch(...) { CErr("Insufficient memory for TOTAL_ANGULAR_MOMENTUM_VECTOR temporary", std::cout); }

      std::vector<double*> tmpS = {Sx,Sy,Sz};
      OnePInts<double>::OnePDriverLocal<3,false>(
          [&](libint2::ShellPair& pair, libint2::Shell& sh1, libint2::Shell& sh2){
            auto ov = RealGTOIntEngine::computeOverlapS(pair, sh1, sh2);
            std::vector<std::vector<double>> out(3);
            out[0].resize(ov[0].size()); out[1].resize(ov[0].size()); out[2].resize(ov[0].size());
            for(size_t i=0;i<ov[0].size();++i){ out[0][i]=ov[0][i]*0.5; out[1][i]=ov[0][i]*0.5; out[2][i]=ov[0][i]*0.5; }
            return out;
          }, basis.shells, tmpS);

      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> SxM(Sx,NB,NB);
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> SyM(Sy,NB,NB);
      Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> SzM(Sz,NB,NB);

      // Add S components into existing pointers (which currently hold L)
      for (size_t comp = 0; comp < 3; ++comp) {
        Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>>
          outM(pointers()[comp], NB, NB);
        if (comp == 0) outM += SxM;
        if (comp == 1) outM += SyM;
        if (comp == 2) outM += SzM;
      }

      CQMemManager::get().free(Sx,Sy,Sz);
    }
      break;
    default:
      CErr("Requested operator is not implemented in VectorInts.");
      break;
    }

  };

  /**
   *  \brief Computes relativistic dipole integrals,
   *         using int1e_sprsp by libcint
   */
  template <>
  void MultipoleInts<dcomplex>::MultipoleRelDriverLibcint(const Molecule&,
      const BasisSet&, const HamiltonianOptions&) {
    CErr("Only real GTOs are allowed",std::cout);
  };
  template <>
  void MultipoleInts<double>::MultipoleRelDriverLibcint(const Molecule &molecule_,
      const BasisSet &originalBasisSet, const HamiltonianOptions &options) {

    if (originalBasisSet.forceCart)
      CErr("Libcint + cartesian GTO NYI.");

    if (this->highOrder_ > 1)
      CErr("MultipoleRelDriverLibcint only handles 1st order Multipole");

    BasisSet basisSet_ = originalBasisSet.groupGeneralContractionBasis();

    size_t buffSize = std::max_element(basisSet_.shells.begin(),
                                       basisSet_.shells.end(),
                                       [](libint2::Shell &a, libint2::Shell &b) {
                                         return a.size() < b.size();
                                       })->size();
    buffSize *= buffSize;
    
    int nAtoms = molecule_.nAtoms;
    int nShells = basisSet_.nShell;

    // ATM_SLOTS = 6; BAS_SLOTS = 8;
    int *atm = CQMemManager::get().template malloc<int>(nAtoms * ATM_SLOTS);
    int *bas = CQMemManager::get().template malloc<int>(nShells * BAS_SLOTS);
    double *env = CQMemManager::get().template malloc<double>(basisSet_.getLibcintEnvLength(molecule_));
    
    basisSet_.setLibcintEnv(molecule_, atm, bas, env, options.finiteWidthNuc);

    clear();
    
    //Container for LL dipoles
    std::vector< Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> > LLMaps;
    LLMaps.reserve(3);
    for(double* ptr: this->dipolePointers()) LLMaps.emplace_back(ptr, NB, NB);

    //Container for SS dipoles
    std::vector< Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>> > SSMaps;
    buffSize *= 12;
    SSMaps.reserve(12);
    for (size_t i = 0; i < 3; i++){
      if(auto di = std::dynamic_pointer_cast<OnePRelInts<double>>((*this)[i])){
        for (double *ptr : di->SOXYZPointers())   SSMaps.emplace_back(ptr, NB, NB);
        SSMaps.emplace_back(di->scalar().pointer(), NB, NB);
      }else{
        CErr("OnePInts Stored in MultipoleInts Not Converted to OnePRelInts");
      }
    }
    

    size_t cache_size = 0;
    for (int i = 0; i < nShells; i++) {
      size_t n;
      int shls[2]{i,i};
      n = int1e_r_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
      cache_size = std::max(cache_size, n);
      n = int1e_sprsp_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
      cache_size = std::max(cache_size, n);
    }
    
    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();

    double *buffAll = CQMemManager::get().template malloc<double>(buffSize*nthreads);
    double *cacheAll = CQMemManager::get().template malloc<double>(cache_size*nthreads);

    #pragma omp parallel
    {
      int thread_id = GetThreadID();
      size_t n1,n2;
      int shls[2];
      double *buff = buffAll + buffSize * thread_id;
      double *cache = cacheAll + cache_size * thread_id;

      // Loop over unique shell pairs
      for(size_t s1(0), bf1_s(0), s12(0); s1 < basisSet_.nShell; bf1_s+=n1, s1++){
        n1 = basisSet_.shells[s1].size(); // Size of Shell 1
      for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++, s12++) {
        n2 = basisSet_.shells[s2].size(); // Size of Shell 2

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s12 % nthreads != thread_id ) continue;
        #endif

        // Assign shells, note row-major in libcint
        shls[0] = int(s2);
        shls[1] = int(s1);

        // Place integral blocks into their respective matricies
        Eigen::Map<
          const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
            Eigen::RowMajor>
        > bufMat(buff, n1, n2);

        // Compute LL dipole integrals
        if(int1e_r_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)) {
          size_t n1n2 = n1*n2;
          // Place x,y,z integral blocks into their respective matricies
          for(auto iMat = 0; iMat < 3; iMat++){
            Eigen::Map<
              const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
                Eigen::RowMajor>>
              bufMat(buff + iMat * n1n2,n1,n2);

            LLMaps[iMat].block(bf1_s,bf2_s,n1,n2) = bufMat.template cast<double>();
            // Symmetrize
            LLMaps[iMat].block(bf2_s,bf1_s,n2,n1) = bufMat.transpose().template cast<double>();
          } // Loop over x y z direction
        }



        // Compute SS dipole integrals
        if (int1e_sprsp_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)) {
          size_t n1n2 = n1*n2;
          // Place 12 integral blocks into their respective matricies
          // Order: X_x, X_y, X_z, X_s, Y_x, Y_y, Y_z, Y_s, Z_x, Z_y, Z_z, Z_s
          for(auto iMat = 0; iMat < SSMaps.size(); iMat++){
            Eigen::Map<
              const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
                Eigen::RowMajor>>
              bufMat(buff + iMat * n1n2,n1,n2);

            SSMaps[iMat].block(bf1_s,bf2_s,n1,n2) = bufMat.template cast<double>();
            // Symmetrize
            SSMaps[iMat].block(bf2_s,bf1_s,n2,n1) = bufMat.transpose().template cast<double>();
          } // Loop over integral blocks
        }
        
        } // Loop over s2 <= s1
      } // Loop over s1

    } // end OpenMP context

    CQMemManager::get().free(cacheAll, buffAll, env, bas, atm);

    //Currently, did not scale dipole by particle charge
  };

  template <>
  void MultipoleInts<double>::computeAOInts(BasisSet &basis, Molecule &mol,
      EMPerturbation &emPert, OPERATOR op, const HamiltonianOptions &options) {
    if (options.basisType != REAL_GTO)
      CErr("Only Real GTOs are allowed in MultipoleInts<double>",std::cout);
    // For 4C, use Libcint to compute multipole
    if (options.OneEScalarRelativity or options.OneESpinOrbit) {
      if (options.Libcint) {
        MultipoleRelDriverLibcint(mol, basis, options);
        return;
      }    
      CErr("Relativistic multipole integrals are implemented with Libint2",std::cout);
    }

    std::vector<double*> _multipole(1, nullptr);
    libint2::Operator libOp;

    switch (op) {
    case OVERLAP:
    case KINETIC:
    case NUCLEAR_POTENTIAL:
      CErr("Requested operator is not implemented in MultipoleInts,"
           " it is implemented in OnePInts",std::cout);
      break;
    case ELECTRON_REPULSION:
      CErr("Electron repulsion integrals are not implemented in MultipoleInts,"
           " they are implemented in TwoEInts",std::cout);
      break;
    case LEN_ELECTRIC_MULTIPOLE:
      try { _multipole[0] = CQMemManager::get().malloc<double>(nBasis()*nBasis()); }
      catch(...) {
        std::cout << std::fixed;
        std::cout << "Insufficient memory for the full INTS tensor ("
                  << (nBasis()*nBasis()/1e9) * sizeof(double) << " GB)" << std::endl;
        std::cout << std::endl << CQMemManager::get() << std::endl;
        CErr();
      }
      if (highOrder() >= 1) {
        std::copy_n(dipolePointers().begin(), 3, std::back_inserter(_multipole));
        libOp = libint2::Operator::emultipole1;
        if (highOrder() >= 2) {
          std::copy_n(quadrupolePointers().begin(), 6, std::back_inserter(_multipole));
          libOp = libint2::Operator::emultipole2;
          if (highOrder() == 3) {
            std::copy_n(octupolePointers().begin(), 10, std::back_inserter(_multipole));
            libOp = libint2::Operator::emultipole3;
          } else
            CErr("Requested operator is NYI in MultipoleInts.",std::cout);
        }
      }
      OnePInts<double>::OnePDriverLibint(libOp,mol,basis,_multipole,options.particle);
      CQMemManager::get().free(_multipole[0]);
      break;
    case VEL_ELECTRIC_MULTIPOLE:
    case MAGNETIC_MULTIPOLE:
      for (VectorInts<double> &vInts: components_) {
        vInts.computeAOInts(basis, mol, emPert, op, options);
      }
      break;
    default:
      CErr("Requested operator is not implemented in MultipoleInts.");
      break;
    }

  };

  template <>
  void OnePRelInts<double>::computeAOInts(BasisSet &basis, Molecule &mol,
      EMPerturbation&, OPERATOR op, const HamiltonianOptions &options) {
    if (options.basisType != REAL_GTO)
      CErr("Only Real GTOs are allowed in OnePRelInts<double>",std::cout);
    if (not options.OneEScalarRelativity or op != NUCLEAR_POTENTIAL)
      CErr("Only relativistic nuclear potential is implemented in OnePRelInts.",std::cout);


    if (options.Libcint) {
      OnePRelDriverLibcint(mol, basis, options);
      return;
#ifdef LIBCINT_FOR_FINITE_NUC
    } else {
      std::cout << "Using Libcint for relativistic integrals." << std::endl;
      OnePRelDriverLibcint(mol, basis, options);
      return;
#endif
    }
    
    // Retain only the classical nuclei for NEO V integrals
    Molecule cmol = mol.retainCNuc();

    std::vector<double*> _potential(1, pointer());
    if (options.finiteWidthNuc) {
      OnePDriverLocal<1,true>(
          [&](libint2::ShellPair& pair, libint2::Shell& sh1,
              libint2::Shell& sh2) -> std::vector<std::vector<double>> {
            return RealGTOIntEngine::computePotentialV(cmol.chargeDist,
                pair,sh1,sh2,cmol);
            }, basis.shells,_potential);
    }
    else
      OnePDriverLibint(libint2::Operator::nuclear,cmol,basis,_potential,options.particle);

    // Point nuclei is used when chargeDist is empty
    const std::vector<libint2::Shell> &chargeDist = options.finiteWidthNuc ?
        cmol.chargeDist : std::vector<libint2::Shell>();

    std::vector<double*> _PVdP(1, scalar().pointer());
    OnePInts<double>::OnePDriverLocal<1,true>(
          [&](libint2::ShellPair& pair, libint2::Shell& sh1,
              libint2::Shell& sh2) -> std::vector<std::vector<double>> {
            return RealGTOIntEngine::computepVdotp(chargeDist,
                pair,sh1,sh2,cmol);
            }, basis.shells, _PVdP);

    if (options.OneESpinOrbit) {
      if (not hasSpinOrbit())
        CErr("computeAOInts: Requested spin-orbit integrals, "
             "but the OnePRelInts object does not contain spin-orbit components");

      OnePInts<double>::OnePDriverLocal<3,false>(
            [&](libint2::ShellPair& pair, libint2::Shell& sh1,
                libint2::Shell& sh2) -> std::vector<std::vector<double>> {
              return RealGTOIntEngine::computeSL(chargeDist,
                  pair,sh1,sh2,cmol);
              }, basis.shells, SOXYZPointers());
    }

  };

  /**
   *  \brief Build 4-component spinor angular momentum operators for Spin
   *
  *  Computes 4C spin components with LL and SS contributions only
   *  following derivations.tex equations 1115-1173 (Spin Pauli spinor representation)
   *
   *  Returns vector of 3 PauliSpinorMatrices (x, y, z components)
  *  with structure: S() = 0
  *                  X/Y/Z() = Pauli-channel coefficients for each Cartesian
  *                            spin component.
  *  
  *  Scaling: in Pauli storage, LL uses the full overlap in the matching
  *           Pauli channel because spinGather contributes the final 1/2 block
  *           prefactor. SS uses 1/2*(1/2mc)^2 = 1/(8 m^2 c^2).
   *  
   *  Ref: derivations.tex lines 1115-1173 (Spin operator simplifications)
   */
  std::vector<cqmatrix::PauliSpinorMatrices<dcomplex>>
  build4CSpinVectorOperator(const BasisSet& basis, const Molecule& mol,
                            const HamiltonianOptions&) {
      BasisSet basisSet_ = basis.groupGeneralContractionBasis();
      size_t NB = basisSet_.nBasis;

      const double factor_1_2mc = 1.0 / (2.0 * SpeedOfLight()); // 1/(2mc) in atomic units
      const double factor_1_8mc2 = 0.5 * factor_1_2mc * factor_1_2mc; // 1/2*(1/2mc)^2

      std::vector<cqmatrix::PauliSpinorMatrices<dcomplex>> result;
      result.reserve(3);

      // SS from int1e_spsigmasp_sph: <SIGMA·P i|SIGMA|SIGMA·P j>
      // For each Cartesian spin component i = x,y,z, store Pauli coefficients X,Y,Z.
      double *SS[3][3];
      for (int i = 0; i < 3; ++i)
        for (int p = 0; p < 3; ++p) {
          SS[i][p] = CQMemManager::get().malloc<double>(NB*NB);
          std::fill_n(SS[i][p], NB*NB, 0.0);
        }

      // LL component stored in Pauli form. This uses the full overlap because
      // spinGather applies the final 1/2 factor when forming spin blocks.
      double *LL = CQMemManager::get().malloc<double>(NB*NB);
      std::fill_n(LL, NB*NB, 0.0);

      int nAtoms = static_cast<int>(basisSet_.centers.size());
      int nShells = static_cast<int>(basisSet_.nShell);

      int *atm = CQMemManager::get().template malloc<int>(nAtoms * ATM_SLOTS);
      int *bas = CQMemManager::get().template malloc<int>(nShells * BAS_SLOTS);
      double *env = CQMemManager::get().template malloc<double>(basisSet_.getLibcintEnvLength(mol));
      basisSet_.setLibcintEnv(mol, atm, bas, env, false);

      // Compute LL and SS components via libcint kernels.
      // LL from overlap, SS from int1e_spsigmasp_sph (12-component sigma-resolved output).

      // Allocate work buffers for libcint
      size_t cache_size = 0;
      for (int i = 0; i < nShells; i++) {
        size_t n;
        int shls[2]{i, i};
        n = int1e_ovlp_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
        n = int1e_spsigmasp_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
      }

      size_t buffSize = 12 * NB * NB;
      int nthreads = 1;
      #ifdef _OPENMP
      nthreads = GetNumThreads();
      #endif

      double* buffAll = CQMemManager::get().malloc<double>(buffSize * nthreads);
      double* cacheAll = CQMemManager::get().malloc<double>(cache_size * nthreads);

      #pragma omp parallel
      {
        int thread_id = 0;
        #ifdef _OPENMP
        thread_id = GetThreadID();
        #endif
        
        FINT shls[2];
        double* buff = buffAll + buffSize * thread_id;
        double* cache = cacheAll + cache_size * thread_id;

        // Loop over unique shell pairs
        for(size_t s1(0), bf1_s(0), s12(0); s1 < basisSet_.nShell; bf1_s += basisSet_.shells[s1].size(), s1++){
          size_t n1 = basisSet_.shells[s1].size();
          for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s += basisSet_.shells[s2].size(), s2++, s12++) {
            size_t n2 = basisSet_.shells[s2].size();

            // Round-robin work distribution
            #ifdef _OPENMP
            if(s12 % nthreads != thread_id) continue;
            #endif

            // Setup shell pair (libcint uses row-major ordering)
            shls[0] = s2;
            shls[1] = s1;

            if (int1e_ovlp_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)) {
              Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
                Eigen::RowMajor>> bufMat(buff, n1, n2);
              Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
                Eigen::ColMajor>> LL_mat(LL, NB, NB);
              LL_mat.block(bf1_s, bf2_s, n1, n2) = bufMat.template cast<double>();
              LL_mat.block(bf2_s, bf1_s, n2, n1) = bufMat.transpose().template cast<double>();
            }

            // Compute int1e_spsigmasp_sph: 12 components per basis pair.
            // Layout per pair n: [xx,xy,xz,0, yx,yy,yz,0, zx,zy,zz,0]
            if (int1e_spsigmasp_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)) {
              size_t n1n2 = n1 * n2;

              const int mapIdx[3][3] = {
                {0, 1, 2},
                {4, 5, 6},
                {8, 9, 10}
              };

              for (int iCart = 0; iCart < 3; ++iCart) {
                for (int pCart = 0; pCart < 3; ++pCart) {
                  const int iBuf = mapIdx[iCart][pCart];
                  Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
                    Eigen::RowMajor>> bufMat(buff + iBuf * n1n2, n1, n2);
                  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
                    Eigen::ColMajor>> targetMat(SS[iCart][pCart], NB, NB);

                  targetMat.block(bf1_s, bf2_s, n1, n2) = bufMat.template cast<double>();
                  if (s1 != s2)
                    targetMat.block(bf2_s, bf1_s, n2, n1) = bufMat.transpose().template cast<double>();
                }
              }
            }
          } // s2 <= s1
        } // s1
      } // omp parallel

      CQMemManager::get().free(cacheAll, buffAll);

      // Scale SS by 1/2*(1/(2mc))^2 for relativistic correction
      // This matches S_i^SS = (1/8m^2c^2) * (sigma·p sigma_i sigma·p).
      for (int i = 0; i < 3; ++i)
        for (int p = 0; p < 3; ++p)
          blas::scal(NB*NB, factor_1_8mc2, SS[i][p], 1);

      // Debug dump for small systems: inspect raw 4C spin-integral blocks.
      // This keeps output manageable while making H/STO-3G type cases transparent.
      // if (NB <= 4) {
      //   auto print_real_block = [&](const std::string &label, const double *ptr) {
      //     std::cout << "\n  [4C spin integral debug] " << label << "\n";
      //     std::cout << std::scientific << std::setprecision(12);
      //     for (size_t r = 0; r < NB; ++r) {
      //       std::cout << "    ";
      //       for (size_t c = 0; c < NB; ++c)
      //         std::cout << std::setw(20) << ptr[r + c * NB];
      //       std::cout << '\n';
      //     }
      //   };

      //   print_real_block("LL (overlap Pauli coefficient)", LL);
      //   for (int i = 0; i < 3; ++i) {
      //     for (int p = 0; p < 3; ++p) {
      //       std::string lbl = "SS[" + std::to_string(i) + "][" + std::to_string(p) + "] (scaled)";
      //       print_real_block(lbl, SS[i][p]);
      //     }
      //   }
      //   std::cout << std::flush;
      // }

      // Assemble 3 PauliSpinorMatrices (x, y, z components)
      // Each represents one Cartesian component of the spin operator
      for (size_t i = 0; i < 3; ++i) {
        cqmatrix::PauliSpinorMatrices<dcomplex> P(NB, true, true);

        // Start from zero scalar component. Cartesian spin operators are carried
        // by the matching Pauli channel (X/Y/Z), not by a common scalar part.
        Eigen::Map<Eigen::Matrix<dcomplex, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
          P_S(P.S().pointer(), NB, NB);
        P_S.setZero();

        // Pauli components for S_i. LL contributes to the matching channel,
        // and SS contributes as a small relativistic correction.
        Eigen::Map<Eigen::Matrix<dcomplex, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
          P_X(P.X().pointer(), NB, NB);
        Eigen::Map<Eigen::Matrix<dcomplex, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
          P_Y(P.Y().pointer(), NB, NB);
        Eigen::Map<Eigen::Matrix<dcomplex, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
          P_Z(P.Z().pointer(), NB, NB);

        P_X.setZero();
        P_Y.setZero();
        P_Z.setZero();

        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
          LL_map(LL, NB, NB);

        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
          SS_X_map(SS[i][0], NB, NB);
        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
          SS_Y_map(SS[i][1], NB, NB);
        Eigen::Map<const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>>
          SS_Z_map(SS[i][2], NB, NB);

        // Keep the full sigma-channel structure from SS, and place LL only
        // in the matching Cartesian Pauli channel (Sx->sigma_x, etc.).
        P_X = SS_X_map.template cast<dcomplex>();
        P_Y = SS_Y_map.template cast<dcomplex>();
        P_Z = SS_Z_map.template cast<dcomplex>();

        if (i == 0) P_X += LL_map.template cast<dcomplex>();
        if (i == 1) P_Y += LL_map.template cast<dcomplex>();
        if (i == 2) P_Z += LL_map.template cast<dcomplex>();

        result.emplace_back(std::move(P));
      }

      // Clean up temporary storage
      CQMemManager::get().free(LL);
      for (int i = 0; i < 3; ++i)
        for (int p = 0; p < 3; ++p)
          CQMemManager::get().free(SS[i][p]);
      CQMemManager::get().free(env, bas, atm);

      return result;
  }

  void build4CSpinVectorOperatorsFull(const BasisSet& basis,
      const Molecule& mol, const HamiltonianOptions& options,
      std::vector<cqmatrix::Matrix<dcomplex>>& total,
      std::vector<cqmatrix::Matrix<dcomplex>>& llOnly,
      std::vector<cqmatrix::Matrix<dcomplex>>& ssOnly) {

      (void)options;
      const size_t nComp = 3;

      total.clear();
      llOnly.clear();
      ssOnly.clear();
      total.reserve(nComp);
      llOnly.reserve(nComp);
      ssOnly.reserve(nComp);

      // Build LL/SS directly from raw kernels so the 4C block decomposition
      // is unambiguous: LL from overlap, SS from int1e_spsigmasp_sph.
      BasisSet basisSet_ = basis.groupGeneralContractionBasis();
      const size_t NB = basisSet_.nBasis;
 
      const double factor_1_2mc = 1.0 / (2.0 * SpeedOfLight());
      const double factor_1_8mc2 = 0.5 * factor_1_2mc * factor_1_2mc;

      double *LL = CQMemManager::get().malloc<double>(NB * NB);
      std::fill_n(LL, NB * NB, 0.0);

      double *SS[3][3];
      for (int i = 0; i < 3; ++i)
        for (int p = 0; p < 3; ++p) {
          SS[i][p] = CQMemManager::get().malloc<double>(NB * NB);
          std::fill_n(SS[i][p], NB * NB, 0.0);
        }

      int nAtoms = static_cast<int>(basisSet_.centers.size());
      int nShells = static_cast<int>(basisSet_.nShell);
      int *atm = CQMemManager::get().template malloc<int>(nAtoms * ATM_SLOTS);
      int *bas = CQMemManager::get().template malloc<int>(nShells * BAS_SLOTS);
      double *env = CQMemManager::get().template malloc<double>(basisSet_.getLibcintEnvLength(mol));
      basisSet_.setLibcintEnv(mol, atm, bas, env, false);

      size_t cache_size = 0;
      for (int i = 0; i < nShells; i++) {
        int shls[2]{i, i};
        size_t n = int1e_ovlp_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
        n = int1e_spsigmasp_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
      }

      size_t buffSize = 12 * NB * NB;
      int nthreads = 1;
      #ifdef _OPENMP
      nthreads = GetNumThreads();
      #endif

      double* buffAll = CQMemManager::get().malloc<double>(buffSize * nthreads);
      double* cacheAll = CQMemManager::get().malloc<double>(cache_size * nthreads);

      #pragma omp parallel
      {
        int thread_id = 0;
        #ifdef _OPENMP
        thread_id = GetThreadID();
        #endif

        FINT shls[2];
        double* buff = buffAll + buffSize * thread_id;
        double* cache = cacheAll + cache_size * thread_id;

        for(size_t s1(0), bf1_s(0), s12(0); s1 < basisSet_.nShell; bf1_s += basisSet_.shells[s1].size(), s1++){
          size_t n1 = basisSet_.shells[s1].size();
          for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s += basisSet_.shells[s2].size(), s2++, s12++) {
            size_t n2 = basisSet_.shells[s2].size();

            #ifdef _OPENMP
            if(s12 % nthreads != thread_id) continue;
            #endif

            shls[0] = s2;
            shls[1] = s1;

            if (int1e_ovlp_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)) {
              Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
                Eigen::RowMajor>> bufMat(buff, n1, n2);
              Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
                Eigen::ColMajor>> LL_mat(LL, NB, NB);

              LL_mat.block(bf1_s, bf2_s, n1, n2) = bufMat.template cast<double>();
              if (s1 != s2)
                LL_mat.block(bf2_s, bf1_s, n2, n1) = bufMat.transpose().template cast<double>();
            }

            if (int1e_spsigmasp_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)) {
              size_t n1n2 = n1 * n2;
              const int mapIdx[3][3] = {
                {0, 1, 2},
                {4, 5, 6},
                {8, 9, 10}
              };

              for (int iCart = 0; iCart < 3; ++iCart) {
                for (int pCart = 0; pCart < 3; ++pCart) {
                  const int iBuf = mapIdx[iCart][pCart];
                  Eigen::Map<const Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
                    Eigen::RowMajor>> bufMat(buff + iBuf * n1n2, n1, n2);
                  Eigen::Map<Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,
                    Eigen::ColMajor>> targetMat(SS[iCart][pCart], NB, NB);

                  targetMat.block(bf1_s, bf2_s, n1, n2) = bufMat.template cast<double>();
                  if (s1 != s2)
                    targetMat.block(bf2_s, bf1_s, n2, n1) = bufMat.transpose().template cast<double>();
                }
              }
            }
          }
        }
      }

      CQMemManager::get().free(cacheAll, buffAll);
      CQMemManager::get().free(env, bas, atm);

      for (int i = 0; i < 3; ++i)
        for (int p = 0; p < 3; ++p)
          blas::scal(NB * NB, factor_1_8mc2, SS[i][p], 1);

      for (size_t iComp = 0; iComp < nComp; ++iComp) {
        const size_t NBp = NB;

        cqmatrix::PauliSpinorMatrices<dcomplex> pll(NBp, true, true);
        pll.S().clear();
        pll.X().clear();
        pll.Y().clear();
        pll.Z().clear();

        cqmatrix::PauliSpinorMatrices<dcomplex> pss(NBp, true, true);
        pss.S().clear();
        pss.X().clear();
        pss.Y().clear();
        pss.Z().clear();

        cqmatrix::Matrix<dcomplex> LLc(NBp);
        for (size_t r = 0; r < NBp; ++r)
          for (size_t c = 0; c < NBp; ++c)
            LLc(r, c) = dcomplex(LL[r + c * NBp], 0.0);

        if (iComp == 0) pll.X() = LLc;
        if (iComp == 1) pll.Y() = LLc;
        if (iComp == 2) pll.Z() = LLc;

        for (size_t r = 0; r < NBp; ++r)
          for (size_t c = 0; c < NBp; ++c) {
            pss.X()(r, c) = dcomplex(SS[iComp][0][r + c * NBp], 0.0);
            pss.Y()(r, c) = dcomplex(SS[iComp][1][r + c * NBp], 0.0);
            pss.Z()(r, c) = dcomplex(SS[iComp][2][r + c * NBp], 0.0);
          }

        auto llSpin = pll.template spinGather<dcomplex>();
        auto ssSpin = pss.template spinGather<dcomplex>();

        cqmatrix::Matrix<dcomplex> zeroBlock(llSpin.nRows(), llSpin.nColumns());
        zeroBlock.clear();

        cqmatrix::Matrix<dcomplex> fullLL(2 * llSpin.nRows(), 2 * llSpin.nColumns());
        fullLL.componentGather(llSpin, zeroBlock, zeroBlock, zeroBlock, false);

        cqmatrix::Matrix<dcomplex> fullSS(2 * llSpin.nRows(), 2 * llSpin.nColumns());
        fullSS.componentGather(zeroBlock, zeroBlock, zeroBlock, ssSpin, false);

        cqmatrix::Matrix<dcomplex> fullTotal(fullLL);
        fullTotal += fullSS;

        total.emplace_back(std::move(fullTotal));
        llOnly.emplace_back(std::move(fullLL));
        ssOnly.emplace_back(std::move(fullSS));
      }

      CQMemManager::get().free(LL);
      for (int i = 0; i < 3; ++i)
        for (int p = 0; p < 3; ++p)
          CQMemManager::get().free(SS[i][p]);
  }

  template<>
  void GradInts<OnePInts,double>::computeAOInts(BasisSet& basis,
    Molecule& mol, EMPerturbation&, OPERATOR op, const HamiltonianOptions& options)
  {

    std::vector<double*> gradPtrs(3*nAtoms_, nullptr);

    for (auto i = 0; i < 3*nAtoms_; i++) {
      gradPtrs[i] = components_[i]->pointer();
    }

    switch (op) {
    case OVERLAP: {
      OnePInts<double>::OnePDriverLibint(
        libint2::Operator::overlap, mol, basis, gradPtrs, options.particle, 1
      );
      break;
    }
    case KINETIC:
      OnePInts<double>::OnePDriverLibint(
        libint2::Operator::kinetic, mol, basis, gradPtrs, options.particle, 1
      );
      break;
    case NUCLEAR_POTENTIAL:
      if (options.finiteWidthNuc)
        CErr("Finite width nuclei potential gradients not yet implemented!");
      else
        OnePInts<double>::OnePDriverLibint(
          libint2::Operator::nuclear, mol, basis, gradPtrs, options.particle, 1
        );
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
        break;
    }


  };

  template<>
  void GradInts<VectorInts, double>::computeAOInts(BasisSet&,
    Molecule&, EMPerturbation&, OPERATOR, const HamiltonianOptions&) {

    CErr("Gradient integrals for VectorInts operators not yet implemented!");

  };

  template<>
  void GradInts<MultipoleInts, double>::computeAOInts(BasisSet&,
    Molecule&, EMPerturbation&, OPERATOR, const HamiltonianOptions&) {

    CErr("Gradient integrals for multipole operators not yet implemented!");

  };

  template<>
  void GradInts<OnePRelInts, double>::computeAOInts(BasisSet&,
    Molecule&, EMPerturbation&, OPERATOR, const HamiltonianOptions&) {

    CErr("Gradient integrals for relativistic operators not yet implemented!");

  };

  template <>
  void GradInts<OnePInts, double>::computeAOInts(BasisSet&, BasisSet&,
    Molecule&, EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("Two basis gradients not implemented for OnePInts");
  };
  template <>
  void GradInts<VectorInts, double>::computeAOInts(BasisSet&, BasisSet&,
    Molecule&, EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("Two basis gradients not implemented for VectorInts");
  };
  template <>
  void GradInts<MultipoleInts, double>::computeAOInts(BasisSet&, BasisSet&,
    Molecule&, EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("Two basis gradients not implemented for MultipoleInts");
  };
  template <>
  void GradInts<OnePRelInts, double>::computeAOInts(BasisSet&, BasisSet&,
    Molecule&, EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("Two basis gradients not implemented for OnePRelInts");
  };


  template void Integrals<double>::computeAOOneP(
      Molecule&, BasisSet&, EMPerturbation&,
      const std::vector<std::pair<OPERATOR,size_t>>&,
      const HamiltonianOptions&);

  template void Integrals<double>::computeGradInts(
      Molecule&, BasisSet&, EMPerturbation&,
      const std::vector<std::pair<OPERATOR,size_t>>&,
      const HamiltonianOptions&);

}; // namespace ChronusQ
