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
      Molecule &mol, BasisSet& basis, std::vector<dcomplex*> mats,
      Particle p, size_t deriv, size_t S0a) {
    CErr("Only real GTOs are allowed",std::cout);
  };

  template <>
  void OnePInts<double>::OnePDriverLibint(libint2::Operator op,
      Molecule &mol, BasisSet& basis, std::vector<double*> mats, 
      Particle p, size_t deriv, size_t S0a) {

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
    for (OnePInts<double> &opi : components_)
      opi.matrix() *= -1.0 * options.particle.charge;

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
      switch (op) {
      case OVERLAP:
      case KINETIC:
      case NUCLEAR_POTENTIAL:
        OnePDriverLibcint(op, mol, basis, options);
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
      return;
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
    if (options.OneEScalarRelativity or options.OneESpinOrbit)
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

    std::vector<double*> _potential(1, pointer());
    if (options.finiteWidthNuc)
      OnePDriverLocal<1,true>(
          [&](libint2::ShellPair& pair, libint2::Shell& sh1,
              libint2::Shell& sh2) -> std::vector<std::vector<double>> {
            return RealGTOIntEngine::computePotentialV(mol.chargeDist,
                pair,sh1,sh2,mol);
            }, basis.shells,_potential);
    else
      OnePDriverLibint(libint2::Operator::nuclear,mol,basis,_potential,options.particle);

    // Point nuclei is used when chargeDist is empty
    const std::vector<libint2::Shell> &chargeDist = options.finiteWidthNuc ?
        mol.chargeDist : std::vector<libint2::Shell>();

    std::vector<double*> _PVdP(1, scalar().pointer());
    OnePInts<double>::OnePDriverLocal<1,true>(
          [&](libint2::ShellPair& pair, libint2::Shell& sh1,
              libint2::Shell& sh2) -> std::vector<std::vector<double>> {
            return RealGTOIntEngine::computepVdotp(chargeDist,
                pair,sh1,sh2,mol);
            }, basis.shells, _PVdP);

    if (options.OneESpinOrbit) {
      if (not hasSpinOrbit())
        CErr("computeAOInts: Requested spin-orbit integrals, "
             "but the OnePRelInts object does not contain spin-orbit components");

      OnePInts<double>::OnePDriverLocal<3,false>(
            [&](libint2::ShellPair& pair, libint2::Shell& sh1,
                libint2::Shell& sh2) -> std::vector<std::vector<double>> {
              return RealGTOIntEngine::computeSL(chargeDist,
                  pair,sh1,sh2,mol);
              }, basis.shells, SOXYZPointers());
    }

  };

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
