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
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <cqlinalg.hpp>
#include <cqlinalg/blasutil.hpp>
#include <util/timer.hpp>
#include <util/matout.hpp>

#include <util/threads.hpp>
#include <chrono>

#include <libcint.hpp>

//#define __DEBUGERI__
#define __INHOUSEGAUGE__
// _REPORT_INCORE_INTEGRAL_TIMINGS 

namespace ChronusQ {

  template <>
  void InCore4indexTPI<dcomplex>::computeERINRCINT(BasisSet&, Molecule&,
      EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("Only real GTOs are allowed",std::cout);
  };

  template <>
  void InCore4indexTPI<double>::computeERINRCINT(BasisSet &basisSet_, Molecule &molecule_,
      EMPerturbation&, OPERATOR, const HamiltonianOptions &hamiltonianOptions) {

    if (basisSet_.forceCart)
      CErr("Libcint + cartesian GTO NYI.");

    int nAtoms = molecule_.nAtoms;
    int nShells = basisSet_.nShell;

    // ATM_SLOTS = 6; BAS_SLOTS = 8;
    int *atm = CQMemManager::get().malloc<int>(nAtoms * ATM_SLOTS);
    int *bas = CQMemManager::get().malloc<int>(nShells * BAS_SLOTS);
    double *env = CQMemManager::get().malloc<double>(basisSet_.getLibcintEnvLength(molecule_));


    basisSet_.setLibcintEnv(molecule_, atm, bas, env);


    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();
 
    // Allocate and zero out ERIs
    size_t NB  = basisSet_.nBasis;
    size_t NB2 = NB*NB;
    size_t NB3 = NB2*NB;
    size_t NB4 = NB2*NB2;

    InCore4indexTPI<double>::clear();


    // Get threads result buffer
    int buffSize = (basisSet_.maxL+1)*(basisSet_.maxL+2)/2;
    size_t buffN4 = buffSize*buffSize*buffSize*buffSize;
    double *buffAll = CQMemManager::get().malloc<double>(buffN4*nthreads);

    std::cout<<"Using Libcint "<<std::endl;

#ifdef _REPORT_INCORE_CONTRACTION_TIMINGS 
    auto topERI4 = tick();
#endif

    #pragma omp parallel
    {
      int thread_id = GetThreadID();

      size_t n1,n2,n3,n4,i,j,k,l,ijkl,bf1,bf2,bf3,bf4;
      size_t s4_max;
      int shls[4];
      double *buff = buffAll+buffN4*thread_id;

      for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShells; 
          bf1_s+=n1, s1++) { 

        n1 = basisSet_.shells[s1].size(); // Size of Shell 1

      for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {

        n2 = basisSet_.shells[s2].size(); // Size of Shell 2

      for(size_t s3(0), bf3_s(0); s3 <= s1; bf3_s+=n3, s3++) {

        n3 = basisSet_.shells[s3].size(); // Size of Shell 3
        s4_max = (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4

      for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {

        n4 = basisSet_.shells[s4].size(); // Size of Shell 4

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s1234 % nthreads != thread_id ) continue;
        #endif

        shls[0] = int(s1);
        shls[1] = int(s2);
        shls[2] = int(s3);
        shls[3] = int(s4);

        if (basisSet_.forceCart) {
          if(int2e_cart(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr)==0) continue;
        } else {
          if(int2e_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr)==0) continue;
        }

        // permutational symmetry
	ijkl = 0ul;
        for(l = 0ul, bf4 = bf4_s ; l < n4; ++l, bf4++)
        for(k = 0ul, bf3 = bf3_s ; k < n3; ++k, bf3++) 
        for(j = 0ul, bf2 = bf2_s ; j < n2; ++j, bf2++) 
        for(i = 0ul, bf1 = bf1_s ; i < n1; ++i, bf1++) 
	{

            // (12 | 34)
            (*this)(bf1, bf2, bf3, bf4) = buff[ijkl];
            // (12 | 43)
            (*this)(bf1, bf2, bf4, bf3) = buff[ijkl];
            // (21 | 34)
            (*this)(bf2, bf1, bf3, bf4) = buff[ijkl];
            // (21 | 43)
            (*this)(bf2, bf1, bf4, bf3) = buff[ijkl];
            // (34 | 12)
            (*this)(bf3, bf4, bf1, bf2) = buff[ijkl];
            // (43 | 12)
            (*this)(bf4, bf3, bf1, bf2) = buff[ijkl];
            // (34 | 21)
            (*this)(bf3, bf4, bf2, bf1) = buff[ijkl];
            // (43 | 21)
            (*this)(bf4, bf3, bf2, bf1) = buff[ijkl];

	    ijkl++;

        }; // ijkl loop
      }; // s4
      }; // s3
      }; // s2
      }; // s1

    }; // omp region

#ifdef  _REPORT_INCORE_INTEGRAL_TIMINGS
    auto durERI4 = tock(topERI4);
    std::cout << "Libcint-ERI4 duration   = " << durERI4 << std::endl;
#endif

    CQMemManager::get().free(buffAll, atm, bas, env);

#ifdef __DEBUGERI__
    // Debug output of the ERIs
    std::cout << std::scientific << std::setprecision(16);
    std::cout << "Libcint ERI (ab|cd)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)(i, j, k, l) << std::endl;
    };
#endif // __DEBUGERI__ 


  } // computeERINRCINT




  template <>
  void InCore4indexRelERI<dcomplex>::computeERICINT(BasisSet&, Molecule&,
      EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("Only real GTOs are allowed",std::cout);
  };

  template <>
  void InCore4indexRelERI<double>::computeERICINT(BasisSet &originalBasisSet, Molecule &molecule_,
      EMPerturbation&, OPERATOR, const HamiltonianOptions &hamiltonianOptions) {

    if (originalBasisSet.forceCart)
      CErr("Libcint + cartesian GTO NYI.");

    BasisSet basisSet_ = originalBasisSet.groupGeneralContractionBasis();

    size_t buffSize = std::max_element(basisSet_.shells.begin(),
                                       basisSet_.shells.end(),
                                       [](libint2::Shell &a, libint2::Shell &b) {
                                         return a.size() < b.size();
                                       })->size();

    int nAtoms = molecule_.nAtoms;
    int nShells = basisSet_.nShell;

    // ATM_SLOTS = 6; BAS_SLOTS = 8;
    int *atm = CQMemManager::get().malloc<int>(nAtoms * ATM_SLOTS);
    int *bas = CQMemManager::get().malloc<int>(nShells * BAS_SLOTS);
    double *env = CQMemManager::get().malloc<double>(basisSet_.getLibcintEnvLength(molecule_));


    basisSet_.setLibcintEnv(molecule_, atm, bas, env);


    size_t cache_size = 0;
    for (int i = 0; i < nShells; i++) {
      size_t n;
      int shls[4]{i,i,i,i};
      if (basisSet_.forceCart) {
        n = int2e_cart(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
        if(hamiltonianOptions.DiracCoulomb) {
          n = int2e_ipvip1_cart(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
        }
        if(hamiltonianOptions.Gaunt) {
          n = int2e_ip1ip2_cart(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
        }
        if(hamiltonianOptions.DiracCoulombSSSS) {
          n = int2e_ipvip1ipvip2_cart(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
        }
        if(hamiltonianOptions.Gauge) {
          n = int2e_gauge_r1_sps1sps2_cart(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
        }
      } else {
        n = int2e_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
        cache_size = std::max(cache_size, n);
        if(hamiltonianOptions.DiracCoulomb) {
          n = int2e_ipvip1_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
        }
        if(hamiltonianOptions.Gaunt) {
          n = int2e_ip1ip2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
        }
        if(hamiltonianOptions.DiracCoulombSSSS) {
          n = int2e_ipvip1ipvip2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
        }
        if(hamiltonianOptions.Gauge) {
          n = int2e_gauge_r1_ssp1sps2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
          n = int2e_gauge_r2_ssp1sps2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
          n = int2e_gauge_r1_ssp1ssp2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
          n = int2e_gauge_r2_ssp1ssp2_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
          cache_size = std::max(cache_size, n);
        }
      }
    }


    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();
 
    // Allocate and zero out ERIs
    size_t NB  = basisSet_.nBasis;
    size_t NB2 = NB*NB;
    size_t NB3 = NB2*NB;
    size_t NB4 = NB2*NB2;

    InCore4indexTPI<double>::clear();


    // Get threads result buffer
    size_t buffN4 = buffSize*buffSize*buffSize*buffSize;
    if(hamiltonianOptions.DiracCoulombSSSS)
      buffN4 *= 81;
    else if (hamiltonianOptions.Gauge)
      buffN4 *= 16;
    else if (hamiltonianOptions.DiracCoulomb or hamiltonianOptions.Gaunt)
      buffN4 *= 9;

    double *buffAll;
    if(hamiltonianOptions.Gauge) buffAll = CQMemManager::get().malloc<double>(buffN4*nthreads*2);
    else buffAll = CQMemManager::get().malloc<double>(buffN4*nthreads);

    double *cacheAll = CQMemManager::get().malloc<double>(cache_size*nthreads);

    std::cout<<"Using Libcint "<<std::endl;

#if 1 // (ij|kl)

#ifdef _REPORT_INCORE_INTEGRAL_TIMINGS 
    auto topERI4 = tick();
#endif

    #pragma omp parallel
    {
      int thread_id = GetThreadID();

      size_t n1,n2,n3,n4,i,j,k,l,ijkl,bf1,bf2,bf3,bf4;
      size_t s4_max;
      int shls[4];
      double *buff = buffAll+buffN4*thread_id;
      double *cache = cacheAll+cache_size*thread_id;

      for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShells; 
          bf1_s+=n1, s1++) { 

        n1 = basisSet_.shells[s1].size(); // Size of Shell 1

      for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {

        n2 = basisSet_.shells[s2].size(); // Size of Shell 2

      for(size_t s3(0), bf3_s(0); s3 <= s1; bf3_s+=n3, s3++) {

        n3 = basisSet_.shells[s3].size(); // Size of Shell 3
        s4_max = (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4

      for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {

        n4 = basisSet_.shells[s4].size(); // Size of Shell 4

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s1234 % nthreads != thread_id ) continue;
        #endif

        shls[0] = int(s1);
        shls[1] = int(s2);
        shls[2] = int(s3);
        shls[3] = int(s4);

        if (basisSet_.forceCart) {
          if(int2e_cart(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
        } else {
          if(int2e_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
        }

        // permutational symmetry
	ijkl = 0ul;
        for(l = 0ul, bf4 = bf4_s ; l < n4; ++l, bf4++)
        for(k = 0ul, bf3 = bf3_s ; k < n3; ++k, bf3++) 
        for(j = 0ul, bf2 = bf2_s ; j < n2; ++j, bf2++) 
        for(i = 0ul, bf1 = bf1_s ; i < n1; ++i, bf1++) 
	{

            // (12 | 34)
            (*this)(bf1, bf2, bf3, bf4) = buff[ijkl];
            // (12 | 43)
            (*this)(bf1, bf2, bf4, bf3) = buff[ijkl];
            // (21 | 34)
            (*this)(bf2, bf1, bf3, bf4) = buff[ijkl];
            // (21 | 43)
            (*this)(bf2, bf1, bf4, bf3) = buff[ijkl];
            // (34 | 12)
            (*this)(bf3, bf4, bf1, bf2) = buff[ijkl];
            // (43 | 12)
            (*this)(bf4, bf3, bf1, bf2) = buff[ijkl];
            // (34 | 21)
            (*this)(bf3, bf4, bf2, bf1) = buff[ijkl];
            // (43 | 21)
            (*this)(bf4, bf3, bf2, bf1) = buff[ijkl];

	    ijkl++;

        }; // ijkl loop
      }; // s4
      }; // s3
      }; // s2
      }; // s1

    }; // omp region

#ifdef _REPORT_INCORE_CONTRACTION_TIMINGS 
    auto durERI4 = tock(topERI4);
    //std::cout << "L = "<< basisSet_.shells[s1].contr[0].l<<" "<<basisSet_.shells[s2].contr[0].l<<" "
    //	               << basisSet_.shells[s3].contr[0].l<<" "<<basisSet_.shells[s4].contr[0].l<<std::endl;
    std::cout << "Libcint-ERI4 duration   = " << durERI4 << std::endl;
#endif

#ifdef __DEBUGERI__
    // Debug output of the ERIs
    std::cout << std::scientific << std::setprecision(16);
    std::cout << "Libcint ERI (ab|cd)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)(i, j, k, l) << std::endl;
    };
#endif // __DEBUGERI__ 

#endif // (ijkl)





    /* Dirac-Coulomb Integrals */
    if(hamiltonianOptions.DiracCoulomb) { // Dirac-Coulomb ∇_i∇_j(ij|kl)

      for (InCore4indexTPI<double>& c : components_) c.clear();
  
      int AxBx = 0;
      int AxBy = 1;
      int AxBz = 2;
      int AyBx = 3;
      int AyBy = 4;
      int AyBz = 5;
      int AzBx = 6;
      int AzBy = 7;
      int AzBz = 8;
  
#ifdef _REPORT_INCORE_INTEGRAL_TIMINGS 
      auto topERIDC = tick();
#endif  
      #pragma omp parallel
      {
        int thread_id = GetThreadID();
  
        size_t n1,n2,n3,n4,i,j,k,l,ijkl,bf1,bf2,bf3,bf4;
        size_t s4_max;
        int shls[4];
        double *buff = buffAll + buffN4*thread_id;
        double *cache = cacheAll+cache_size*thread_id;
  
        for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShells; bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2
  
        for(size_t s3(0), bf3_s(0); s3 < nShells ; bf3_s+=n3, s3++) {
  
          n3 = basisSet_.shells[s3].size(); // Size of Shell 3
  
        for(size_t s4(0), bf4_s(0); s4 <= s3; bf4_s+=n4, s4++, s1234++) {
  
          n4 = basisSet_.shells[s4].size(); // Size of Shell 4
  
          // Round Robbin work distribution
          #ifdef _OPENMP
          if( s1234 % nthreads != thread_id ) continue;
          #endif
  
          shls[0] = int(s1);
          shls[1] = int(s2);
          shls[2] = int(s3);
          shls[3] = int(s4);
  
          if (basisSet_.forceCart) {
            if(int2e_ipvip1_cart(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
          } else {
            if(int2e_ipvip1_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
          }

          ijkl = 0ul;
  	  auto nQuad = n1*n2*n3*n4;
          for(l = 0ul, bf4 = bf4_s ; l < n4; ++l, bf4++)
          for(k = 0ul, bf3 = bf3_s ; k < n3; ++k, bf3++) 
          for(j = 0ul, bf2 = bf2_s ; j < n2; ++j, bf2++) 
          for(i = 0ul, bf1 = bf1_s ; i < n1; ++i, bf1++) {
  
#ifdef __DEBUGERI__
  
          std::cout << std::scientific << std::setprecision(16);
  	  std::cout <<"Libcint ∇A∙∇B(ij|kl)"<<std::endl;
  	  std::cout<<buff[AxBx*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AxBy*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AxBz*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AyBx*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AyBy*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AyBz*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AzBx*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AzBy*nQuad+ijkl]<<std::endl;
  	  std::cout<<buff[AzBz*nQuad+ijkl]<<std::endl;
  
#endif
  
          // ∇A∙∇B(ij|kl)
          auto dAdotdB = buff[AxBx*nQuad+ijkl] + buff[AyBy*nQuad+ijkl] + buff[AzBz*nQuad+ijkl];
          // ∇Ax∇B(ijkl)
          auto dAcrossdB_x =  buff[AyBz*nQuad+ijkl] - buff[AzBy*nQuad+ijkl];
          auto dAcrossdB_y = -buff[AxBz*nQuad+ijkl] + buff[AzBx*nQuad+ijkl];
          auto dAcrossdB_z =  buff[AxBy*nQuad+ijkl] - buff[AyBx*nQuad+ijkl];
  
          auto IJKL = bf1 + bf2*NB + bf3*NB2 + bf4*NB3;
          auto IJLK = bf1 + bf2*NB + bf4*NB2 + bf3*NB3;
          auto JIKL = bf2 + bf1*NB + bf3*NB2 + bf4*NB3;
          auto JILK = bf2 + bf1*NB + bf4*NB2 + bf3*NB3;
  
          // ∇A∙∇B(ij|kl) followed by ∇Ax∇B(ij|kl) X, Y, and Z
          // (ij|kl)
          (*this)[0].pointer()[IJKL] = dAdotdB;
          (*this)[1].pointer()[IJKL] = dAcrossdB_x;
          (*this)[2].pointer()[IJKL] = dAcrossdB_y;
          (*this)[3].pointer()[IJKL] = dAcrossdB_z;
          // (ij|lk)
          (*this)[0].pointer()[IJLK] = dAdotdB;
          (*this)[1].pointer()[IJLK] = dAcrossdB_x;
          (*this)[2].pointer()[IJLK] = dAcrossdB_y;
          (*this)[3].pointer()[IJLK] = dAcrossdB_z;
          // (ji|kl)
          (*this)[0].pointer()[JIKL] = dAdotdB;
          (*this)[1].pointer()[JIKL] = -dAcrossdB_x;
          (*this)[2].pointer()[JIKL] = -dAcrossdB_y;
          (*this)[3].pointer()[JIKL] = -dAcrossdB_z;
          // (ji|lk)
          (*this)[0].pointer()[JILK] = dAdotdB;
          (*this)[1].pointer()[JILK] = -dAcrossdB_x;
          (*this)[2].pointer()[JILK] = -dAcrossdB_y;
          (*this)[3].pointer()[JILK] = -dAcrossdB_z;
  
  	  ijkl++;
  
          }; // ijkl loop
        }; // s4
        }; // s3
        }; // s2
        }; // s1
  
      }; // omp region
 
#ifdef  _REPORT_INCORE_INTEGRAL_TIMINGS
      auto durERIDC = tock(topERIDC);
      std::cout << "Libcint-ERI-Dirac-Coulomb duration   = " << durERIDC << std::endl;
#endif   


#if 0
      std::cout << std::scientific << std::setprecision(16);
      std::cout << "ERI00-03: ∇A∙∇B(ab|cd)  ∇Ax∇B(ab|cd)-X  ∇Ax∇B(ab|cd)-Y  ∇Ax∇B(ab|cd)-Z" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[0](i, j, k, l);
        std::cout << "   ";
        std::cout << (*this)[1](i, j, k, l);
        std::cout << "   ";
        std::cout << (*this)[2](i, j, k, l);
        std::cout << "   ";
        std::cout << (*this)[3](i, j, k, l) << std::endl;
      };
#endif

    } // Dirac-Coulomb ∇_i∇_j(ij|kl)












    /* Gaunt Integrals */
    // ∇_j∇_k(ij|kl)
    if(hamiltonianOptions.Gaunt) {

      auto nERIRef = 0;
      if (hamiltonianOptions.DiracCoulomb) nERIRef +=4; // Dirac-Coulomb

      int AxCx = 0;
      int AxCy = 1;
      int AxCz = 2;
      int AyCx = 3;
      int AyCy = 4;
      int AyCz = 5;
      int AzCx = 6;
      int AzCy = 7;
      int AzCz = 8;
  
      int BxCx = 0;
      int BxCy = 1;
      int BxCz = 2;
      int ByCx = 3;
      int ByCy = 4;
      int ByCz = 5;
      int BzCx = 6;
      int BzCy = 7;
      int BzCz = 8;
  
#ifdef _REPORT_INCORE_INTEGRAL_TIMINGS 
      auto topERIGaunt = tick();
#endif

      #pragma omp parallel
      {
        int thread_id = GetThreadID();

        size_t n1,n2,n3,n4,i,j,k,l,ijkl,bf1,bf2,bf3,bf4;
        size_t s4_max;
        int shls[4];
        double *buff = buffAll + buffN4*thread_id;
        double *cache = cacheAll+cache_size*thread_id;
  
        for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShells; 
            bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 < nShells; bf2_s+=n2, s2++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2
  
        for(size_t s3(0), bf3_s(0); s3 <= s1 ; bf3_s+=n3, s3++) {
  
          n3 = basisSet_.shells[s3].size(); // Size of Shell 3
  
        for(size_t s4(0), bf4_s(0); s4 < nShells ; bf4_s+=n4, s4++, s1234++) {
  
          n4 = basisSet_.shells[s4].size(); // Size of Shell 4
  
          // Round Robbin work distribution
          #ifdef _OPENMP
          if( s1234 % nthreads != thread_id ) continue;
          #endif
  
          shls[0] = int(s1);
          shls[1] = int(s2);
          shls[2] = int(s3);
          shls[3] = int(s4);

          if (basisSet_.forceCart) {
            if(int2e_ip1ip2_cart(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
          } else {
            if(int2e_ip1ip2_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
          }

  
  	  ijkl = 0ul;
  	  auto nQuad = n1*n2*n3*n4;
          for(l = 0ul, bf4 = bf4_s ; l < n4; ++l, bf4++)
          for(k = 0ul, bf3 = bf3_s ; k < n3; ++k, bf3++) 
          for(j = 0ul, bf2 = bf2_s ; j < n2; ++j, bf2++) 
          for(i = 0ul, bf1 = bf1_s ; i < n1; ++i, bf1++) {
  
#ifdef __DEBUGERI__
            std::cout << std::scientific << std::setprecision(16);
            std::cout <<"Libcint ∇A∙∇C(ij|kl)"<<std::endl;
  	    std::cout<<buff[AxCx*nQuad+ijkl]<<std::endl;
  	    std::cout<<buff[AxCy*nQuad+ijkl]<<std::endl;
  	    std::cout<<buff[AxCz*nQuad+ijkl]<<std::endl;
  	    std::cout<<buff[AyCx*nQuad+ijkl]<<std::endl;
  	    std::cout<<buff[AyCy*nQuad+ijkl]<<std::endl;
  	    std::cout<<buff[AyCz*nQuad+ijkl]<<std::endl;
  	    std::cout<<buff[AzCx*nQuad+ijkl]<<std::endl;
  	    std::cout<<buff[AzCy*nQuad+ijkl]<<std::endl;
  	    std::cout<<buff[AzCz*nQuad+ijkl]<<std::endl;
#endif
  
            // ∇A∙∇C(ij|kl)
            auto dAdotdC = buff[AxCx*nQuad+ijkl] + buff[AyCy*nQuad+ijkl] + buff[AzCz*nQuad+ijkl];
            // ∇Ax∇C(ijkl)
            auto dAcrossdC_x =  buff[AyCz*nQuad+ijkl] - buff[AzCy*nQuad+ijkl];
            auto dAcrossdC_y = -buff[AxCz*nQuad+ijkl] + buff[AzCx*nQuad+ijkl];
            auto dAcrossdC_z =  buff[AxCy*nQuad+ijkl] - buff[AyCx*nQuad+ijkl];
  
  	    // Change the index so that we do ∇B∙∇C(ij|kl) using the ∇A∙∇C engine
            auto IJKL = bf2 + bf1*NB + bf3*NB2 + bf4*NB3;
            auto LKJI = bf4 + bf3*NB + bf1*NB2 + bf2*NB3;
  
            // ∇B∙∇C(ij|kl) followed by ∇Bx∇C(ij|kl) X, Y, and Z
            // (ij|kl)
            (*this)[nERIRef].pointer()[IJKL] = dAdotdC;
            (*this)[nERIRef+1].pointer()[IJKL] = dAcrossdC_x;
            (*this)[nERIRef+2].pointer()[IJKL] = dAcrossdC_y;
            (*this)[nERIRef+3].pointer()[IJKL] = dAcrossdC_z;
  	    // (lk|ji)
            (*this)[nERIRef].pointer()[LKJI] = dAdotdC;
            (*this)[nERIRef+1].pointer()[LKJI] = -dAcrossdC_x;
            (*this)[nERIRef+2].pointer()[LKJI] = -dAcrossdC_y;
            (*this)[nERIRef+3].pointer()[LKJI] = -dAcrossdC_z;
  
  
            // ∇B_x∇C_y(ij|kl) + ∇B_y∇C_x(ij|kl)
            // (ij|kl)
            (*this)[nERIRef+4].pointer()[IJKL] = buff[BxCy*nQuad+ijkl] + buff[ByCx*nQuad+ijkl];
            // (lk|ji)
            (*this)[nERIRef+4].pointer()[LKJI] = buff[BxCy*nQuad+ijkl] + buff[ByCx*nQuad+ijkl];
  
  
            // ∇B_y∇C_x(ij|kl)
            // (ij|kl)
            (*this)[nERIRef+5].pointer()[IJKL] = buff[ByCx*nQuad+ijkl];
            // (lk|ji)
            (*this)[nERIRef+5].pointer()[LKJI] = buff[BxCy*nQuad+ijkl];
  
  
            // ∇B_x∇C_z(ij|kl) + ∇B_z∇C_x(ij|kl)
            // (ij|kl)
            (*this)[nERIRef+6].pointer()[IJKL] = buff[BxCz*nQuad+ijkl] + buff[BzCx*nQuad+ijkl];
            // (lk|ji)
            (*this)[nERIRef+6].pointer()[LKJI] = buff[BxCz*nQuad+ijkl] + buff[BzCx*nQuad+ijkl];
  
  
            // ∇B_z∇C_x(ij|kl)
            // (ij|kl)
            (*this)[nERIRef+7].pointer()[IJKL] = buff[BzCx*nQuad+ijkl];
            // (lk|ji)
            (*this)[nERIRef+7].pointer()[LKJI] = buff[BxCz*nQuad+ijkl];
  
  
            // ∇B_y∇C_z(ij|kl) + ∇B_z∇C_y(ij|kl)
            // (ij|kl)
            (*this)[nERIRef+8].pointer()[IJKL] = buff[ByCz*nQuad+ijkl] + buff[BzCy*nQuad+ijkl];
            // (lk|ji)
            (*this)[nERIRef+8].pointer()[LKJI] = buff[ByCz*nQuad+ijkl] + buff[BzCy*nQuad+ijkl];
  
  
            // ∇B_z∇C_y(ij|kl)
            // (ij|kl)
            (*this)[nERIRef+9].pointer()[IJKL] = buff[BzCy*nQuad+ijkl];
            // (lk|ji)
            (*this)[nERIRef+9].pointer()[LKJI] = buff[ByCz*nQuad+ijkl];
  
  
            // - ∇B_x∇C_x(ij|kl) - ∇B_y∇C_y(ij|kl) + ∇B_z∇C_z(ij|kl)
            // (ij|kl)
            (*this)[nERIRef+10].pointer()[IJKL] = - buff[BxCx*nQuad+ijkl] - buff[ByCy*nQuad+ijkl] + buff[BzCz*nQuad+ijkl];
            // (lk|ji)
            (*this)[nERIRef+10].pointer()[LKJI] = - buff[BxCx*nQuad+ijkl] - buff[ByCy*nQuad+ijkl] + buff[BzCz*nQuad+ijkl];
  
  
            // ∇B_x∇C_x(ij|kl) - ∇B_y∇C_y(ij|kl) - ∇B_z∇C_z(ij|kl)
            // (ij|kl)
            (*this)[nERIRef+11].pointer()[IJKL] = buff[BxCx*nQuad+ijkl] - buff[ByCy*nQuad+ijkl] - buff[BzCz*nQuad+ijkl];
            // (lk|ji)
            (*this)[nERIRef+11].pointer()[LKJI] = buff[BxCx*nQuad+ijkl] - buff[ByCy*nQuad+ijkl] - buff[BzCz*nQuad+ijkl];
  
  
            // - ∇B_x∇C_x(ij|kl) + ∇B_y∇C_y(ij|kl) - ∇B_z∇C_z(ij|kl)
            // (ij|kl)
            (*this)[nERIRef+12].pointer()[IJKL] = - buff[BxCx*nQuad+ijkl] + buff[ByCy*nQuad+ijkl] - buff[BzCz*nQuad+ijkl];
            // (lk|ji)
            (*this)[nERIRef+12].pointer()[LKJI] = - buff[BxCx*nQuad+ijkl] + buff[ByCy*nQuad+ijkl] - buff[BzCz*nQuad+ijkl];
  
  
            // ∇B_x∇C_x(ij|kl)
            // (ij|kl)
            (*this)[nERIRef+13].pointer()[IJKL] = buff[BxCx*nQuad+ijkl];
            // (lk|ji)
            (*this)[nERIRef+13].pointer()[LKJI] = buff[BxCx*nQuad+ijkl];
  
  
            // ∇B_x∇C_y(ij|kl)
            // (ij|kl)
            (*this)[nERIRef+14].pointer()[IJKL] = buff[BxCy*nQuad+ijkl];
            // (lk|ji)
            (*this)[nERIRef+14].pointer()[LKJI] = buff[ByCx*nQuad+ijkl];
  
  
            // ∇B_x∇C_z(ij|kl)
            // (ij|kl)
            (*this)[nERIRef+15].pointer()[IJKL] = buff[BxCz*nQuad+ijkl];
            // (lk|ji)
            (*this)[nERIRef+15].pointer()[LKJI] = buff[BzCx*nQuad+ijkl];
  

            // ∇B_y∇C_y(ij|kl)
            // (ij|kl)
            (*this)[nERIRef+16].pointer()[IJKL] = buff[ByCy*nQuad+ijkl];
            // (lk|ji)
            (*this)[nERIRef+16].pointer()[LKJI] = buff[ByCy*nQuad+ijkl];
  

            // ∇B_y∇C_z(ij|kl)
            // (ij|kl)
            (*this)[nERIRef+17].pointer()[IJKL] = buff[ByCz*nQuad+ijkl];
            // (lk|ji)
            (*this)[nERIRef+17].pointer()[LKJI] = buff[BzCy*nQuad+ijkl];
  
  
            // ∇B_z∇C_z(ij|kl)
            // (ij|kl)
            (*this)[nERIRef+18].pointer()[IJKL] = buff[BzCz*nQuad+ijkl];
            // (lk|ji)
            (*this)[nERIRef+18].pointer()[LKJI] = buff[BzCz*nQuad+ijkl];
  
  
  	        ijkl++;
          }; // ijkl loop
        }; // s4
        }; // s3
        }; // s2
        }; // s1
  
      }; // omp region
 
#ifdef  _REPORT_INCORE_INTEGRAL_TIMINGS
      auto durERIGaunt = tock(topERIGaunt);
      std::cout << "Libcint-ERI-Gaunt duration   = " << durERIGaunt << std::endl;
#endif

  
#ifdef __DEBUGERI__

      std::cout << std::scientific << std::setprecision(16);
  
      std::cout << "ERI04-07: ∇B∙∇C(ab|cd)  ∇Bx∇C(ab|cd)-X  ∇Bx∇C(ab|cd)-Y  ∇Bx∇C(ab|cd)-Z" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef](i, j, k, l);
        std::cout << "   ";
        std::cout << (*this)[nERIRef+1](i, j, k, l);
        std::cout << "   ";
        std::cout << (*this)[nERIRef+2](i, j, k, l);
        std::cout << "   ";
        std::cout << (*this)[nERIRef+3](i, j, k, l) << std::endl;
      };
  
      std::cout << "ERI08: ∇B_x∇C_y(ij|kl) + ∇B_y∇C_x(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+4](i, j, k, l) << std::endl;
      };
  
      std::cout << "ERI09: ∇B_y∇C_x(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+5](i, j, k, l) << std::endl;
      };
  
      std::cout << "ERI10: ∇B_x∇C_z(ij|kl) + ∇B_z∇C_x(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+6](i, j, k, l) << std::endl;
      };
  
      std::cout << "ERI11: ∇B_z∇C_x(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+7](i, j, k, l) << std::endl;
      };
  
      std::cout << "ERI12: ∇B_y∇C_z(ij|kl) + ∇B_z∇C_y(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+8](i, j, k, l) << std::endl;
      };
  
      std::cout << "ERI13: ∇B_z∇C_y(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+9](i, j, k, l) << std::endl;
      };
  
      std::cout << "ERI14: - ∇B_x∇C_x(ij|kl) - ∇B_y∇C_y(ij|kl) + ∇B_z∇C_z(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+10](i, j, k, l) << std::endl;
      };
  
      std::cout << "ERI15: ∇B_x∇C_x(ij|kl) - ∇B_y∇C_y(ij|kl) - ∇B_z∇C_z(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+11](i, j, k, l) << std::endl;
      };
  
      std::cout << "ERI16: - ∇B_x∇C_x(ij|kl) + ∇B_y∇C_y(ij|kl) - ∇B_z∇C_z(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+12](i, j, k, l) << std::endl;
      };
  
      std::cout << "ERI17: ∇B_x∇C_x(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+13](i, j, k, l) << std::endl;
      };
  
      std::cout << "ERI18: ∇B_x∇C_y(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+14](i, j, k, l) << std::endl;
      };
  
      std::cout << "ERI19: ∇B_x∇C_z(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+15](i, j, k, l) << std::endl;
      };
  
      std::cout << "ERI20: ∇B_y∇C_y(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+16](i, j, k, l) << std::endl;
      };
  
      std::cout << "ERI21: ∇B_y∇C_z(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+17](i, j, k, l) << std::endl;
      };
  
      std::cout << "ERI22: ∇B_z∇C_z(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+18](i, j, k, l) << std::endl;
      };

#endif
  
    } // Gaunt






    /* Dirac-Coulomb (SSSS) Integrals */
    if(hamiltonianOptions.DiracCoulombSSSS) { // Dirac-Coulomb ∇_i∇_j∇_k∇_l(ij|kl)

      auto nERIRef = 0;
      if(hamiltonianOptions.DiracCoulomb) nERIRef +=4;
      if(hamiltonianOptions.Gaunt) nERIRef += 19;

      int AxBxCxDx = 0;
      int AxBxCxDy = 1;
      int AxBxCxDz = 2;
      int AxBxCyDx = 3;
      int AxBxCyDy = 4;
      int AxBxCyDz = 5;
      int AxBxCzDx = 6;
      int AxBxCzDy = 7;
      int AxBxCzDz = 8;

      int AxByCxDx = 9;
      int AxByCxDy = 10;
      int AxByCxDz = 11;
      int AxByCyDx = 12;
      int AxByCyDy = 13;
      int AxByCyDz = 14;
      int AxByCzDx = 15;
      int AxByCzDy = 16;
      int AxByCzDz = 17;
  
      int AxBzCxDx = 18;
      int AxBzCxDy = 19;
      int AxBzCxDz = 20;
      int AxBzCyDx = 21;
      int AxBzCyDy = 22;
      int AxBzCyDz = 23;
      int AxBzCzDx = 24;
      int AxBzCzDy = 25;
      int AxBzCzDz = 26;



      int AyBxCxDx = 27;
      int AyBxCxDy = 28;
      int AyBxCxDz = 29;
      int AyBxCyDx = 30;
      int AyBxCyDy = 31;
      int AyBxCyDz = 32;
      int AyBxCzDx = 33;
      int AyBxCzDy = 34;
      int AyBxCzDz = 35;

      int AyByCxDx = 36;
      int AyByCxDy = 37;
      int AyByCxDz = 38;
      int AyByCyDx = 39;
      int AyByCyDy = 40;
      int AyByCyDz = 41;
      int AyByCzDx = 42;
      int AyByCzDy = 43;
      int AyByCzDz = 44;
  
      int AyBzCxDx = 45;
      int AyBzCxDy = 46;
      int AyBzCxDz = 47;
      int AyBzCyDx = 48;
      int AyBzCyDy = 49;
      int AyBzCyDz = 50;
      int AyBzCzDx = 51;
      int AyBzCzDy = 52;
      int AyBzCzDz = 53;



      int AzBxCxDx = 54;
      int AzBxCxDy = 55;
      int AzBxCxDz = 56;
      int AzBxCyDx = 57;
      int AzBxCyDy = 58;
      int AzBxCyDz = 59;
      int AzBxCzDx = 60;
      int AzBxCzDy = 61;
      int AzBxCzDz = 62;

      int AzByCxDx = 63;
      int AzByCxDy = 64;
      int AzByCxDz = 65;
      int AzByCyDx = 66;
      int AzByCyDy = 67;
      int AzByCyDz = 68;
      int AzByCzDx = 69;
      int AzByCzDy = 70;
      int AzByCzDz = 71;
  
      int AzBzCxDx = 72;
      int AzBzCxDy = 73;
      int AzBzCxDz = 74;
      int AzBzCyDx = 75;
      int AzBzCyDy = 76;
      int AzBzCyDz = 77;
      int AzBzCzDx = 78;
      int AzBzCzDy = 79;
      int AzBzCzDz = 80;

#ifdef _REPORT_INCORE_INTEGRAL_TIMINGS 
      auto topERIDCSSSS = tick();
#endif

      #pragma omp parallel
      {
        int thread_id = GetThreadID();
  
        size_t n1,n2,n3,n4,i,j,k,l,ijkl,bf1,bf2,bf3,bf4;
        size_t s4_max;
        int shls[4];
        double *buff = buffAll + buffN4*thread_id;
        double *cache = cacheAll+cache_size*thread_id;

        for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShells; bf1_s+=n1, s1++) { 
  
          n1 = basisSet_.shells[s1].size(); // Size of Shell 1
  
        for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {
        //for(size_t s2(0), bf2_s(0); s2 < nShells; bf2_s+=n2, s2++) {
  
          n2 = basisSet_.shells[s2].size(); // Size of Shell 2

        for(size_t s3(0), bf3_s(0); s3 <= s1; bf3_s+=n3, s3++) {
        //for(size_t s3(0), bf3_s(0); s3 < nShells; bf3_s+=n3, s3++) {

          n3 = basisSet_.shells[s3].size(); // Size of Shell 3
          s4_max = (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4

        for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {
        //for(size_t s4(0), bf4_s(0); s4 < nShells; bf4_s+=n4, s4++, s1234++) {
  
          n4 = basisSet_.shells[s4].size(); // Size of Shell 4
  
          // Round Robbin work distribution
          #ifdef _OPENMP
          if( s1234 % nthreads != thread_id ) continue;
          #endif
  
          shls[0] = int(s1);
          shls[1] = int(s2);
          shls[2] = int(s3);
          shls[3] = int(s4);

          if (basisSet_.forceCart) {
            if(int2e_ipvip1ipvip2_cart(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
          } else {
            if(int2e_ipvip1ipvip2_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
          }

          ijkl = 0ul;
          auto nQuad = n1 * n2 * n3 * n4;
          for (l = 0ul, bf4 = bf4_s; l < n4; ++l, bf4++)
          for (k = 0ul, bf3 = bf3_s; k < n3; ++k, bf3++)
          for (j = 0ul, bf2 = bf2_s; j < n2; ++j, bf2++)
          for (i = 0ul, bf1 = bf1_s; i < n1; ++i, bf1++) {

            // (∇A∙∇B)(∇C∙∇D)(ij|kl)
            auto dAdotdBdCdotdD =
              buff[AxBxCxDx * nQuad + ijkl] + buff[AxBxCyDy * nQuad + ijkl] + buff[AxBxCzDz * nQuad + ijkl]
              + buff[AyByCxDx * nQuad + ijkl] + buff[AyByCyDy * nQuad + ijkl] + buff[AyByCzDz * nQuad + ijkl]
              + buff[AzBzCxDx * nQuad + ijkl] + buff[AzBzCyDy * nQuad + ijkl] + buff[AzBzCzDz * nQuad + ijkl];

            // (∇Ax∇B)(∇C∙∇D)(ijkl)
            auto dAcrossdB_xdCdotdD = buff[AyBzCxDx * nQuad + ijkl] - buff[AzByCxDx * nQuad + ijkl]
                                      + buff[AyBzCyDy * nQuad + ijkl] - buff[AzByCyDy * nQuad + ijkl]
                                      + buff[AyBzCzDz * nQuad + ijkl] - buff[AzByCzDz * nQuad + ijkl];

            auto dAcrossdB_ydCdotdD = -buff[AxBzCxDx * nQuad + ijkl] + buff[AzBxCxDx * nQuad + ijkl]
                                      - buff[AxBzCyDy * nQuad + ijkl] + buff[AzBxCyDy * nQuad + ijkl]
                                      - buff[AxBzCzDz * nQuad + ijkl] + buff[AzBxCzDz * nQuad + ijkl];

            auto dAcrossdB_zdCdotdD = buff[AxByCxDx * nQuad + ijkl] - buff[AyBxCxDx * nQuad + ijkl]
                                      + buff[AxByCyDy * nQuad + ijkl] - buff[AyBxCyDy * nQuad + ijkl]
                                      + buff[AxByCzDz * nQuad + ijkl] - buff[AyBxCzDz * nQuad + ijkl];

            // (∇A∙∇B)(∇Cx∇D)(ijkl)
            auto dAdotdBdCcrossdD_x = buff[AxBxCyDz * nQuad + ijkl] - buff[AxBxCzDy * nQuad + ijkl]
                                      + buff[AyByCyDz * nQuad + ijkl] - buff[AyByCzDy * nQuad + ijkl]
                                      + buff[AzBzCyDz * nQuad + ijkl] - buff[AzBzCzDy * nQuad + ijkl];

            auto dAdotdBdCcrossdD_y = -buff[AxBxCxDz * nQuad + ijkl] + buff[AxBxCzDx * nQuad + ijkl]
                                      - buff[AyByCxDz * nQuad + ijkl] + buff[AyByCzDx * nQuad + ijkl]
                                      - buff[AzBzCxDz * nQuad + ijkl] + buff[AzBzCzDx * nQuad + ijkl];

            auto dAdotdBdCcrossdD_z = buff[AxBxCxDy * nQuad + ijkl] - buff[AxBxCyDx * nQuad + ijkl]
                                      + buff[AyByCxDy * nQuad + ijkl] - buff[AyByCyDx * nQuad + ijkl]
                                      + buff[AzBzCxDy * nQuad + ijkl] - buff[AzBzCyDx * nQuad + ijkl];

            // (∇Ax∇B)(∇Cx∇D)(ijkl)
            auto dAcrossdB_xdCcrossdD_x = buff[AyBzCyDz * nQuad + ijkl] - buff[AzByCyDz * nQuad + ijkl]
                                          - buff[AyBzCzDy * nQuad + ijkl] + buff[AzByCzDy * nQuad + ijkl];

            auto dAcrossdB_xdCcrossdD_y = buff[AyBzCzDx * nQuad + ijkl] - buff[AzByCzDx * nQuad + ijkl]
                                          - buff[AyBzCxDz * nQuad + ijkl] + buff[AzByCxDz * nQuad + ijkl];

            auto dAcrossdB_xdCcrossdD_z = buff[AyBzCxDy * nQuad + ijkl] - buff[AzByCxDy * nQuad + ijkl]
                                          - buff[AyBzCyDx * nQuad + ijkl] + buff[AzByCyDx * nQuad + ijkl];

            auto dAcrossdB_ydCcrossdD_x = buff[AzBxCyDz * nQuad + ijkl] - buff[AxBzCyDz * nQuad + ijkl]
                                          - buff[AzBxCzDy * nQuad + ijkl] + buff[AxBzCzDy * nQuad + ijkl];

            auto dAcrossdB_ydCcrossdD_y = buff[AzBxCzDx * nQuad + ijkl] - buff[AxBzCzDx * nQuad + ijkl]
                                          - buff[AzBxCxDz * nQuad + ijkl] + buff[AxBzCxDz * nQuad + ijkl];

            auto dAcrossdB_ydCcrossdD_z = buff[AzBxCxDy * nQuad + ijkl] - buff[AxBzCxDy * nQuad + ijkl]
                                          - buff[AzBxCyDx * nQuad + ijkl] + buff[AxBzCyDx * nQuad + ijkl];

            auto dAcrossdB_zdCcrossdD_x = buff[AxByCyDz * nQuad + ijkl] - buff[AyBxCyDz * nQuad + ijkl]
                                          - buff[AxByCzDy * nQuad + ijkl] + buff[AyBxCzDy * nQuad + ijkl];

            auto dAcrossdB_zdCcrossdD_y = buff[AxByCzDx * nQuad + ijkl] - buff[AyBxCzDx * nQuad + ijkl]
                                          - buff[AxByCxDz * nQuad + ijkl] + buff[AyBxCxDz * nQuad + ijkl];

            auto dAcrossdB_zdCcrossdD_z = buff[AxByCxDy * nQuad + ijkl] - buff[AyBxCxDy * nQuad + ijkl]
                                          - buff[AxByCyDx * nQuad + ijkl] + buff[AyBxCyDx * nQuad + ijkl];

            //std::cout << std::scientific << std::setprecision(16);
            //std::cout << "(" << bf1 << "," << bf2 << "|" << bf3 << "," << bf4 << ")  ";
            //std::cout << dAcrossdB_xdCdotdD << std::endl;


            auto IJKL = bf1 + bf2 * NB + bf3 * NB2 + bf4 * NB3;
            auto IJLK = bf1 + bf2 * NB + bf4 * NB2 + bf3 * NB3;
            auto JIKL = bf2 + bf1 * NB + bf3 * NB2 + bf4 * NB3;
            auto JILK = bf2 + bf1 * NB + bf4 * NB2 + bf3 * NB3;
            auto KLIJ = bf3 + bf4 * NB + bf1 * NB2 + bf2 * NB3;
            auto LKIJ = bf4 + bf3 * NB + bf1 * NB2 + bf2 * NB3;
            auto KLJI = bf3 + bf4 * NB + bf2 * NB2 + bf1 * NB3;
            auto LKJI = bf4 + bf3 * NB + bf2 * NB2 + bf1 * NB3;

            //auto KLIJ = bf1 + bf2*NB + bf3*NB2 + bf4*NB3;
            //auto LKIJ = bf1 + bf2*NB + bf4*NB2 + bf3*NB3;
            //auto KLJI = bf2 + bf1*NB + bf3*NB2 + bf4*NB3;
            //auto LKJI = bf2 + bf1*NB + bf4*NB2 + bf3*NB3;
            //auto IJKL = bf3 + bf4*NB + bf1*NB2 + bf2*NB3;
            //auto IJLK = bf4 + bf3*NB + bf1*NB2 + bf2*NB3;
            //auto JIKL = bf3 + bf4*NB + bf2*NB2 + bf1*NB3;
            //auto JILK = bf4 + bf3*NB + bf2*NB2 + bf1*NB3;


            // (∇A∙∇B)(∇C∙∇D)(ij|kl) 
            // (ij|kl)
            (*this)[nERIRef].pointer()[IJKL] = dAdotdBdCdotdD;
            (*this)[nERIRef].pointer()[IJLK] = dAdotdBdCdotdD;
            (*this)[nERIRef].pointer()[JIKL] = dAdotdBdCdotdD;
            (*this)[nERIRef].pointer()[JILK] = dAdotdBdCdotdD;
            (*this)[nERIRef].pointer()[KLIJ] = dAdotdBdCdotdD;
            (*this)[nERIRef].pointer()[LKIJ] = dAdotdBdCdotdD;
            (*this)[nERIRef].pointer()[KLJI] = dAdotdBdCdotdD;
            (*this)[nERIRef].pointer()[LKJI] = dAdotdBdCdotdD;

            // (∇Ax∇B)_x(∇C∙∇D)(ijkl)
            (*this)[nERIRef + 1].pointer()[IJKL] = dAcrossdB_xdCdotdD;
            (*this)[nERIRef + 1].pointer()[IJLK] = dAcrossdB_xdCdotdD;
            (*this)[nERIRef + 1].pointer()[JIKL] = -dAcrossdB_xdCdotdD;
            (*this)[nERIRef + 1].pointer()[JILK] = -dAcrossdB_xdCdotdD;
            (*this)[nERIRef + 1].pointer()[KLIJ] = dAdotdBdCcrossdD_x;
            (*this)[nERIRef + 1].pointer()[LKIJ] = -dAdotdBdCcrossdD_x;
            (*this)[nERIRef + 1].pointer()[KLJI] = dAdotdBdCcrossdD_x;
            (*this)[nERIRef + 1].pointer()[LKJI] = -dAdotdBdCcrossdD_x;

            // (∇Ax∇B)_y(∇C∙∇D)(ijkl)
            (*this)[nERIRef + 2].pointer()[IJKL] = dAcrossdB_ydCdotdD;
            (*this)[nERIRef + 2].pointer()[IJLK] = dAcrossdB_ydCdotdD;
            (*this)[nERIRef + 2].pointer()[JIKL] = -dAcrossdB_ydCdotdD;
            (*this)[nERIRef + 2].pointer()[JILK] = -dAcrossdB_ydCdotdD;
            (*this)[nERIRef + 2].pointer()[KLIJ] = dAdotdBdCcrossdD_y;
            (*this)[nERIRef + 2].pointer()[LKIJ] = -dAdotdBdCcrossdD_y;
            (*this)[nERIRef + 2].pointer()[KLJI] = dAdotdBdCcrossdD_y;
            (*this)[nERIRef + 2].pointer()[LKJI] = -dAdotdBdCcrossdD_y;

            // (∇Ax∇B)_z(∇C∙∇D)(ijkl)
            (*this)[nERIRef + 3].pointer()[IJKL] = dAcrossdB_zdCdotdD;
            (*this)[nERIRef + 3].pointer()[IJLK] = dAcrossdB_zdCdotdD;
            (*this)[nERIRef + 3].pointer()[JIKL] = -dAcrossdB_zdCdotdD;
            (*this)[nERIRef + 3].pointer()[JILK] = -dAcrossdB_zdCdotdD;
            (*this)[nERIRef + 3].pointer()[KLIJ] = dAdotdBdCcrossdD_z;
            (*this)[nERIRef + 3].pointer()[LKIJ] = -dAdotdBdCcrossdD_z;
            (*this)[nERIRef + 3].pointer()[KLJI] = dAdotdBdCcrossdD_z;
            (*this)[nERIRef + 3].pointer()[LKJI] = -dAdotdBdCcrossdD_z;



            // (∇A∙∇B)(∇Cx∇D)_x(ijkl)
            (*this)[nERIRef + 4].pointer()[IJKL] = dAdotdBdCcrossdD_x;
            (*this)[nERIRef + 4].pointer()[IJLK] = -dAdotdBdCcrossdD_x;
            (*this)[nERIRef + 4].pointer()[JIKL] = dAdotdBdCcrossdD_x;
            (*this)[nERIRef + 4].pointer()[JILK] = -dAdotdBdCcrossdD_x;
            (*this)[nERIRef + 4].pointer()[KLIJ] = dAcrossdB_xdCdotdD;
            (*this)[nERIRef + 4].pointer()[LKIJ] = dAcrossdB_xdCdotdD;
            (*this)[nERIRef + 4].pointer()[KLJI] = -dAcrossdB_xdCdotdD;
            (*this)[nERIRef + 4].pointer()[LKJI] = -dAcrossdB_xdCdotdD;

            // (∇A∙∇B)(∇Cx∇D)_y(ijkl)
            (*this)[nERIRef + 5].pointer()[IJKL] = dAdotdBdCcrossdD_y;
            (*this)[nERIRef + 5].pointer()[IJLK] = -dAdotdBdCcrossdD_y;
            (*this)[nERIRef + 5].pointer()[JIKL] = dAdotdBdCcrossdD_y;
            (*this)[nERIRef + 5].pointer()[JILK] = -dAdotdBdCcrossdD_y;
            (*this)[nERIRef + 5].pointer()[KLIJ] = dAcrossdB_ydCdotdD;
            (*this)[nERIRef + 5].pointer()[LKIJ] = dAcrossdB_ydCdotdD;
            (*this)[nERIRef + 5].pointer()[KLJI] = -dAcrossdB_ydCdotdD;
            (*this)[nERIRef + 5].pointer()[LKJI] = -dAcrossdB_ydCdotdD;

            // (∇A∙∇B)(∇Cx∇D)_z(ijkl)
            (*this)[nERIRef + 6].pointer()[IJKL] = dAdotdBdCcrossdD_z;
            (*this)[nERIRef + 6].pointer()[IJLK] = -dAdotdBdCcrossdD_z;
            (*this)[nERIRef + 6].pointer()[JIKL] = dAdotdBdCcrossdD_z;
            (*this)[nERIRef + 6].pointer()[JILK] = -dAdotdBdCcrossdD_z;
            (*this)[nERIRef + 6].pointer()[KLIJ] = dAcrossdB_zdCdotdD;
            (*this)[nERIRef + 6].pointer()[LKIJ] = dAcrossdB_zdCdotdD;
            (*this)[nERIRef + 6].pointer()[KLJI] = -dAcrossdB_zdCdotdD;
            (*this)[nERIRef + 6].pointer()[LKJI] = -dAcrossdB_zdCdotdD;



            // (∇Ax∇B)_x(∇Cx∇D)_x(ijkl)
            (*this)[nERIRef + 7].pointer()[IJKL] = dAcrossdB_xdCcrossdD_x;
            (*this)[nERIRef + 7].pointer()[IJLK] = -dAcrossdB_xdCcrossdD_x;
            (*this)[nERIRef + 7].pointer()[JIKL] = -dAcrossdB_xdCcrossdD_x;
            (*this)[nERIRef + 7].pointer()[JILK] = dAcrossdB_xdCcrossdD_x;
            (*this)[nERIRef + 7].pointer()[KLIJ] = dAcrossdB_xdCcrossdD_x;
            (*this)[nERIRef + 7].pointer()[LKIJ] = -dAcrossdB_xdCcrossdD_x;
            (*this)[nERIRef + 7].pointer()[KLJI] = -dAcrossdB_xdCcrossdD_x;
            (*this)[nERIRef + 7].pointer()[LKJI] = dAcrossdB_xdCcrossdD_x;

            // (∇Ax∇B)_x(∇Cx∇D)_y(ijkl)
            (*this)[nERIRef + 8].pointer()[IJKL] = dAcrossdB_xdCcrossdD_y;
            (*this)[nERIRef + 8].pointer()[IJLK] = -dAcrossdB_xdCcrossdD_y;
            (*this)[nERIRef + 8].pointer()[JIKL] = -dAcrossdB_xdCcrossdD_y;
            (*this)[nERIRef + 8].pointer()[JILK] = dAcrossdB_xdCcrossdD_y;
            (*this)[nERIRef + 8].pointer()[KLIJ] = dAcrossdB_ydCcrossdD_x;
            (*this)[nERIRef + 8].pointer()[LKIJ] = -dAcrossdB_ydCcrossdD_x;
            (*this)[nERIRef + 8].pointer()[KLJI] = -dAcrossdB_ydCcrossdD_x;
            (*this)[nERIRef + 8].pointer()[LKJI] = dAcrossdB_ydCcrossdD_x;

            // (∇Ax∇B)_x(∇Cx∇D)_z(ijkl)
            (*this)[nERIRef + 9].pointer()[IJKL] = dAcrossdB_xdCcrossdD_z;
            (*this)[nERIRef + 9].pointer()[IJLK] = -dAcrossdB_xdCcrossdD_z;
            (*this)[nERIRef + 9].pointer()[JIKL] = -dAcrossdB_xdCcrossdD_z;
            (*this)[nERIRef + 9].pointer()[JILK] = dAcrossdB_xdCcrossdD_z;
            (*this)[nERIRef + 9].pointer()[KLIJ] = dAcrossdB_zdCcrossdD_x;
            (*this)[nERIRef + 9].pointer()[LKIJ] = -dAcrossdB_zdCcrossdD_x;
            (*this)[nERIRef + 9].pointer()[KLJI] = -dAcrossdB_zdCcrossdD_x;
            (*this)[nERIRef + 9].pointer()[LKJI] = dAcrossdB_zdCcrossdD_x;



            // (∇Ax∇B)_y(∇Cx∇D)_x(ijkl)
            (*this)[nERIRef + 10].pointer()[IJKL] = dAcrossdB_ydCcrossdD_x;
            (*this)[nERIRef + 10].pointer()[IJLK] = -dAcrossdB_ydCcrossdD_x;
            (*this)[nERIRef + 10].pointer()[JIKL] = -dAcrossdB_ydCcrossdD_x;
            (*this)[nERIRef + 10].pointer()[JILK] = dAcrossdB_ydCcrossdD_x;
            (*this)[nERIRef + 10].pointer()[KLIJ] = dAcrossdB_xdCcrossdD_y;
            (*this)[nERIRef + 10].pointer()[LKIJ] = -dAcrossdB_xdCcrossdD_y;
            (*this)[nERIRef + 10].pointer()[KLJI] = -dAcrossdB_xdCcrossdD_y;
            (*this)[nERIRef + 10].pointer()[LKJI] = dAcrossdB_xdCcrossdD_y;

            // (∇Ax∇B)_y(∇Cx∇D)_y(ijkl)
            (*this)[nERIRef + 11].pointer()[IJKL] = dAcrossdB_ydCcrossdD_y;
            (*this)[nERIRef + 11].pointer()[IJLK] = -dAcrossdB_ydCcrossdD_y;
            (*this)[nERIRef + 11].pointer()[JIKL] = -dAcrossdB_ydCcrossdD_y;
            (*this)[nERIRef + 11].pointer()[JILK] = dAcrossdB_ydCcrossdD_y;
            (*this)[nERIRef + 11].pointer()[KLIJ] = dAcrossdB_ydCcrossdD_y;
            (*this)[nERIRef + 11].pointer()[LKIJ] = -dAcrossdB_ydCcrossdD_y;
            (*this)[nERIRef + 11].pointer()[KLJI] = -dAcrossdB_ydCcrossdD_y;
            (*this)[nERIRef + 11].pointer()[LKJI] = dAcrossdB_ydCcrossdD_y;

            // (∇Ax∇B)_y(∇Cx∇D)_z(ijkl)
            (*this)[nERIRef + 12].pointer()[IJKL] = dAcrossdB_ydCcrossdD_z;
            (*this)[nERIRef + 12].pointer()[IJLK] = -dAcrossdB_ydCcrossdD_z;
            (*this)[nERIRef + 12].pointer()[JIKL] = -dAcrossdB_ydCcrossdD_z;
            (*this)[nERIRef + 12].pointer()[JILK] = dAcrossdB_ydCcrossdD_z;
            (*this)[nERIRef + 12].pointer()[KLIJ] = dAcrossdB_zdCcrossdD_y;
            (*this)[nERIRef + 12].pointer()[LKIJ] = -dAcrossdB_zdCcrossdD_y;
            (*this)[nERIRef + 12].pointer()[KLJI] = -dAcrossdB_zdCcrossdD_y;
            (*this)[nERIRef + 12].pointer()[LKJI] = dAcrossdB_zdCcrossdD_y;



            // (∇Ax∇B)_z(∇Cx∇D)_x(ijkl)
            (*this)[nERIRef + 13].pointer()[IJKL] = dAcrossdB_zdCcrossdD_x;
            (*this)[nERIRef + 13].pointer()[IJLK] = -dAcrossdB_zdCcrossdD_x;
            (*this)[nERIRef + 13].pointer()[JIKL] = -dAcrossdB_zdCcrossdD_x;
            (*this)[nERIRef + 13].pointer()[JILK] = dAcrossdB_zdCcrossdD_x;
            (*this)[nERIRef + 13].pointer()[KLIJ] = dAcrossdB_xdCcrossdD_z;
            (*this)[nERIRef + 13].pointer()[LKIJ] = -dAcrossdB_xdCcrossdD_z;
            (*this)[nERIRef + 13].pointer()[KLJI] = -dAcrossdB_xdCcrossdD_z;
            (*this)[nERIRef + 13].pointer()[LKJI] = dAcrossdB_xdCcrossdD_z;

            // (∇Ax∇B)_z(∇Cx∇D)_y(ijkl)
            (*this)[nERIRef + 14].pointer()[IJKL] = dAcrossdB_zdCcrossdD_y;
            (*this)[nERIRef + 14].pointer()[IJLK] = -dAcrossdB_zdCcrossdD_y;
            (*this)[nERIRef + 14].pointer()[JIKL] = -dAcrossdB_zdCcrossdD_y;
            (*this)[nERIRef + 14].pointer()[JILK] = dAcrossdB_zdCcrossdD_y;
            (*this)[nERIRef + 14].pointer()[KLIJ] = dAcrossdB_ydCcrossdD_z;
            (*this)[nERIRef + 14].pointer()[LKIJ] = -dAcrossdB_ydCcrossdD_z;
            (*this)[nERIRef + 14].pointer()[KLJI] = -dAcrossdB_ydCcrossdD_z;
            (*this)[nERIRef + 14].pointer()[LKJI] = dAcrossdB_ydCcrossdD_z;

            // (∇Ax∇B)_z(∇Cx∇D)_z(ijkl)
            (*this)[nERIRef + 15].pointer()[IJKL] = dAcrossdB_zdCcrossdD_z;
            (*this)[nERIRef + 15].pointer()[IJLK] = -dAcrossdB_zdCcrossdD_z;
            (*this)[nERIRef + 15].pointer()[JIKL] = -dAcrossdB_zdCcrossdD_z;
            (*this)[nERIRef + 15].pointer()[JILK] = dAcrossdB_zdCcrossdD_z;
            (*this)[nERIRef + 15].pointer()[KLIJ] = dAcrossdB_zdCcrossdD_z;
            (*this)[nERIRef + 15].pointer()[LKIJ] = -dAcrossdB_zdCcrossdD_z;
            (*this)[nERIRef + 15].pointer()[KLJI] = -dAcrossdB_zdCcrossdD_z;
            (*this)[nERIRef + 15].pointer()[LKJI] = dAcrossdB_zdCcrossdD_z;


            ijkl++;
  
          }; // ijkl loop
        }; // s4
        }; // s3
        }; // s2
        }; // s1
  
      }; // omp region
 

#ifdef _REPORT_INCORE_CONTRACTION_TIMINGS 
      auto durERIDCSSSS = tock(topERIDCSSSS);
      std::cout << "Libcint-ERI-Dirac-Coulomb-SSSS duration   = " << durERIDCSSSS << std::endl;
#endif

#if 0 
      std::cout << std::scientific << std::setprecision(16);

      std::cout << "(∇A∙∇B)(∇C∙∇D)(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_x(∇C∙∇D)(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+1](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_y(∇C∙∇D)(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+2](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_z(∇C∙∇D)(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+3](i, j, k, l) << std::endl;
      };

      std::cout << "(∇A∙∇B)(∇Cx∇D)_x(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+4](i, j, k, l) << std::endl;
      };

      std::cout << "(∇A∙∇B)(∇Cx∇D)_y(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+5](i, j, k, l) << std::endl;
      };

      std::cout << "(∇A∙∇B)(∇Cx∇D)_z(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+6](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_x(∇Cx∇D)_x(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+7](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_x(∇Cx∇D)_y(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+8](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_x(∇Cx∇D)_z(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+9](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_y(∇Cx∇D)_x(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+10](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_y(∇Cx∇D)_y(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+11](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_y(∇Cx∇D)_z(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+12](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_z(∇Cx∇D)_x(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+13](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_z(∇Cx∇D)_y(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+14](i, j, k, l) << std::endl;
      };

      std::cout << "(∇Ax∇B)_z(∇Cx∇D)_z(ij|kl)" << std::endl;
      for(auto i = 0ul; i < NB; i++)
      for(auto j = 0ul; j < NB; j++)
      for(auto k = 0ul; k < NB; k++)
      for(auto l = 0ul; l < NB; l++){
        std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
        std::cout << (*this)[nERIRef+15](i, j, k, l) << std::endl;
      };



#endif


    } // Dirac-Coulomb (SSSS) ∇_i∇_j∇_k∇_l(ij|kl)





    /* Gauge Integrals */
    if (hamiltonianOptions.Gauge) {

      auto nERIRef = 0;
      if (hamiltonianOptions.DiracCoulomb) nERIRef +=4; // Dirac-Coulomb
      if (hamiltonianOptions.Gaunt) nERIRef += 19; // Gaunt
      if (hamiltonianOptions.DiracCoulombSSSS) nERIRef += 16; // Dirac-Coulomb-SSSS

#ifdef _REPORT_INCORE_INTEGRAL_TIMINGS 
      auto topERIGauge = tick();
#endif

#pragma omp parallel
      {
        int thread_id = GetThreadID();

        size_t n1, n2, n3, n4, i, j, k, l, ijkl, bf1, bf2, bf3, bf4;
        size_t s4_max;
        int shls[4];
        double *buffr1 = buffAll+buffN4*thread_id;
        double *buffr2 = buffAll+nthreads*buffN4+buffN4*thread_id;
        double *cache = cacheAll+cache_size*thread_id;


        for (size_t s1(0), bf1_s(0), s1234(0); s1 < nShells; bf1_s += n1, s1++) {

          n1 = basisSet_.shells[s1].size(); // Size of Shell 1

        for (size_t s2(0), bf2_s(0); s2 < nShells; bf2_s += n2, s2++) {

          n2 = basisSet_.shells[s2].size(); // Size of Shell 2

        //for (size_t s3(0), bf3_s(0); s3 <= s2; bf3_s += n3, s3++) {
        for (size_t s3(0), bf3_s(0); s3 < nShells; bf3_s += n3, s3++) {

          n3 = basisSet_.shells[s3].size(); // Size of Shell 3
          s4_max = (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4

        for (size_t s4(0), bf4_s(0); s4 < nShells; bf4_s += n4, s4++, s1234++) {

          n4 = basisSet_.shells[s4].size(); // Size of Shell 4

          // Round Robbin work distribution
#ifdef _OPENMP
          if (s1234 % nthreads != thread_id) continue;
#endif


          shls[0] = int(s1);
          shls[1] = int(s2);
          shls[2] = int(s3);
          shls[3] = int(s4);

          if (basisSet_.forceCart) {
            if(int2e_gauge_r1_sps1sps2_cart(buffr1, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0 and
               int2e_gauge_r2_sps1sps2_cart(buffr2, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0)
              continue;
          } else {
            //∇A∇C
            //int2e_gauge_r1_sps1sps2_sph(buffr1, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);
            //int2e_gauge_r2_sps1sps2_sph(buffr2, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache);
	    
            //∇A∇D
            //int2e_gauge_r1_sps1ssp2_sph(buffr1, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache); 
            //int2e_gauge_r2_sps1ssp2_sph(buffr2, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache); 

            //∇B∇C
            int2e_gauge_r1_ssp1sps2_sph(buffr1, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache); 
            int2e_gauge_r2_ssp1sps2_sph(buffr2, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache); 
          }
          
          ijkl = 0ul;
          auto nQuad = n1 * n2 * n3 * n4;
          for (l = 0ul, bf4 = bf4_s; l < n4; ++l, bf4++) 
          for (k = 0ul, bf3 = bf3_s; k < n3; ++k, bf3++)
          for (j = 0ul, bf2 = bf2_s; j < n2; ++j, bf2++)
          for (i = 0ul, bf1 = bf1_s; i < n1; ++i, bf1++) {
          //for (i = 0ul, bf1 = bf1_s; i < n1; ++i, bf1++)
          //for (j = 0ul, bf2 = bf2_s; j < n2; ++j, bf2++)
          //for (k = 0ul, bf3 = bf3_s; k < n3; ++k, bf3++)
          //for (l = 0ul, bf4 = bf4_s; l < n4; ++l, bf4++) {

            auto IJKL = bf1 + bf2 * NB + bf3 * NB2 + bf4 * NB3;
            auto LKJI = bf4 + bf3 * NB + bf2 * NB2 + bf1 * NB3;


            // σ_x * σ_x     0
            // σ_y * σ_x     1
            // σ_z * σ_x     2
            // I   * σ_x     3
            // σ_x * σ_y     4
            // σ_y * σ_y     5
            // σ_z * σ_y     6
            // I   * σ_y     7
            // σ_x * σ_z     8
            // σ_y * σ_z     9
            // σ_z * σ_z     10
            // I   * σ_z     11
            // σ_x * I       12
            // σ_y * I       13
            // σ_z * I       14
            // I   * I       15


            // (ss)
            (*this)[nERIRef + 0].pointer()[IJKL] = buffr1[15*nQuad+ijkl] - buffr2[15*nQuad+ijkl];
            (*this)[nERIRef + 0].pointer()[LKJI] = buffr1[15*nQuad+ijkl] - buffr2[15*nQuad+ijkl];

            // (sσ)_x
            (*this)[nERIRef + 1].pointer()[IJKL] = buffr1[3*nQuad+ijkl] - buffr2[3*nQuad+ijkl];
            (*this)[nERIRef + 1].pointer()[LKJI] = buffr2[12*nQuad+ijkl] - buffr1[12*nQuad+ijkl];

            // (sσ)_y
            (*this)[nERIRef + 2].pointer()[IJKL] = buffr1[7*nQuad+ijkl] - buffr2[7*nQuad+ijkl];
            (*this)[nERIRef + 2].pointer()[LKJI] = buffr2[13*nQuad+ijkl] - buffr1[13*nQuad+ijkl];

            // (sσ)_z
            (*this)[nERIRef + 3].pointer()[IJKL] = buffr1[11*nQuad+ijkl] - buffr2[11*nQuad+ijkl];
            (*this)[nERIRef + 3].pointer()[LKJI] = buffr2[14*nQuad+ijkl] - buffr1[14*nQuad+ijkl];

            // (σs)_x
            (*this)[nERIRef + 4].pointer()[IJKL] = buffr1[12*nQuad+ijkl] - buffr2[12*nQuad+ijkl];
            (*this)[nERIRef + 4].pointer()[LKJI] = buffr2[3*nQuad+ijkl] - buffr1[3*nQuad+ijkl];

            // (σs)_y
            (*this)[nERIRef + 5].pointer()[IJKL] = buffr1[13*nQuad+ijkl] - buffr2[13*nQuad+ijkl];
            (*this)[nERIRef + 5].pointer()[LKJI] = buffr2[7*nQuad+ijkl] - buffr1[7*nQuad+ijkl];

            // (σs)_z
            (*this)[nERIRef + 6].pointer()[IJKL] = buffr1[14*nQuad+ijkl] - buffr2[14*nQuad+ijkl];
            (*this)[nERIRef + 6].pointer()[LKJI] = buffr2[11*nQuad+ijkl] - buffr1[11*nQuad+ijkl];

            // (σ∙σ)
            (*this)[nERIRef + 7].pointer()[IJKL] =-(buffr1[ijkl] - buffr2[ijkl])
                                                  -(buffr1[5*nQuad+ijkl] - buffr2[5*nQuad+ijkl])
                                                  -(buffr1[10*nQuad+ijkl] - buffr2[10*nQuad+ijkl]);
            (*this)[nERIRef + 7].pointer()[LKJI] =-(buffr1[ijkl] - buffr2[ijkl])
                                                  -(buffr1[5*nQuad+ijkl] - buffr2[5*nQuad+ijkl])
                                                  -(buffr1[10*nQuad+ijkl] - buffr2[10*nQuad+ijkl]);

            // (σxσ)_x = σ_y*σ_z - σ_z*σ_y
            (*this)[nERIRef + 8].pointer()[IJKL] =-(buffr1[9*nQuad+ijkl] - buffr2[9*nQuad+ijkl]) 
                                                  +(buffr1[6*nQuad+ijkl] - buffr2[6*nQuad+ijkl]);
            (*this)[nERIRef + 8].pointer()[LKJI] = (buffr1[9*nQuad+ijkl] - buffr2[9*nQuad+ijkl]) 
                                                  -(buffr1[6*nQuad+ijkl] - buffr2[6*nQuad+ijkl]);

            // (σxσ)_y = σ_z*σ_x - σ_x*σ_z
            (*this)[nERIRef + 9].pointer()[IJKL] =-(buffr1[2*nQuad+ijkl] - buffr2[2*nQuad+ijkl])
                                                  +(buffr1[8*nQuad+ijkl] - buffr2[8*nQuad+ijkl]);
            (*this)[nERIRef + 9].pointer()[LKJI] = (buffr1[2*nQuad+ijkl] - buffr2[2*nQuad+ijkl])
                                                  -(buffr1[8*nQuad+ijkl] - buffr2[8*nQuad+ijkl]);

            // (σxσ)_z = σ_x*σ_y - σ_y*σ_x
            (*this)[nERIRef + 10].pointer()[IJKL] =-(buffr1[4*nQuad+ijkl] - buffr2[4*nQuad+ijkl]) 
                                                   +(buffr1[1*nQuad+ijkl] - buffr2[1*nQuad+ijkl]);
            (*this)[nERIRef + 10].pointer()[LKJI] = (buffr1[4*nQuad+ijkl] - buffr2[4*nQuad+ijkl]) 
                                                   -(buffr1[1*nQuad+ijkl] - buffr2[1*nQuad+ijkl]);

            // (σσ) (xx - yy - zz)
            (*this)[nERIRef + 11].pointer()[IJKL] =-(buffr1[ijkl] - buffr2[ijkl])
                                                   +(buffr1[5*nQuad+ijkl] - buffr2[5*nQuad+ijkl])
                                                   +(buffr1[10*nQuad+ijkl] - buffr2[10*nQuad+ijkl]);
            (*this)[nERIRef + 11].pointer()[LKJI] =-(buffr1[ijkl] - buffr2[ijkl])
                                                   +(buffr1[5*nQuad+ijkl] - buffr2[5*nQuad+ijkl])
                                                   +(buffr1[10*nQuad+ijkl] - buffr2[10*nQuad+ijkl]);

            // (σσ) (-xx + yy - zz)
            (*this)[nERIRef + 12].pointer()[IJKL] = (buffr1[ijkl] - buffr2[ijkl])
                                                   -(buffr1[5*nQuad+ijkl] - buffr2[5*nQuad+ijkl])
                                                   +(buffr1[10*nQuad+ijkl] - buffr2[10*nQuad+ijkl]);
            (*this)[nERIRef + 12].pointer()[LKJI] = (buffr1[ijkl] - buffr2[ijkl])
                                                   -(buffr1[5*nQuad+ijkl] - buffr2[5*nQuad+ijkl])
                                                   +(buffr1[10*nQuad+ijkl] - buffr2[10*nQuad+ijkl]);

            // (σσ) (-xx - yy + zz)
            (*this)[nERIRef + 13].pointer()[IJKL] = (buffr1[ijkl] - buffr2[ijkl])
                                                   +(buffr1[5*nQuad+ijkl] - buffr2[5*nQuad+ijkl])
                                                   -(buffr1[10*nQuad+ijkl] - buffr2[10*nQuad+ijkl]);
            (*this)[nERIRef + 13].pointer()[LKJI] = (buffr1[ijkl] - buffr2[ijkl])
                                                   +(buffr1[5*nQuad+ijkl] - buffr2[5*nQuad+ijkl])
                                                   -(buffr1[10*nQuad+ijkl] - buffr2[10*nQuad+ijkl]);

            //  σ_x*σ_y + σ_y*σ_x
            (*this)[nERIRef + 14].pointer()[IJKL] =-(buffr1[4*nQuad+ijkl] - buffr2[4*nQuad+ijkl])
                                                   -(buffr1[1*nQuad+ijkl] - buffr2[1*nQuad+ijkl]);
            (*this)[nERIRef + 14].pointer()[LKJI] =-(buffr1[4*nQuad+ijkl] - buffr2[4*nQuad+ijkl])
                                                   -(buffr1[1*nQuad+ijkl] - buffr2[1*nQuad+ijkl]);

            //  σ_z*σ_x + σ_x*σ_z
            (*this)[nERIRef + 15].pointer()[IJKL] =-(buffr1[2*nQuad+ijkl] - buffr2[2*nQuad+ijkl])
                                                   -(buffr1[8*nQuad+ijkl] - buffr2[8*nQuad+ijkl]);
            (*this)[nERIRef + 15].pointer()[LKJI] =-(buffr1[2*nQuad+ijkl] - buffr2[2*nQuad+ijkl])
                                                   -(buffr1[8*nQuad+ijkl] - buffr2[8*nQuad+ijkl]);

            //  σ_y*σ_z + σ_z*σ_y
            (*this)[nERIRef + 16].pointer()[IJKL] =-(buffr1[9*nQuad+ijkl] - buffr2[9*nQuad+ijkl])
                                                   -(buffr1[6*nQuad+ijkl] - buffr2[6*nQuad+ijkl]);
            (*this)[nERIRef + 16].pointer()[LKJI] =-(buffr1[9*nQuad+ijkl] - buffr2[9*nQuad+ijkl])
                                                   -(buffr1[6*nQuad+ijkl] - buffr2[6*nQuad+ijkl]);

            // σ_x*σ_x
            (*this)[nERIRef + 17].pointer()[IJKL] = -(buffr1[ijkl] - buffr2[ijkl]);
            (*this)[nERIRef + 17].pointer()[LKJI] = -(buffr1[ijkl] - buffr2[ijkl]);

            // σ_x*σ_y
            (*this)[nERIRef + 18].pointer()[IJKL] = -(buffr1[4*nQuad+ijkl] - buffr2[4*nQuad+ijkl]);
            (*this)[nERIRef + 18].pointer()[LKJI] = -(buffr1[1*nQuad+ijkl] - buffr2[1*nQuad+ijkl]);

            // σ_x*σ_z
            (*this)[nERIRef + 19].pointer()[IJKL] = -(buffr1[8*nQuad+ijkl] - buffr2[8*nQuad+ijkl]);
            (*this)[nERIRef + 19].pointer()[LKJI] = -(buffr1[2*nQuad+ijkl] - buffr2[2*nQuad+ijkl]);

            // σ_y*σ_x
            (*this)[nERIRef + 20].pointer()[IJKL] = -(buffr1[1*nQuad+ijkl] - buffr2[1*nQuad+ijkl]);
            (*this)[nERIRef + 20].pointer()[LKJI] = -(buffr1[4*nQuad+ijkl] - buffr2[4*nQuad+ijkl]);

            // σ_y*σ_y
            (*this)[nERIRef + 21].pointer()[IJKL] = -(buffr1[5*nQuad+ijkl] - buffr2[5*nQuad+ijkl]);
            (*this)[nERIRef + 21].pointer()[LKJI] = -(buffr1[5*nQuad+ijkl] - buffr2[5*nQuad+ijkl]);

            // σ_y*σ_z
            (*this)[nERIRef + 22].pointer()[IJKL] = -(buffr1[9*nQuad+ijkl] - buffr2[9*nQuad+ijkl]);
            (*this)[nERIRef + 22].pointer()[LKJI] = -(buffr1[6*nQuad+ijkl] - buffr2[6*nQuad+ijkl]);

            // σ_z*σ_x
            (*this)[nERIRef + 23].pointer()[IJKL] = -(buffr1[2*nQuad+ijkl] - buffr2[2*nQuad+ijkl]);
            (*this)[nERIRef + 23].pointer()[LKJI] = -(buffr1[8*nQuad+ijkl] - buffr2[8*nQuad+ijkl]);

            // σ_z*σ_y
            (*this)[nERIRef + 24].pointer()[IJKL] = -(buffr1[6*nQuad+ijkl] - buffr2[6*nQuad+ijkl]);
            (*this)[nERIRef + 24].pointer()[LKJI] = -(buffr1[9*nQuad+ijkl] - buffr2[9*nQuad+ijkl]);

            // σ_z*σ_z
            (*this)[nERIRef + 25].pointer()[IJKL] = -(buffr1[10*nQuad+ijkl] - buffr2[10*nQuad+ijkl]);
            (*this)[nERIRef + 25].pointer()[LKJI] = -(buffr1[10*nQuad+ijkl] - buffr2[10*nQuad+ijkl]);

            ijkl++;
          } //ijkl loop

        }
        }
        }
        }
      }


#ifdef _REPORT_INCORE_CONTRACTION_TIMINGS 
      auto durERIGauge = tock(topERIGauge);
      std::cout << "Libcint-ERI-Gauge duration   = " << durERIGauge << std::endl;
#endif

#if 0
      // LSSL
      for(auto index=0; index<26; index++){
        std::cout << "Libcint ("<<index<<")(ij|kl)" << std::endl;
        for(auto i = 0ul; i < NB; i++)
        for(auto j = 0ul; j < NB; j++)
        for(auto k = 0ul; k < NB; k++)
        for(auto l = 0ul; l < NB; l++){
          std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
          std::cout << (*this)[nERIRef+index](i, j, k, l) << std::endl;
        };
      };
#endif



    } // Gauge



    CQMemManager::get().free(cacheAll, buffAll, env, bas, atm);



#if 0
    prettyPrintSmart(std::cout,"Rank-2 ERI00 ∇∇(ab|cd)",(*this)[0].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI01 ∇∇(ab|cd)",(*this)[1].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI02 ∇∇(ab|cd)",(*this)[2].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI03 ∇∇(ab|cd)",(*this)[3].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI04 ∇∇(ab|cd)",(*this)[4].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI05 ∇∇(ab|cd)",(*this)[5].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI06 ∇∇(ab|cd)",(*this)[6].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI07 ∇∇(ab|cd)",(*this)[7].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI08 ∇∇(ab|cd)",(*this)[8].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI09 ∇∇(ab|cd)",(*this)[9].pointer(),  NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI10 ∇∇(ab|cd)",(*this)[10].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI11 ∇∇(ab|cd)",(*this)[11].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI12 ∇∇(ab|cd)",(*this)[12].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI13 ∇∇(ab|cd)",(*this)[13].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI14 ∇∇(ab|cd)",(*this)[14].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI15 ∇∇(ab|cd)",(*this)[15].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI16 ∇∇(ab|cd)",(*this)[16].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI17 ∇∇(ab|cd)",(*this)[17].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI18 ∇∇(ab|cd)",(*this)[18].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI19 ∇∇(ab|cd)",(*this)[19].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI20 ∇∇(ab|cd)",(*this)[20].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI21 ∇∇(ab|cd)",(*this)[21].pointer(), NB*NB,NB*NB,NB*NB);
    prettyPrintSmart(std::cout,"Rank-2 ERI22 ∇∇(ab|cd)",(*this)[22].pointer(), NB*NB,NB*NB,NB*NB);
#endif

  }; // InCore4indexRelERI<double>::computeERICINT





  template <>
  void InCore4indexTPI<dcomplex>::computeERIGCCINT(BasisSet&, Molecule&,
      EMPerturbation&, OPERATOR, const HamiltonianOptions&) {
    CErr("Only real GTOs are allowed",std::cout);
  };

  template <>
  void InCore4indexTPI<double>::computeERIGCCINT(BasisSet &originalBasisSet, Molecule &molecule_,
      EMPerturbation&, OPERATOR, const HamiltonianOptions&) {

    if (originalBasisSet.forceCart)
      CErr("Libcint + cartesian GTO NYI.");

    BasisSet basisSet_ = originalBasisSet.groupGeneralContractionBasis();

    size_t buffSize = std::max_element(basisSet_.shells.begin(),
                                       basisSet_.shells.end(),
                                       [](libint2::Shell &a, libint2::Shell &b) {
                                         return a.size() < b.size();
                                       })->size();

    int nAtoms = molecule_.nAtoms;
    int nShells = basisSet_.nShell;

    // ATM_SLOTS = 6; BAS_SLOTS = 8;
    int *atm = CQMemManager::get().malloc<int>(nAtoms * ATM_SLOTS);
    int *bas = CQMemManager::get().malloc<int>(nShells * BAS_SLOTS);
    double *env = CQMemManager::get().malloc<double>(basisSet_.getLibcintEnvLength(molecule_));


    basisSet_.setLibcintEnv(molecule_, atm, bas, env);


    size_t cache_size = 0;
    for (int i = 0; i < nShells; i++) {
      size_t n;
      int shls[4]{i,i,i,i};
      if (basisSet_.forceCart) {
        n = int2e_cart(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
      } else {
        n = int2e_sph(nullptr, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, nullptr);
      }
      cache_size = std::max(cache_size, n);
    }

    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();

    // Allocate and zero out ERIs
    size_t NB  = basisSet_.nBasis;
    size_t NB2 = NB*NB;
    size_t NB3 = NB2*NB;
    size_t NB4 = NB2*NB2;

    InCore4indexTPI<double>::clear();


    // Get threads result buffer
    size_t buffN4 = buffSize*buffSize*buffSize*buffSize;
    double *buffAll = CQMemManager::get().malloc<double>(buffN4*nthreads);
    double *cacheAll = CQMemManager::get().malloc<double>(cache_size*nthreads);

    std::cout<<"Using Libcint "<<std::endl;

#ifdef _REPORT_INCORE_INTEGRAL_TIMINGS 
    auto topERI4 = tick();
#endif

    #pragma omp parallel
    {
      int thread_id = GetThreadID();

      size_t n1,n2,n3,n4,i,j,k,l,ijkl,bf1,bf2,bf3,bf4;
      size_t s4_max;
      int shls[4];
      double *buff = buffAll+buffN4*thread_id;
      double *cache = cacheAll+cache_size*thread_id;

      for(size_t s1(0), bf1_s(0), s1234(0); s1 < nShells;
          bf1_s+=n1, s1++) {

        n1 = basisSet_.shells[s1].size(); // Size of Shell 1

      for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {

        n2 = basisSet_.shells[s2].size(); // Size of Shell 2

      for(size_t s3(0), bf3_s(0); s3 <= s1; bf3_s+=n3, s3++) {

        n3 = basisSet_.shells[s3].size(); // Size of Shell 3
        s4_max = (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4

      for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {

        n4 = basisSet_.shells[s4].size(); // Size of Shell 4

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s1234 % nthreads != thread_id ) continue;
        #endif

        shls[0] = int(s1);
        shls[1] = int(s2);
        shls[2] = int(s3);
        shls[3] = int(s4);

        if (basisSet_.forceCart) {
          if(int2e_cart(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
        } else {
          if(int2e_sph(buff, nullptr, shls, atm, nAtoms, bas, nShells, env, nullptr, cache)==0) continue;
        }

        // permutational symmetry
        ijkl = 0ul;
        for(l = 0ul, bf4 = bf4_s ; l < n4; ++l, bf4++)
        for(k = 0ul, bf3 = bf3_s ; k < n3; ++k, bf3++)
        for(j = 0ul, bf2 = bf2_s ; j < n2; ++j, bf2++)
        for(i = 0ul, bf1 = bf1_s ; i < n1; ++i, bf1++) {

          // (12 | 34)
          (*this)(bf1, bf2, bf3, bf4) = buff[ijkl];
          // (12 | 43)
          (*this)(bf1, bf2, bf4, bf3) = buff[ijkl];
          // (21 | 34)
          (*this)(bf2, bf1, bf3, bf4) = buff[ijkl];
          // (21 | 43)
          (*this)(bf2, bf1, bf4, bf3) = buff[ijkl];
          // (34 | 12)
          (*this)(bf3, bf4, bf1, bf2) = buff[ijkl];
          // (43 | 12)
          (*this)(bf4, bf3, bf1, bf2) = buff[ijkl];
          // (34 | 21)
          (*this)(bf3, bf4, bf2, bf1) = buff[ijkl];
          // (43 | 21)
          (*this)(bf4, bf3, bf2, bf1) = buff[ijkl];

          ijkl++;

        }; // ijkl loop

      }; // s4
      }; // s3
      }; // s2
      }; // s1

    }; // omp region


#ifdef _REPORT_INCORE_INTEGRAL_TIMINGS 
    auto durERI4 = tock(topERI4);
    std::cout << "Libcint-ERI4 duration   = " << durERI4 << std::endl;
#endif

    CQMemManager::get().free(cacheAll, buffAll, env, bas, atm);

#ifdef __DEBUGERI__
    // Debug output of the ERIs
    std::cout << std::scientific << std::setprecision(16);
    std::cout << "Libcint ERI (ab|cd)" << std::endl;
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++)
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout << (*this)(i, j, k, l) << std::endl;
    };
#endif // __DEBUGERI__

  }; // InCore4indexRelERI<double>::computeERIGCCINT


}; // namespace ChronusQ

//#endif
