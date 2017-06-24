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

#include <particleintegrals/inhouseaointegral.hpp>
#include <cqlinalg.hpp>
#include <cqlinalg/blasutil.hpp>
#include <util/matout.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>

#include <util/threads.hpp>
#include <chrono>

//#define _DEBUGORTHO
//#define _DEBUGERI
//#define _DEBUGGIAOERI //SS
//#define _DEBUGGIAOONEE //SS 
#define bottomupGIAO //SS
// Debug directives

namespace ChronusQ {
 
  /**
   *  \brief Allocate, compute and store the full rank-4 complex ERI tensor using
   *  in house GIAO code over the CGTO basis.
   */ 
  template <>
  void InCore4indexTPI<dcomplex>::computeERINR(BasisSet &basisSet, BasisSet &basisSet2, 
      Molecule&, EMPerturbation &emPert, OPERATOR op, const HamiltonianOptions &options) {

    if (&basisSet != &basisSet2)
      CErr("Only same basis is allowed in InCore4indexTPI<dcomplex>",std::cout);

    if (op != ELECTRON_REPULSION)
      CErr("Only Electron repulsion integrals in InCore4indexTPI<dcomplex>",std::cout);
    if (options.basisType == REAL_GTO)
      CErr("Real GTOs are not allowed in InCore4indexTPI<dcomplex>",std::cout);
    if (options.basisType == COMPLEX_GTO)
      CErr("Complex GTOs NYI in InCore4indexTPI<dcomplex>",std::cout);

    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();
    
//SS: for debug    
/*
    // Create a vector of libint2::Engines for possible threading
      std::vector<libint2::Engine> engines(1);

    // Initialize the first engine for the integral evaluation
    
    engines[0] = libint2::Engine(libint2::Operator::coulomb,
      basisSet.maxPrim,basisSet.maxL,0);
    engines[0].set_precision(0.);
*/    
// SS: end

    // Copy over the engines to other threads if need be
    // for(size_t i = 1; i < nthreads; i++) engines[i] = engines[0];


    // define magnetic field

    // Allocate and zero out ERIs

    InCore4indexTPI<dcomplex> &eri4I = *this;
    std::fill_n(eri4I.pointer(),NB2*NB2,0.);

    #pragma omp parallel
    {
      int thread_id = GetThreadID();
/*
      // Get threads result buffer
      const auto& buf_vec = engines[thread_id].results();
*/

// SS: for debug
//      const auto& buf_vec = engines[0].results();
// SS: end


      auto magAmp = emPert.getDipoleAmp(Magnetic);
      // std::cout<<"magAmp 2e 0: "<<magAmp[0]<<" 1: "<<magAmp[1]<<" 2: "<<magAmp[2]<<std::endl;   


      size_t n1,n2,n3,n4,i,j,k,l,ijkl,bf1,bf2,bf3,bf4;
      size_t s4_max;
      for(size_t s1(0), bf1_s(0), s1234(0); s1 < basisSet.nShell;
          bf1_s+=n1, s1++) { 

        n1 = basisSet.shells[s1].size(); // Size of Shell 1

      for(size_t s2(0), bf2_s(0); s2 <= s1; bf2_s+=n2, s2++) {

        n2 = basisSet.shells[s2].size(); // Size of Shell 2

        //SS Start generate shellpair1 

        libint2::ShellPair pair1_to_use;
        pair1_to_use.init( basisSet.shells[s1],basisSet.shells[s2],-1000);

        libint2::ShellPair pair1_to_use_switch;
        // switch s1 and s2
        pair1_to_use_switch.init( basisSet.shells[s2],basisSet.shells[s1],-1000);


      for(size_t s3(0), bf3_s(0); s3 <= s1; bf3_s+=n3, s3++) {

        n3 = basisSet.shells[s3].size(); // Size of Shell 3
        s4_max = (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4

      for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {

        n4 = basisSet.shells[s4].size(); // Size of Shell 4

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s1234 % nthreads != thread_id ) continue;
        #endif

        //SS start generate shellpair2 and calculate GIAO ERI

        libint2::ShellPair pair2_to_use;
        
        pair2_to_use.init( basisSet.shells[s3],basisSet.shells[s4],-1000);

#ifdef _DEBUGGIAOERI
 std::cout<<" s1 "<<s1<<" s2 "<<s2<<" s3 "<<s3<<" s4 "<<s4<<std::endl;
#endif

// SS: for debug

// std::cout<<" s1 "<<s1<<" s2 "<<s2<<" s3 "<<s3<<" s4 "<<s4<<std::endl;

// SS: end

// SS bottom up start
#ifdef bottomupGIAO
        auto two2buff = ComplexGIAOIntEngine::bottomupcomplexERI(pair1_to_use,pair2_to_use,
          basisSet.shells[s1],basisSet.shells[s2],
          basisSet.shells[s3],basisSet.shells[s4],&magAmp[0]);
//std::cout<<"calculate bottom up GIAO ERI fuck!!!"<<std::endl;
        auto two2buff_switch = ComplexGIAOIntEngine::bottomupcomplexERI(pair1_to_use_switch,pair2_to_use,
          basisSet.shells[s2],basisSet.shells[s1],
          basisSet.shells[s3],basisSet.shells[s4],&magAmp[0]);
// SS bottom up end
#else
        // calculate integral (s1,s2|s3,s4)
        auto two2buff = ComplexGIAOIntEngine::computeGIAOERIabcd(pair1_to_use,pair2_to_use,
          basisSet.shells[s1],basisSet.shells[s2],
          basisSet.shells[s3],basisSet.shells[s4],&magAmp[0]);


        // calculate integral (s2,s1|s3,s4)
        auto two2buff_switch = ComplexGIAOIntEngine::computeGIAOERIabcd(pair1_to_use_switch,pair2_to_use,
          basisSet.shells[s2],basisSet.shells[s1],
          basisSet.shells[s3],basisSet.shells[s4],&magAmp[0]);

#endif
        
/*
        auto realbuff = RealGTOIntEngine::computeERIabcd(pair1_to_use,pair2_to_use,
          basisSet_.shells[s2],basisSet_.shells[s1],basisSet_.shells[s3],basisSet_.shells[s4]); 
*/

/*
        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s1234 % nthreads != thread_id ) continue;
        #endif
*/
// SS: for debug 
/*
        // Evaluate ERI for shell quartet
        engines[0].compute2<
          libint2::Operator::coulomb, libint2::BraKet::xx_xx, 0>(
          basisSet_.shells[s1],
          basisSet_.shells[s2],
          basisSet_.shells[s3],
          basisSet_.shells[s4]
        );
        // Libint2 internal screening
        const double *buff = buf_vec[0];
        if(buff == nullptr) continue;
*/
// SS:  end

        // Place shell quartet into persistent storage with
        // permutational symmetry
        for(i = 0ul, bf1 = bf1_s, ijkl = 0ul ; i < n1; ++i, bf1++) 
        for(j = 0ul, bf2 = bf2_s             ; j < n2; ++j, bf2++) 
        for(k = 0ul, bf3 = bf3_s             ; k < n3; ++k, bf3++) 
        for(l = 0ul, bf4 = bf4_s             ; l < n4; ++l, bf4++, ++ijkl) {

          int jikl;
          jikl = j*n1*n3*n4 + i*n3*n4 + k*n4 + l;
          
// SS: for debug
/*  
if ( std::abs(two2buff[ijkl]-buff[ijkl]) > 1.0e-10  ) {
//if ( std::abs(realbuff[ijkl]-buff[ijkl]) > 1.0e-7  ) {
  std::cout<<"LA "<<basisSet_.shells[s1].contr[0].l
  <<" LB "<<basisSet_.shells[s2].contr[0].l
  <<" LC "<<basisSet_.shells[s3].contr[0].l 
  <<" LD "<<basisSet_.shells[s4].contr[0].l<<std::endl;
  std::cout<<"  GIAO integral "<<std::setprecision(12)<<two2buff[ijkl];
  std::cout<<"  libint integral  "<<std::setprecision(12)<<buff[ijkl]<<"i "<<i<<" j "<<j<<" k "<<k<<" l "<<l<<std::endl;
}
*/
//SS: end 

// SS start compare the difference between non switch GIAO and switched-GIAO integrals

/*
if ( std::abs(two2buff[ijkl]-two2buff_switch[ijkl]) > 1.0e-11  ) {
  std::cout<<"LA "<<basisSet_.shells[s1].contr[0].l
  <<" LB "<<basisSet_.shells[s2].contr[0].l
  <<" LC "<<basisSet_.shells[s3].contr[0].l 
  <<" LD "<<basisSet_.shells[s4].contr[0].l<<std::endl;
  std::cout<<"  GIAO integral "<<std::setprecision(12)<<two2buff[ijkl];
  std::cout<<"  GIAO switch integral  "<<std::setprecision(12)<<two2buff_switch[jikl]<<std::endl;
}
*/
// SS end

/*
            // (12 | 34)
            ERI[bf1 + bf2*NB + bf3*NB2 + bf4*NB3] = two2nonbuff[ijkl];
            // (12 | 43)
            ERI[bf1 + bf2*NB + bf4*NB2 + bf3*NB3] = two2nonbuff[ijkl];
            // (21 | 34)
            ERI[bf2 + bf1*NB + bf3*NB2 + bf4*NB3] = two2nonbuff[ijkl];
            // (21 | 43)
            ERI[bf2 + bf1*NB + bf4*NB2 + bf3*NB3] = two2nonbuff[ijkl];
            // (34 | 12)
            ERI[bf3 + bf4*NB + bf1*NB2 + bf2*NB3] = two2nonbuff[ijkl];
            // (43 | 12)
            ERI[bf4 + bf3*NB + bf1*NB2 + bf2*NB3] = two2nonbuff[ijkl];
            // (34 | 21)
            ERI[bf3 + bf4*NB + bf2*NB2 + bf1*NB3] = two2nonbuff[ijkl];
            // (43 | 21)
            ERI[bf4 + bf3*NB + bf2*NB2 + bf1*NB3] = two2nonbuff[ijkl];
*/
            // (12 | 34)
            eri4I(bf1, bf2, bf3, bf4) = two2buff[ijkl];
            // (34 | 12)
            eri4I(bf3, bf4, bf1, bf2) = two2buff[ijkl];
            // (21 | 43)
            eri4I(bf2, bf1, bf4, bf3) = std::conj(two2buff[ijkl]);
            // (43 | 21)
            eri4I(bf4, bf3, bf2, bf1) = std::conj(two2buff[ijkl]);


            // (21 | 34)
            eri4I(bf2, bf1, bf3, bf4) = two2buff_switch[jikl];
            // (34 | 21)
            eri4I(bf3, bf4, bf2, bf1) = two2buff_switch[jikl];
            // (12 | 43)
            eri4I(bf1, bf2, bf4, bf3) = std::conj(two2buff_switch[jikl]);
            // (43 | 12)
            eri4I(bf4, bf3, bf1, bf2) = std::conj(two2buff_switch[jikl]);



        }; // ijkl loop
      }; // s4
      }; // s3
      }; // s2
      }; // s1
    }; // omp region

    // Debug output of the ERIs
#if _DEBUGGIAOERI
    std::cout << "Two-Electron GIAO Integrals (GIAO ERIs)" << std::endl;
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++)
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout <<std::setprecision(12)<< ERI[i + j*NB  + k*NB2 + l*NB3] << std::endl;
    };
#endif
  }; // InCore4indexERI<dcomplex>::computeAOInts


}; // namespace ChronusQ

//#endif
