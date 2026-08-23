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
#include <particleintegrals/gradints.hpp>
#include <particleintegrals/twopints/gtodirectreleri.hpp>
#include <util/timer.hpp>
#include <util/threads.hpp>
#include <chrono>

//#define _DEBUGORTHO
//#define _DEBUGERI
//#define _DEBUGGIAOERI //SS
//#define _DEBUGGIAOERIDERIV //TDD
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

    ProgramTimer::tick("Form ERI");

    bool sameBasis = (&basisSet == &basisSet2);
    if (!sameBasis and op != EP_ATTRACTION)
      CErr("(ee|pp) needs op==EP_ATTRACTION in InCore4indexTPI<dcomplex>)",std::cout);
    if (op != ELECTRON_REPULSION and op != EP_ATTRACTION)
      CErr("Only e-p attraction/e-e/p-p repulsion integrals in InCore4indexTPI<dcomplex>",std::cout);
    if (options.basisType == REAL_GTO)
      CErr("Real GTOs are not allowed in InCore4indexTPI<dcomplex>",std::cout);
    if (options.basisType == COMPLEX_GTO)
      CErr("Complex GTOs NYI in InCore4indexTPI<dcomplex>",std::cout);

    // GIAO London phase sign is per-side particle charge (electron(-1) is default, proton(+1) flips)
    double braCharge = options.particle.charge;
    double ketCharge = options.particle2.charge;


    //if (op == EP_ATTRACTION) std::cout << "ATTENTION: DEBUG START" << std::endl;

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
    size_t NB2 = this->NB2;
    size_t MB2 = this->sNB2;
    InCore4indexTPI<dcomplex> &eri4I = *this;
    std::fill_n(eri4I.pointer(),NB2*MB2,0.);

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
      size_t s3_max, s4_max;
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

      s3_max = sameBasis ? s1 : basisSet2.nShell - 1;

      for(size_t s3(0), bf3_s(0); s3 <= s3_max; bf3_s+=n3, s3++) {

        n3 = basisSet2.shells[s3].size(); // Size of Shell 3
        s4_max = sameBasis && (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4

      for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {

        n4 = basisSet2.shells[s4].size(); // Size of Shell 4

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s1234 % nthreads != thread_id ) continue;
        #endif

        //SS start generate shellpair2 and calculate GIAO ERI

        libint2::ShellPair pair2_to_use;
        
        pair2_to_use.init( basisSet2.shells[s3],basisSet2.shells[s4],-1000);

#ifdef _DEBUGGIAOERI
 std::cout<<"current shell number:  s1 "<<s1<<" s2 "<<s2<<" s3 "<<s3<<" s4 "<<s4<<std::endl;
#endif

// SS: for debug

// std::cout<<" s1 "<<s1<<" s2 "<<s2<<" s3 "<<s3<<" s4 "<<s4<<std::endl;

// SS: end

// SS bottom up start
#ifdef bottomupGIAO
        auto two2buff = ComplexGIAOIntEngine::bottomupcomplexERI(pair1_to_use,pair2_to_use,
          basisSet.shells[s1],basisSet.shells[s2],
          basisSet2.shells[s3],basisSet2.shells[s4],&magAmp[0],braCharge,ketCharge);
        auto two2buff_switch = ComplexGIAOIntEngine::bottomupcomplexERI(pair1_to_use_switch,pair2_to_use,
          basisSet.shells[s2],basisSet.shells[s1],
          basisSet2.shells[s3],basisSet2.shells[s4],&magAmp[0],braCharge,ketCharge);
// SS bottom up end
#else

#ifdef _DEBUGGIAOERI
 std::cout<<"doing shell gradient ("<<s1<<" "<<s2<<"|"<<s3<<" "<<s4<<") "<<std::endl;
#endif
        // calculate integral (s1,s2|s3,s4)
        auto two2buff = ComplexGIAOIntEngine::computeGIAOERIabcd(pair1_to_use,pair2_to_use,
          basisSet.shells[s1],basisSet.shells[s2],
          basisSet2.shells[s3],basisSet2.shells[s4],&magAmp[0],braCharge,ketCharge);

#ifdef _DEBUGGIAOERI
 std::cout<<"doing shell gradient ("<<s2<<" "<<s1<<"|"<<s3<<" "<<s4<<") "<<std::endl;
#endif

        // calculate integral (s2,s1|s3,s4)
        auto two2buff_switch = ComplexGIAOIntEngine::computeGIAOERIabcd(pair1_to_use_switch,pair2_to_use,
          basisSet.shells[s2],basisSet.shells[s1],
          basisSet2.shells[s3],basisSet2.shells[s4],&magAmp[0],braCharge,ketCharge);

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
            
            // 4-fold symmetry for GIAO
            // 4-fold symmetry only if left basis is the same as right basis

            // (12 | 34)
            eri4I(bf1, bf2, bf3, bf4) = two2buff[ijkl];
            // (21 | 43)
            eri4I(bf2, bf1, bf4, bf3) = std::conj(two2buff[ijkl]);

            if( sameBasis ) {
              // (34 | 12)
              eri4I(bf3, bf4, bf1, bf2) = two2buff[ijkl];
              // (43 | 21)
              eri4I(bf4, bf3, bf2, bf1) = std::conj(two2buff[ijkl]);
            }

            // (21 | 34)
            eri4I(bf2, bf1, bf3, bf4) = two2buff_switch[jikl];
            // (12 | 43)
            eri4I(bf1, bf2, bf4, bf3) = std::conj(two2buff_switch[jikl]);

            if( sameBasis ) {
            // (34 | 21)
            eri4I(bf3, bf4, bf2, bf1) = two2buff_switch[jikl];
            // (43 | 12)
            eri4I(bf4, bf3, bf1, bf2) = std::conj(two2buff_switch[jikl]);
            }


        }; // ijkl loop
      }; // s4
      }; // s3
      }; // s2
      }; // s1
    }; // omp region

    // Debug output of the ERIs
#ifdef _DEBUGGIAOERI
    std::cout << "braCharge=" << braCharge << " ketCharge=" << ketCharge << std::endl;
    auto magAmp = emPert.getDipoleAmp(Magnetic);
    std::cout<<"magAmp 2e 0: "<< magAmp[0]<<" 1: "<< magAmp[1]<<" 2: "<< magAmp[2]<<std::endl; 
    std::cout << "Two-Electron GIAO Integrals (GIAO ERIs)" << std::endl;
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++)
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout <<std::setprecision(12)<< eri4I(i, j, k, l) << std::endl;
    };
#endif

    ProgramTimer::tock("Form ERI"); 

  }; // InCore4indexERI<dcomplex>::computeAOInts


  /**
   *  \brief Allocate, compute and store the full rank-4 complex ERI order 1 gradient tensor 
   *  using in house GIAO code over the CGTO basis.
   * 
   *  TangDD: I'm trying to unify the implementation of gradient of GTO and GIAO, but
   *  these two codes differ too much. In the future, a refactor is expected to make all GIAO ERI
   *  go through 'computeAOInts' and be general for any order of gradients.
   * 
   */ 
  template<>
  void GradInts<TwoPInts,dcomplex>::computeAOInts(BasisSet& basisSet,
    BasisSet& basisSet2, Molecule& mol, EMPerturbation& emPert, OPERATOR op,
    const HamiltonianOptions &options)
  {

    ProgramTimer::tick("Form ERI");

    if (std::dynamic_pointer_cast<DirectTPI<dcomplex>>(components_[0]))
      return;

    // Get vector of internal storages
    std::vector<dcomplex*> eris;
    std::transform(components_.begin(), components_.end(),
      std::back_inserter(eris),
      [](std::shared_ptr<TwoPInts<dcomplex>>& p) {
        return std::dynamic_pointer_cast<InCore4indexTPI<dcomplex>>(p)->pointer();
      }
    );

    size_t NB = basisSet.nBasis;
    size_t MB = basisSet2.nBasis;
    size_t NB2 = NB * NB;
    size_t NB3 = NB2 * MB;
    size_t MB2 = MB * MB;
    size_t NB2MB2 = NB2 * MB2;

    // Clear previous ERIs
    std::for_each(eris.begin(), eris.end(),
      [&](dcomplex* p){std::fill_n(p,NB2MB2,dcomplex(0.0));}
    );

    // NEO Options
    bool sameBasis = (&basisSet == &basisSet2);
    if (!sameBasis and op != EP_ATTRACTION)
      CErr("(ee|pp) needs op==EP_ATTRACTION in InCore4indexTPI<dcomplex>)",std::cout);
    if (op != ELECTRON_REPULSION and op != EP_ATTRACTION)
      CErr("Only e-p attraction/e-e/p-p repulsion integrals in InCore4indexTPI<dcomplex>",std::cout);
    if (options.basisType == REAL_GTO)
      CErr("Real GTOs are not allowed in InCore4indexTPI<dcomplex>",std::cout);
    if (options.basisType == COMPLEX_GTO)
      CErr("Complex GTOs NYI in InCore4indexTPI<dcomplex>",std::cout);

    // GIAO London phase sign is per-side particle charge (electron(-1) is default, proton(+1) flips)
    double braCharge = options.particle.charge;
    double ketCharge = options.particle2.charge;

    // Determine the number of OpenMP threads
    int nthreads = GetNumThreads();

int indx_switch[12]; 
indx_switch[0] = 3;
indx_switch[1] = 4;
indx_switch[2] = 5;
indx_switch[3] = 0;
indx_switch[4] = 1;
indx_switch[5] = 2;
indx_switch[6] = 6;
indx_switch[7] = 7;
indx_switch[8] = 8;
indx_switch[9] = 9;
indx_switch[10] = 10;
indx_switch[11] = 11;

    // Parallel Region Start
    // We will first distribute the ERI gradient by shellpairs, this makes no difference from ERI case.
    #pragma omp parallel
    {
      int thread_id = GetThreadID();

      auto magAmp = emPert.getDipoleAmp(Magnetic); 

      size_t n1,n2,n3,n4,i,j,k,l,ijkl,bf1,bf2,bf3,bf4;
      size_t s3_max, s4_max;
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

      s3_max = sameBasis ? s1 : basisSet2.nShell - 1;

      for(size_t s3(0), bf3_s(0); s3 <= s3_max; bf3_s+=n3, s3++) {

        n3 = basisSet2.shells[s3].size(); // Size of Shell 3
        s4_max = sameBasis && (s1 == s3) ? s2 : s3; // Determine the unique max of Shell 4

      for(size_t s4(0), bf4_s(0); s4 <= s4_max; bf4_s+=n4, s4++, s1234++) {

        n4 = basisSet2.shells[s4].size(); // Size of Shell 4

        // Round Robbin work distribution
        #ifdef _OPENMP
        if( s1234 % nthreads != thread_id ) continue;
        #endif

        //SS start generate shellpair2 and calculate GIAO ERI

        libint2::ShellPair pair2_to_use;
        
        pair2_to_use.init( basisSet2.shells[s3],basisSet2.shells[s4],-1000);

#ifdef _DEBUGGIAOERI
std::cout<<" s1 "<<s1<<" s2 "<<s2<<" s3 "<<s3<<" s4 "<<s4<<std::endl;
#endif

        // Find out which atom
        std::vector<size_t> ac;
        ac.push_back(basisSet.mapSh2Cen[s1]);
        ac.push_back(basisSet.mapSh2Cen[s2]);
        ac.push_back(basisSet2.mapSh2Cen[s3]);
        ac.push_back(basisSet2.mapSh2Cen[s4]);

        // No gradient for (AA|AA)
        //if (ac[0] == ac[1])
        //  if (ac[2] == ac[3])
        //    if (ac[0] == ac[2])
        //      break;

// SS bottom up start
// TangDD: Gradient always use Bottomup

#ifdef _DEBUGGIAOERI
 std::cout<<"doing shell gradient ("<<s1<<" "<<s2<<"|"<<s3<<" "<<s4<<") "<<std::endl;
#endif

        auto two2buff = ComplexGIAOIntEngine::bottomupcomplexERI_deriv1(pair1_to_use,pair2_to_use,
          basisSet.shells[s1],basisSet.shells[s2],
          basisSet2.shells[s3],basisSet2.shells[s4],&magAmp[0],braCharge,ketCharge);

#ifdef _DEBUGGIAOERI
 std::cout<<"doing shell gradient ("<<s2<<" "<<s1<<"|"<<s3<<" "<<s4<<") "<<std::endl;
#endif

        auto two2buff_switch = ComplexGIAOIntEngine::bottomupcomplexERI_deriv1(pair1_to_use_switch,pair2_to_use,
          basisSet.shells[s2],basisSet.shells[s1],
          basisSet2.shells[s3],basisSet2.shells[s4],&magAmp[0],braCharge,ketCharge);
// SS bottom up end

        // TangDD Start
        // The result vector contains information [IC][IXYZ][dERI].
        // We need to assign it to the output.

        // Loop over centers
        for ( auto iC = 0, itot = 0; iC < 4; iC++ )
        for ( auto iXYZ = 0; iXYZ < 3; iXYZ++, itot++) {

          // No internal screening yet;
          // if ( results[itot] == nullptr ) continue;

          // Place shell quartet into persistent storage with
          // permutational symmetry
          for(i = 0ul, bf1 = bf1_s, ijkl = 0ul ; i < n1; ++i, bf1++) 
          for(j = 0ul, bf2 = bf2_s             ; j < n2; ++j, bf2++) 
          for(k = 0ul, bf3 = bf3_s             ; k < n3; ++k, bf3++) 
          for(l = 0ul, bf4 = bf4_s             ; l < n4; ++l, bf4++, ++ijkl) {

            #ifdef _DEBUGGIAOERI
            std::cout << "(" << bf1_s+i << "," << bf2_s+j << "|" << bf3_s+k << "," << bf4_s+l << ")  " << std::endl;
            #endif

            int jikl;
            jikl = j*n1*n3*n4 + i*n3*n4 + k*n4 + l;
            
            // 4-fold symmetry for GIAO
            // 4-fold symmetry only if left basis is the same as right basis

            //std::cout << "d(ERI)[" <<bf1<<bf2<<bf3<<bf4<< "]/d["<<ac[iC]<< iXYZ<< " assigned as called " << itot << std::endl; 
            //std::cout << "d(ERI)[" <<bf2<<bf1<<bf3<<bf4<< "]/d["<<ac[iC]<< iXYZ<< " assigned as called conj " << itot << std::endl; 

            // (12 | 34)
            #ifdef _DEBUGGIAOERI
            std::cout << "d(ERI)[" <<bf1<<bf2<<bf3<<bf4<< "]/d["<<ac[iC]<< "] at xyz [" << iXYZ<< "] assigned as called " << itot << " for cart index " << i<<j<<k<<l << std::endl;
            std::cout << eris[3*ac[iC]+iXYZ][bf1 + bf2*NB + bf3*NB2 + bf4*NB3] << " += " << two2buff[itot][ijkl] << std::endl; 
            #endif
            eris[3*ac[iC]+iXYZ][bf1 + bf2*NB + bf3*NB2 + bf4*NB3] += two2buff[itot][ijkl];
  
            if ( s3 != s4 ) {
              // (12 | 43)
              #ifdef _DEBUGGIAOERI
              std::cout << "d(ERI)[" <<bf1<<bf2<<bf4<<bf3<< "]/d["<<ac[iC]<< "] at xyz [" << iXYZ<< "] assigned as called switch " << indx_switch[itot] << " for cart index " << j<<i<<k<<l << std::endl;
              std::cout << eris[3*ac[iC]+iXYZ][bf1 + bf2*NB + bf4*NB2 + bf3*NB3] << " += " << std::conj(two2buff_switch[indx_switch[itot]][jikl]) << std::endl;
              #endif   
              eris[3*ac[iC]+iXYZ][bf1 + bf2*NB + bf4*NB2 + bf3*NB3] += std::conj(two2buff_switch[indx_switch[itot]][jikl]);
            }
  
            if ( s1 != s2 ) {
              // (21 | 34)
              #ifdef _DEBUGGIAOERI
              std::cout << "d(ERI)[" <<bf2<<bf1<<bf3<<bf4<< "]/d["<<ac[iC]<< "] at xyz [" << iXYZ<< "] assigned as called switch " << indx_switch[itot] << " for cart index " << j<<i<<k<<l << std::endl;
              std::cout << eris[3*ac[iC]+iXYZ][bf2 + bf1*NB + bf3*NB2 + bf4*NB3] << " += " << two2buff_switch[indx_switch[itot]][jikl] << std::endl; 
              #endif
              eris[3*ac[iC]+iXYZ][bf2 + bf1*NB + bf3*NB2 + bf4*NB3] += two2buff_switch[indx_switch[itot]][jikl];
              if ( s3 != s4 ) {
                // (21 | 43)
                #ifdef _DEBUGGIAOERI
                std::cout << "d(ERI)[" <<bf2<<bf1<<bf4<<bf3<< "]/d["<<ac[iC]<< "] at xyz [" << iXYZ<< "] assigned as called " << itot << " for cart index " << i<<j<<k<<l << std::endl;
                std::cout << eris[3*ac[iC]+iXYZ][bf2 + bf1*NB + bf4*NB2 + bf3*NB3] << " += " << std::conj(two2buff[itot][ijkl]) << std::endl;  
                #endif
                eris[3*ac[iC]+iXYZ][bf2 + bf1*NB + bf4*NB2 + bf3*NB3] += std::conj(two2buff[itot][ijkl]);
              }
            } // sh1/2  
  
            if( sameBasis ) {
              if ( s1 != s3 || s2 != s4 ) {
                // (34 | 12)
                #ifdef _DEBUGGIAOERI
                std::cout << "d(ERI)[" <<bf3<<bf4<<bf1<<bf2<< "]/d["<<ac[iC]<< "] at xyz [" << iXYZ<< "] assigned as called " << itot << " for cart index " << i<<j<<k<<l << std::endl;
                std::cout << eris[3*ac[iC]+iXYZ][bf3 + bf4*NB + bf1*NB2 + bf2*NB3] << " += " << two2buff[itot][ijkl] << std::endl;  
                #endif
                eris[3*ac[iC]+iXYZ][bf3 + bf4*NB + bf1*NB2 + bf2*NB3] += two2buff[itot][ijkl];
  
                if ( s3 != s4 ) {
                  // (43 | 12)
                  #ifdef _DEBUGGIAOERI
                  std::cout << "d(ERI)[" <<bf4<<bf3<<bf1<<bf2<< "]/d["<<ac[iC]<< "] at xyz [" << iXYZ<< "] assigned as called switch " << indx_switch[itot] << " for cart index " << j<<i<<k<<l << std::endl;
                  std::cout << eris[3*ac[iC]+iXYZ][bf4 + bf3*NB + bf1*NB2 + bf2*NB3] << " += " << std::conj(two2buff_switch[indx_switch[itot]][jikl]) << std::endl; 
                  #endif 
                  eris[3*ac[iC]+iXYZ][bf4 + bf3*NB + bf1*NB2 + bf2*NB3] += std::conj(two2buff_switch[indx_switch[itot]][jikl]);
                }
  
                if ( s1 != s2 ) {
                  // (34 | 21)
                  #ifdef _DEBUGGIAOERI
                  std::cout << "d(ERI)[" <<bf3<<bf4<<bf2<<bf1<< "]/d["<<ac[iC]<< "] at xyz [" << iXYZ<< "] assigned as called switch " << indx_switch[itot] << " for cart index " << j<<i<<k<<l << std::endl;
                  std::cout << eris[3*ac[iC]+iXYZ][bf3 + bf4*NB + bf2*NB2 + bf1*NB3] << " += " << two2buff_switch[indx_switch[itot]][jikl] << std::endl;  
                  #endif
                  eris[3*ac[iC]+iXYZ][bf3 + bf4*NB + bf2*NB2 + bf1*NB3] += two2buff_switch[indx_switch[itot]][jikl];
  
                  if ( s3 != s4 ) {
                    // (43 | 21)
                    #ifdef _DEBUGGIAOERI
                    std::cout << "d(ERI)[" <<bf4<<bf3<<bf2<<bf1<< "]/d["<<ac[iC]<< "] at xyz [" << iXYZ<< "] assigned as called " << itot << " for cart index " << i<<j<<k<<l << std::endl;
                    std::cout << eris[3*ac[iC]+iXYZ][bf4 + bf3*NB + bf2*NB2 + bf1*NB3] << " += " << std::conj(two2buff[itot][ijkl]) << std::endl;  
                    #endif
                    eris[3*ac[iC]+iXYZ][bf4 + bf3*NB + bf2*NB2 + bf1*NB3] += std::conj(two2buff[itot][ijkl]);
                  }
                  
                } // sh1/2
              } // sh1/3 or sh2/4
            } 

          }; // ijkl loop
        }; // IC/IXYZ loop

      }; // s4
      }; // s3
      }; // s2
      }; // s1
    }; // omp region

#ifdef _DEBUGGIAOERIDERIV
auto magAmp = emPert.getDipoleAmp(Magnetic);
std::cout<<"magAmp 2e 0: "<< magAmp[0]<<" 1: "<< magAmp[1]<<" 2: "<< magAmp[2]<<std::endl; 
std::cout << "eepp GIAO Integrals order 1 gradient (GIAO ERIs) 0x" << std::endl;
for(auto ii = 0ul; ii < 6; ii++)
for(auto k = 0ul; k < MB; k++)
for(auto l = 0ul; l < MB; l++)
for(auto i = 0ul; i < NB; i++)
for(auto j = 0ul; j < NB; j++){
  std::cout << "d(" << i << "," << j << "|" << k << "," << l << ")/dX" << ii << " ";
  std::cout <<std::setprecision(12)<< eris[ii][i + j*NB + k*NB2 + l*NB3] << std::endl;
}
#endif
//CErr("Normal Termination of TDDTest.");

// Debug output of the ERIs
/*
#ifdef _DEBUGGIAOERI
    std::cout << "braCharge=" << braCharge << " ketCharge=" << ketCharge << std::endl;
    auto magAmp = emPert.getDipoleAmp(Magnetic);
    std::cout<<"magAmp 2e 0: "<< magAmp[0]<<" 1: "<< magAmp[1]<<" 2: "<< magAmp[2]<<std::endl; 
    std::cout << "Two-Electron GIAO Integrals order 1 gradient (GIAO ERIs) 0x" << std::endl;
    for(auto k = 0ul; k < NB; k++)
    for(auto l = 0ul; l < NB; l++)
    for(auto i = 0ul; i < NB; i++)
    for(auto j = 0ul; j < NB; j++){
      std::cout << "(" << i << "," << j << "|" << k << "," << l << ")  ";
      std::cout <<std::setprecision(12)<< eris[0](i, j, k, l) << std::endl;
    };
#endif
*/
    ProgramTimer::tick("Form ERI");

  } // void GradInts<TwoPInts,dcomplex>::computeAOInts

  template<>
  void GradInts<TwoPInts,dcomplex>::computeAOInts(BasisSet& basisSet,
    Molecule& mol, EMPerturbation& pert, OPERATOR op,
    const HamiltonianOptions &options)
  {
    computeAOInts(basisSet, basisSet, mol, pert, op, options);
  }

}; // namespace ChronusQ

//#endif
