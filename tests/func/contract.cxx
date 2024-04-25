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
#include <func.hpp>

#include <cxxapi/input.hpp>
#include <cxxapi/options.hpp>
#include <cxxapi/boilerplate.hpp>

#include <util/threads.hpp>
#include <util/mpi.hpp>

#include <memmanager.hpp>
#include <cerr.hpp>
#include <molecule.hpp>
#include <basisset.hpp>
#include <integrals.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>

#include <cqlinalg/blasext.hpp>


using namespace ChronusQ;






/**
 *  \brief Generate a randon number of different
 *  types using the STL random number generator
 */
template <typename T>
T RAND_NUMBER(std::default_random_engine &e, 
  std::uniform_real_distribution<> &dis);

template <>
double RAND_NUMBER(std::default_random_engine &e, 
  std::uniform_real_distribution<> &dis) {

  return dis(e);

}; // RAND_NUMBER<double>

template <>
dcomplex RAND_NUMBER(std::default_random_engine &e, 
  std::uniform_real_distribution<> &dis) {

  return dcomplex(dis(e),dis(e));

}; // RAND_NUMBER<dcomplex>


/**
 *  Structure of matricies for contration tests
 */
enum HER_MAT {
  HERMETIAN,
  NONHERMETIAN
};

/**
 *  \brief Template wrapper around functions to optionally
 *  symmetrize a matrix. Null call for H == NONHERMETIAN
 */
template <typename T, HER_MAT H>
void HerMatLoc(char UPLO, size_t N, T *A) {

  if(H == HERMETIAN) ChronusQ::HerMat(UPLO,N,A,N);
  else { ; }

};


template <typename FIELD, HER_MAT HER>
void CONTRACT_TEST(TWOBODY_CONTRACTION_TYPE type, std::string storage) {

#ifdef _CQ_GENERATE_TESTS
  bool exists = false;
#else
  bool exists = true;
#endif


  // Reference File
  SafeFile refFile(FUNC_REFERENCE "contract.hdf5", exists); 


  // Input file for the constructions of Basis sets and Molecules
  CQInputFile input(FUNC_INPUT "contract_ref.inp");
  input.parse();
  
  // Memory
  CQMiscOptions(std::cout,input); 
  // Dummy scrName since READGEOM=INPUTFILE
  std::string scrName;
  
  // Molecule and BasisSet
  std::cout << input << std::endl;
  Molecule mol(std::move(CQMoleculeOptions(std::cout,input,scrName)));
  std::shared_ptr<BasisSet> basis = CQBasisSetOptions(std::cout,input,mol,"BASIS");

  // AOIntegrals object
  Integrals<double> aoints;
  aoints.TPI =
      std::make_shared<DirectTPI<double>>(*basis,*basis,mol,1e-12);
  
  // Scratch memory
  size_t NB = basis->nBasis;
  FIELD *SX  = CQMemManager::get().malloc<FIELD>(NB*NB); 
  FIELD *SX2 = CQMemManager::get().malloc<FIELD>(NB*NB); 
  FIELD *Rand = CQMemManager::get().malloc<FIELD>(NB*NB); 

  // Set up direct contraction
  std::vector<TwoBodyContraction<FIELD>> cont = 
    { { Rand, SX, (HER == HERMETIAN), type  } };
  std::fill_n(SX,NB*NB,0.); // zero out scratch space
  
  EMPerturbation pert; // Dummy perturbation

#ifdef _CQ_GENERATE_TESTS

  // Generate random "X" matrix
  std::random_device r; 
  std::default_random_engine e(r());
  std::uniform_real_distribution<> dis(-50,68); 
  
  for(auto i = 0; i < NB; i++)
  for(auto j = 0; j < NB; j++)
    Rand[j + i*NB] = RAND_NUMBER<FIELD>(e,dis);
  
  // Optionally symmetrize X
  HerMatLoc<FIELD,HER>('U',NB,Rand); 

  // Write X to disk
  refFile.safeWriteData(storage + "/X",Rand,{NB,NB});
  
  // Perform incore ERI contraction and write result to disk
  GTODirectTPIContraction<FIELD,double> TPI(aoints.TPI);
  TPI.twoBodyContract(MPI_COMM_WORLD,true,cont,pert);
  refFile.safeWriteData(storage + "/AX",SX,{NB,NB});

#else

  // Read in X and G[X] from disk
  refFile.readData(storage + "/X",Rand);
  refFile.readData(storage + "/AX",SX2);
  
  // Form G[X] directly
  GTODirectTPIContraction<FIELD,double> TPI(aoints.TPI);
  TPI.twoBodyContract(MPI_COMM_WORLD,true,cont,pert);
  
  // Compare with reference result
  double maxDiff(0.);
  for(auto i = 0; i < NB*NB; i++) 
    maxDiff = std::max(maxDiff,std::abs(SX[i] - SX2[i]));
  
  EXPECT_TRUE(maxDiff < 1e-9) << maxDiff;

#endif

  CQMemManager::get().free(SX,SX2,Rand);


}







// Real Hermetian "J" contraction test
TEST( DIRECT_CONTRACTION_REAL, HER_J_CONTRACT ) {

  CONTRACT_TEST<double,HERMETIAN>(COULOMB,"CONTRACTION/HER/REAL/J");

}

// Real Non-Hermetian "J" contraction test
TEST( DIRECT_CONTRACTION_REAL, NONHER_J_CONTRACT ) {

  CONTRACT_TEST<double,NONHERMETIAN>(COULOMB,"CONTRACTION/NONHER/REAL/J");

}

// Real Hermetian "K" contraction test
TEST( DIRECT_CONTRACTION_REAL, HER_K_CONTRACT ) {

  CONTRACT_TEST<double,HERMETIAN>(EXCHANGE,"CONTRACTION/HER/REAL/K");

}

// Real Non-Hermetian "K" contraction test
TEST( DIRECT_CONTRACTION_REAL, NONHER_K_CONTRACT ) {

  CONTRACT_TEST<double,NONHERMETIAN>(COULOMB,"CONTRACTION/NONHER/REAL/K");

}

// Parallel Real direct contraction tests
#ifdef _CQ_DO_PARTESTS

// Parallel Real Hermetian "J" contraction test
TEST( DIRECT_CONTRACTION_REAL, PAR_HER_J_CONTRACT ) {

  CONTRACT_TEST<double,HERMETIAN>(COULOMB,"CONTRACTION/HER/REAL/J");

}

// Parallel Real Non-Hermetian "J" contraction test
TEST( DIRECT_CONTRACTION_REAL, PAR_NONHER_J_CONTRACT ) {

  CONTRACT_TEST<double,NONHERMETIAN>(COULOMB,"CONTRACTION/NONHER/REAL/J");

}

// Parallel Real Hermetian "K" contraction test
TEST( DIRECT_CONTRACTION_REAL, PAR_HER_K_CONTRACT ) {

  CONTRACT_TEST<double,HERMETIAN>(EXCHANGE,"CONTRACTION/HER/REAL/K");

}

// Parallel Real Non-Hermetian "K" contraction test
TEST( DIRECT_CONTRACTION_REAL, PAR_NONHER_K_CONTRACT ) {

  CONTRACT_TEST<double,NONHERMETIAN>(COULOMB,"CONTRACTION/NONHER/REAL/K");

}


#endif


// Complex Hermetian "J" contraction test
TEST( DIRECT_CONTRACTION_COMPLEX, HER_J_CONTRACT ) {

  CONTRACT_TEST<dcomplex,HERMETIAN>(COULOMB,"CONTRACTION/HER/COMPLEX/J");

}

// Complex Non-Hermetian "J" contraction test
TEST( DIRECT_CONTRACTION_COMPLEX, NONHER_J_CONTRACT ) {

  CONTRACT_TEST<dcomplex,NONHERMETIAN>(COULOMB,"CONTRACTION/NONHER/COMPLEX/J");

}


// Complex Hermetian "K" contraction test
TEST( DIRECT_CONTRACTION_COMPLEX, HER_K_CONTRACT ) {

  CONTRACT_TEST<dcomplex,HERMETIAN>(EXCHANGE,"CONTRACTION/HER/COMPLEX/K");

}

// Complex Non-Hermetian "K" contraction test
TEST( DIRECT_CONTRACTION_COMPLEX, NONHER_K_CONTRACT ) {

  CONTRACT_TEST<dcomplex,HERMETIAN>(EXCHANGE,"CONTRACTION/NONHER/COMPLEX/K");

}



// Parallel Complex direct contraction tests
#ifdef _CQ_DO_PARTESTS 

// Parallel Complex Hermetian "J" contraction test
TEST( DIRECT_CONTRACTION_COMPLEX, PAR_HER_J_CONTRACT ) {

  CONTRACT_TEST<dcomplex,HERMETIAN>(COULOMB,"CONTRACTION/HER/COMPLEX/J");

}

// Parallel Complex Non-Hermetian "J" contraction test
TEST( DIRECT_CONTRACTION_COMPLEX, PAR_NONHER_J_CONTRACT ) {

  CONTRACT_TEST<dcomplex,NONHERMETIAN>(COULOMB,"CONTRACTION/NONHER/COMPLEX/J");

}


// Parallel Complex Hermetian "K" contraction test
TEST( DIRECT_CONTRACTION_COMPLEX, PAR_HER_K_CONTRACT ) {

  CONTRACT_TEST<dcomplex,HERMETIAN>(EXCHANGE,"CONTRACTION/HER/COMPLEX/K");

}

// Parallel Complex Non-Hermetian "K" contraction test
TEST( DIRECT_CONTRACTION_COMPLEX, PAR_NONHER_K_CONTRACT ) {

  CONTRACT_TEST<dcomplex,HERMETIAN>(EXCHANGE,"CONTRACTION/NONHER/COMPLEX/K");

}

#endif

