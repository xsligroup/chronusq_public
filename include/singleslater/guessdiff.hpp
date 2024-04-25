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

namespace ChronusQ {

  extern std::unordered_map<int,std::string> refMap;

  /**
   *  \brief Reads in MO from bin file
   *  of different type as calculation
   *  and uses it as initial guess.
   *
   **/
  template <typename MatsT, typename IntsT>
  template <typename ScrMatsT>
  void SingleSlater<MatsT,IntsT>::getScrMO(SafeFile& scrBin) {

    // dimension of mo1 and mo2
    auto NB = this->nC * this->nAlphaOrbital();
    auto NB2 = NB*NB;

    this->mo[0].clear();

    std::string prefix = "/SCF/";

    auto MO1dims = scrBin.getDims( prefix + "MO1" );
    auto MO2dims = scrBin.getDims( prefix + "MO2" );

    int scrRefType, binRefType;
    scrBin.readData("REF/REFTYPE",&scrRefType);
    savFile.readData("REF/REFTYPE",&binRefType);

    std::cout << "    * Converting from " << refMap[scrRefType] << " to "
      << refMap[binRefType] << std::endl;

    // Find errors in MO1
    if( MO1dims.size() == 0 )
      CErr(prefix + "MO1 does not exist in " + scrBin.fName(), std::cout);

    if( MO1dims.size() != 2 )
      CErr(prefix + "MO1 not saved as a rank-2 tensor in " + scrBin.fName(),
          std::cout);

    // MOs on scr bin file
    std::vector<cqmatrix::Matrix<ScrMatsT>> motmp;
    motmp.emplace_back(MO1dims[0]);
    if( scrRefType == RefType::isURef or scrRefType == RefType::isRORef ) motmp.emplace_back(MO2dims[0]);

    // Read in MO1
    std::cout << "    * Found SCF/MO1 !" << std::endl;
    scrBin.readData(prefix + "MO1",motmp[0].pointer());

    // Unrestricted calculations
    if( scrRefType == RefType::isURef or scrRefType == RefType::isRORef ) {

      if( MO2dims.size() == 0 )
        std::cout << "    * WARNING: SCF/MO2 does not exist in "
          << scrBin.fName() << " -- Copying SCF/MO1 -> SCF/MO2 " << std::endl;

      if( MO2dims.size() > 2  )

        CErr("SCF/MO2 not saved as a rank-2 tensor in " + scrBin.fName(),
            std::cout);

      // Read in MO2
      std::cout << "    * Found SCF/MO2 !" << std::endl;
      scrBin.readData(prefix + "MO2",motmp[1].pointer());

    }

    // Handle motmp->mo
    // motmp size (assumes alpha and beta same size)
    auto scrMOSize = motmp[0].dimension();

    // Guess mo same size as calculation mo
    if( scrMOSize == NB ){

      // Same size and same type
      if( scrRefType == binRefType ){

        SetMat('N',NB,NB,MatsT(1.),motmp[0].pointer(),NB,this->mo[0].pointer(),NB);
        if( binRefType == RefType::isRORef or binRefType == RefType::isURef )
          SetMat('N',NB,NB,MatsT(1.),motmp[1].pointer(),NB,this->mo[1].pointer(),NB);

      // ROHF guesses
      }else if( binRefType == RefType::isRORef ){

        // RHF->ROHF
        if( scrRefType == RefType::isRRef ){

          SetMat('N',NB,NB,MatsT(1.),motmp[0].pointer(),NB,this->mo[0].pointer(),NB);
          SetMat('N',NB,NB,MatsT(1.),motmp[0].pointer(),NB,this->mo[1].pointer(),NB);

        } else {
          CErr("Same size guess conversion for ROHF failed.");
        }

      // UHF guesses
      }else if( binRefType == RefType::isURef ){

        // RHF/ROHF->UHF
        if( scrRefType == RefType::isRRef or scrRefType == RefType::isRORef ){

          SetMat('N',NB,NB,MatsT(1.),motmp[0].pointer(),NB,this->mo[0].pointer(),NB);
          SetMat('N',NB,NB,MatsT(1.),motmp[0].pointer(),NB,this->mo[1].pointer(),NB);

        } else {
          CErr("Same size guess conversion for UHF failed.");
        }

      } else {
        CErr("This case for guesses from different calcs of same size NYI.");
      }

    // Guess mo different size as calculation mo
    } else if( scrMOSize < NB ){

      // 2c guesses
      if( binRefType == RefType::isTwoCRef ){

        // RHF/ROHF->2c
        if( scrRefType == RefType::isRRef or scrRefType == RefType::isRORef ){

          convert1CRto2CU(motmp,this->mo);

        // UHF->2c
        } else if( scrRefType == RefType::isURef ){

          convert1CUto2CU(motmp,this->mo);

        } else {
          CErr("Initial Guess MO Conversion for 2c Failed");
        } // end 2c guesses

      // 4c guesses
      } else if( binRefType == RefType::isFourCRef ){

        // 4c MOs stored as alpha-large, alpha-small, beta-large, beta-small
        // Negative MOs come before positive MOs

        // RHF/ROHF->4c
        if( scrRefType == RefType::isRRef or scrRefType == RefType::isRORef ){

          std::cout << "    * WARNING: Small component and negative-energy solutions set to zero" << std::endl;

          convert1CRto4CU(motmp,this->mo);

        }else if( scrRefType == RefType::isURef ){ // UHF->4c

          std::cout << "    * WARNING: Small component and negative-energy solutions set to zero" << std::endl;

          convert1CUto4CU(motmp,this->mo);

        // 2c->4c
        }else if( scrRefType == RefType::isTwoCRef ){

          convert2CUto4CU(motmp,this->mo,scrBin);

        } else {
          CErr("Initial Guess MO Conversion for 4c Failed");
        } // end 4c guesses

      } else { // end check on all types of guesses
        CErr("Initial Guess MO Conversion NYI");
      }

    } else CErr("Cannot use a guess of larger size.");

    std::cout << "\n" << std::endl;

    motmp.clear();

  } // SingleSlater<T>::getScrMO()

  /**
   *  \brief Converts 1-component restricted MOs to 2-component unrestricted
   *
   **/
  template <typename MatsT, typename IntsT>
  template <typename ScrMatsT>
  void SingleSlater<MatsT,IntsT>::convert1CRto2CU(std::vector<cqmatrix::Matrix<ScrMatsT>>& inputMO, std::vector<cqmatrix::Matrix<MatsT>>& outputMO) {

    size_t NB = outputMO[0].dimension();
    size_t scrMOSize = inputMO[0].dimension();

    size_t smallMO=0; // RHF/ROHF index

    for( size_t iMO=0; iMO<NB; iMO++ ){

      smallMO = iMO%2==0 ? iMO/2 : (iMO-1)/2;

      // alpha spinor is even and beta is odd
      if( iMO%2 == 0 )
        SetMat('N',scrMOSize,1,MatsT(1.),inputMO[0].pointer()+scrMOSize*smallMO,scrMOSize,outputMO[0].pointer()+NB*iMO,NB);
      else if( iMO%2 != 0 )
        SetMat('N',scrMOSize,1,MatsT(1.),inputMO[0].pointer()+scrMOSize*smallMO,scrMOSize,outputMO[0].pointer()+NB*iMO+NB/2,NB);

    }

  } // SingleSlater<MatsT,IntsT>::convert1CRto2CU

  /**
   *  \brief Converts 1-component unrestricted MOs to 2-component unrestricted
   *
   **/
  template <typename MatsT, typename IntsT>
  template <typename ScrMatsT>
  void SingleSlater<MatsT,IntsT>::convert1CUto2CU(std::vector<cqmatrix::Matrix<ScrMatsT>>& inputMO, std::vector<cqmatrix::Matrix<MatsT>>& outputMO) {

    size_t NB = outputMO[0].dimension();
    size_t scrMOSize = inputMO[0].dimension();

    size_t smallMO=0; // UHF index

    for( size_t iMO=0; iMO<NB; iMO++ ){

      smallMO = iMO%2==0 ? iMO/2 : (iMO-1)/2;

      // alpha spinor is even and beta is odd
      if( iMO%2 == 0 )
        SetMat('N',scrMOSize,1,MatsT(1.),inputMO[0].pointer()+scrMOSize*smallMO,scrMOSize,outputMO[0].pointer()+NB*iMO,NB);
      else if( iMO%2 != 0 )
        SetMat('N',scrMOSize,1,MatsT(1.),inputMO[1].pointer()+scrMOSize*smallMO,scrMOSize,outputMO[0].pointer()+NB*iMO+NB/this->nC,NB);

    }

  } // SingleSlater<MatsT,IntsT>::convert1CUto2CU

  /**
   *  \brief Converts 1-component restricted MOs to 4-component unrestricted
   *
   **/
  template <typename MatsT, typename IntsT>
  template <typename ScrMatsT>
  void SingleSlater<MatsT,IntsT>::convert1CRto4CU(std::vector<cqmatrix::Matrix<ScrMatsT>>& inputMO, std::vector<cqmatrix::Matrix<MatsT>>& outputMO) {

    size_t NB = outputMO[0].dimension();
    size_t scrMOSize = inputMO[0].dimension();

    size_t smallMO=0; // RHF/ROHF index

    // start from negative energy spinors for now
    for( size_t iMO=0; iMO<NB; iMO++ ){

      // negative spinors
      if( iMO < NB/2 ){
        //skip for now
      } else {  // positive spinors

        smallMO = iMO%2==0 ? (iMO-NB/2)/2 : (iMO-NB/2-1)/2;

        // alpha spinor is even and beta is odd
        if( iMO%2 == 0 )
          SetMat('N',scrMOSize,1,MatsT(1.),inputMO[0].pointer()+scrMOSize*smallMO,scrMOSize,outputMO[0].pointer()+NB*iMO,NB);
        else if( iMO%2 != 0 )
          SetMat('N',scrMOSize,1,MatsT(1.),inputMO[0].pointer()+scrMOSize*smallMO,scrMOSize,outputMO[0].pointer()+NB*iMO+NB/2,NB);

      }

    }

  } // SingleSlater<MatsT,IntsT>::convert1CRto4CU

  /**
   *  \brief Converts 1-component unrestricted MOs to 4-component unrestricted
   *
   **/
  template <typename MatsT, typename IntsT>
  template <typename ScrMatsT>
  void SingleSlater<MatsT,IntsT>::convert1CUto4CU(std::vector<cqmatrix::Matrix<ScrMatsT>>& inputMO, std::vector<cqmatrix::Matrix<MatsT>>& outputMO) {

    size_t NB = outputMO[0].dimension();
    size_t scrMOSize = inputMO[0].dimension();

    size_t smallMO=0; // UHF index

    // start from negative energy spinors for now
    for( size_t iMO=0; iMO<NB; iMO++ ){

      // negative spinors
      if( iMO < NB/2 ){
        //skip for now
      } else {  // positive spinors

        smallMO = iMO%2==0 ? (iMO-NB/2)/2 : (iMO-NB/2-1)/2;

        // alpha spinor is even and beta is odd
        if( iMO%2 == 0 )
          SetMat('N',scrMOSize,1,MatsT(1.),inputMO[0].pointer()+scrMOSize*smallMO,scrMOSize,outputMO[0].pointer()+NB*iMO,NB);
        else if( iMO%2 != 0 )
          SetMat('N',scrMOSize,1,MatsT(1.),inputMO[1].pointer()+scrMOSize*smallMO,scrMOSize,outputMO[0].pointer()+NB*iMO+NB/2,NB);

      }

    }

  } // SingleSlater<MatsT,IntsT>::convert1CUto4CU

  /**
   *  \brief Converts 2-component unrestricted MOs to 4-component unrestricted
   *
   **/
  template <typename MatsT, typename IntsT>
  template <typename ScrMatsT>
  void SingleSlater<MatsT,IntsT>::convert2CUto4CU(std::vector<cqmatrix::Matrix<ScrMatsT>>& inputMO, std::vector<cqmatrix::Matrix<MatsT>>& outputMO, SafeFile& scrBin) {


    std::cout << "    * Looking for U matrices on scratch file: " << scrBin.fName() << std::endl;

    bool doSmallCGuess = true;

    size_t NB = this->basisSet().nBasis;
    size_t NP = this->basisSet().nPrimitive;
    size_t NBC = outputMO[0].dimension();
    size_t scrMOSize = inputMO[0].dimension();
    size_t NBC2 = NBC*NBC;

    // dimensions: row 2*NP, column 2*NB
    size_t Urow = 2*NP;
    size_t Ucol = 2*NB;

    if( Urow != Ucol ) CErr("Only implemented for uncontracted basis set");

    MatsT *readUL = CQMemManager::get().malloc<MatsT>(Urow*Ucol);
    MatsT *readUS = CQMemManager::get().malloc<MatsT>(Urow*Ucol);

    // Read in U matrices
    std::string prefix = "X2C/";

    try{
      scrBin.readData(prefix + "UL", readUL);
    } catch (...) {
      std::cout << "    * Cannot find " + prefix + "UL on scratch file!" << std::endl;
      doSmallCGuess = false;
    }

    try{
      scrBin.readData(prefix + "US", readUS);
    } catch (...) {
      std::cout << "    * Cannot find " + prefix + "US on scratch file!" << std::endl;
      doSmallCGuess = false;
    }

    if( doSmallCGuess ){

      // Note on picture change. We use U for back transformation instead of U^{dagger} because
      // CQ has convention of U^{dagger}HU instead of UHU^{dagger}
//    prettyPrintSmart(std::cout,"UL after read",readUL,Urow,Ucol,Urow);
//    prettyPrintSmart(std::cout,"US after read",readUS,Urow,Ucol,Urow);

      if( scrMOSize != Urow or scrMOSize != Ucol ) CErr("2c MO and U matrix need to have same dimensions!");

      // Make a temporary copy for 2c transformed MO
      cqmatrix::Matrix<MatsT> tmpMO(scrMOSize);

//    prettyPrintSmart(std::cout, "phi^2c in 2CUto4CU", tmpMO.pointer(),scrMOSize,scrMOSize,scrMOSize);

      // Large component: UL phi
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,scrMOSize,scrMOSize,scrMOSize,MatsT(1.),readUL,scrMOSize,
       inputMO[0].pointer(),scrMOSize,MatsT(0.),tmpMO.pointer(),scrMOSize);

//    prettyPrintSmart(std::cout, "UL phi^2c in 2CUto4CU", tmpMO.pointer(),scrMOSize,scrMOSize,scrMOSize);

      // initialize plus large alpha
      SetMat('N',scrMOSize/2,scrMOSize,MatsT(1.),tmpMO.pointer(),scrMOSize,outputMO[0].pointer()+NBC2/2,NBC);
      // initialize plus large beta
      SetMat('N',scrMOSize/2,scrMOSize,MatsT(1.),tmpMO.pointer()+scrMOSize/2,scrMOSize,outputMO[0].pointer()+NBC2/2+NBC/2,NBC);

//    prettyPrintSmart(std::cout, "4c initial guess with just large component in 2CUto4CU", outputMO[0].pointer(),NBC,NBC,NBC);

      // Small component: US phi
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,scrMOSize,scrMOSize,scrMOSize,MatsT(1.),readUS,scrMOSize,
       inputMO[0].pointer(),scrMOSize,MatsT(0.),tmpMO.pointer(),scrMOSize);

//    prettyPrintSmart(std::cout, "US phi^2c in 2CUto4CU", tmpMO.pointer(),scrMOSize,scrMOSize,scrMOSize);

      // initialize plus small alpha
      SetMat('N',scrMOSize/2,scrMOSize,MatsT(1.),tmpMO.pointer(),scrMOSize,outputMO[0].pointer()+NBC2/2+NBC/4,NBC);
      // initialize plus small beta
      SetMat('N',scrMOSize/2,scrMOSize,MatsT(1.),tmpMO.pointer()+scrMOSize/2,scrMOSize,outputMO[0].pointer()+NBC2/2+3*NBC/4,NBC);

//    prettyPrintSmart(std::cout, "4c initial guess with small component in 2CUto4CU", outputMO[0].pointer(),NBC,NBC,NBC);

      tmpMO.clear();

    } else {

      std::cout << "    * WARNING: Small component guess set to zero" << std::endl;

      // 2c for plus large alpha
      SetMat('N',scrMOSize/2,scrMOSize,MatsT(1.),inputMO[0].pointer(),scrMOSize,outputMO[0].pointer()+NBC2/2,NBC);
      // 2c for plus large beta
      SetMat('N',scrMOSize/2,scrMOSize,MatsT(1.),inputMO[0].pointer()+scrMOSize/2,scrMOSize,outputMO[0].pointer()+NBC2/2+NBC/2,NBC);

    }

    CQMemManager::get().free(readUL,readUS);

  } // SingleSlater<MatsT,IntsT>::convert2CUto4CU

  template <>
  template <>
  void SingleSlater<double,double>::getScrMO<dcomplex>(SafeFile& scrBin) {

    CErr("Cannot do complex guess MOs for real calculation.");

  }

  template <>
  template <>
  void SingleSlater<double,dcomplex>::getScrMO<dcomplex>(SafeFile& scrBin) {

    CErr("Cannot do complex guess MOs for real calculation.");

  }

}
