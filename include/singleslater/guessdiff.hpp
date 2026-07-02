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

  // Does double->double, double->dcomplex, dcomplex->double, and dcomplex->dcomplex
  template<typename Out, typename In>
  Out rough_cast(const In& x){
    if constexpr (std::is_convertible_v<In, Out>) {
        return static_cast<Out>(x);     // valid only when convertible
    } else if constexpr (
        std::is_same_v<Out,double> &&
        std::is_same_v<In,std::complex<double>>
    ) {
        return x.real();                // complex->real
    } else if constexpr (
        std::is_same_v<Out,std::complex<double>> &&
        std::is_arithmetic_v<In>
    ) {
        return Out(x,0.0);              // real->complex
    } else {
        static_assert(sizeof(Out)==0, "Unsupported type conversion in rough_cast");
    }
  }

  /*
   *  \brief Driver for basis set projection.
   *  Returns both the projected matrix and projection matrix (for re-use)
   **/
  template <typename MatsT, typename IntsT>
  std::tuple<std::shared_ptr<cqmatrix::Matrix<MatsT>>, std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>>> SingleSlater<MatsT,IntsT>::projectMatrix( 
    const Molecule& mol, const BasisSet& fromBasis, const BasisSet& toBasis, const std::shared_ptr<cqmatrix::Matrix<MatsT>>& fromMatrix, bool doFull ){

    const size_t tNB = toBasis.nBasis;
    const size_t fNB  = fromBasis.nBasis;

    // Generate overlap matrix between fromBasis and and toBasis
    // IntsT = dcomplex NYI
    auto overlapSTF = cqmatrix::Matrix<double>(tNB, fNB);
    std::vector<double*> MatVecS21;
    MatVecS21.emplace_back(overlapSTF.pointer());
    OnePInts<double>::OnePDriverLibint( libint2::Operator::overlap, mol, toBasis, fromBasis, MatVecS21, Particle {-1.,1.}, 0);

    // Generate overlap between toBasis and itself
    auto overlapSTT = cqmatrix::Matrix<double>(tNB, tNB);
    std::vector<double*> MatVecS22;
    MatVecS22.emplace_back(overlapSTT.pointer());
    OnePInts<double>::OnePDriverLibint( libint2::Operator::overlap, mol, toBasis, MatVecS22, Particle {-1.,1.}, 0);

    // Get projection matrix from overlaps
    std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>> projectionMat;
    projectionMat.emplace_back(getProjectionMatrix( overlapSTF, overlapSTT ));

    if( this-> nC == 4 ){

      // Small component overlaps do not need 1/2c^2 since cancels out upon projectionMatrix construction
      // Mixed-basis overlap
      auto overlapSTFss = cqmatrix::Matrix<double>(tNB, fNB);
      std::vector<double*> MatVecS21ss;
      MatVecS21ss.emplace_back(overlapSTFss.pointer());
      OnePInts<double>::OnePDriverLibint( libint2::Operator::kinetic, mol, toBasis, fromBasis, MatVecS21ss, Particle {-1.,1.}, 0);

      // Generate overlap between toBasis and itself
      auto overlapSTTss = cqmatrix::Matrix<double>(tNB, tNB);
      std::vector<double*> MatVecS22ss;
      MatVecS22ss.emplace_back(overlapSTTss.pointer());
      OnePInts<double>::OnePDriverLibint( libint2::Operator::kinetic, mol, toBasis, MatVecS22ss, Particle {-1.,1.}, 0);

      projectionMat.emplace_back(getProjectionMatrix( overlapSTFss, overlapSTTss ));

    }

    auto projectedMat = doFull ? projectFullMatrix(projectionMat, fromMatrix)
                               : projectHalfMatrix(projectionMat, fromMatrix);

    return std::make_tuple(
        projectedMat, projectionMat);
  } // Matrix<T> projectMatrix()

  /*
   * \brief Returns the projection matrix that maps matrices in basis 1 to basis 2
   *        (matrix in basis 2) = (Proj).(matrix in basis 1).(Proj)^T
   **/
  template <typename MatsT, typename IntsT>
  std::shared_ptr<cqmatrix::Matrix<MatsT>> SingleSlater<MatsT,IntsT>::getProjectionMatrix( const cqmatrix::Matrix<MatsT>& overlap21, const cqmatrix::Matrix<MatsT>& overlap22) {
    const size_t NB_1 = overlap21.nColumns();
    const size_t NB_2 = overlap22.nColumns();

        if( overlap21.nRows() != NB_2 )
          CErr("Bad dimensions in getOrthoProjection");

        // Allocate temporaries for eigendecomposition
        auto Vmat     = overlap22;                         //<<< (will be overwritten) Orthogonal V matrix (transposed)
        auto Diag     = std::vector<double>(NB_2);         //<<< Holds eigenvalues
        auto DiagMat  = cqmatrix::Matrix<MatsT>(NB_2, NB_2);  //<<< Holds eigenvalues on its diagonal
        auto tmp      = cqmatrix::Matrix<MatsT>(NB_2, NB_2);
        auto overlap22_inv  = cqmatrix::Matrix<MatsT>(NB_2, NB_2);
        std::fill_n(DiagMat.pointer(), NB_2*NB_2, 0.0);    // Only DiagMat needs to be zeroed out since the others are overwritten entirely

        // Compute inverse of S_22 using eigendecomposition
        HermitianEigen('V', 'L', NB_2, Vmat.pointer(), NB_2, Diag.data());
        for(size_t it(0); it<NB_2; ++it) {
          DiagMat(it,it) = 1.0/Diag[it];
        }

        // DiagMat^inv * V^T --> tmp
        // XXX: We could eliminate this blas call because DiagMat is... well... diagonal
        // (aka just rescale the matrix and multiply each col by diag element)
        blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::Trans, NB_2, NB_2, NB_2, MatsT(1.), DiagMat.pointer(), NB_2, Vmat.pointer(), NB_2, MatsT(0.), tmp.pointer(), NB_2);

        // V * tmp --> overlap22_inv
        blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
            NB_2, NB_2, NB_2, MatsT(1.), Vmat.pointer(), NB_2, tmp.pointer(), NB_2, MatsT(0.), overlap22_inv.pointer(), NB_2);

        // Form left-projection matrix: S22_inv * S21
        auto LeftProj = cqmatrix::Matrix<MatsT>(NB_2, NB_1); //<<< Final projection matrix

        blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans,
            NB_2, NB_1, NB_2, MatsT(1.), overlap22_inv.pointer(), NB_2, overlap21.pointer(), NB_2, MatsT(0.), LeftProj.pointer(), NB_2);

        auto LeftProjptr = std::make_shared<cqmatrix::Matrix<MatsT>>(LeftProj);

        return LeftProjptr;

  } // Matrix<MatsT> getProjectionMatrix

  /*
   * \brief Driver for basis set projection with a given projection matrix
   **/
  template <typename MatsT, typename IntsT>
  std::shared_ptr<cqmatrix::Matrix<MatsT>> SingleSlater<MatsT,IntsT>::projectFullMatrix( 
    const std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>>& projMat, const std::shared_ptr<cqmatrix::Matrix<MatsT>>&fromMatrix ) {

    // Allocate return and temporary matrices
    const size_t fourCFact = (this->nC == 4 ? 2 : 1);
    const size_t tNB = projMat[0]->nRows();
    const size_t fNB = projMat[0]->nColumns();
    const size_t tNBF = projMat[0]->nRows()*fourCFact;
    const size_t fNBF = projMat[0]->nColumns()*fourCFact;
    auto toMatrix = std::make_shared<cqmatrix::Matrix<MatsT>>(tNBF, tNBF);
    auto Intermediate = std::make_shared<cqmatrix::Matrix<MatsT>>(fNBF, tNBF);
    std::fill(toMatrix->pointer(),toMatrix->pointer()+tNBF*tNBF,MatsT(0.0));

    // **********************************************************************************
    // Form density in "to" basis: projMat * fromMatrix * projMat^T --> toMatrix
    // **********************************************************************************

    if( this->nC == 1 or this->nC == 2 ){

      // Do fromMatrix  * ( projMat  )^T --> Intermediate
      //     (fNB, fNB) * (tNB, fNB)^T     --> (fNB, tNB)
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, fNB, tNB, fNB, MatsT(1.), fromMatrix->pointer(), fNB, projMat[0]->pointer(), tNB, MatsT(0.), Intermediate->pointer(), fNB);

      // Do     projMat * Intermediate  --> toMatrix
      //  (tNB, fNB) * (fNB, tNB)    --> (tNB, tNB)
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, tNB, fNB, MatsT(1.), projMat[0]->pointer(), tNB, Intermediate->pointer(), fNB, MatsT(0.), toMatrix->pointer(), tNB);

    }else if( this->nC == 4 ){

      //prettyPrintSmart(std::cout,"Starting 1-PDM",fromMatrix->pointer(),fNBF,fNBF,fNBF);
      //for( size_t i=0; i<fNBF*fNBF; i++) std::cout << "starting 1-PDM[" << i << "]= " << fromMatrix->pointer()[i] << std::endl;

      // LL block
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, fNB, tNB, fNB, MatsT(1.), fromMatrix->pointer(), fNBF, projMat[0]->pointer(), tNB, MatsT(0.), Intermediate->pointer(), fNBF);
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, tNB, fNB, MatsT(1.), projMat[0]->pointer(), tNB, Intermediate->pointer(), fNBF, MatsT(0.), toMatrix->pointer(), tNBF);
      //prettyPrintSmart(std::cout,"1-PDM after LL proj",toMatrix->pointer(),tNBF,tNBF,tNBF);

      // LS block
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, fNB, tNB, fNB, MatsT(1.), fromMatrix->pointer()+fNB, fNBF, projMat[1]->pointer(), tNB, MatsT(0.), Intermediate->pointer()+fNB, fNBF);
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, tNB, fNB, MatsT(1.), projMat[0]->pointer(), tNB, Intermediate->pointer()+fNB, fNBF, MatsT(1.), toMatrix->pointer()+tNB, tNBF);
      //prettyPrintSmart(std::cout,"1-PDM after LS proj",toMatrix->pointer(),tNBF,tNBF,tNBF);

      // SL block
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, fNB, tNB, fNB, MatsT(1.), fromMatrix->pointer()+2*fNB*fNB, fNBF, projMat[0]->pointer(), tNB, MatsT(0.), Intermediate->pointer()+2*fNB*tNB, fNBF);
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, tNB, fNB, MatsT(1.), projMat[1]->pointer(), tNB, Intermediate->pointer()+2*fNB*tNB, fNBF, MatsT(1.), toMatrix->pointer()+2*tNB*tNB, tNBF);
      //prettyPrintSmart(std::cout,"1-PDM after SL proj",toMatrix->pointer(),tNBF,tNBF,tNBF);

      // SS block
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, fNB, tNB, fNB, MatsT(1.), fromMatrix->pointer()+2*fNB*fNB+fNB, fNBF, projMat[1]->pointer(), tNB, MatsT(0.), Intermediate->pointer()+2*fNB*tNB+fNB, fNBF);
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, tNB, fNB, MatsT(1.), projMat[1]->pointer(), tNB, Intermediate->pointer()+2*fNB*tNB+fNB, fNBF, MatsT(1.), toMatrix->pointer()+2*tNB*tNB+tNB, tNBF);
      //prettyPrintSmart(std::cout,"1-PDM after SS proj",toMatrix->pointer(),tNBF,tNBF,tNBF);

    }else{ CErr("Do not know how to do PDM projection.");}


    return toMatrix;

  }

  /*
   * \brief Driver for basis set projection of AOxMO matrix with a given projection matrix
   **/
  template <typename MatsT, typename IntsT>
  std::shared_ptr<cqmatrix::Matrix<MatsT>> SingleSlater<MatsT,IntsT>::projectHalfMatrix( 
    const std::vector<std::shared_ptr<cqmatrix::Matrix<MatsT>>>& projMat, const std::shared_ptr<cqmatrix::Matrix<MatsT>>&fromMatrix ) {

    // Allocate return and temporary matrices
    const size_t tNB = projMat[0]->nRows();
    const size_t fNB = projMat[0]->nColumns();
    const size_t fNBC = fromMatrix->nColumns(); // only creats fNBC MOs in toMatrix
    const size_t tNBC = tNB * this->nC;
    auto toMatrix = std::make_shared<cqmatrix::Matrix<MatsT>>(tNBC, tNBC);
    std::fill(toMatrix->pointer(),toMatrix->pointer()+tNBC*tNBC,MatsT(0.0));

    // **********************************************************************************
    // Form MO in "to" basis: projMat * fromMatrix --> toMatrix
    // **********************************************************************************

    //  (tNB, fNB) * (fNB, fNB*nC)    --> (tNB, fNB*nC)

    if( this->nC == 1 ){

      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, tNB, fNB, MatsT(1.), projMat[0]->pointer(), tNB, fromMatrix->pointer(), fNB, MatsT(0.), toMatrix->pointer(), tNB);

    }else if( this->nC == 2 ){

      // project alpha basis functions
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, fNBC, fNB, MatsT(1.), projMat[0]->pointer(), tNB, fromMatrix->pointer(), fNBC, MatsT(0.), toMatrix->pointer(), tNBC);

      // project beta basis functions
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, fNBC, fNB, MatsT(1.), projMat[0]->pointer(), tNB, fromMatrix->pointer()+fNB, fNBC, MatsT(0.), toMatrix->pointer()+tNB, tNBC);

    }else if( this->nC == 4 ){

      // Performing separate projections for negative and positive-energy spinors
      // If not, would populate all negative-energy spinors first.
      const size_t fNBCh = fNBC/2;
      const size_t tNBCh = tNBC/2;
      const size_t fromOffSet = fNBCh*fNBC;
      const size_t toOffSet = tNBCh*tNBC;

      //for( size_t i=0; i<fNBC*fNBC; i++) std::cout << "starting mo[" << i << "]= " << fromMatrix->pointer()[i] << std::endl;

      // Negative-energy spinors
      // project alpha large basis functions
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, fNBCh, fNB, MatsT(1.), projMat[0]->pointer(), tNB, fromMatrix->pointer(), fNBC, MatsT(0.), toMatrix->pointer(), tNBC);

      // project alpha small basis functions
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, fNBCh, fNB, MatsT(1.), projMat[1]->pointer(), tNB, fromMatrix->pointer()+fNB, fNBC, MatsT(0.), toMatrix->pointer()+tNB, tNBC);

      // project beta large basis functions
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, fNBCh, fNB, MatsT(1.), projMat[0]->pointer(), tNB, fromMatrix->pointer()+2*fNB, fNBC, MatsT(0.), toMatrix->pointer()+2*tNB, tNBC);

      // project beta small basis functions
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, fNBCh, fNB, MatsT(1.), projMat[1]->pointer(), tNB, fromMatrix->pointer()+3*fNB, fNBC, MatsT(0.), toMatrix->pointer()+3*tNB, tNBC);

      // Positive-energy spinors
      // project alpha large basis functions
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, fNBCh, fNB, MatsT(1.), projMat[0]->pointer(), tNB, fromMatrix->pointer()+fromOffSet, fNBC, MatsT(0.), toMatrix->pointer()+toOffSet, tNBC);

      // project alpha small basis functions
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, fNBCh, fNB, MatsT(1.), projMat[1]->pointer(), tNB, fromMatrix->pointer()+fromOffSet+fNB, fNBC, MatsT(0.), toMatrix->pointer()+toOffSet+tNB, tNBC);

      // project beta large basis functions
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, fNBCh, fNB, MatsT(1.), projMat[0]->pointer(), tNB, fromMatrix->pointer()+fromOffSet+2*fNB, fNBC, MatsT(0.), toMatrix->pointer()+toOffSet+2*tNB, tNBC);

      // project beta small basis functions
      blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::NoTrans, tNB, fNBCh, fNB, MatsT(1.), projMat[1]->pointer(), tNB, fromMatrix->pointer()+fromOffSet+3*fNB, fNBC, MatsT(0.), toMatrix->pointer()+toOffSet+3*tNB, tNBC);

      //for( size_t i=0; i<tNBC*tNBC; i++) std::cout << "projected mo[" << i << "]= " << toMatrix->pointer()[i] << std::endl;


    }else{ CErr("Do not know how to do MO projection.");}

    return toMatrix;

  }

   /**
   *  \brief Reads in 1PDM from bin file
   *  of different type as calculation
   *  and uses it as initial guess.
   *
   **/
  template <typename MatsT, typename IntsT>
  template <typename ScrMatsT>
  void SingleSlater<MatsT,IntsT>::getScr1PDM(SafeFile& scrBin) {
    getScr1PDM<ScrMatsT>(scrBin, nullptr);
  }


  template <typename MatsT, typename IntsT>
  template <typename ScrMatsT>
  void SingleSlater<MatsT,IntsT>::getScr1PDM(SafeFile& scrBin, const std::shared_ptr<BasisSet> guessBasisSet ) {

    if( MPIRank(comm) == 0 ) {

      // dimension of 1PDM
      auto NB = basisSet().nBasis;
      if( this->nC == 4 ) NB=2*NB;
      auto NB2 = NB*NB;

      std::string prefix = "/SCF/";
      if (this->particle.charge == 1.0)
          prefix = "/PROT_SCF/";

      auto DSdims = scrBin.getDims( prefix + "1PDM_SCALAR" );
      auto DZdims = scrBin.getDims( prefix + "1PDM_MZ" );
      auto DYdims = scrBin.getDims( prefix + "1PDM_MY" );
      auto DXdims = scrBin.getDims( prefix + "1PDM_MX" );

      bool hasDS = DSdims.size() != 0;
      bool hasDZ = DZdims.size() != 0;
      bool hasDY = DYdims.size() != 0;
      bool hasDX = DXdims.size() != 0;

      bool r2DS = DSdims.size() == 2;
      bool r2DZ = DZdims.size() == 2;
      bool r2DY = DYdims.size() == 2;
      bool r2DX = DXdims.size() == 2;

      int scrRefType, binRefType;
      scrBin.readData("REF/REFTYPE",&scrRefType);
      savFile.readData("REF/REFTYPE",&binRefType);

      std::cout << "    * Converting from " << refMap[scrRefType] << " to "
        << refMap[binRefType] << std::endl;


      // onePDM on scr bin file
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<ScrMatsT>> onePDMtmp;
      onePDMtmp = std::make_shared<cqmatrix::PauliSpinorMatrices<ScrMatsT>>(DSdims[0],hasDY,hasDZ);


      // Errors in 1PDM SCALAR
      if( not hasDS )
        CErr(prefix + "1PDM_SCALAR does not exist in " + scrBin.fName(), std::cout);

      else if( not r2DS )
        CErr(prefix + "1PDM_SCALAR not saved as a rank-2 tensor in " +
            scrBin.fName(), std::cout);

      // Error out if any dimensions don't line up
      //size_t NBCheck = guessBasisSet ? guessBasisSet->nBasis : NB;
      //std::cout << "DSdims[0] = " << DSdims[0] << std::endl;
      //std::cout << "nC = " << this->nC << std::endl;
      //std::cout << "NBCheck = " << NBCheck << std::endl;
      //if( NBCheck != DSdims[0] or (hasDZ and (NBCheck != DZdims[0])) or (hasDY and (NBCheck != DYdims[0])) or (hasDX and (NBCheck != DXdims[0])) )
      //  CErr("Scratch file 1PDM dimensions do not match basis dimensions!");

      // Let the user know guessbasis is being used
      if( guessBasisSet ){
        std::cout << "    * GUESSBASIS section specified, projecting basis set " <<
         guessBasisSet->basisName << " -> " << this->basisSet().basisName << std::endl;
        std::cout << "      WARNING: If not projecting for intitial guess," << std::endl;
        std::cout << "               check agreement between READMO and READDEN." << std::endl;
      }

      // Read in 1PDM SCALAR
      std::cout << "    * Looking for " << prefix << "1PDM_SCALAR... ";
      scrBin.readData(prefix + "1PDM_SCALAR", onePDMtmp->S().pointer());
      std::cout << "Found." << std::endl;

      // MZ
      if( onePDMtmp->hasZ() ){

        std::cout << "    * Looking for " << prefix << "1PDM_MZ... " << std::endl;
        if( not r2DZ )
          CErr(prefix + "1PDM_MZ not saved as a rank-2 tensor in " +
            scrBin.fName(), std::cout);
        scrBin.readData(prefix + "1PDM_MZ", onePDMtmp->Z().pointer());
        std::cout << "Found." << std::endl;
      }

      // MY
      if( onePDMtmp->hasXY() ){

        std::cout << "    * Looking for " << prefix << "1PDM_MX... " << std::endl;
        if( not r2DX )
          CErr(prefix + "1PDM_MX not saved as a rank-2 tensor in " +
            scrBin.fName(), std::cout);
        scrBin.readData(prefix + "1PDM_MX",onePDMtmp->X().pointer());
        std::cout << "Found." << std::endl;

        std::cout << "    * Looking for " << prefix << "1PDM_MY... " << std::endl;
        if( not r2DY )
          CErr(prefix + "1PDM_MY not saved as a rank-2 tensor in " +
            scrBin.fName(), std::cout);
        scrBin.readData(prefix + "1PDM_MY",onePDMtmp->Y().pointer());
        std::cout << "Found." << std::endl;

      }

      // Initialize onePDM
      auto scr1PDMSize = onePDMtmp->nRows();
      if( not guessBasisSet ) {
        // Guess 1PDM same size as calculation 1PDM
        if( scr1PDMSize == NB ) *this->onePDM = *onePDMtmp;
        // Guess 1PDM smaller than 1PDM
        else if( scr1PDMSize < NB ){
            auto p1Comps = this->onePDM->SZYXPointers();
            auto p2Comps = onePDMtmp->SZYXPointers();
            auto nComp = p1Comps.size();
            auto n2Comp = p2Comps.size();
            for( auto iComp=0; iComp<nComp; iComp++ ){
              if( iComp < n2Comp )
                SetMat('N',scr1PDMSize,scr1PDMSize,MatsT(1.),
                   p2Comps[iComp],scr1PDMSize,p1Comps[iComp],NB);
            }
          } else CErr("Cannot use a guess of larger size. Specify GUESSBASIS section if guess is in a different basis.");
      } else {

        // If GUESSBASIS section is specified, project that basis!
        std::cout << "    * Projecting 1PDM_SCALAR" << std::endl;
        auto [onePDMS, projMat] = projectMatrix( this->molecule(), *guessBasisSet, this->basisSet(), std::make_shared<cqmatrix::Matrix<MatsT>>(onePDMtmp->S()) );

        this->onePDM->S() = *onePDMS;
        if( onePDMtmp->hasZ() ) {
          if( this->onePDM->hasZ() ) {
            std::cout << "    * Projecting 1PDM_MZ" << std::endl;
            this->onePDM->Z() = *(projectFullMatrix( projMat, std::make_shared<cqmatrix::Matrix<MatsT>>(onePDMtmp->Z()) ));
          }
          else std::cout << "    * WARNING: Guess has 1PDM_MZ but this reference doesn't! Zeroing out guess MZ..." << std::endl;
        }
        if( onePDMtmp->hasXY() ) {
          if( this->onePDM->hasXY() ) {
            std::cout << "    * Projecting 1PDM_MY" << std::endl;
            this->onePDM->Y() = *(projectFullMatrix( projMat, std::make_shared<cqmatrix::Matrix<MatsT>>(onePDMtmp->Y()) ));
            std::cout << "    * Projecting 1PDM_MX" << std::endl;
            this->onePDM->X() = *(projectFullMatrix( projMat, std::make_shared<cqmatrix::Matrix<MatsT>>(onePDMtmp->X()) ));
          }
          else std::cout << "    * WARNING: Guess has 1PDM_MY/MX but this reference doesn't! Zeroing out guess MY/MX..." << std::endl;
        }
      }

      std::cout << "\n" << std::endl;
      onePDMtmp = nullptr;

    }

  } // SingleSlater<T>::getScr1PDM()

  template <>
  template <>
  void SingleSlater<double,double>::getScr1PDM<dcomplex>(SafeFile& scrBin, const std::shared_ptr<BasisSet> guessBasisSet) {

    CErr("Cannot do complex guess density for real calculation.");

  }

  template <>
  template <>
  void SingleSlater<double,dcomplex>::getScr1PDM<dcomplex>(SafeFile& scrBin, const std::shared_ptr<BasisSet> guessBasisSet) {

    CErr("Cannot do complex guess density for real calculation.");

  }

  /**
   *  \brief Reads in MO from bin file
   *  of different type as calculation
   *  and uses it as initial guess.
   *
   **/
  template <typename MatsT, typename IntsT>
  template <typename ScrMatsT>
  void SingleSlater<MatsT,IntsT>::getScrMO(SafeFile& scrBin) {
    getScrMO<ScrMatsT>(scrBin,nullptr);
  }


  template <typename MatsT, typename IntsT>
  template <typename ScrMatsT>
  void SingleSlater<MatsT,IntsT>::getScrMO(SafeFile& scrBin, const std::shared_ptr<BasisSet> guessBasisSet ) {

    // dimension of mo1 and mo2
    auto NB = this->nC * this->nAlphaOrbital();
    auto NB2 = NB*NB;

    this->mo[0].clear();

    std::string prefix = "/SCF/";
    if (this->particle.charge == 1.0)
        prefix = "/PROT_SCF/";

    auto MO1dims = scrBin.getDims( prefix + "MO1" );
    auto MO2dims = scrBin.getDims( prefix + "MO2" );

    int scrRefType, binRefType;
    scrBin.readData("REF/REFTYPE",&scrRefType);
    savFile.readData("REF/REFTYPE",&binRefType);

    std::cout << "    * Converting from " << refMap[scrRefType] << " to "
      << refMap[binRefType] << std::endl;

    // Let the user know guessbasis is being used
    if( guessBasisSet ){
      std::cout << "    * GUESSBASIS section specified, projecting basis set " <<
        guessBasisSet->basisName << " -> " << this->basisSet().basisName << std::endl;
      if( refMap[scrRefType] != refMap[binRefType] ) CErr("Change of reference type and change of basis set at the same time NYI. Do one at a time.");
      std::cout << "      WARNING: If not projecting for intitial guess," << std::endl;
      std::cout << "               check agreement between READMO and READDEN." << std::endl;
    }

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
    std::cout << "    * Found " << prefix << "MO1 !" << std::endl;
    scrBin.readData(prefix + "MO1",motmp[0].pointer());

    // Unrestricted calculations
    if( scrRefType == RefType::isURef or scrRefType == RefType::isRORef ) {

      if( MO2dims.size() == 0 )
        std::cout << "    * WARNING: " << prefix << "MO2 does not exist in "
          << scrBin.fName() << " -- Copying " << prefix << "MO1 -> " << prefix << "MO2 " << std::endl;

      if( MO2dims.size() > 2  )

        CErr(prefix + "MO2 not saved as a rank-2 tensor in " + scrBin.fName(),
            std::cout);

      // Read in MO2
      std::cout << "    * Found " << prefix << "MO2 !" << std::endl;
      scrBin.readData(prefix + "MO2",motmp[1].pointer());

    }

    // Handle motmp->mo
    // motmp size (assumes alpha and beta same size)
    auto scrMOSize = motmp[0].nRows();

    // Code for different basis set of same ref type
    if( guessBasisSet ){

      // If GUESSBASIS section is specified, project that basis!
      auto [motmp1, projMat] = projectMatrix( this->molecule(), *guessBasisSet, this->basisSet(), std::make_shared<cqmatrix::Matrix<MatsT>>(motmp[0]), false );
      this->mo[0] = *motmp1;

      // project MO2
      if( scrRefType == RefType::isRORef or scrRefType == RefType::isURef ){

         // If GUESSBASIS section is specified, project that basis!
         auto [motmp2, projMat] = projectMatrix( this->molecule(), *guessBasisSet, this->basisSet(), std::make_shared<cqmatrix::Matrix<MatsT>>(motmp[1]), false );
         this->mo[1] = *motmp2;

      }

    }else{     //Code for different ref type of same basis set
    
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

         // RHF->2c
         if( scrRefType == RefType::isRRef ){

           convert1CRto2CU(motmp,this->mo);

         // UHF/ROHF->2c
         } else if( scrRefType == RefType::isURef  or scrRefType == RefType::isRORef ){

          //More information is needed for open-shell systems
          size_t nOccA,nOccB;
          scrBin.readData("REF/NOCCA",&nOccA);
          scrBin.readData("REF/NOCCB",&nOccB); 

          convert1CUto2CU(motmp,this->mo,nOccA,nOccB);

         } else {
            CErr("Initial Guess MO Conversion for 2c Failed");
         } // end 2c guesses

        // 4c guesses
        } else if( binRefType == RefType::isFourCRef ){

          // 4c MOs stored as alpha-large, alpha-small, beta-large, beta-small
          // Negative MOs come before positive MOs

          // RHF->4c
          if( scrRefType == RefType::isRRef ){

           convert1CRto4CU(motmp,this->mo,scrBin);

          }else if( scrRefType == RefType::isURef or scrRefType == RefType::isRORef ){ // UHF/ROHF->4c

           convert1CUto4CU(motmp,this->mo,scrBin);

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

    }

    motmp.clear();

  } // SingleSlater<T>::getScrMO()

  /**
   *  \brief Converts 1-component restricted MOs to 2-component unrestricted
   *
   **/
  template <typename MatsT, typename IntsT>
  template <typename ScrMatsT>
  void SingleSlater<MatsT,IntsT>::convert1CRto2CU(std::vector<cqmatrix::Matrix<ScrMatsT>>& inputMO, std::vector<cqmatrix::Matrix<MatsT>>& outputMO) {

    size_t NB = outputMO[0].nRows();
    size_t scrMOSize = inputMO[0].nRows();

    size_t smallMO=0; // RHF index

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
  void SingleSlater<MatsT,IntsT>::convert1CUto2CU(std::vector<cqmatrix::Matrix<ScrMatsT>>& inputMO,
      std::vector<cqmatrix::Matrix<MatsT>>& outputMO, size_t nA, size_t nB) {

    size_t NB = outputMO[0].nRows();
    size_t scrMOSize = inputMO[0].nRows();
    size_t numC = 2; // Need to hardcode for function so can be used in tandem with other functions

    size_t numUnpaired = nA-nB;
    size_t numPairedMOs = 2*nB;
    size_t numOccMOs = numPairedMOs + numUnpaired;
    size_t plusDisplacedBetas = numOccMOs + numUnpaired;
//  std::cout << "NB, nC, scrMOSize, numUnpaired: " << NB << ", " << numC << ", " << scrMOSize << ", " << numUnpaired << std::endl;

    for( size_t iMO=0; iMO<NB; iMO++ ){

      size_t smallMO = 0; // Which UHF MO to copy
      bool isAlpha = true; // Copy alpha or beta MO from UHF?

      // Fills in 2c MOs as pairs: one for alpha and one for beta occupieds
      if( iMO < numPairedMOs ){

        smallMO = iMO/2; // Does floor function here for beta orbitals
//      std::cout << "inPaired: iMO, smallMO = " << iMO << ", " << smallMO << std::endl;
        isAlpha = (iMO % 2 == 0); // even alpha and odd beta

      // Handles unpaired electrons (will not reach for closed shell)    
      } else if( iMO < numOccMOs ){

        // Assume each 2c MO is alpha for unpaireds
        smallMO = nB + (iMO - numPairedMOs);
//      std::cout << "inUnpaired: iMO, smallMO = " << iMO << ", " << smallMO << std::endl;
        isAlpha = true;

      // Handles those beta virtuals that pair to the occupied alphas
      } else if( iMO < plusDisplacedBetas ){

        // Assume each 2c MO is alpha for unpaireds
        smallMO = nB + (iMO - numOccMOs);
//      std::cout << "inUnpairedBetas: iMO, smallMO = " << iMO << ", " << smallMO << std::endl;
        isAlpha = false;

      // Handles virtuals
      } else {

        // alpha starts at nA and beta starts at nB
        size_t shiftedVirt = iMO - plusDisplacedBetas;

        if( shiftedVirt % 2 == 0 ){
          // alpha virtual
          isAlpha = true;
          smallMO = nA + (shiftedVirt/2);
        } else {
          // beta virtual
          isAlpha = false;
          smallMO = nA + (shiftedVirt/2);
        }
//      std::cout << "inVirtuals: iMO, smallMO = " << iMO << ", " << smallMO << std::endl;

      }
      
      // Alpha block is upper part of 2c spinor and beta block is lower
      if( isAlpha )
        SetMat('N',scrMOSize,1,MatsT(1.),inputMO[0].pointer()+scrMOSize*smallMO,scrMOSize,outputMO[0].pointer()+NB*iMO,NB);
      else{
        SetMat('N',scrMOSize,1,MatsT(1.),inputMO[1].pointer()+scrMOSize*smallMO,scrMOSize,outputMO[0].pointer()+NB*iMO+NB/numC,NB);
      }

    }

  } // SingleSlater<MatsT,IntsT>::convert1CUto2CU

  /**
   *  \brief Converts 1-component unrestricted MOs to 2-component unrestricted with ScrMatsT
   *
   **/
  template <typename MatsT, typename IntsT>
  template <typename ScrMatsT>
  void SingleSlater<MatsT,IntsT>::convert1CUto2CU_sameType(std::vector<cqmatrix::Matrix<ScrMatsT>>& inputMO,
      std::vector<cqmatrix::Matrix<ScrMatsT>>& outputMO, size_t nA, size_t nB) {

    size_t NB = inputMO[0].nRows();
    size_t twoCMOSize = NB*2;

    std::vector<cqmatrix::Matrix<MatsT>> tmp2CMO;
    tmp2CMO.emplace_back(twoCMOSize);
    std::fill_n(tmp2CMO[0].pointer(), twoCMOSize*twoCMOSize, MatsT(0.0));


//  prettyPrintSmart(std::cout, "phi_a start same_type", inputMO[0].pointer(),NB,NB,NB);
//  prettyPrintSmart(std::cout, "phi_b start same_type", inputMO[1].pointer(),NB,NB,NB);

    // Goes from 1CU(ScrMatsT)->2CU(MatsT)
    convert1CUto2CU(inputMO,tmp2CMO,nA,nB);
//  prettyPrintSmart(std::cout, "phi after 1CUto2CU", tmp2CMO[0].pointer(),twoCMOSize,twoCMOSize,twoCMOSize);

    // Element copy to handle unusual type conversions    
    for( size_t i=0; i<twoCMOSize*twoCMOSize; i++ ){
      outputMO[0].pointer()[i] = rough_cast<ScrMatsT>(tmp2CMO[0].pointer()[i]);
    }

  } // SingleSlater<MatsT,IntsT>::convert1CUto2CU_sameType

  /**
   *  \brief Converts 1-component restricted MOs to 4-component unrestricted
   *
   **/
  template <typename MatsT, typename IntsT>
  template <typename ScrMatsT>
  void SingleSlater<MatsT,IntsT>::convert1CRto4CU(std::vector<cqmatrix::Matrix<ScrMatsT>>& inputMO, std::vector<cqmatrix::Matrix<MatsT>>& outputMO, SafeFile& scrBin) {

      size_t NB = inputMO[0].nRows();
      size_t twoCMOSize = inputMO[0].nRows()*2;

//    prettyPrintSmart(std::cout, "phi before 1CUto2CU", inputMO[0].pointer(),NB,NB,NB);

      // Make a temporary copy for 1c transformed MO
      std::vector<cqmatrix::Matrix<ScrMatsT>> tmp1CMO;
      tmp1CMO.emplace_back(inputMO[0].nRows());
      tmp1CMO.emplace_back(inputMO[0].nRows());
      std::fill_n(tmp1CMO[0].pointer(), NB*NB, ScrMatsT(0.0));
      std::fill_n(tmp1CMO[1].pointer(), NB*NB, ScrMatsT(0.0));

      // Make a temporary copy for 2c transformed MO
      std::vector<cqmatrix::Matrix<ScrMatsT>> tmp2CMO;
      tmp2CMO.emplace_back(twoCMOSize);
      std::fill_n(tmp2CMO[0].pointer(), twoCMOSize*twoCMOSize, ScrMatsT(0.0));

      // Converting from 1CR to 1CU
      SetMat('N',NB,NB,ScrMatsT(1.),inputMO[0].pointer(),NB,tmp1CMO[0].pointer(),NB);
      SetMat('N',NB,NB,ScrMatsT(1.),inputMO[0].pointer(),NB,tmp1CMO[1].pointer(),NB);

      //More information is needed for open-shell systems
      size_t nOccA,nOccB;
      scrBin.readData("REF/NOCCA",&nOccA);
      scrBin.readData("REF/NOCCB",&nOccB);

      convert1CUto2CU_sameType(tmp1CMO,tmp2CMO,nOccA,nOccB);
//    prettyPrintSmart(std::cout, "phi after 1CUto2CU-same", tmp2CMO[0].pointer(),twoCMOSize,twoCMOSize,twoCMOSize);

      convert2CUto4CU(tmp2CMO,outputMO,scrBin);
//    prettyPrintSmart(std::cout, "phi after 2CUto4CU", outputMO[0].pointer(),outputMO[0].nRows(),outputMO[0].nRows(),outputMO[0].nRows());


  } // SingleSlater<MatsT,IntsT>::convert1CRto4CU

  /**
   *  \brief Converts 1-component unrestricted MOs to 4-component unrestricted
   *
   **/
  template <typename MatsT, typename IntsT>
  template <typename ScrMatsT>
  void SingleSlater<MatsT,IntsT>::convert1CUto4CU(std::vector<cqmatrix::Matrix<ScrMatsT>>& inputMO, std::vector<cqmatrix::Matrix<MatsT>>& outputMO, SafeFile& scrBin) {

      size_t twoCMOSize = inputMO[0].nRows()*2;

      // Make a temporary copy for 2c transformed MO
      std::vector<cqmatrix::Matrix<ScrMatsT>> tmp2CMO;
      tmp2CMO.emplace_back(twoCMOSize);
      std::fill_n(tmp2CMO[0].pointer(), twoCMOSize*twoCMOSize, ScrMatsT(0.0));

//    prettyPrintSmart(std::cout, "phi at beginning of 1CUto4CU", inputMO[0].pointer(),inputMO[0].nRows(),inputMO[0].nRows(),inputMO[0].nRows());
      
      //More information is needed for open-shell systems
      size_t nOccA,nOccB;
      scrBin.readData("REF/NOCCA",&nOccA);
      scrBin.readData("REF/NOCCB",&nOccB);

      convert1CUto2CU_sameType(inputMO,tmp2CMO,nOccA,nOccB);
//    prettyPrintSmart(std::cout, "phi after 1CUto2CU", tmpMO[0].pointer(),twoCMOSize,twoCMOSize,twoCMOSize);

      convert2CUto4CU(tmp2CMO,outputMO,scrBin);
//    prettyPrintSmart(std::cout, "phi after 2CUto4CU", outputMO[0].pointer(),outputMO[0].nRows(),outputMO[0].nRows(),outputMO[0].nRows());

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
    size_t NBC = outputMO[0].nRows();
    size_t scrMOSize = inputMO[0].nRows();
    size_t NBC2 = NBC*NBC;

    // dimensions: row 2*NP, column 2*NB
    size_t Urow = 2*NP;
    size_t Ucol = 2*NB;

    if( Urow != Ucol ) CErr("Only implemented for uncontracted basis set");

    ScrMatsT *readUL = CQMemManager::get().malloc<ScrMatsT>(Urow*Ucol);
    ScrMatsT *readUS = CQMemManager::get().malloc<ScrMatsT>(Urow*Ucol);

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
      std::fill_n(tmpMO.pointer(), scrMOSize*scrMOSize, MatsT(0.0));

//    prettyPrintSmart(std::cout, "phi^2c in 2CUto4CU", inputMO[0].pointer(),scrMOSize,scrMOSize,scrMOSize);

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
  void SingleSlater<double,double>::getScrMO<dcomplex>(SafeFile& scrBin, const std::shared_ptr<BasisSet> guessBasisSet ) {

    CErr("Cannot do complex guess MOs for real calculation.");

  }

  template <>
  template <>
  void SingleSlater<double,dcomplex>::getScrMO<dcomplex>(SafeFile& scrBin, const std::shared_ptr<BasisSet> guessBasisSet ) {

    CErr("Cannot do complex guess MOs for real calculation.");

  }

}
