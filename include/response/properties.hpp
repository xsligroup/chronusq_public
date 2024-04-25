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

#include <response/tbase.hpp>

#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasext.hpp>
#include <cqlinalg/blasutil.hpp>

#include <physcon.hpp>

#include <util/matout.hpp>

namespace ChronusQ {


  template <typename T>
  template <typename U>
  std::map<ResponseOperator, U*> 
    ResponseTBase<T>::evalProperties(
      std::vector<ResponseOperator> ops, size_t nVec, U* V) {

    std::map<ResponseOperator, U*> opMap; 


    size_t N      = nSingleDim_;

    T* g; size_t nProp;

    for(auto &op : ops) {

      // Compute property gradient
      std::tie(nProp,g) = formPropGrad(op);


      // Allocate space for the property
      opMap[op] = CQMemManager::get().malloc<U>(nProp*nVec);
      std::fill_n(opMap[op],nProp*nVec,U(0.));

      // Evaluate the property (ensures proper behaviour for mixed type)
      blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,nProp,nVec,N,U(1.),g,N,V,N,U(0.),opMap[op],nProp);
      IMatCopy('C',nProp,nVec,U(1.),opMap[op],nProp,nVec);


      // Free scratch space
      CQMemManager::get().free(g);

    }

    return opMap;

  }



  template <typename T>
  void ResponseTBase<T>::residueTMoments() {

    if( genSettings.aOps.size() == 0 ) return;


    std::vector<ResponseOperator> herOps, antiHerOps;


    std::copy_if(genSettings.aOps.begin(),genSettings.aOps.end(),
        std::back_inserter(herOps),isHerOp);
    std::copy_if(genSettings.aOps.begin(),genSettings.aOps.end(),
        std::back_inserter(antiHerOps),isAntiHerOp);
      

    T* VR = resResults.VR; 
    T* VL = resSettings.needVL ? resResults.VL : VR;

    auto herOpMap     = evalProperties(herOps,resSettings.nRoots,VR);
    auto antiHerOpMap = evalProperties(antiHerOps,resSettings.nRoots,VL);

    auto opMap = herOpMap;
    opMap.insert(antiHerOpMap.begin(),antiHerOpMap.end());

    // Place the property evaluations
      
    // Electric Dipole (Length)
    if( opMap.find(LenElectricDipole) != opMap.end() ) {

      resResults.tLenElecDipole_ge = opMap[LenElectricDipole];

      IMatCopy('T',resSettings.nRoots,3,T(1.),resResults.tLenElecDipole_ge,
        resSettings.nRoots,3);

      if( savFile.exists() )
        savFile.safeWriteData("/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_LENGTH",
          resResults.tLenElecDipole_ge, {resSettings.nRoots,3});
    }




    // Electric Quadrupole (Length)
    if( opMap.find(LenElectricQuadrupole) != opMap.end() ) {

      resResults.tLenElecQuadrupole_ge = opMap[LenElectricQuadrupole];

      IMatCopy('T',resSettings.nRoots,6,T(1.),resResults.tLenElecQuadrupole_ge,
        resSettings.nRoots,6);

      if( savFile.exists() )
        savFile.safeWriteData(
          "/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_LENGTH",
          resResults.tLenElecQuadrupole_ge, {resSettings.nRoots, 6});
    }




    // Electric Octupole (Length)
    if( opMap.find(LenElectricOctupole) != opMap.end() ) {
      resResults.tLenElecOctupole_ge = opMap[LenElectricOctupole];

      IMatCopy('T',resSettings.nRoots,10,T(1.),resResults.tLenElecOctupole_ge,
        resSettings.nRoots,10);

      if( savFile.exists() )
        savFile.safeWriteData(
          "/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_LENGTH",
          resResults.tLenElecOctupole_ge, {resSettings.nRoots, 10});
    }



    // Electric Dipole (Velocity)
    if( opMap.find(VelElectricDipole) != opMap.end() ) {

      resResults.tVelElecDipole_ge = opMap[VelElectricDipole];

      IMatCopy('T',resSettings.nRoots,3,T(1.),resResults.tVelElecDipole_ge,
        resSettings.nRoots,3);

      if( savFile.exists() )
        savFile.safeWriteData(
          "/RESP/RESIDUE/TRANSITION_ELECTRIC_DIPOLE_VELOCITY",
          resResults.tVelElecDipole_ge, {resSettings.nRoots,3});
    }




    // Electric Quadrupole (Velocity)
    if( opMap.find(VelElectricQuadrupole) != opMap.end() ) {
      resResults.tVelElecQuadrupole_ge = opMap[VelElectricQuadrupole];

      IMatCopy('T',resSettings.nRoots,6,T(1.),resResults.tVelElecQuadrupole_ge,
        resSettings.nRoots,6);

      if( savFile.exists() )
        savFile.safeWriteData(
          "/RESP/RESIDUE/TRANSITION_ELECTRIC_QUADRUPOLE_VELOCITY",
          resResults.tVelElecQuadrupole_ge, {resSettings.nRoots, 6});
    }




    // Electric Octupole (Velocity)
    if( opMap.find(VelElectricOctupole) != opMap.end() ) {
      resResults.tVelElecOctupole_ge = opMap[VelElectricOctupole];

      IMatCopy('T',resSettings.nRoots,10,T(1.),resResults.tVelElecOctupole_ge,
        resSettings.nRoots,10);

      if( savFile.exists() )
        savFile.safeWriteData(
          "/RESP/RESIDUE/TRANSITION_ELECTRIC_OCTUPOLE_VELOCITY",
          resResults.tVelElecOctupole_ge, {resSettings.nRoots, 10});
    }



    // Magnetic Dipole
    if( opMap.find(MagneticDipole) != opMap.end() ) {
      resResults.tMagDipole_ge = opMap[MagneticDipole];

      IMatCopy('T',resSettings.nRoots,3,T(1.),resResults.tMagDipole_ge,
        resSettings.nRoots,3);

      if( savFile.exists() )
        savFile.safeWriteData("/RESP/RESIDUE/TRANSITION_MAGNETIC_DIPOLE",
          resResults.tMagDipole_ge, {resSettings.nRoots,3});
    }




    // Magnetic Quadrupole
    if( opMap.find(MagneticQuadrupole) != opMap.end() ) {
      resResults.tMagQuadrupole_ge = opMap[MagneticQuadrupole];

      IMatCopy('T',resSettings.nRoots,6,T(1.),resResults.tMagQuadrupole_ge,
        resSettings.nRoots,6);

      if( savFile.exists() )
        savFile.safeWriteData("/RESP/RESIDUE/TRANSITION_MAGNETIC_QUADRUPOLE",
          resResults.tMagQuadrupole_ge, {resSettings.nRoots,6});
    }
      
  }

  template <typename T>
  void ResponseTBase<T>::residueObservables() {

    size_t nRoots = resSettings.nRoots;
    double * W = resResults.W;

    if( resResults.tLenElecDipole_ge ) {

      resObs.oscStrength = 
        CQMemManager::get().malloc<double>(nRoots);

      for(auto iO = 0; iO < nRoots; iO++) {

        double omega = W[iO];
        T* tDipole = resResults.tLenElecDipole_ge + 3*iO;

        double tDipoleX = std::abs(tDipole[0]);
        double tDipoleY = std::abs(tDipole[1]);
        double tDipoleZ = std::abs(tDipole[2]);


        resObs.oscStrength[iO] = 2. * omega / 3. * 
          (tDipoleX * tDipoleX + tDipoleY * tDipoleY + tDipoleZ * tDipoleZ);

      }

      if( savFile.exists() )
        savFile.safeWriteData("/RESP/RESIDUE/OSC_STRENGTH",resObs.oscStrength,
          {nRoots});

    }






    if( resResults.tLenElecDipole_ge and resResults.tMagDipole_ge ) {

      resObs.rotatory_len_RM = 
        CQMemManager::get().malloc<double>(nRoots);

      for(auto iO = 0; iO < nRoots; iO++) {

        double omega = W[iO];

        T* tElecDipole = resResults.tLenElecDipole_ge + 3*iO;
        T* tMagDipole  = resResults.tMagDipole_ge     + 3*iO;

        double rXX = std::real( tElecDipole[0] * tMagDipole[0] );
        double rYY = std::real( tElecDipole[1] * tMagDipole[1] );
        double rZZ = std::real( tElecDipole[2] * tMagDipole[2] );

        resObs.rotatory_len_RM[iO] = -0.5 * ( rXX + rYY + rZZ );

      }

      if( savFile.exists() )
        savFile.safeWriteData("/RESP/RESIDUE/ROT_STRENGTH_LEN_EDA",
          resObs.rotatory_len_RM, {nRoots});

    }





  };


  template <typename T>
  template <typename U>
  void ResponseTBase<T>::fdrRF(FDResponseResults<T,U> &results) {

    if( genSettings.aOps.size() == 0 ) return;


    std::vector<ResponseOperator> ops = genSettings.aOps;
      

    size_t nOmega = fdrSettings.bFreq.size();
    size_t nRHS   = fdrSettings.nRHS;


    auto opMap = evalProperties(ops,nOmega*nRHS,results.SOL);



    auto getOpOffset = [&]( ResponseOperator op ) -> int {

      auto opPos =
        std::find(genSettings.bOps.begin(),genSettings.bOps.end(),op);

      return std::accumulate(genSettings.bOps.begin(),opPos,0,
               [](int prev, ResponseOperator op2){ 
                 return prev + OperatorSize[op2];
               });

    };

    auto hasOp = 
      []( std::vector<ResponseOperator> &tops, ResponseOperator op ) -> bool {

      return std::find(tops.begin(),tops.end(),op) != tops.end();

    };


    int edlOff = getOpOffset(LenElectricDipole);
    int mdOff  = getOpOffset(MagneticDipole);


    bool AHasEDL = hasOp(genSettings.aOps,LenElectricDipole);
    bool AHasEQL = hasOp(genSettings.aOps,LenElectricQuadrupole);
    bool AHasEOL = hasOp(genSettings.aOps,LenElectricOctupole);
    bool AHasMD  = hasOp(genSettings.aOps,MagneticDipole);
    bool AHasMQ  = hasOp(genSettings.aOps,MagneticQuadrupole);


    bool BHasEDL = hasOp(genSettings.bOps,LenElectricDipole);
    bool BHasEQL = hasOp(genSettings.bOps,LenElectricQuadrupole);
    bool BHasEOL = hasOp(genSettings.bOps,LenElectricOctupole);
    bool BHasMD  = hasOp(genSettings.bOps,MagneticDipole);
    bool BHasMQ  = hasOp(genSettings.bOps,MagneticQuadrupole);

    // Electric Dipole - Electric Dipole Polarizability (Len)
    if( AHasEDL and BHasEDL ) {


      results.ed_ed_Polar = 
        CQMemManager::get().malloc<U>(9*nOmega);

      for(auto iOmega = 0; iOmega < nOmega; iOmega++) 

        SetMat('N',3,3,U(1.),opMap[LenElectricDipole] + edlOff + nRHS*iOmega,
          nRHS*nOmega, results.ed_ed_Polar + iOmega*9 ,3);

      if( savFile.exists() )
        savFile.safeWriteData("/RESP/FDR/ED_ED_POLARIZABILITY_LENGTH",
          results.ed_ed_Polar, {nOmega,3,3});


    }

    // Electric Quadrupole - Electric Dipole Polarizability (Len)
    if( AHasEQL and BHasEDL ) {


      results.eq_ed_Polar = 
        CQMemManager::get().malloc<U>(3*6*nOmega);

      for(auto iOmega = 0; iOmega < nOmega; iOmega++) {

        SetMat('N',3,6,U(1.),opMap[LenElectricQuadrupole] +edlOff+ nRHS*iOmega,
          nRHS*nOmega, results.eq_ed_Polar + iOmega*3*6 ,3);



        // Make traceless
        U* qdStart = results.eq_ed_Polar + iOmega*3*6; 

        for(auto i = 0; i < 3; i++) {
  
          U qTrace = qdStart[i + 0*3] + 
                     qdStart[i + 3*3] + 
                     qdStart[i + 5*3];

          qdStart[i + 0*3] -= qTrace / 3.; 
          qdStart[i + 3*3] -= qTrace / 3.; 
          qdStart[i + 5*3] -= qTrace / 3.;

        }
                   
        // Transpose in place
        IMatCopy('T',3,6,U(1.),qdStart,3,6);

      }

      if( savFile.exists() )
        savFile.safeWriteData("/RESP/FDR/EQ_ED_POLARIZABILITY_LENGTH",
          results.eq_ed_Polar, {nOmega,3,6});

    }





    // Magnetic Dipole - Electric Dipole Polarizability (Len)
    if( AHasMD and BHasEDL ) {


      results.md_ed_Polar = 
        CQMemManager::get().malloc<U>(9*nOmega);

      for(auto iOmega = 0; iOmega < nOmega; iOmega++){ 

        SetMat('N',3,3,U(1.),opMap[MagneticDipole] + edlOff + nRHS*iOmega,
          nRHS*nOmega, results.md_ed_Polar + iOmega*9 ,3);

        // Transpose in place
        IMatCopy('T',3,3,U(1.),results.md_ed_Polar + iOmega*9,3,3);

      }

      if( savFile.exists() )
        savFile.safeWriteData("/RESP/FDR/MD_ED_POLARIZABILITY_LENGTH",
          results.md_ed_Polar, {nOmega,3,3});

    }


    if( AHasMD and BHasMD ) {


      results.md_md_Polar = 
        CQMemManager::get().malloc<U>(9*nOmega);

      for(auto iOmega = 0; iOmega < nOmega; iOmega++){ 

        SetMat('N',3,3,U(1.),opMap[MagneticDipole] + mdOff + nRHS*iOmega,
          nRHS*nOmega, results.md_md_Polar + iOmega*9 ,3);

        // Transpose in place
        IMatCopy('T',3,3,U(1.),results.md_md_Polar + iOmega*9,3,3);

      }

      if( savFile.exists() )
        savFile.safeWriteData("/RESP/FDR/MD_MD_POLARIZABILITY",
          results.md_md_Polar, {nOmega,3,3});


    }


    // Free up the full memory
    for(auto &op : opMap) CQMemManager::get().free(op.second);

    
  };



  template <typename T>
  template <typename U>
  void ResponseTBase<T>::fdrObservables(FDResponseResults<T,U> &results) {

    

    size_t nOmega = fdrSettings.bFreq.size();

    if( results.ed_ed_Polar ) {

      fdObs.edStrength = 
        CQMemManager::get().malloc<double>(nOmega);

      fdObs.opaCross_eda = 
        CQMemManager::get().malloc<double>(nOmega);

      for(auto iOmega = 0; iOmega < nOmega; iOmega++) {

        double omega = fdrSettings.bFreq[iOmega];
        U*     eded  = results.ed_ed_Polar + iOmega*9;

        U trace = eded[0] + eded[4] + eded[8];

        fdObs.edStrength[iOmega]   = std::imag( trace );

        fdObs.opaCross_eda[iOmega] = 4. * M_PI / SpeedOfLight *
          omega * fdObs.edStrength[iOmega];

      }

      if( savFile.exists() )
        savFile.safeWriteData("/RESP/FDR/OPA_CROSS_SECTION_EDA",
          fdObs.opaCross_eda,{nOmega});

    }

  }

  


}; // namespace ChronusQ

