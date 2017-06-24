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

#include <mcwavefunction.hpp>
#include <mointstransformer/impl.hpp>
#include <particleintegrals/twopints/incore4indexreleri.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <cqlinalg/blas1.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/matfunc.hpp>
#include <cxxapi/output.hpp>

#include <util/matout.hpp>

// #define DEBUG_MCWFN_MOINTS

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::setMORanges() {
    mointsTF->setMORanges(dynamic_cast<const MCWaveFunctionBase &>(*this));
  }; // MCWaveFunction::setMORanges

  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::transformInts(EMPerturbation & pert) {

    // Precompute the field - nuclear moment contributions
    precompute_NucEField(pert);

    // If the field has changed, we need to rebuild the AOHCore cache
    const double FIELD_DIFF_EPSILON = 1e-15;
    if( pert_has_type(pert,Electric) ) {
      std::array<double, 3> dip_field = pert.getDipoleAmp(Electric);
      double diff = 0.0;
      for(auto iXYZ = 0;    iXYZ < 3;     iXYZ++){
        diff += std::pow(dip_field[iXYZ] - old_dip_field[iXYZ], 2.0);
        old_dip_field[iXYZ] = dip_field[iXYZ];
      }
      diff = std::sqrt(diff);
      this->field_changed = (diff > FIELD_DIFF_EPSILON);
    }
    if(field_changed)
      mointsTF->clearAllCache();

    // TODO: expand to RI
    // build object and allocate memory
    size_t nTOrb  = this->MOPartition.nMO;
    size_t nCorrO = this->MOPartition.nCorrO;
    size_t nInact = this->MOPartition.nInact;
    size_t nInact2 = nInact * nInact;
    
    /*
     * compute inactive core energy
     */
    MatsT * h1e_ii  = CQMemManager::get().malloc<MatsT>(nInact);
    MatsT * GDjj_ii = CQMemManager::get().malloc<MatsT>(nInact);

    double fc1C = (this->reference().nC == 1) ? 2.0 : 1.0;
    
    mointsTF->transformHCore(pert, h1e_ii, "ii", true);
    mointsTF->transformGD(pert, 'i', GDjj_ii, "jj", true, true, "WithInactive-i");  
    
    MatsT ECore = 0.;
    // compute core enenrgy
    for (auto i = 0ul; i < nInact; i++) {
      ECore += h1e_ii[i] + 0.5 * GDjj_ii[i];
    }
    
    CQMemManager::get().free(h1e_ii, GDjj_ii);

    this->InactEnergy = std::real(ECore) * fc1C;
    
    /*
     * compute hCore and ERI in correlated space
     */ 
    OnePInts<MatsT> hCore_tu(nCorrO);
    OnePInts<MatsT> hCoreP_tu(nCorrO);
    
    mointsTF->transformHCore(pert, hCore_tu.pointer(), "tu", false, 'i');
    
    // For RT, the ERI don't need retransformed
    std::shared_ptr<InCore4indexTPI<MatsT>> ERI_tuvw = this->moints->template getIntegral<InCore4indexTPI,MatsT>("ERI_Correlated_Space");
    if(!ERI_tuvw)
    {
      ERI_tuvw = std::make_shared<InCore4indexTPI<MatsT>>(nCorrO);
      mointsTF->transformTPI(pert, ERI_tuvw->pointer(), "tuvw", this->cacheHalfTransTPI_);
      this->moints->addIntegral("ERI_Correlated_Space",ERI_tuvw);
    }
    
    /*
     * compute hCoreP
     */
    hCoreP_tu.clear();
    if (this->MOPartition.scheme == RAS) { // form hCoreP for RAS case
      std::cout<<"RAS hcore prime forming starts:"<<std::endl;
      std::vector<size_t> nActO = this->MOPartition.nActOs;
      std::vector<std::pair<size_t, size_t>> rasloop;
      rasloop.push_back(std::make_pair(0, nActO[0]));
      rasloop.push_back(std::make_pair(nActO[0], nActO[0]+nActO[1]));
      rasloop.push_back(std::make_pair(nActO[0]+nActO[1], nCorrO));

      for (auto i = 0; i < nCorrO; i++)
      for (auto j = 0; j < nCorrO; j++)
        hCoreP_tu(i, j) = hCore_tu(i, j);
      double symmFc;
      for (auto pblk = 0; pblk < 3; pblk++)
      for (auto sblk = 0; sblk < 3; sblk++)
      for (auto qrblk = 0; qrblk < 3; qrblk++) {
        if ((pblk + 3*qrblk) < (qrblk + 3*sblk)) symmFc = 1.0;
        else if ((pblk + 3*qrblk) == (qrblk + 3*sblk)) symmFc = 0.5;
        else continue;

        for (auto p = rasloop[pblk].first; p < rasloop[pblk].second; p++)
        for (auto s = rasloop[sblk].first; s < rasloop[sblk].second; s++)
        for (auto qr = rasloop[qrblk].first; qr < rasloop[qrblk].second; qr++)
          if (pblk == qrblk or sblk == qrblk)
            hCoreP_tu(p, s) -= symmFc * (*ERI_tuvw)(p, qr, qr, s);
          else if (pblk < qrblk)
            hCoreP_tu(p, s) = hCoreP_tu(p, s) - (*ERI_tuvw)(p, qr, qr, s) + (*ERI_tuvw)(qr, qr, p, s);
      }
    } else {

#pragma omp parallel for schedule(static) collapse(2) default(shared)       
      for (auto u = 0ul; u < nCorrO; u++) 
      for (auto t = 0ul; t < nCorrO; t++) {

        MatsT tmp = 0.;
        for (auto v = 0ul; v < nCorrO; v++)
          tmp += 0.5 * (*ERI_tuvw)(t, v, v, u);

        hCoreP_tu(t, u) = hCore_tu(t, u) - tmp;
      }
    }

    this->moints->addIntegral("hCore_Correlated_Space", 
      std::make_shared<OnePInts<MatsT>>(hCore_tu));
    this->moints->addIntegral("hCoreP_Correlated_Space", 
      std::make_shared<OnePInts<MatsT>>(hCoreP_tu));

  }; // MCWaveFunction::transformInts

  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::precompute_NucEField(EMPerturbation & pert)
  {
    // Zero out in case this has already been calculated & stored
    // (for example from RT-CI)
    this->EFieldNuc = 0.0;

    if(!pert_has_type(pert,Electric))
      return;

    std::array<double,3> nucmoment = {0.0,0.0,0.0};
    for(auto & atom : this->reference().molecule().atoms)
    {
      if(atom.quantum) continue;
      MatAdd('N','N',3,1,1.,&nucmoment[0],3,atom.nucCharge,&atom.coord[0],3,&nucmoment[0],3);
    }
      
    auto elecDipoleField = pert.getDipoleAmp(Electric);
    this->EFieldNuc+=blas::dot(3,&nucmoment[0],1,&elecDipoleField[0],1);
  
  // The code below takes an alternate approach:  Instead of folding the field
  // contributions into the one electron integrals, it calculates the dipole 
  // moment for each individual basis function (in this case, slater determinants)
  // as well as the cross terms between them (ie, E\dot<SD_i|x,y,z|SD_j>).
  // This is far less computationally efficient than folding the field contribution
  // into the integrals, but might have future use if say the dipole moment of a 
  // single non-aufbau slater determinant is required for some reason.

  // To use the below code, you must also make corresponding changes in the 
  // hcore portion of mointstransformer in order to prevent the field being
  // double counted by both folding into the integrals and calculating each 
  // slater determinants dipole moment

//    if(this->reference().nC!=1)
//      CErr("CI Dipole Diagonal EField for nC!=1 NYI");
//
//    SingleSlater<MatsT,IntsT> * ss_ptr = &(this->reference());
//
//    size_t NDet = this->NDet;
//    size_t nAO = this->reference().nAlphaOrbital();
//    size_t nI = this->MOPartition.nFCore + this->MOPartition.nInact;
//    size_t nCorrO = this->MOPartition.nCorrO;
//
//    this->EFieldDiag = CQMemManager::get().malloc<MatsT>(NDet*NDet);
//    std::fill_n(this->EFieldDiag,NDet*NDet,MatsT(0.0));
//
//    MatsT * dummyCIi = CQMemManager::get().malloc<MatsT>(NDet);
//    MatsT * dummyCIj = CQMemManager::get().malloc<MatsT>(NDet);
//    SquareMatrix<MatsT> dummyRDM(NDet);
//    SquareMatrix<MatsT> dummyRDM2(NDet);
//    SquareMatrix<MatsT> dummyPDM(nAO);
//
//    auto elecDipoleField = pert.getDipoleAmp(Electric);
//    
//    // Loop over determinants
//    for(size_t i = 0; i < NDet; i++)
//    {
//      // Zero out the CI vector and RDM
//      std::fill_n(dummyCIi,NDet,MatsT(0.0));
//      dummyCIi[i]=1.0;
//      for(size_t j = 0; j < NDet; j++)
//      {
//        if(i==j)
//        {
//          std::fill_n(dummyCIj,NDet,MatsT(0.0));
//          dummyCIj[j]=1.0;
//          dummyRDM.clear();
//          dummyPDM.clear();
//          for(size_t k = 0; k < 3; k++)
//            ss_ptr->elecDipole[k] = 0.0;
//
//          // Fill in the dummy CI Vector and turn it into the RDM
//          this->ciBuilder->computeTDM(*this,dummyCIi,dummyCIj,dummyRDM);
//
//          // Convert the RDM to a PDM
//          // Fold in core orbitals
//          for(size_t k = 0; k < nI; k++)
//            dummyPDM(k,k) = 2.0;
//          // Add RDM contributions
//          SetMat('R',nCorrO,nCorrO,1.0,dummyRDM.pointer(),nCorrO,dummyPDM.pointer()+nI*(nAO+1),nAO);
//          // Transform
//          dummyPDM = dummyPDM.transform('C',this->reference().mo[0].pointer(),nAO,nAO);
//          // Set the reference PDM to this PDM
//          this->reference().onePDM->S() = dummyPDM;
//
//          ss_ptr->computeMultipole(pert);
//
//          this->EFieldDiag[i+j*NDet] = ss_ptr->elecDipole[0]*elecDipoleField[0]+
//                                      ss_ptr->elecDipole[1]*elecDipoleField[1]+
//                                      ss_ptr->elecDipole[2]*elecDipoleField[2];
//        }
//        else
//        {
//          std::fill_n(dummyCIj,NDet,MatsT(0.0));
//          dummyCIj[j]=1.0;
//          dummyRDM.clear();
//          dummyRDM2.clear();
//          dummyPDM.clear();
//          // Fill in the dummy CI Vector and turn it into the RDM
//          this->ciBuilder->computeTDM(*this,dummyCIi,dummyCIj,dummyRDM);
//          this->ciBuilder->computeTDM(*this,dummyCIj,dummyCIi,dummyRDM2);
//
//          auto MOdipole = this->moints.template getIntegral<VectorInts,MatsT>("MOdipole");
//          if (!MOdipole) {
//            std::shared_ptr<VectorInts<IntsT>> AOdipole =
//                      std::make_shared<VectorInts<IntsT>>( nAO, 1, true);
//            std::shared_ptr<VectorInts<MatsT>> MOdipole_scr =
//                      std::make_shared<VectorInts<MatsT>>( nCorrO, 1, true);
//
//            std::vector<std::pair<size_t, size_t>> active(2, {this->MOPartition.nFCore+this->MOPartition.nInact, nCorrO});
//
//            for(auto iXYZ = 0; iXYZ < 3; iXYZ++) {
//              (*AOdipole)[iXYZ] = std::make_shared<OnePInts<IntsT>>( *((*this->reference().aoints_->lenElectric)[iXYZ]) );
//              
//              (*AOdipole)[iXYZ]->subsetTransform('N',this->reference().mo[0].pointer(),
//                  nAO, active, (*MOdipole_scr)[iXYZ]->pointer(), false);
//            }
//
//            this->moints.addIntegral("MOdipole", MOdipole_scr);
//          }
//
//          MOdipole = this->moints.template getIntegral<VectorInts,MatsT>("MOdipole");
//
//          MatsT Etemp = 0.0;
//          for(size_t ixyz = 0; ixyz < 3; ixyz++)
//          {
//            Etemp += elecDipoleField[ixyz]*blas::dotu(nCorrO*nCorrO,dummyRDM.pointer(),1,(*MOdipole)[ixyz]->pointer(),1);
//          }
//          this->EFieldDiag[i+j*NDet]=-Etemp;
//          this->EFieldDiag[j+i*NDet]=-Etemp;
// 
//        }
//      }
//    }
//    // prettyPrintSmart(std::cout,"Matrix-Dipole additions",this->EFieldDiag,NDet,NDet,NDet);
//
  }; //MCWaveFunction::precompute_NucEField

}; // namespace ChronusQ
