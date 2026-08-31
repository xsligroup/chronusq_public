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
#pragma once

#include <mointstransformer.hpp>

namespace ChronusQ {

  /**
   *  \brief form inactive core density
   */
  template <typename MatsT, typename IntsT>
  std::shared_ptr<cqmatrix::Matrix<MatsT>> MOIntsTransformer<MatsT,IntsT>::formInactDen(
    const char coreIndex) {

      size_t nAO  = ss_.nAlphaOrbital() * ss_.nC;
      auto Den = std::make_shared<cqmatrix::Matrix<MatsT>>(nAO);
      
      auto off_size = parseMOType(coreIndex);
      size_t ioff = off_size.first;
      size_t ni   = off_size.second;
      
      if (ni != 0) {
        blas::gemm(blas::Layout::ColMajor, blas::Op::NoTrans, blas::Op::ConjTrans, 
          nAO, nAO, ni, MatsT(1.), ss_.mo[0].pointer() + ioff*nAO, nAO,
          ss_.mo[0].pointer() + ioff*nAO, nAO, MatsT(0.), Den->pointer(), nAO);
      } else {
        Den->clear();
      }

      return Den;

  }; // form inactive core density
  
  /**
   *  \brief transform MO G(D)  
   */
  template <typename MatsT, typename IntsT>
  std::shared_ptr<OnePInts<MatsT>> MOIntsTransformer<MatsT,IntsT>::formAOGD(
    EMPerturbation & pert, const cqmatrix::Matrix<MatsT> & Den, bool HerDen, 
    bool cacheAOGD, const std::string & cacheId) {  
     
      std::string cacheAOGDStr = "AOGD-" + cacheId; 
      auto AOGD = ints_cache_.template getIntegral<OnePInts, MatsT>(cacheAOGDStr);
     
#ifdef _DEBUG_MOINTSTRANSFORMER_CACHE
      std::cout << "AOCache = " << std::setw(25) << cacheAOGDStr;
      if (AOGD) std::cout << "----Found cache!!!" << std::endl;
      else {
        std::cout << "----Didn't find cached ints, doing transformation" << std::endl;
#else       
      if (not AOGD) { 
#endif
       
        size_t nAO  = ss_.nAlphaOrbital() * ss_.nC;
        // hack thru ss_.forkbuilder->formGD
        if (ss_.nC == 1) {
          *ss_.onePDM = cqmatrix::PauliSpinorMatrices<MatsT>::spinBlockScatterBuild(Den);
        } else{
          *ss_.onePDM = Den.template spinScatter<MatsT>();  
        }
           
        ss_.fockBuilder->formGD(ss_, pert, false, 1.0, HerDen);
        
         
        if (ss_.nC == 1) {
          AOGD = std::make_shared<OnePInts<MatsT>>(0.5 * ss_.twoeH->S());
        } else {
          AOGD = std::make_shared<OnePInts<MatsT>>(
            ss_.twoeH->template spinGather<MatsT>()); 
        }  
        
        if (cacheAOGD) ints_cache_.addIntegral(cacheAOGDStr, AOGD);
      }
       
      return AOGD; 
  }; // MOIntsTransformer::getAOHCoreInCore
  
  /**
   *  \brief transform a subset of HCore 
   */
  template <typename MatsT, typename IntsT>
  std::shared_ptr<OnePInts<MatsT>> MOIntsTransformer<MatsT,IntsT>::formAOHCore(
    EMPerturbation & pert, bool cacheAOHCore, const char coreIndex) {
    
      bool withInactiveCore = coreIndex != '\0';

      std::string cacheAOHCoreStr = "AOHCore";
      if (withInactiveCore) {
        cacheAOHCoreStr += "-WithInactive-";
        cacheAOHCoreStr += coreIndex;
      }
      
      auto AOHCore = ints_cache_.template getIntegral<OnePInts, MatsT>(cacheAOHCoreStr);
       
#ifdef _DEBUG_MOINTSTRANSFORMER_CACHE
      std::cout << "AOCache = " << std::setw(25) << cacheAOHCoreStr;
      if (AOHCore) std::cout << "----Found cache!!!" << std::endl;
      else {
        std::cout << "----Didn't find cached ints, doing transformation" << std::endl;
#else       
      if (not AOHCore) {
#endif
        
        AOHCore = ints_cache_.template getIntegral<OnePInts, MatsT>("AOHCore");
        if (not AOHCore) {
            std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> pertContributions;
            size_t NB = ss_.basisSet().nBasis;
            if (ss_.nC == 4 ) NB = 2 * NB;
            if(ss_.nC > 1)
              pertContributions = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, true);
            else if (not ss_.iCS)
              pertContributions = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, false);
            else
              pertContributions = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(NB, false, false);
            pertContributions->clear();
          // FIXME: the magnetic field contribution is currently absorbed into HCore but that is planned to change
          // https://github.com/xsligroup/chronusq_dev/blob/56c94dfac7a1c28e7fdb4e217052f38ea7d9013a/include/fockbuilder/impl.hpp#L296
          // For RT manipulation, fields will not be added to HCore so we will add the Electric field back in here.
          if( pert_has_type(pert,Electric) ) {
            auto dipAmp = pert.getDipoleAmp(Electric);
            for(auto i = 0;    i < 3;     i++){
                if (dipAmp[i] != 0.0) {
                  pertContributions->S() -= 2. * dipAmp[i] * (*ss_.aoints_->lenElectric)[i]->matrix();
                }
            }
          }
          if(ss_.nC == 1) {
            AOHCore = std::make_shared<OnePInts<MatsT>>(0.5* (ss_.coreH->S() + pertContributions->S()));
          } else { 
            AOHCore = std::make_shared<OnePInts<MatsT>>(
              ss_.coreH->template spinGather<MatsT>() +
              pertContributions-> template spinGather<MatsT>());
          }
          if (cacheAOHCore) ints_cache_.addIntegral("AOHCore", AOHCore);
        }
        
        if (withInactiveCore) {

          auto AOH1e = ints_cache_.template getIntegral<OnePInts, MatsT>("AOHCore");

          auto Den = formInactDen(coreIndex);
          std::string cacheAOGDStr = "WithInactive-";
          cacheAOGDStr += coreIndex;
          auto AOGD = formAOGD(pert, *Den, true, true, cacheAOGDStr); 
          size_t nAO  = ss_.nAlphaOrbital() * ss_.nC;
          
          AOHCore = std::make_shared<OnePInts<MatsT>>(nAO);
          AOHCore->matrix() = AOGD->matrix() + AOH1e->matrix(); 
          
          if (cacheAOHCore) ints_cache_.addIntegral(cacheAOHCoreStr, AOHCore);
        }
      }
  
      return AOHCore; 
  }; // MOIntsTransformer::getAOHCoreInCore
  
  /**
   *  \brief transform a subset of one-partical intgrals 
   */
  template <typename MatsT, typename IntsT>
  void MOIntsTransformer<MatsT,IntsT>::subsetTransformOPI(
    const std::vector<std::pair<size_t,size_t>> &off_sizes, 
    const OnePInts<MatsT> & AOOPI, MatsT* MOOPI, bool deltaPQ) {
    
      size_t nAO = ss_.nAlphaOrbital() * ss_.nC;
      if (not deltaPQ) {
        AOOPI.subsetTransform('N', ss_.mo[0].pointer(), nAO, off_sizes, MOOPI, false); 
      } else {
        
        size_t poff = off_sizes[0].first; 
        size_t np   = off_sizes[0].second; 
        MatsT * SCR = CQMemManager::get().malloc<MatsT>(nAO); 

        for (auto p = 0ul; p < np; p++) {

          auto pMO = ss_.mo[0].pointer() + (p + poff) * nAO;      
          
          // SCR(nu) = MO(mu, p)^H  AOHCOre(mu, nu) 
          blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans, 
            1, nAO, nAO, MatsT(1.), pMO, nAO, AOOPI.pointer(), nAO,
            MatsT(0.), SCR, 1); 
            
          // MOTPI(p,p) = SCR(nu) MO(nu, p)
          MOOPI[p] = blas::dotu(nAO, SCR, 1, pMO, 1);
        }
        CQMemManager::get().free(SCR);
      } 
  
  }; // MOIntsTransformer::subsetTransformOPI
  
  
  /**
   *  \brief transform a subset of HCore 
   */
  template <typename MatsT, typename IntsT>
  void MOIntsTransformer<MatsT,IntsT>::subsetTransformHCore(EMPerturbation & pert,
    const std::vector<std::pair<size_t,size_t>> &off_sizes, MatsT* MOHCore,
    bool deltaPQ, const char coreIndex) {
    
      auto AOHCore = formAOHCore(pert, true, coreIndex);
      subsetTransformOPI(off_sizes, *AOHCore, MOHCore, deltaPQ); 
      
  }; // MOIntsTransformer::subsetTransformHCore 

  /**
   *  \brief transform a subset of GD 
   */
  template <typename MatsT, typename IntsT>
  void MOIntsTransformer<MatsT,IntsT>::subsetTransformGD(EMPerturbation & pert,
    const cqmatrix::Matrix<MatsT> & Den, bool HerDen, MatsT* MOGD,
    const std::vector<std::pair<size_t,size_t>> &off_sizes, 
    bool deltaPQ, bool cacheAOGD, const std::string & cacheId) {
    
      auto AOGD = formAOGD(pert, Den, HerDen, cacheAOGD, cacheId);
      subsetTransformOPI(off_sizes, *AOGD, MOGD, deltaPQ); 
  
  }; // MOIntsTransformer::subsetTransformAOGD 

  template <typename MatsT, typename IntsT>
  double MixedMOIntsTransformer<MatsT,IntsT>::handleCrossCoreInts(std::string label1, 
                                                                  std::string label2,
                                                                  std::shared_ptr<MOIntsTransformer<MatsT,IntsT>> p1tf,
                                                                  std::shared_ptr<MOIntsTransformer<MatsT,IntsT>> p2tf,
                                                                  std::shared_ptr<IntegralsCollection> intscache)
  {
    auto interfockp1p2 = mpss_->getInterFockBuilder(label1,label2);
    auto interfockp2p1 = mpss_->getInterFockBuilder(label2,label1);

    auto hasCore = [&](std::shared_ptr<MOIntsTransformer<MatsT,IntsT>> tf)
    {
      std::vector<std::pair<size_t,size_t>> off_sizes = tf->parseMOType("ijkl");
      bool hasCore = false;
      for(auto os : off_sizes)
        if(os.second) hasCore = true;
      return hasCore;
    };
    bool label1_hasCore = hasCore(p1tf);
    bool label2_hasCore = hasCore(p2tf);

    auto getnCorrO = [&](std::shared_ptr<MOIntsTransformer<MatsT,IntsT>> tf)
    {
      std::vector<std::pair<size_t,size_t>> off_sizes = tf->parseMOType("t");
      return off_sizes[0].second - off_sizes[0].first;
    };

    // All field perturbations are taken care of by the substituent transformers, 
    // so just need a dummy perturbation here
    EMPerturbation pert;

    // At this point we assume that all SingleSlater objects have appropriately had their
    // internal densities overwritten to contain the core densities
    if(label1_hasCore)
    {
      // Build the J matrix using the core density from particle 2
      interfockp2p1->formInterParticleCoulomb(*p2_,false);
      // The outmat is the raw P_i,j (i,j|K,L) -> J_K,L contraction with no additional factors of charge
      auto mat = interfockp2p1->getOutMat();
      OnePInts<MatsT> AOCrossCore(mat->nRows());
      AOCrossCore.matrix() = *mat;
      
      // Preparing to transform J_K,L into the MO basis
      size_t p2nCorrO = getnCorrO(p2tf);
      OnePInts<MatsT> CrossCoreInts(p2nCorrO);
      // Do the transformation
      p2tf->subsetTransformOPI(p2tf->parseMOType("tu"),AOCrossCore,CrossCoreInts.pointer(),false);

      // CrossCoreInts now contains the appropriately sized matrix in the MO basis
      // Grab the respective MOInts from the moints cache
      std::string existingMOIntsString = label2 + "hCoreP_Correlated_Space";
      std::shared_ptr<OnePInts<MatsT>> existingMOInts = intscache->template getIntegral<OnePInts, MatsT>(existingMOIntsString);
      if(!existingMOInts) CErr("Error in grabbing the proper integrals for handling cross core particles!");

      double factor = p1_->particle.charge * p2_->particle.charge;
      existingMOInts->matrix() += factor * CrossCoreInts.matrix();
      intscache->addIntegral(existingMOIntsString,existingMOInts);

    }

    if(label2_hasCore)
    {
      CErr("Currently we assume only core electrons, which will always be label1!");
      // SMG 08/17/26
      // This logic can directly follow the above logic, in fact it probably should
      // be extracted to a lambda or something of the sort to handle both cases
    }

    double coreEnergy = 0.0;
    // If they both have core integrals, we need to add that energy to the inactive
    // energy total (which is what this function ultimately returns)
    if(label1_hasCore && label2_hasCore)
    {
      CErr("Currently multiple core particles not supported!");      
    }
    return coreEnergy;

  }

}; // namespace ChronusQ
