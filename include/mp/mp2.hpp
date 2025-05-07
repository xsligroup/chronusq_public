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

#include <singleslater.hpp>
#include <mointstransformer.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>

namespace ChronusQ {

   void MP2Header()
   {
      std::cout << BannerTop << std::endl;
      std::cout << "Begin MP2 calculation" << std::endl;
      std::cout << BannerEnd << std::endl;
   }

   void MP2Footer()
   {
      std::cout << BannerTop << std::endl;
      std::cout << "End MP2 calculation" << std::endl;
      std::cout << BannerEnd << std::endl;
   }

   void summarizeMP2(double eos, double ess, double ehf)
   {
      std::cout << std::endl;
      std::cout << bannerTop << std::endl;
      std::cout << "MP2 Summary " << std::endl;
      std::cout << bannerTop << std::endl;

      std::cout << "MP2 Opposite Spin:     " << std::setprecision(15) << eos << std::endl;
      std::cout << "MP2 Same Spin:         " << std::setprecision(15) << ess << std::endl;
      std::cout << "MP2 Total Correlation: " << std::setprecision(15) << eos+ess << std::endl;
      std::cout << std::endl;
      std::cout << "MP2 Total Energy:      " << std::setprecision(15) << ehf + eos + ess << std::endl;
      std::cout << bannerEnd << std::endl;
      std::cout << std::endl;

   };

   template <typename MatsT, typename IntsT>
   void MP2<MatsT,IntsT>::saveState()
   {
      ROOT_ONLY(this->ref_->comm);
      // PostHartreeFockBase holds onto savFile as a public member
      if(this->savFile.exists())
      {
        this->savFile.safeWriteData("MP2/HF_ENERGY",&this->ref_->totalEnergy,{1});
        this->savFile.safeWriteData("MP2/ECORR_OS",&eos,{1});
        this->savFile.safeWriteData("MP2/ECORR_SS",&ess,{1});
        double ecorr = eos+ess;
        this->savFile.safeWriteData("MP2/ECORR",&ecorr,{1});
        double emp2 = eos+ess+this->ref_->totalEnergy;
        this->savFile.safeWriteData("MP2/MP2_ENERGY",&emp2,{1});
      }
      return;
   }

   template <typename MatsT, typename IntsT>
   void MP2<MatsT,IntsT>::runMP2(EMPerturbation & pert)
   {
      MP2Header();

      ProgramTimer::tick("MP2 Total");

      MP2EnergyEval(pert);
      if(this->makeMP2NOs)
      {
         make1RDM();
         generateMP2NOs();
      }

      summarizeMP2(eos,ess,this->ref_->totalEnergy);

      saveState();

      ProgramTimer::tock("MP2 Total");

      MP2Footer();
   }

   template <typename MatsT, typename IntsT>
   void MP2<MatsT,IntsT>::MP2EnergyEval(EMPerturbation & pert)
   {
      ProgramTimer::tick("MP2 Energy Eval");

      // Because the Direct_N6 integral transform for MP2 does NOT create the full
      // (ia|jb) i,j -> occ, a,b -> virt
      // tensor, instead doing it piece wise, we have separate code to actually do
      // the energy evaluation
      if(this->ref_->aoints_->TPITransAlg == TPI_TRANSFORMATION_ALG::DIRECT_N6)
      {
        if(this->makeMP2NOs)
        {
          CErr("MP2 Natural Orbitals with Direct Integrals NYI!");
        }  
        MP2DirectEnergyEval(pert);
      }
      else
      {
        MP2InCoreEnergyEval(pert);
      }
        
    // Handle factor of 2 for NEO Protons
    if(this->ref_->particle.charge==1.0)
    {
        eos = 0.0;
        ess/= 2.0;
    }

    ProgramTimer::tock("MP2 Energy Eval");

     return;
   } // MP2::MP2EnergyEval

  template <typename MatsT, typename IntsT>
  void MP2<MatsT,IntsT>::MP2DirectEnergyEval(EMPerturbation & pert)
  {
    // Generate MOIntsTransformer
    this->mointsTF = this->ref_->generateMOIntsTransformer(TPI_TRANSFORMATION_ALG::DIRECT_N6);

    // Determine how big our blocks should be
    size_t occvirpairsize = nocc * nvir;

    // Divide the (ia|jb) into blocks
    // SMG 08/12/24 TODO:
    // There's almost certainly some optimal blocking analysis which can be done
    // here to figure out the optimal sizes
    // For now, take the limit of i,j always contained within the block, and a,b 
    // split into blocks of 10 (remember, we need all unique i,j pairs)
    size_t ab_block_size = 10;
    // Check if we have enough available memory for the moERI tensor
    size_t moeri_size = nocc * nocc * ab_block_size * ab_block_size;
    if(moeri_size != CQMemManager::get().max_avail_allocatable<MatsT>(1,moeri_size))
    {
      CErr("Not enough memory for incore evaluation of MP2 block size!");
    }
    // Allocate the temporary moERI tensor
    moERI = std::make_shared<InCore4indexTPI<MatsT>>(nocc,ab_block_size);

    // Generate the unique pairs of (a,b) 
    std::vector<std::pair<size_t,size_t>> ab_pairs;
    for(size_t a = 0; a < nvir; a+=ab_block_size)
    {
      // Only need unique a,b pairs since we will calculate both (ia|jb) and (ib|ja) 
      // at the same time
      for(size_t b = a; b < nvir; b+=ab_block_size)
      {
        ab_pairs.push_back({a,b});
      }
    }

    // If only one unique a,b pairing, no need to allocate the second moERI block
    if(ab_pairs.size() > 1)
    {
      exmoERI = std::make_shared<InCore4indexTPI<MatsT>>(nocc,ab_block_size);
    }

    // Loop over a,b pairs
    for(auto abpair : ab_pairs)
    {
      ProgramTimer::tick("Integral Transform");

      // Required to zero out since the underlying integral transform
      // may not calculate integrals
      std::fill_n(moERI->pointer(),moeri_size,MatsT(0.0));
      if(abpair.first != abpair.second)
      {
        std::fill_n(exmoERI->pointer(),moeri_size,MatsT(0.0));
      }

      // Prepare the transformer for the new ranges
      this->mointsTF->resetMORanges();
      this->mointsTF->addMORanges({'i','j'},{0,nocc});
      size_t num_a = abpair.first + ab_block_size < nvir ? ab_block_size : nvir-abpair.first;
      size_t num_b = abpair.second + ab_block_size < nvir ? ab_block_size : nvir-abpair.second;
      this->mointsTF->addMORanges({'a'},{abpair.first+nocc,num_a});
      this->mointsTF->addMORanges({'b'},{abpair.second+nocc,num_b});
      this->mointsTF->printMORangesSummary();

      // Transform the integrals
      this->mointsTF->directTransformTPI(pert,moERI->pointer(),"iajb");
      if(abpair.first!=abpair.second)
      {
        this->mointsTF->directTransformTPI(pert,exmoERI->pointer(),"ibja");
      }

      // Scale the integrals by 0.5
      for(size_t i = 0; i < nocc * nocc * num_a * num_b; i++)
      {
        moERI->pointer()[i]*=0.5;
        if(abpair.first!=abpair.second)
        {
          exmoERI->pointer()[i]*=0.5;
        }
      } 
    
      ProgramTimer::tock("Integral Transform");

      // Add this block's contribution to the total energy
      if(abpair.first == abpair.second)
      {
        MP2EnergyLooper(moERI->pointer(),moERI->pointer(),0,nocc,abpair.first,num_a,0,nocc,abpair.second,num_b);
      }
      else
      {
        MP2EnergyLooper(moERI->pointer(),exmoERI->pointer(),0,nocc,abpair.first,num_a,0,nocc,abpair.second,num_b);
      }
    } // Loop of a,b pairs
  } // MP2::MP2DirectEnergyEval

  template <typename MatsT, typename IntsT>
  void MP2<MatsT,IntsT>::MP2InCoreEnergyEval(EMPerturbation & pert)
  {
    ProgramTimer::tick("Integral Transform");
    std::cout << "Transforming integrals..." << std::endl;

    // moERI is size nocc*nocc*nvir*nvir, although we need it as nocc*nvir*nocc*nvir
    // This will mean we will access the underlying integrals in a custom way assuming
    // the storage is p + q * NB1 + r * NB1*NB1 + s * NB1 * NB1 * NB2
    // IF THE INCORE4INDEXTPI STORAGE CHANGES, THIS WILL BREAK
    moERI = std::make_shared<InCore4indexTPI<MatsT>>(nocc,nvir);
    // Only allocate storage for T2 integrals if we're making the 1 RDM
    if(this->makeMP2NOs)
    {
      std::cout << std::endl;
      std::cout << "Natural orbitals requested" << std::endl;
      std::cout << "Generating and storing T2 amplitudes!" << std::endl;
      std::cout << std::endl;
      // T2 storage is <ij|ab> in physisicts notation
      T2 = std::make_shared<InCore4indexTPI<MatsT>>(nocc,nvir); 
    }
    
    // Transform integrals
    this->mointsTF = this->ref_->generateMOIntsTransformer(TPI_TRANSFORMATION_ALG::INCORE_N5);
    // Convince mointsTF that there are nocc active and NB-nocc virtual
    this->mointsTF->resetMORanges();
    this->mointsTF->addMORanges({'t','u','v','w'},{0,nocc});
    this->mointsTF->addMORanges({'a','b','c','d'},{nocc,nvir});
    this->mointsTF->printMORangesSummary();
    this->mointsTF->transformTPI(pert,moERI->pointer(),"taub",false);

    std::cout << "Done transforming integrals!" << std::endl;
    
    ProgramTimer::tock("Integral Transform");
      
    // Evaluate the energy contribution
    MP2EnergyLooper(moERI->pointer(),moERI->pointer(),0,nocc,0,nvir,0,nocc,0,nvir);

    return;
  } // MP2::MP2IncoreEnergyEval

  // Loops over (ia|jb) integral and adds their opposite and same spin contributions
  // to the MP2 energy expression
  // Note that astart and bstart are offset by nocc
  // that is to say, astart = 0 implies the first virtual orbital
  template <typename MatsT, typename IntsT>
  void MP2<MatsT,IntsT>::MP2EnergyLooper(MatsT * ERIbase, MatsT * exERIbase,
                                         size_t istart, size_t istride,
                                         size_t astart, size_t astride,
                                         size_t jstart, size_t jstride,
                                         size_t bstart, size_t bstride)
  {
    #define MP2_TPI_GET(base_,i_,j_,a_,b_) \
      base_[(i_)+(a_)*istride+(j_)*istride*astride+(b_)*istride*astride*jstride]

    // Assumes storage of complement integrals is symmetric 
    #define MP2_TPI_exGET(base_,i_,j_,a_,b_) \
      base_[(i_)+(b_)*istride+(j_)*istride*bstride+(a_)*istride*bstride*jstride]

    #define MP2_TPI_SET(base_,i_,j_,a_,b_,val_) \
      base_[(i_)+(j_)*istride+(a_)*istride*jstride+(b_)*istride*jstride*astride] = val_

    MatsT * T2base;
    if(this->makeMP2NOs) T2base  = T2->pointer();
    MatsT eps;
    MatsT val;

    // Temporary reduction variables (wouldn't compile without?)
    double eos_temp = 0.0;
    double ess_temp = 0.0;

#pragma omp parallel for reduction(+:eos_temp,ess_temp) schedule(dynamic,1) private(eps,val)
    for(size_t i = istart; i < istride; i++)
    {
      for(size_t j = jstart; j < jstride; j++)
      {
        for(size_t a = 0; a < astride; a++)
        {
          for(size_t b = 0; b < bstride; b++)
          {
              // (ia|jb) contribution
              eps = this->ref_->eps1[i]+this->ref_->eps1[j]-this->ref_->eps1[a+nocc+astart]-this->ref_->eps1[b+nocc+bstart];
              val = MP2_TPI_GET(ERIbase,i,j,a,b)/eps;
              if(this->makeMP2NOs)
              {
                MP2_TPI_SET(T2base,i,j,a,b,val);
              }
              eos_temp+=std::real(val*MP2_TPI_GET(ERIbase,i,j,a,b));
              ess_temp+=std::real(val*(MP2_TPI_GET(ERIbase,i,j,a,b)-MP2_TPI_exGET(exERIbase,i,j,a,b)));

              // (ib|ja) contribution, if a,b different classes of orbitals
              if(ERIbase!=exERIbase)
              {
                val = MP2_TPI_exGET(exERIbase,i,j,a,b)/eps;
                if(this->makeMP2NOs)
                {
                  // SMG 08/13/24
                  // NOT DEBUGGED
                  // SHOULD NEVER REACH THIS POINT, CERR MUCH EARLIER
                  CErr("MP2 Natural Orbitals with Direct Integrals NYI!");
                }
                eos_temp+=std::real(val*MP2_TPI_exGET(exERIbase,i,j,a,b));
                ess_temp+=std::real(val*(MP2_TPI_exGET(exERIbase,i,j,a,b)-MP2_TPI_GET(ERIbase,i,j,a,b)));
              }
          }
        }
      }
    }

    eos += eos_temp;
    ess += ess_temp;

    return;
  } // MP2::EnergyLooper

  template <typename MatsT, typename IntsT>
  void MP2<MatsT,IntsT>::make1RDM()
  {

    ProgramTimer::tick("MP2 RDM");

    // Use PostHartreeFock::oneRDM
    this->oneRDM.push_back(std::make_shared<cqmatrix::Matrix<MatsT>>(NB));
    // For convenience, pointer to oneRDM base
    MatsT * base = this->oneRDM[0]->pointer();
    std::fill_n(base,NB*NB,MatsT(0.0));

    // The diagonal is occupied with the standard occupation,
    // which is 2 for RHF and 1 for NEO UHF
    MatsT docc = this->ref_->particle.charge == 1.0 ? 1.0 : 2.0;

    #define MP2_1RDM_MAT(i_,j_) \
      base[(i_)+(j_)*NB]

    for(size_t ii = 0; ii < nocc; ii++)
    {
        MP2_1RDM_MAT(ii,ii) = docc;
    }

    // 

    for(size_t i = 0; i < nocc; i++)
    {
      for(size_t j = 0; j < nocc; j++)
      {
        for(size_t k = 0; k < nocc; k++)
        {
          for(size_t a = 0; a < nvir; a++)
          {
            for(size_t b = 0; b < nvir; b++)
            {
               MP2_1RDM_MAT(i,j)-=docc*(*T2)(i,k,a,b)*(docc*(*T2)(j,k,a,b)-(*T2)(j,k,b,a));
            }
          }
        }
      }
    }


    for(size_t a = 0; a < nvir; a++)
    {
      for(size_t b = 0; b < nvir; b++)
      {
        for(size_t c = 0; c < nvir; c++)
        {
          for(size_t i = 0; i < nocc; i++)
          {
            for(size_t j = 0; j < nocc; j++)
            {
              MP2_1RDM_MAT(a+nocc,b+nocc)+=docc*(*T2)(i,j,b,c)*(docc*(*T2)(i,j,a,c)-(*T2)(i,j,c,a));
            }
          }
        }
      }
    }

    ProgramTimer::tock("MP2 RDM");

    return;
  } // MP2::make1RDM

  template <typename MatsT, typename IntsT>
  void MP2<MatsT,IntsT>::generateMP2NOs()
  {
     ProgramTimer::tick("MP2 RDM Diag");

     // Save the 1RDM here because it is possibly modified by multiple functions in a NEO context 
     if(this->savFile.exists())
     {
       std::string loc = "MP2/1RDM";
       if(this->ref_->particle.charge > 0)
        loc = "MP2/P_1RDM";
       this->savFile.safeWriteData(loc,this->oneRDM[0]->pointer(),{NB,NB});
     }

     std::shared_ptr<cqmatrix::Matrix<MatsT>> MP2NOs = std::make_shared<cqmatrix::Matrix<MatsT>>(NB);

     std::cout << std::endl;
     std::cout << "Diagonalizing the MP2 1RDM" << std::endl;
     std::cout << std::endl;

     // Diagonalize the MP2 1RDM
     MP2_NO_occs.resize(NB);
     HermitianEigen('V', 'U', NB, this->oneRDM[0]->pointer(), NB, MP2_NO_occs.data());

     // Need to reverse the order of the 1RDM
     std::reverse(MP2_NO_occs.begin(),MP2_NO_occs.end()); 
     cqmatrix::Matrix<MatsT> Unitary(NB);
     MatsT * Eigvec = this->oneRDM[0]->pointer() + (NB-1)*NB;
     MatsT * SCR = Unitary.pointer();
     for(size_t n = 0; n < NB; n++, SCR+=NB, Eigvec-=NB)
     {
        std::copy_n(Eigvec,NB,SCR);
     }

     std::cout << std::endl;
     std::cout << bannerTop << std::endl;

     std::cout << " Natural Orbital Occupation Numbers" << std::endl;
     std::cout << bannerTop << std::endl;
     std::cout << "   Occupied Orbitals" << std::endl;
     std::cout << bannerMid << std::endl;
     for(size_t i = 0; i < this->ref_->nOA; i++)
     {
        std::cout << "   " << std::setprecision(6) << std::setw(10) << MP2_NO_occs[i] << std::endl;
     }
     
     std::cout << bannerMid << std::endl;
     std::cout << "   Virtual Orbitals" << std::endl;
     std::cout << bannerMid << std::endl;
     for(size_t i = this->ref_->nOA; i < NB; i++)
     {
        std::cout << "   " << std::setprecision(6) << std::setw(10) << MP2_NO_occs[i] << std::endl;
     }

     std::cout << bannerEnd << std::endl;

     blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,NB,NB,NB,
        MatsT(1.0),this->ref_->mo[0].pointer(),NB,Unitary.pointer(),NB,MatsT(0.0),
        MP2NOs->pointer(),NB);
     std::copy_n(MP2NOs->pointer(),NB*NB,this->ref_->mo[0].pointer());

     // Write to bin the updated orbitals
     // TODO: Build a cubegen interface which generates the MP2NO cubes
     this->ref_->saveCurrentState();

     // Write to bin the natural orbital occupations
     ROOT_ONLY(this->ref_->comm);
     if(this->savFile.exists())
     {
       std::string loc = "MP2/NO_OCCS";
       if(this->ref_->particle.charge > 0)
        loc = "MP2/P_NO_OCCS";
       this->savFile.safeWriteData(loc,&this->MP2_NO_occs[0],{NB});
     }

     ProgramTimer::tock("MP2 RDM Diag");

  } // MP2::generateMP2NOs

}; // namespace ChronusQ
