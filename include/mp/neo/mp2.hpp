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
#include <mp/mp2.hpp>


// Numbers fof SOS` scaling
// From https://doi.org/10.1021/acs.jpclett.2c01357
#define NEOSOSPrime_E_SS 0.0
#define NEOSOSPrime_E_OS 1.3
#define NEOSOSPrime_EP   1.3
#define NEOSOSPrime_PP   1.0


namespace ChronusQ {
      
    void summarizeNEOMP2(double eos, double ess, double ep, double pp, double ehf, bool ignorevpp)
    {
      std::cout << std::endl;
      std::cout << bannerTop << std::endl;
      std::cout << "MP2 Summary" << std::endl;
      std::cout << bannerTop << std::endl;

      // WARN THE USER OF MULTIPLE REFERENCE HAMILTONIAN DEFINITIONS
      std::cout << std::endl;
      std::cout << bannerTop << std::endl; 
      std::cout << " ** WARNING **" << std::endl;
      std::cout << bannerTop << std::endl; 
      std::cout << "Your NEO-MP2 Energy is reference hamiltonian" << std::endl;
      std::cout << "dependent! See:" << std::endl;
      std::cout << "https://doi.org/10.1016/j.cplett.2005.01.115" << std::endl;
      std::cout << "You should ensure your MP2 energy is" << std::endl;
      std::cout << "consistent with your choice for \"IGNOREPROTONTWOBODY\"" << std::endl;
      std::cout << "in the [PROTQM] section" << std::endl;
      std::cout << bannerMid << std::endl;
      std::string ignorevppstring = ignorevpp ? "TRUE" : "FALSE";
      std::cout << "For this calculation, IGNOREPROTONTWOBODY= " << ignorevppstring << std::endl;
      std::cout << bannerEnd << std::endl; 

      std::cout << std::endl << std::endl;

      std::cout << bannerTop << std::endl;
      std::cout << "MP2 Opposite Spin:     " << std::setprecision(15) << eos << std::endl;
      std::cout << "MP2 Same Spin:         " << std::setprecision(15) << ess << std::endl;
      std::cout << "MP2 Proton-Proton:     " << std::setprecision(15) << pp << std::endl;
      std::cout << "MP2 Electron-Proton:   " << std::setprecision(15) << ep << std::endl;
      std::cout << "MP2 Total Correlation: " << std::setprecision(15) << eos+ess+pp+ep << std::endl;
      std::cout << bannerMid << std::endl;
      std::cout << "MP2 Total Energy:      " << std::setprecision(15) << ehf + eos + ess + ep + pp << std::endl;
      std::cout << bannerEnd << std::endl;
      std::cout << std::endl;

      std::cout << std::endl;

      double SOSPrimeCorrE = NEOSOSPrime_E_OS*eos + NEOSOSPrime_E_SS*ess + NEOSOSPrime_PP*pp + NEOSOSPrime_EP*ep;

      std::cout << bannerTop << std::endl;
      std::cout << "MP2 SOS' Opposite Spin:      " << std::setprecision(15) << NEOSOSPrime_E_OS * eos << std::endl;
      std::cout << "MP2 SOS' Same SPin:          " << std::setprecision(15) << NEOSOSPrime_E_SS*ess << std::endl;
      std::cout << "MP2 SOS' Proton-Proton:      " << std::setprecision(15) << NEOSOSPrime_PP*pp << std::endl;
      std::cout << "MP2 SOS' Electron-Proton:    " << std::setprecision(15) << NEOSOSPrime_EP*ep << std::endl;
      std::cout << "MP2 SOS' Total Correlation:  " << std::setprecision(15) << SOSPrimeCorrE << std::endl;
      std::cout << bannerMid << std::endl;
      std::cout << "MP2 SOS' Total Energy:       " << std::setprecision(15) << ehf + SOSPrimeCorrE << std::endl;
      std::cout << bannerEnd << std::endl;
      std::cout << std::endl;

    }

    template <typename MatsT, typename IntsT>
    void NEOMP2<MatsT,IntsT>::saveState()
    {
      ROOT_ONLY(this->ref_->comm);
      // PostHartreeFockBase holds onto savFile as a public member
      if(this->savFile.exists())
      {
        this->savFile.safeWriteData("MP2/HF_ENERGY",&this->ref_->totalEnergy,{1});
        this->savFile.safeWriteData("MP2/ECORR_OS",&refs[0]->eos,{1});
        this->savFile.safeWriteData("MP2/ECORR_SS",&refs[0]->ess,{1});
        this->savFile.safeWriteData("MP2/ECORR_EP",&ep,{1});
        this->savFile.safeWriteData("MP2/ECORR_PP",&refs[1]->ess,{1});
        double ecorr = refs[0]->eos+refs[0]->ess+refs[1]->ess+ep;
        this->savFile.safeWriteData("MP2/ECORR",&ecorr,{1});
        double emp2 = ecorr+this->ref_->totalEnergy;
        this->savFile.safeWriteData("MP2/MP2_ENERGY",&emp2,{1});
      }
      return;

    }

    template <typename MatsT, typename IntsT>
    void NEOMP2<MatsT,IntsT>::runMP2(EMPerturbation & pert)
    {

        ProgramTimer::tick("MP2 Total");

        MP2Header();

        // Each subsystem adds its own contributions
        for(auto ssmp : subsystems)
        {
          std::cout << std::endl;
          std::cout << "Running MP2 on a subsystem..." << std::endl;

          auto mp = std::dynamic_pointer_cast<MP2<MatsT,IntsT>>(ssmp);
          refs.push_back(mp);
          NB.push_back(mp->NB);
          nocc.push_back(mp->nocc);
          nvir.push_back(mp->nvir);
          // Set NO request
          mp->makeMP2NOs = this->makeMP2NOs;
          // Calculate each Subsystems individual MP2 energies
          mp->MP2EnergyEval(pert);
          if(this->makeMP2NOs)
          {
              mp->make1RDM();
          }
          // Free the memory for reuse
          mp->dealloc();

          std::cout << "Subsystem MP2 complete!" << std::endl;
          std::cout << std::endl;
        }

        ProgramTimer::tick("Integral Transform");

        std::cout << std::endl;
        std::cout << "Calculating cross terms..." << std::endl;

        std::cout << "Transforming integrals..." << std::endl;

        // Make the cross term transformer
        std::shared_ptr<MixedMOIntsTransformer<MatsT,IntsT>> epTF = nullptr;
        //auto epTF = std::make_shared<MixedMOIntsTransformer<MatsT,IntsT>>
        //    (*neoref_->template getSubsystem<SingleSlater>(std::string("Electronic")),
        //     *neoref_->template getSubsystem<SingleSlater>(std::string("Protonic")),
        //     neoref_->getCrossTPIs(std::string("Electronic"),std::string("Protonic")).second,
        //     neoref_->getCrossTPIs(std::string("Electronic"),std::string("Protonic")).first);

        // Necessary components from each subsystem
        size_t index1 = std::max(nocc[0],nvir[0]);
        size_t index2 = std::max(nocc[1],nvir[1]);
        epERI = std::make_shared<InCore4indexTPI<MatsT>>(index1,index2);
        std::fill_n(epERI->pointer(),index1*index1*index2*index2,MatsT(0.0));
        auto moints1 = std::dynamic_pointer_cast<MP2<MatsT,IntsT>>(subsystems[0])->mointsTF;
        auto moints2 = std::dynamic_pointer_cast<MP2<MatsT,IntsT>>(subsystems[1])->mointsTF;

        // SMG 08/14/24
        // Need to reset the MORanges in the individual transformers because those 
        // could be direct
        moints1->resetMORanges();
        moints1->addMORanges({'t','u','v','w'},{0,nocc[0]});
        moints1->addMORanges({'a','b','c','d'},{nocc[0],nvir[0]});
        moints2->resetMORanges();
        moints2->addMORanges({'t','u','v','w'},{0,nocc[1]});
        moints2->addMORanges({'a','b','c','d'},{nocc[1],nvir[1]});

        // Store T2 is making NOs;
        MatsT * T2base;
        if(this->makeMP2NOs)
        {
            MixedT2 = std::make_shared<InCore4indexTPI<MatsT>>(index1,index2);
            T2base = MixedT2->pointer();
        }

        epTF->transformAsymmTPI(pert,
                                epERI->pointer(),
                                moints1,
                                moints2,
                                "tata",
                                "tata",
                                false,
                                false);
        
        std::cout << "Done transforming integrals!" << std::endl;
        
        ProgramTimer::tock("Integral Transform");

        ProgramTimer::tick("MP2 Energy Eval");

        MatsT * base = epERI->pointer();

        #define NEOMP2_TPI_GET(i_,a_,J_,B_) \
          base[(i_)+(a_)*nocc[0]+(J_)*nocc[0]*nvir[0]+(B_)*nocc[0]*nvir[0]*nocc[1]]

        #define NEOMP2_TPI_SET(i_,a_,J_,B_,val_) \
          T2base[(i_)+(a_)*nocc[0]+(J_)*nocc[0]*nvir[0]+(B_)*nocc[0]*nvir[0]*nocc[1]] = (val_)


        MatsT eps;
        MatsT val;

        double ep_temp = 0.0;

#pragma omp parallel for reduction(+:ep_temp) schedule(dynamic,1) private(eps,val)
      for(size_t i = 0; i < nocc[0]; i++)
      {
        for(size_t a = 0; a < nvir[0]; a++)
        {
          for(size_t J = 0; J < nocc[1]; J++)
          {
            for(size_t B = 0; B < nvir[1]; B++)
            {
              eps = refs[0]->ref_->eps1[i]-refs[0]->ref_->eps1[a+nocc[0]]+ refs[1]->ref_->eps1[J]-refs[1]->ref_->eps1[B+nocc[1]];
              val = NEOMP2_TPI_GET(i,a,J,B)/eps;
              if(this->makeMP2NOs)
              {
                  NEOMP2_TPI_SET(i,a,J,B,val);
              }
              ep_temp+=std::real(2.0*val*NEOMP2_TPI_GET(i,a,J,B));
            }
          }
        }
      }

      ep = ep_temp;

      std::cout << "Done calculating cross terms!" << std::endl;
      std::cout << std::endl;
        
      ProgramTimer::tick("MP2 Energy Eval");

      if(this->makeMP2NOs)
      {
        addepto1RDM();
        for(auto mp : refs)
        {
          mp->generateMP2NOs();
        }
      }

      // Because we know this is NEO, we need to get one of the subsystem fockBuilders
      // (either works since the hamiltonianOptions are the same) 
      auto ppss = this->neoref_->template getSubsystem<SingleSlater>("Protonic");
      bool ignorevpp = ppss->fockBuilder->hamiltonianOptions_.ignoreProtonTwoBody;
      summarizeNEOMP2(refs[0]->eos, refs[0]->ess, ep, refs[1]->ess, this->ref_->totalEnergy,ignorevpp);

      MP2Footer();

      saveState();

      ProgramTimer::tock("MP2 Total");

    } // NEOMP2::runMP2

    template <typename MatsT, typename IntsT>
    void NEOMP2<MatsT,IntsT>::addepto1RDM()
    {
      MatsT * ebase = refs[0]->oneRDM[0]->pointer();
      MatsT * pbase = refs[1]->oneRDM[0]->pointer();
      MatsT * T2Base = MixedT2->pointer();

      size_t eNB = NB[0];
      size_t pNB = NB[1];


      #define E1RDM(i_,j_) \
        ebase[(i_)+(j_)*eNB] 
      #define P1RDM(i_,j_) \
        pbase[(i_)+(j_)*pNB] 
      #undef NEOMP2_TPI_GET
      #define NEOMP2_TPI_GET(i_,a_,J_,B_) \
        T2Base[(i_)+(a_)*nocc[0]+(J_)*nocc[0]*nvir[0]+(B_)*nocc[0]*nvir[0]*nocc[1]]

      // Electron-Electron block
      // Occupied-Occupied Block
      for(size_t i = 0; i < nocc[0]; i++)
      {
        for(size_t j = 0; j < nocc[0]; j++)
        {
          for(size_t I = 0; I < nocc[1]; I++)
          {
            for(size_t a = 0; a < nvir[0]; a++)
            {
              for(size_t A = 0; A < nvir[1]; A++)
              {
                E1RDM(i,j)-=2.0*NEOMP2_TPI_GET(i,a,I,A)*NEOMP2_TPI_GET(j,a,I,A);
              }
            } 
          }
        }
      }
      // Virtual-Virtual Block
      for(size_t a = 0; a < nvir[0]; a++)
      {
        for(size_t b = 0; b < nvir[0]; b++)
        {
          for(size_t A = 0; A < nvir[1]; A++)
          {
            for(size_t i = 0; i < nocc[0]; i++)
            {
              for(size_t I = 0; I < nocc[1]; I++)
              {
                E1RDM(a+nocc[0],b+nocc[0])+=2.0*NEOMP2_TPI_GET(i,a,I,A)*NEOMP2_TPI_GET(i,b,I,A);
              }
            }
          }
        }
      }

      // Proton 1 RDM
      // Occupied-Occupied Block
      for(size_t I = 0; I < nocc[1]; I++)
      {
        for(size_t J = 0; J < nocc[1]; J++)
        {
          for(size_t i = 0; i < nocc[0]; i++)
          {
            for(size_t A = 0; A < nvir[1]; A++)
            {
              for(size_t a = 0; a < nvir[0]; a++)
              {
                P1RDM(I,J)-=2.0*NEOMP2_TPI_GET(i,a,I,A)*NEOMP2_TPI_GET(i,a,J,A);
              }
            } 
          }
        }
      }
      // Virtual-Virtual Block
      for(size_t A = 0; A < nvir[1]; A++)
      {
        for(size_t B = 0; B < nvir[1]; B++)
        {
          for(size_t a = 0; a < nvir[0]; a++)
          {
            for(size_t I = 0; I < nocc[1]; I++)
            {
              for(size_t i = 0; i < nocc[0]; i++)
              {
                P1RDM(A+nocc[1],B+nocc[1])+=2.0*NEOMP2_TPI_GET(i,a,I,A)*NEOMP2_TPI_GET(i,a,I,B);
              }
            }
          }
        }
      }

      return;
    } // addepto1RDM


}; // namespace ChronusQ
