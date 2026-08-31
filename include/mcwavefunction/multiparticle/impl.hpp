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

#include <mcwavefunction.hpp>
#include <quantum/preprocessor.hpp>
#include <util/preprocessor.hpp>
#include <util/print.hpp>
 
namespace ChronusQ {

    template <typename MatsT, typename IntsT>
    std::shared_ptr<MCWaveFunctionBase> MultiParticleMCWaveFunction<MatsT,IntsT>::getSubMCWaveFunctionBase(std::string label)
    {
        return std::dynamic_pointer_cast<MCWaveFunctionBase>(subsystems[label]);
    }
    template <typename MatsT, typename IntsT>
    std::shared_ptr<MCWaveFunction<MatsT,IntsT>> MultiParticleMCWaveFunction<MatsT,IntsT>::getSubMCWaveFunction(std::string label)
    {
        return subsystems[label];
    }

    template <typename MatsT, typename IntsT>
    std::vector<std::string> MultiParticleMCWaveFunction<MatsT,IntsT>::getLabels()
    {
        return order_;
    }

    template <typename MatsT, typename IntsT>
    void MultiParticleMCWaveFunction<MatsT,IntsT>::addMCWaveFunction(std::string label, std::shared_ptr<MCWaveFunctionBase> mcwfnbase)
    {

        if(label.empty())           CErr("MCWaveFunction Subsystem label cannot be empty");
        if(!mcwfnbase)              CErr("Cannot add a null mcwavefunction subsystem");
        if(subsystems.count(label)) CErr("Duplicate quantum subsystem label: " + label);

        order_.push_back(label);
        std::shared_ptr<MCWaveFunction<MatsT,IntsT>> mcwfn = std::dynamic_pointer_cast<MCWaveFunction<MatsT,IntsT>>(mcwfnbase);
        if(!mcwfn) CErr("Error in attempting to add MCWaveFunction to MultiParticleMCWavefunction with Label: " + label);
        subsystems.emplace(label,mcwfn);

        mcwfn->StateEnergy = this->StateEnergy;
        mcwfn->moints = this->moints;
        this->NDet *= mcwfn->NDet;

        interIntegrals.try_emplace(label);
        interparticleMOTransformers.try_emplace(label);

    }

    template <typename MatsT, typename IntsT>
    void MultiParticleMCWaveFunction<MatsT,IntsT>::addInteraction(const std::string label1, const std::string label2, const std::shared_ptr<IntegralsBase> & ints)
    {
        if(label1 == label2)          CErr("Interparticle interaction requires two different subsystems.");
        if(!ints)                     CErr("Cannot add null interparticle integrals for " + label1 + "-" + label2);
        if(!subsystems.count(label1)) CErr("Unknown subsystem label: " + label1);
        if(!subsystems.count(label2)) CErr("Unknown subsystem label: " + label2);
        if(interIntegrals.at(label1).count(label2)) 
            CErr("Interaction already exists between " + label1 + " and " + label2);

        auto pairInts = std::dynamic_pointer_cast<Integrals<IntsT>>(ints);
        std::shared_ptr<InCore4indexTPI<IntsT>> tpi = std::dynamic_pointer_cast<InCore4indexTPI<IntsT>>(pairInts->TPI);
        if(!tpi) CErr("MultiparticleMCWavefunction only implemented with Incore TPI currently");

        interIntegrals.at(label1).emplace(label2,tpi);

        // Build the cross-particle MOIntsTransformer
        std::shared_ptr<MixedMOIntsTransformer<MatsT,IntsT>> mixedmointsTF = 
            std::make_shared<MixedMOIntsTransformer<MatsT,IntsT>>(mcSSref_,
                                                                  mcSSref_->getSubSS(label1),
                                                                  mcSSref_->getSubSS(label2),
                                                                  interIntegrals.at(label1).at(label2),
                                                                  true);
        interparticleMOTransformers.at(label1).emplace(label2,mixedmointsTF);
        
    }

    template <typename MatsT, typename IntsT>
    void MultiParticleMCWaveFunction<MatsT,IntsT>::alloc(bool)
    {
        this->CIVecs = std::vector<MatsT*>(this->NStates);
        this->ciBuilder = std::make_shared<MultiParticleCASCI<MatsT,IntsT>>();
        this->multiparticleciBuilder = std::dynamic_pointer_cast<MultiParticleCASCI<MatsT,IntsT>>(this->ciBuilder);
        for(size_t i = 0; i < this->NStates; i++)
        {
            this->CIVecs[i] = CQMemManager::get().template malloc<MatsT>(this->NDet);
        }
        if(this->readCI) this->ReadGuessCIVector();
        ApplyToEach([&](SubMCWfnPtr & mcwfn){mcwfn->alloc(false);});

        // Build the structure for which the CIBuilder will iterate
        for(size_t i = 0; i < order_.size(); i++)
        {
            std::string istring = order_[i];
            std::string istringA = istring;
            std::shared_ptr<MCWaveFunction<MatsT,IntsT>> isubwfn = subsystems.at(istring);
            // Even if both alpha and beta, they will share strings here!
            // Although this should allow unrestricted in the future if that would
            // be interesting (I wouldn't want to do it though)
            std::string ioneBodyIntsStr = istring + "hCoreP_Correlated_Space";
            std::string itwoBodyIntsStr = istring + "ERI_Correlated_Space";

            // For Beta on i
            std::string istringB;
            std::shared_ptr<const ExcitationList> iexList_b;
            size_t iNDetB;
            size_t inCorrA, inCorrB;
            size_t jnCorrA, jnCorrB;

            // SMG 08/03/26
            // Assumes CASStringManager, will break if not!
            std::shared_ptr<const ExcitationList> iexList_a = std::dynamic_pointer_cast<CASStringManager>(isubwfn->detStr)->excitationList();
            if(!iexList_a) CErr("Error in grabbing Excitation List!");
            size_t iNDetA = iexList_a->nString();

            bool ihasBeta = isubwfn->detStrBeta ? true : false;
            // Register the alpha part
            istringA += ihasBeta ? "A" : "";
            inCorrA = isubwfn->MOPartition.nCorrEA;
            OneParticleCIBuilderHelper icibuildera(istringA,iNDetA,ioneBodyIntsStr,itwoBodyIntsStr,iexList_a,isubwfn);
            oneparticlebuilders.insert({istringA,icibuildera});
            //oneparticlebuilders.emplace(istringA,istringA,iNDetA,ioneBodyIntsStr,itwoBodyIntsStr,iexList_a,isubwfn);
            twoparticlebuilders.try_emplace(istringA);
            ciBuilderorder_.push_back(istringA);
            rdmorder_.insert({istringA,i});

            if(ihasBeta)
            {
                istringB = istring + "B";
                iexList_b = std::dynamic_pointer_cast<CASStringManager>(isubwfn->detStrBeta)->excitationList();
                iNDetB = iexList_b->nString();
                inCorrB = isubwfn->MOPartition.nCorrEB;
                OneParticleCIBuilderHelper icibuilderb(istringB,iNDetB,ioneBodyIntsStr,itwoBodyIntsStr,iexList_b,isubwfn);
                oneparticlebuilders.insert({istringB,icibuilderb});
                // Add the two body contribution to the front of the vector so it's done first
                TwoParticleCIBuilderHelper twoparticlehelper(iNDetA,iNDetB,inCorrA,inCorrB,istringA,istringB,itwoBodyIntsStr,iexList_a,iexList_b,isubwfn,isubwfn);
                twoparticlebuilders.at(istringA).push_back(twoparticlehelper);
                twoparticlebuilders.try_emplace(istringB);
                ciBuilderorder_.push_back(istringB);
                rdmorder_.insert({istringB,i});

            }



            // Loop through other particles adding the pairwise interactions
            for(size_t j = i+1; j < order_.size(); j++)
            {
                std::string jstring = order_[j];
                std::string jstringA = jstring;
                std::shared_ptr<MCWaveFunction<MatsT,IntsT>> jsubwfn = subsystems.at(jstring);
                jnCorrA = jsubwfn->MOPartition.nCorrEA;
                jnCorrB = jsubwfn->MOPartition.nCorrEB;
                std::string twoParticleTwoBodyIntsStr = istring + jstring + "ERI_Correlated_Space";
                std::shared_ptr<const ExcitationList> jexList_a = std::dynamic_pointer_cast<CASStringManager>(jsubwfn->detStr)->excitationList();
                if(!jexList_a) CErr("Error in grabbing Excitation List!");

                // For Beta on j
                std::string jstringB;
                std::shared_ptr<const ExcitationList> jexList_b;
                size_t jNDetB;

                bool jhasBeta = jsubwfn->detStrBeta ? true : false;

                size_t jNDetA = jexList_a->nString();
                TwoParticleCIBuilderHelper ijaatwoparticlehelper(iNDetA,jNDetA,inCorrA,jnCorrA,istringA,jstringA,twoParticleTwoBodyIntsStr,iexList_a,jexList_a,isubwfn,jsubwfn);
                twoparticlebuilders.at(istringA).push_back(ijaatwoparticlehelper);

                if(jhasBeta)
                {
                    jstringB = jstring + "B";
                    jexList_b = std::dynamic_pointer_cast<CASStringManager>(jsubwfn->detStrBeta)->excitationList();
                    jNDetB = jexList_b->nString();
                    TwoParticleCIBuilderHelper ijabtwoparticlehelper(iNDetA,jNDetB,inCorrA,jnCorrB,istringA,jstringB,twoParticleTwoBodyIntsStr,iexList_a,jexList_b,isubwfn,jsubwfn);

                }
                if(ihasBeta)
                {
                    TwoParticleCIBuilderHelper ijbatwoparticlehelper(iNDetB,jNDetA,inCorrB,jnCorrA,istringB,jstringA,twoParticleTwoBodyIntsStr,iexList_b,jexList_a,isubwfn,jsubwfn);
                    twoparticlebuilders.at(istringB).push_back(ijbatwoparticlehelper);
                    if(jhasBeta)
                    {
                        TwoParticleCIBuilderHelper ijbbtwoparticlehelper(iNDetB,jNDetB,inCorrB,jnCorrB,istringB,jstringB,twoParticleTwoBodyIntsStr,iexList_b,jexList_b,isubwfn,jsubwfn);
                        twoparticlebuilders.at(istringB).push_back(ijbbtwoparticlehelper);
                    }
                }
            }
        }
        std::cout << BannerTop << std::endl;
        std::cout << "MultiComponent CI Builder Plan: " << std::endl;
        for(const auto & label : ciBuilderorder_)
        {
            oneparticlebuilders.at(label).print(std::cout);
            for(const auto & twobody : twoparticlebuilders.at(label))
            {
                twobody.print(std::cout);
            }
        }
        std::cout << std::endl;
        std::cout << "RDM Builder Plan:" << std::endl;
        for(const auto & label : ciBuilderorder_)
        {
            std::cout << "Subsystem " + label + " contributes to RDM indexed: " << rdmorder_.at(label) << std::endl;
        }
        std::cout << std::endl;
        std::cout << BannerEnd << std::endl;
    }

    void OneParticleCIBuilderHelper::print(std::ostream & out) const
    {
        out << bannerTop << std::endl;
        out << "One Particle CI Builder" << std::endl;
        out << "Particle Label: " << label << std::endl;
        out << "NDet: " << NDet << std::endl;
        out << "One Body Integrals: " << oneBodyIntsStr << std::endl;
        out << "Two Body Integrals: " << twoBodyIntsStr << std::endl;
        out << bannerEnd << std::endl;
    };
    void TwoParticleCIBuilderHelper::print(std::ostream & out) const
    {
        out << bannerTop << std::endl;
        out << "Two Particle CI Builder" << std::endl;
        out << "Interaction between: " << labelp1 << " and " << labelp2 << std::endl;
        out << "NDet each: " << NDetp1 << " and " << NDetp2 << std::endl;
        out << "Total Dets: " << NDet << std::endl;
        out << "Two Body Integrals: " << twoBodyIntsString << std::endl;
        out << bannerEnd << std::endl;
    };

    template <typename MatsT, typename IntsT>
    void MultiParticleMCWaveFunction<MatsT,IntsT>::transformInts(EMPerturbation & pert, std::string)
    {
        this->InactEnergy = 0.0;
        // Each subsystem transform it's own integrals
        ApplyToEachLabeled([&](SubMCWfnPtr & mcwfn, std::string label){
            mcwfn->saveMOInts = this->saveMOInts;
            mcwfn->transformInts(pert,label);
            this->InactEnergy += mcwfn->InactEnergy;});

        // We now transform all the cross particle integrals
        for(size_t i = 0; i < order_.size(); i++)
        {
          for(size_t j = i+1; j < order_.size(); j++)
          {
            std::string l1 = order_[i];
            std::string l2 = order_[j];

            std::string thisERI = l1+l2+"ERI_Correlated_Space";
            std::shared_ptr<InCore4indexTPI<MatsT>> tempERI = this->moints->template getIntegral<InCore4indexTPI,MatsT>(thisERI);
            // If we already have this integral done, move on
            if(tempERI) continue;


            // Get p1 transformer
            auto p1tf = subsystems.at(l1)->mointsTF;
            size_t nCorrP1 = subsystems.at(l1)->MOPartition.nCorrO;
            // Get p2 transformer
            auto p2tf = subsystems.at(l2)->mointsTF;
            size_t nCorrP2 = subsystems.at(l2)->MOPartition.nCorrO;
            // Allocate memory for the transformed integrals
            tempERI = std::make_shared<InCore4indexTPI<MatsT>>(nCorrP1,nCorrP2);
            auto crosstransformer = interparticleMOTransformers.at(l1).at(l2);
            crosstransformer->transformAsymmTPI(pert,tempERI->pointer(),p1tf,p2tf,"tuvw","tuvw");
            this->moints->addIntegral(thisERI,tempERI);

            // Handles core particles
            this->InactEnergy += crosstransformer->handleCrossCoreInts(l1,l2,p1tf,p2tf,this->moints);

          }
        }

        if (this->saveMOInts and MPIRank(this->comm) == 0 and this->savFile.exists()) {
          const auto dataSetLabel = [](const std::string &label) {
            return label == "E" ? std::string("e") : label;
          };

          for (const auto &label : order_) {
            const auto nCorrO = subsystems.at(label)->MOPartition.nCorrO;
            const auto binLabel = dataSetLabel(label);
            auto hCore = this->moints->template getIntegral<OnePInts, MatsT>(
                label + "hCore_Correlated_Space");
            auto hCoreP = this->moints->template getIntegral<OnePInts, MatsT>(
                label + "hCoreP_Correlated_Space");
            auto eri = this->moints->template getIntegral<InCore4indexTPI, MatsT>(
                label + "ERI_Correlated_Space");

            if (hCore)
              this->savFile.safeWriteData("MOINTS/" + binLabel + "_oneelec",
                                          hCore->pointer(), {nCorrO, nCorrO});
            if (hCoreP)
              this->savFile.safeWriteData("MOINTS/" + binLabel + "FOLDED_oneelec",
                                          hCoreP->pointer(), {nCorrO, nCorrO});
            if (eri)
              this->savFile.safeWriteData("MOINTS/" + binLabel + binLabel + "_ERI",
                                          eri->pointer(),
                                          {nCorrO, nCorrO, nCorrO, nCorrO});
          }

          for (size_t i = 0; i < order_.size(); i++) {
            for (size_t j = i + 1; j < order_.size(); j++) {
              const auto &label1 = order_[i];
              const auto &label2 = order_[j];
              const auto nCorrO1 = subsystems.at(label1)->MOPartition.nCorrO;
              const auto nCorrO2 = subsystems.at(label2)->MOPartition.nCorrO;
              auto eri = this->moints->template getIntegral<InCore4indexTPI, MatsT>(
                  label1 + label2 + "ERI_Correlated_Space");
              if (eri)
                this->savFile.safeWriteData(
                    "MOINTS/" + dataSetLabel(label1) + dataSetLabel(label2) + "_ERI",
                    eri->pointer(), {nCorrO1, nCorrO1, nCorrO2, nCorrO2});
            }
          }
        }

    }

}; // namespace ChronusQ

#include <mcwavefunction/multiparticle/print.hpp>
#include <mcwavefunction/multiparticle/rdm.hpp>
#include <mcwavefunction/multiparticle/property.hpp>
