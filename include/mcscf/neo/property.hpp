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

#include <mcscf.hpp>
#include <mcwavefunction.hpp>
#include <util/matout.hpp>
#include <util/print.hpp>
#include <cxxapi/output.hpp>
#include <mointstransformer/moranges.hpp>
#include <cqlinalg/blas3.hpp>
#include <cqlinalg/blasutil.hpp>
#include <cqlinalg/matfunc.hpp>


namespace ChronusQ {

    template<typename MatsT, typename IntsT>
    void NEOMCSCF<MatsT,IntsT>::computeMultipole()
    {
        SingleSlater<MatsT,IntsT> * protSS = dynamic_cast<SingleSlater<MatsT,IntsT>*>(&this->pwfn_->reference());
        size_t nTotalP = protSS->molecule_.nTotalP;
        for(size_t i = 0; i < this->NStates; i++)
        {
            this->computeMultipole(i);
            if(nTotalP==1)
            {
                std::cout << std::endl;
                std::cout << " *---------------------------------------------*" << std::endl;
                std::cout << " *Additional Proton Properties for state:" + std::to_string(i+1) + "             *" << std::endl;
                std::cout << " *---------------------------------------------*" << std::endl;
                this->ProtonExpectationValue(i);
                this->ProtonVariance(i);
                this->ProtonKE(i);
                std::cout << " *---------------------------------------------*" << std::endl;
            }
        }
    }

    template <typename MatsT, typename IntsT>
    void NEOMCSCF<MatsT,IntsT>::computeMultipole(size_t s1)
    {
        // Grab a reference to the NEOSS object
        NEOSS<MatsT,IntsT> * neo_ss = &neoref_;

        this->ewfn_->rdm2pdm(this->oneRDM[s1]);
        this->pwfn_->rdm2pdm(this->PoneRDM[s1]);

        EMPerturbation emPert;
        neo_ss->computeMultipole(emPert);

        std::cout << std::endl;
        std::cout << " *---------------------------------------------*" << std::endl;
        std::cout << " *          Multipole for State " + std::to_string(s1+1) + "             *" << std::endl;
        std::cout << " *---------------------------------------------*" << std::endl;
        neo_ss->printMultipoles(std::cout);

        // Adding the multipole in a vector over states
        this->elecDipoles.push_back(neo_ss->elecDipole);
        this->elecQuadrupoles.push_back(neo_ss->elecQuadrupole);
        this->elecOctupoles.push_back(neo_ss->elecOctupole);

        return;
    }

    template<typename MatsT, typename IntsT>
    double NEOMCSCF<MatsT,IntsT>::oscillator_strength(size_t s2, size_t s1)
    {
        // Need to do both the electronic and nuclear dipole calculations
        // result is the sum of the two
        MatsT D = MatsT(0.0);

        // Electronic component
        size_t nAO = this->ewfn_->reference().nAlphaOrbital();
        size_t nCorrO = this->ewfn_->MOPartition.nCorrO;
        size_t nInact = this->ewfn_->MOPartition.nInact;
        cqmatrix::Matrix<MatsT> tmpTDM1(nCorrO);
        cqmatrix::Matrix<MatsT> tmpTDM2(nCorrO);
        this->ciBuilder->computeTDM(*this,this->CIVecs[s1],this->CIVecs[s2],tmpTDM1);
        this->ciBuilder->computeTDM(*this,this->CIVecs[s2],this->CIVecs[s1],tmpTDM2);

        auto MOdipole = this->moints->template getIntegral<VectorInts,MatsT>("MOdipole");

        if(not MOdipole)
        {
            std::shared_ptr<VectorInts<IntsT>> AOdipole = std::make_shared<VectorInts<IntsT>>(nAO,1,true);
            std::shared_ptr<VectorInts<MatsT>> MOdipole_scr = std::make_shared<VectorInts<MatsT>>(nCorrO,1,true);
            // This is just offsets for the subset transform call
            std::vector<std::pair<size_t, size_t>> active(2, {this->ewfn_->MOPartition.nFCore+nInact, nCorrO});
            for(size_t iXYZ = 0; iXYZ < 3; iXYZ++)
            {
                (*AOdipole)[iXYZ] = (*(this->ewfn_->reference()).aoints_->lenElectric)[iXYZ];
//                (*AOdipole)[iXYZ].subsetTransform('N',this->ewfn_->reference().mo[0].pointer(),nAO,active,(*MOdipole_scr)[iXYZ].pointer(),false);
                (*AOdipole)[iXYZ]->subsetTransform('N',this->ewfn_->reference().mo[0].pointer(),nAO, active, (*MOdipole_scr)[iXYZ]->pointer(), false);
            }
            this->moints->addIntegral("MOdipole",MOdipole_scr);
        }
        
        MOdipole = this->moints->template getIntegral<VectorInts,MatsT>("MOdipole");
        for(auto iXYZ = 0; iXYZ < 3; iXYZ++)
        {
            D+=blas::dotu(nCorrO*nCorrO,tmpTDM1.pointer(),1,(*MOdipole)[iXYZ]->pointer(),1)
              *blas::dotu(nCorrO*nCorrO,tmpTDM2.pointer(),1,(*MOdipole)[iXYZ]->pointer(),1);
        }

        // Protonic component
        nAO = this->pwfn_->reference().nAlphaOrbital();
        nCorrO = this->pwfn_->MOPartition.nCorrO;
        nInact = this->pwfn_->MOPartition.nInact;
        cqmatrix::Matrix<MatsT> ptmpTDM1(nCorrO);
        cqmatrix::Matrix<MatsT> ptmpTDM2(nCorrO);
        NEOCIBuilder->computePTDM(*this,this->CIVecs[s1],this->CIVecs[s2],ptmpTDM1);
        NEOCIBuilder->computePTDM(*this,this->CIVecs[s2],this->CIVecs[s1],ptmpTDM2);

        auto PMOdipole = this->moints->template getIntegral<VectorInts,MatsT>("PMOdipole");

        if(not PMOdipole)
        {
            std::shared_ptr<VectorInts<IntsT>> PAOdipole = std::make_shared<VectorInts<IntsT>>(nAO,1,true);
            std::shared_ptr<VectorInts<MatsT>> PMOdipole_scr = std::make_shared<VectorInts<MatsT>>(nCorrO,1,true);
            // This is just offsets for the subset transform call
            std::vector<std::pair<size_t, size_t>> active(2, {this->pwfn_->MOPartition.nFCore+nInact, nCorrO});
            for(size_t iXYZ = 0; iXYZ < 3; iXYZ++)
            {
                (*PAOdipole)[iXYZ] = (*(this->pwfn_->reference()).aoints_->lenElectric)[iXYZ];
                (*PAOdipole)[iXYZ]->subsetTransform('N',this->pwfn_->reference().mo[0].pointer(),nAO, active, (*PMOdipole_scr)[iXYZ]->pointer(), false);
            }
            this->moints->addIntegral("PMOdipole",PMOdipole_scr);
        }
        PMOdipole = this->moints->template getIntegral<VectorInts,MatsT>("PMOdipole");
        for(auto iXYZ = 0; iXYZ < 3; iXYZ++)
        {
            D+=blas::dotu(nCorrO*nCorrO,ptmpTDM1.pointer(),1,(*PMOdipole)[iXYZ]->pointer(),1)
              *blas::dotu(nCorrO*nCorrO,ptmpTDM2.pointer(),1,(*PMOdipole)[iXYZ]->pointer(),1);
        } 

        double f = (2./3.) * (this->StateEnergy[s2]-this->StateEnergy[s1]) * std::real(D);

        std::cout << "Excited State: " << std::setw(3) << std::right << s2+1
                << " to state: " << std::setw(3) << std::right << s1+1 << ":";
        std::cout << std::setw(15) << std::right << "E(Eh) = "
                << std::setprecision(8) << std::fixed << (this->StateEnergy[s2] - this->StateEnergy[s1])*EVPerHartree;
        std::cout << std::setw(15) << std::right << "f = "
                << std::setprecision(12) << std::fixed << f << std::endl;

        return f;
    }

    template<typename MatsT, typename IntsT>
    void NEOMCSCF<MatsT,IntsT>::ProtonExpectationValue(size_t s1)
    {
        SingleSlater<MatsT,IntsT> * protSS = dynamic_cast<SingleSlater<MatsT,IntsT>*>(&this->pwfn_->reference());
        std::shared_ptr<cqmatrix::Matrix<double>> PDM = this->getPOnePDM()[s1];
        std::array<cqmatrix::Matrix<MatsT>,3> xyz = {(*protSS->aoints_->lenElectric)[0]->matrix(),
                                                     (*protSS->aoints_->lenElectric)[1]->matrix(),
                                                     (*protSS->aoints_->lenElectric)[2]->matrix()};
        std::vector<std::string> xyzstrs({"<X> = ","<Y> = ","<Z> = "});
        size_t NB = protSS->basisSet_.nBasis;
        std::array<double,3> pos = ProtonExpectationValue(*PDM);
        std::cout << "Proton expectation value (angstrom): " << std::endl; 
        for(size_t iXYZ = 0; iXYZ < 3; iXYZ++)
            std::cout << "      " << xyzstrs[iXYZ] << std::fixed << std::setw(10) << std::left << std::setprecision(6) << pos[iXYZ]*AngPerBohr << std::endl;

    } // Proton Expectation value

    template<typename MatsT, typename IntsT>
    std::array<double,3> NEOMCSCF<MatsT,IntsT>::ProtonExpectationValue(cqmatrix::Matrix<double> PDM)
    {
        SingleSlater<MatsT,IntsT> * protSS = dynamic_cast<SingleSlater<MatsT,IntsT>*>(&this->pwfn_->reference());
        std::array<cqmatrix::Matrix<MatsT>,3> xyz = {(*protSS->aoints_->lenElectric)[0]->matrix(),
                                                     (*protSS->aoints_->lenElectric)[1]->matrix(),
                                                     (*protSS->aoints_->lenElectric)[2]->matrix()};
        size_t NB = protSS->basisSet_.nBasis;
        std::array<double,3> pos;
        for(size_t iXYZ = 0; iXYZ < 3; iXYZ++)
        {
            MatsT r = -blas::dot(NB*NB,PDM.pointer(),1,xyz[iXYZ].pointer(),1);
            pos[iXYZ] = std::real(r);
        }
        return pos;
    } // Proton Expectation value

    template<typename MatsT, typename IntsT>
    void NEOMCSCF<MatsT,IntsT>::ProtonVariance(size_t s1)
    {
        std::shared_ptr<cqmatrix::Matrix<double>> PDM = this->getPOnePDM()[s1];
        std::vector<std::string> xyzstrs({"<X^2>-<X>^2 = ","<Y^2>-<Y>^2 = ","<Z^2>-<Z>^2 = "});
        std::cout << "Poton variance (bohr^2): " << std::endl; 
        std::array<double,3> var = ProtonVariance(*PDM);
        for(size_t iXYZ = 0; iXYZ < 3; iXYZ++)
            std::cout << "      " << xyzstrs[iXYZ] << std::fixed << std::setw(10) << std::left << std::setprecision(6) << var[iXYZ] << std::endl;
    } // Proton Variance

    template<typename MatsT, typename IntsT>
    std::array<double,3> NEOMCSCF<MatsT,IntsT>::ProtonVariance(cqmatrix::Matrix<double> PDM)
    {
        SingleSlater<MatsT,IntsT> * protSS = dynamic_cast<SingleSlater<MatsT,IntsT>*>(&this->pwfn_->reference());
        std::array<cqmatrix::Matrix<MatsT>,3> xyz = {(*protSS->aoints_->lenElectric)[0]->matrix(),
                                                     (*protSS->aoints_->lenElectric)[1]->matrix(),
                                                     (*protSS->aoints_->lenElectric)[2]->matrix()};
        std::array<cqmatrix::Matrix<MatsT>,3> xxyyzz = {(*protSS->aoints_->lenElectric)[3]->matrix(),
                                                        (*protSS->aoints_->lenElectric)[6]->matrix(),
                                                        (*protSS->aoints_->lenElectric)[8]->matrix()};
        size_t NB = protSS->basisSet_.nBasis;
        std::array<double,3> var;
        for(size_t iXYZ = 0; iXYZ < 3; iXYZ++)
        {
            MatsT r = -blas::dot(NB*NB,PDM.pointer(),1,xyz[iXYZ].pointer(),1);
            MatsT rr = -blas::dot(NB*NB,PDM.pointer(),1,xxyyzz[iXYZ].pointer(),1);
            var[iXYZ] = std::real(rr - r * r);
        }
        return var;
    } // Proton Variance

    template<typename MatsT, typename IntsT>
    void NEOMCSCF<MatsT,IntsT>::ProtonKE(size_t s1)
    {
        std::shared_ptr<cqmatrix::Matrix<double>> PDM = this->getPOnePDM()[s1];
        std::cout << "Kinetc Energy Expectation value (au): ";
        std::cout << std::setw(10) << std::left << std::setprecision(6) << ProtonKE(*PDM) << std::endl;
    } // Proton KE

    template<typename MatsT, typename IntsT>
    double NEOMCSCF<MatsT,IntsT>::ProtonKE(cqmatrix::Matrix<double> PDM)
    {
        SingleSlater<MatsT,IntsT> * protSS = dynamic_cast<SingleSlater<MatsT,IntsT>*>(&this->pwfn_->reference());
        cqmatrix::Matrix<MatsT> AOKE = protSS->aoints_->kinetic->matrix();
        size_t NB = protSS->basisSet_.nBasis;
        return std::real(blas::dot(NB*NB,PDM.pointer(),1,AOKE.pointer(),1));

    } // Proton KE


}; // namespace ChronusQ

