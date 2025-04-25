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
#include <util/matout.hpp>
#include <util/print.hpp>
#include <cxxapi/output.hpp>
#include <mointstransformer/moranges.hpp>

namespace ChronusQ {

    template <typename MatsT, typename IntsT>
    void NEOMCSCF<MatsT,IntsT>::printMCSCFHeader(EMPerturbation & pert)
    {
        std::cout << "Electronic space paritioning: " << std::endl;
        this->ewfn_->printMOSpacePatition();
        std::cout << "Protonic space paritioning: " << std::endl;
        this->pwfn_->printMOSpacePatition();

        // Field print
        if( pert.fields.size() != 0 ) {

            std::cout << "\n\n  * MCSCF will be performed in the presence of an EM "
                << "perturbation:\n\n";

            for(auto &field : pert.fields) {

                auto amp = field->getAmp();

                std::cout << "     * ";
                if( field->emFieldTyp == Electric ) std::cout << "Electric";
                else                                std::cout << "Magnetic";
                
                std::cout << " ";
            
                if( field->size == 3 )        std::cout << "Dipole";
                else if ( field->size == 6 )  std::cout << "Quadrupole";
                else if ( field->size == 10 ) std::cout << "Octupole";
                
                std::cout << " Field: ";
                std::cout << "{ ";
                for(auto i = 0; i < amp.size(); i++) {
                std::cout << amp[i]; if(i != amp.size() - 1) std::cout << ", ";
                }
                std::cout << " }\n";

                }
        }

        // SMG 04/24/24
        // Removing this behavior as it was mostly useful in debugging
        // The printDetOrder functionality prints out how determinants are occupied in order
        // i.e. |0,1>, |0,2>, |0,3>, |1,2>
        // Retaining this code in case it's helpful in the future
        //if(this->printProtonDets)
        //{
        //    printDetOrder(std::cout,this->ewfn_->detStr);
        //    printDetOrder(std::cout,this->ewfn_->detStrBeta);
        //    printDetOrder(std::cout,this->pwfn_->detStr);
        //}

// Code for generating signs of excitation lists
/*
        std::shared_ptr<const ExcitationList> exList_a = std::dynamic_pointer_cast<CASStringManager>(this->ewfn_->detStr)->excitationList();
        size_t nStr_a = exList_a->nString();
        size_t nNZa = exList_a->nNonZero();

        int k,l,Ka,signkl;
        int i,j,Ja,signij;
        // Goal is the sign convention is all relative to the 0'th (if HF) det
        const int * exList_La = exList_a->pointerAtDet(0);
        for(size_t Ekl = 0; Ekl < nNZa; Ekl++, exList_La+=4)
        {
            UNPACK_EXCITATIONLIST_4(exList_La,k,l,Ka,signkl);
            std::cout << "Single excitation: " << l << " " << k << " to index: " << Ka << " with sign: " << signkl << std::endl;
            const int * exList_Ka = exList_a->pointerAtDet(Ka);
            for(size_t Eij= 0; Eij < nNZa; Eij++, exList_Ka+=4)
            {
                UNPACK_EXCITATIONLIST_4(exList_Ka,i,j,Ja,signij);
                std::cout << "Double excitation: " << j << " " << i << " to index: " << Ja << " with sign: " << signij * signkl << std::endl;
            }
        }
*/

    }

    void printOccKet(std::ostream & out, std::vector<size_t> occ)
    {
        // If no particles of the type, print empty ket
        // This can happen if say doing a high spin radical calculation
        // where there are no beta electrons in the active space
        if(!occ.size())
        {
            out << "| > ";
            return;
        }
        out << "|";
        for(size_t i = 0; i < occ.size()-1; i++)
            out << occ[i] << " ";
        out << occ[occ.size()-1] << "> ";
    }
    
    template <typename MatsT, typename IntsT>
    void NEOMCSCF<MatsT,IntsT>::printNEOMCSCFState(std::ostream & out, 
    size_t i, double energy, MatsT * C, std::vector<size_t> & sorted_CAddr,
    size_t N)
    {
        out << std::fixed << std::right << std::setprecision(10);
        out.fill(' ');

        out << std::endl <<  "State:" << std::setw(4) << i + 1 << "  Energy (Hartree):" 
        << std::setw(16) << energy <<  std::endl;
        
        out << std::fixed << std::right<< std::setprecision(7);

        size_t CAddr;
        size_t C_length = 10;

        size_t a_index;
        size_t b_index;
        size_t p_index;

        std::shared_ptr<CASStringManager> astr = std::dynamic_pointer_cast<CASStringManager>(this->ewfn_->detStr); 
        std::shared_ptr<CASStringManager> bstr = std::dynamic_pointer_cast<CASStringManager>(this->ewfn_->detStrBeta); 
        std::shared_ptr<CASStringManager> pstr = std::dynamic_pointer_cast<CASStringManager>(this->pwfn_->detStr); 

        size_t nStr_a = astr->excitationList()->nString();
        size_t nStr_b = bstr->excitationList()->nString();
        size_t nStr_p = pstr->excitationList()->nString();

        std::vector<std::vector<int>> adarray = astr->buildDeAddressingArray(astr->get_addrArray());
        std::vector<std::vector<int>> bdarray = bstr->buildDeAddressingArray(bstr->get_addrArray());
        std::vector<std::vector<int>> pdarray = pstr->buildDeAddressingArray(pstr->get_addrArray());
            
        std::vector<size_t> a_occ;
        std::vector<size_t> b_occ;
        std::vector<size_t> p_occ;

        for(size_t j = 0; j < N; j++)
        {
            CAddr = sorted_CAddr[j];
            out << std::setw(C_length) << std::real(C[CAddr]) << " ";
            a_index = CAddr % nStr_a;
            b_index = (CAddr % (nStr_a * nStr_b)) / nStr_a;
            p_index = CAddr / (nStr_a * nStr_b);
            a_occ = astr->address2DetString(a_index,astr->get_addrArray(),adarray);
            b_occ = bstr->address2DetString(b_index,bstr->get_addrArray(),bdarray);
            p_occ = pstr->address2DetString(p_index,pstr->get_addrArray(),pdarray);
            printOccKet(out,a_occ);
            printOccKet(out,b_occ);
            printOccKet(out,p_occ);
            out << std::endl;
        }

        return;
    }

    template <typename MatsT, typename IntsT>
    void NEOMCSCF<MatsT,IntsT>::printDetOrder(std::ostream & out, std::shared_ptr<DetStringManager> & detstr)
    {
        std::shared_ptr<CASStringManager> subsys= std::dynamic_pointer_cast<CASStringManager>(detstr);
        size_t nStr = subsys->excitationList()->nString();
        std::vector<std::vector<int>> darray = subsys->buildDeAddressingArray(subsys->get_addrArray());
        std::vector<size_t> occ;

        out << std::endl;
        std::cout << " *------------------------------------------------*" << std::endl;      
        out <<       " * Printing Determinant Basis (In Order)          *" << std::endl;
        std::cout << " *------------------------------------------------*" << std::endl;      
        for(size_t i = 0; i < nStr; i++)
        {
            out << i << "\t";
            occ = subsys->address2DetString(i,subsys->get_addrArray(),darray);
            printOccKet(out,occ);
            out << std::endl;
        }
        std::cout << " *------------------------------------------------*" << std::endl;      
        out <<       " * End Determinants                               *" << std::endl;
        std::cout << " *------------------------------------------------*" << std::endl;      
        out << std::endl << std::endl;
    }


/*
    template <typename MatsT, typename IntsT>
    void NEOMCSCF<MatsT,IntsT>::printDetOrder(std::ostream & out)
    {
        std::shared_ptr<CASStringManager> pstr = std::dynamic_pointer_cast<CASStringManager>(this->pwfn_->detStr);
        size_t nStr_p = pstr->excitationList()->nString();
        std::vector<std::vector<int>> pdarray = pstr->buildDeAddressingArray(pstr->addrArray_);
        std::vector<size_t> p_occ;

        out << std::endl;
        std::cout << " *------------------------------------------------*" << std::endl;      
        out <<       " * Printing Protonic Determinant Basis (In Order) *" << std::endl;
        std::cout << " *------------------------------------------------*" << std::endl;      
        for(size_t i = 0; i < nStr_p; i++)
        {
            out << i << "\t";
            p_occ = pstr->address2DetString(i,pstr->addrArray_,pdarray);
            printOccKet(out,p_occ);
            out << std::endl;
        }
        std::cout << " *------------------------------------------------*" << std::endl;      
        out <<       " * End Protonic Determinants                      *" << std::endl;
        std::cout << " *------------------------------------------------*" << std::endl;      
        out << std::endl << std::endl;
    }
*/
    template <typename MatsT, typename IntsT>
    void NEOMCSCF<MatsT,IntsT>::printMCSCFFooter()
    {
        std::cout << std::endl << "MCSCF Results:" << std::endl;
        std::cout << BannerTop << std::endl;

        this->printMOInfo(std::cout);

        std::cout << " *---------------------------------------------*" << std::endl;      
        std::cout << " * Configuration Interaction (CI) Eigen States *" << std::endl;      
        std::cout << " *---------------------------------------------*" << std::endl;      
        
        std::cout << std::endl << BannerTop << std::endl;

        size_t NDet  = this->NDet;
        size_t nS    = this->NStates;
        size_t NPrintC = std::min(NDet, size_t(25));   
        if(this->NDetPrint)
        {
            if(this->NDetPrint == DetPrint::ALLDET)
            {
                NPrintC = NDet;
            }
            else
            {
                NPrintC = this->NDetPrint;
            }
        }

        for (auto i = 0ul; i < nS; i++) { 
        
            // sort coeffients based on the norm 
            auto C = this->CIVecs[i];
            std::vector<size_t> Cindx(NDet);        
            std::iota(Cindx.begin(), Cindx.end(), 0);
            
            std::stable_sort(Cindx.begin(), Cindx.end(), 
                [&] (size_t i , size_t j) {
                return std::norm(C[i]) > std::norm(C[j]);
                }
            );
            
            // only print k Largerst coefficient
            this->printNEOMCSCFState(std::cout, i, this->StateEnergy[i], C, Cindx, NPrintC);  
        }


        this->print1RDMs();
        
        std::cout << BannerTop << std::endl;
    
    
    }; // NEOMCSCF::printMCSCFFooter

}; // namespace ChronusQ