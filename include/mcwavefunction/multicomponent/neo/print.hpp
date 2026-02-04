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
#include <mcwavefunction.hpp>

namespace ChronusQ {

    template <typename MatsT, typename IntsT>
    void NEOMCWaveFunction<MatsT,IntsT>::printMOSpacePartition()
    {
        std::cout << "Electronic space paritioning: " << std::endl;
        this->ewfn_->printMOSpacePartition();
        std::cout << "Protonic space paritioning: " << std::endl;
        this->pwfn_->printMOSpacePartition();
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


    static void printOccKet(std::ostream & out, std::vector<size_t> occ)
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
    void NEOMCWaveFunction<MatsT,IntsT>::printMCState(std::ostream &out, size_t i, double energy, 
                                                      MatsT* C, std::vector<size_t> & sorted_CAddr, 
                                                      size_t N, const size_t n_item_per_row)
    {
        printNEOMCState(out,i,energy,C,sorted_CAddr,N);
    }
 
    template <typename MatsT, typename IntsT>
    void NEOMCWaveFunction<MatsT,IntsT>::printNEOMCState(std::ostream & out, 
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
    void NEOMCWaveFunction<MatsT,IntsT>::printDetOrder(std::ostream & out, std::shared_ptr<DetStringManager> & detstr)
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


    template<typename MatsT, typename IntsT>
    void NEOMCWaveFunction<MatsT,IntsT>::printMOInfo(std::ostream & out,
                                            size_t printMOLevel)
    {
        std::cout << BannerTop << std::endl;
        std::cout << "NEO-MCSCF Print MO Info NYI" << std::endl;
        std::cout << std::endl << BannerTop << std::endl;
    }

}; // namespace ChronusQ