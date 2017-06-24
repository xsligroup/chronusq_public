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

  void MCSCFSettings::print(bool fourComp, size_t nStates) {
    
    std::cout << std::endl;
    FormattedLine(std::cout,"* Computation Parameters:");
    if(ciAlg == CIDiagonalizationAlgorithm::CI_FULL_MATRIX) {
      FormattedLine(std::cout,"  CI Algorithm:",  "Full Matrix");
    } else if (ciAlg == CIDiagonalizationAlgorithm::CI_DAVIDSON) {
      FormattedLine(std::cout,"  CI Algorithm:",  "Davidson");
      FormattedLine(std::cout,"  CI Maxmium Number of Iteration:",  maxCIIter);
      FormattedLine(std::cout,"  CI Vector Convergence Threshold:", ciVectorConv);
      FormattedLine(std::cout,"  Max Len of Davidson Subspace (x NRoots):", maxDavidsonSpace);
      FormattedLine(std::cout,"  Number of Davidson Guess(x NRoots):", nDavidsonGuess);

      if (!energyRefs.empty()) {
        FormattedLine(std::cout,"  Energy specific settings:");
        size_t nLowRoots = nStates;
        for (auto & pair: energyRefs) {
          FormattedLine(std::cout, "  Energy threshold:", pair.first, " #Roots:", pair.second);        
          nLowRoots -= pair.second;
        }
        FormattedLine(std::cout,"  Number of low energy roots:", nLowRoots);
      }

    } else CErr("NYI CI Algorithm");
    
    if(this->doSCF) {
      std::cout << std::endl;
      ORSettings.print(fourComp);
      
      std::cout << std::endl;
      FormattedLine(std::cout,"  SCF Maxmium Number of Iteration:",    maxSCFIter);
      FormattedLine(std::cout,"  SCF Energy Convergence Threshold:",   scfEnergyConv);
      FormattedLine(std::cout,"  SCF Gradient Convergence Threshold:", scfGradientConv);
    }

  }; // MCSCFSettings::print

  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::printMCSCFHeader(EMPerturbation & pert) {
    
    auto & mopart = this->MOPartition;
    auto & ref = this->reference();

    std::cout << std::endl;
    std::cout << bannerTop << std::endl;
    std::cout << "MCSCF Settings:" << std::boolalpha << std::endl << std::endl;
     
    std::cout << std::left << std::setprecision(3) << std::scientific;
    
    std::string job_title = "* Job:  " + std::to_string(ref.nC) + "C-";
    
    if (mopart.scheme == CAS) {
      job_title += "CAS(" + std::to_string(mopart.nCorrE) + "," + std::to_string(mopart.nCorrO)
                  +")-";
    } else if (mopart.scheme == RAS) {
      job_title += "RAS(" + std::to_string(mopart.nCorrE) + "," + std::to_string(mopart.nCorrO)+")-";
    } else {
      CErr("Other than CAS has not been implemented yet "); 
    }
    if (this->settings.doSCF) {
      job_title += "SCF";
    } else {
      job_title += "CI";
    }
    FormattedLine(std::cout,job_title);
    
    std::cout << std::endl;
    if (ref.nC == 4)  {
      FormattedLine(std::cout,"* No-pair Approxmiation:", this->FourCompNoPair);
    }
    
    std::cout << std::endl;
    this->printMOSpacePatition();

    settings.print(ref.nC == 4, this->NStates);
 
    std::cout << std::endl;
    if(this->settings.doSCF and this->StateAverage) {
      FormattedLine(std::cout,"  State Average is ON, with weights:");
      auto & weights = this->SAWeight;
      for (auto i = 0ul; i < this->NStates; i ++)
        std::cout << "        State " << std::setw(5) << std::right << i << ":" 
                  << std::setw(16) <<  weights[i] << std::endl;
      
      std::cout << std::left << std::endl;
    }
   
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
    std::cout << std::endl << bannerTop << std::endl << std::endl;
  
  }; //MCSCF::printMCSCFHeader

  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::printStateEnergy() {
    
    std::cout << std::left << std::endl;
    
    FormattedLine(std::cout, "Energy at this Cycle:");
    std::cout << std::right << std::setprecision(10) << std::scientific;
    
    for (auto i = 0ul; i < this->NStates; i++)
      std::cout << "      State " << std::setw(5) << i + 1 << ":" 
                << std::setw(20) << this->StateEnergy[i] << std::endl;;
  
    std::cout << std::left << std::endl;
  } // printStateEnergy 
  
  
  template <typename MatsT>
  void printMCSCFState(std::ostream &out, size_t i, double energy, 
    MatsT* C, std::vector<size_t> & sorted_CAddr, 
    size_t N, const size_t n_item_per_row = 5) {
    
    out << std::fixed << std::right<< std::setprecision(10);
    out.fill(' ');
    
    out << std::endl <<  "State:" << std::setw(4) << i + 1 << "  Energy (Hartree):" 
    << std::setw(16) << energy <<  std::endl;
    
    out << std::fixed << std::right<< std::setprecision(7);
    
    size_t C_length = 10;
    size_t CAddr;
    
    for (auto j = 0ul; j < (N-1)/n_item_per_row + 1; j++) {
      size_t l = j*n_item_per_row;
      size_t r = std::min((j+1)*n_item_per_row, N);
      
      for (auto p = l; p < r; p++) {
        CAddr = sorted_CAddr[p]; 
        out << "(" << std::setw(5) << CAddr << ") " 
            << std::setw(C_length) << std::real(C[CAddr]); 
        if(std::is_same<MatsT, dcomplex>::value)  
          out << " " << std::setw(C_length) << std::imag(C[CAddr]);
        out << "  ";  
      }
      out << std::endl;
    }
  
  }; // printMCSCFState
  
  
  template <typename MatsT, typename IntsT>
  void MCSCF<MatsT,IntsT>::printMCSCFFooter( ) {
    
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
    std::vector<size_t> kLargestCAddr = std::vector<size_t>(NPrintC, 0);
    std::vector<MatsT> kLargestC = std::vector<MatsT>(NPrintC, MatsT(0.));

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
       printMCSCFState(std::cout, i, this->StateEnergy[i], C, Cindx, NPrintC);  
    }

    this->print1RDMs();
    
    std::cout << BannerTop << std::endl;
  
  }; //MCSCF::printMCSCFFooter
 
}; // namespace ChronusQ


