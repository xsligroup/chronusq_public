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

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::printMOSpacePartition() {
      
    auto & mopart = this->MOPartition;
    auto & ref = this->reference();

    std::string particletype = this->reference().particle.charge < 0 ? "Electron" : "Proton";
  
    std::cout << std::left << std::endl;
    FormattedLine(std::cout,"* MO Space Partition:");
    FormattedLine(std::cout,"  Number of " + particletype + "ic MOs:", mopart.nElecMO);
    if (ref.nC == 4)
      FormattedLine(std::cout,"  Number of Negative MOs:", mopart.nNegMO);
    FormattedLine(std::cout,"  Number of " + particletype + "s in Full Space:",  ref.nO);
    
    std::cout << std::endl;
    FormattedLine(std::cout,"  Number of Inactive Orbitals:",     mopart.nInact);
    FormattedLine(std::cout,"  Number of Frozen Virtual Orbitals:",  mopart.nFVirt);
    
    FormattedLine(std::cout,"  Number of Correlated Orbitals:",      mopart.nCorrO);
    if (mopart.scheme == RAS) {
      FormattedLine(std::cout,"  Number of Orbitals in RAS 1:",      mopart.nActOs[0]);
      FormattedLine(std::cout,"  Number of Orbitals in RAS 2:",      mopart.nActOs[1]);
      FormattedLine(std::cout,"  Number of Orbitals in RAS 3:",      mopart.nActOs[2]);
      FormattedLine(std::cout,"  Maximum Number of Holes in RAS 1:",      mopart.mxHole);
      FormattedLine(std::cout,"  Maximum Number of Electrons in RAS 3:",  mopart.mxElec);
    }

    std::cout << std::endl;
    FormattedLine(std::cout,"  Number of Correlated " + particletype + "s:",     mopart.nCorrE);
    if (ref.nC == 1) {
      FormattedLine(std::cout,"  Number of Correlated Alpha " + particletype + "s:", mopart.nCorrEA);
      FormattedLine(std::cout,"  Number of Correlated Beta " + particletype + "s:",  mopart.nCorrEB);
    }
    
    std::cout << std::endl;
    FormattedLine(std::cout,"  Number of Roots Requested:",this->NStates);
    FormattedLine(std::cout,"  Number of Determinants:",   this->NDet);
    
    std::cout << std::endl;
    std::cout << std::right << std::endl;
    this->mointsTF->printMORangesSummary();
    std::cout << std::left  << std::endl;
  
  }; // MCWaveFunction::printMOSpacePartition
  
  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::print1RDMs() {

    if( printRDMs==0 ) return;

    std::cout << std::endl << "MCWaveFunction 1RDMs:" << std::endl;
    std::cout << bannerTop << std::endl;

    size_t nCorrO = this->MOPartition.nCorrO;
    size_t nS     = this->NStates;
    size_t nRDMEle= 0; // used for line breaker

    //Print full 1RDM
    if( printRDMs==1 ){
      std::cout << "Printing full 1-RDM" << std::endl;
      for (auto i = 0ul; i < nS; i++) {

        prettyPrintSmart(std::cout, "State " + std::to_string(i),
          this->oneRDM[i].pointer(), nCorrO, nCorrO, nCorrO);

        //if (this->SpinAnalysis)
          this->spinAnalysis(i);
      }
    } else if( printRDMs==2 ){
      std::cout.precision(2);
      std::cout << "Printing large (>" << std::fixed << rdmCut << ") diagonal elements of real 1-RDM" << std::endl;

      for (auto i = 0ul; i < nS; i++) {

        std::cout << std::endl << "State " << i+1 << ": ";
        nRDMEle = 0;

        for( auto ipp = 0ul; ipp < nCorrO; ipp++){

          if( std::real(this->oneRDM[i].pointer()[ipp*nCorrO+ipp]) > rdmCut ){
            if( nRDMEle > 6 ){
              std::cout << std::endl << "         ";
              nRDMEle = 0;
            }
            std::cout << std::setw(3) << ipp+1 << "(" << std::real(this->oneRDM[i].pointer()[ipp*nCorrO+ipp]) << ")" << " ";
            nRDMEle++;
          }
        }

        std::cout << std::endl;
        //if (this->SpinAnalysis)
          this->spinAnalysis(i);
      }
    }

    std::cout << bannerTop << std::endl;

  }; //MCWaveFunction::print1RDMs

  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::printMOInfo(std::ostream& out, 
                      size_t printMOLevel) {
    
    SingleSlater<MatsT,IntsT> * ss_ptr = &reference(); 

    std::fill_n(ss_ptr->eps1, ss_ptr->nC * ss_ptr->nAlphaOrbital(),double(0.));
    if( ss_ptr->eps2 != nullptr )
      std::fill_n(ss_ptr->eps2, ss_ptr->nC * ss_ptr->nAlphaOrbital(),double(0.));

    if (not printMOLevel)
      printMOLevel = printMOCoeffs;
    ss_ptr->WaveFunction<MatsT,IntsT>::printMOInfo(out, printMOLevel);

  } //MCWaveFunction::printMOInfo

  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::printMCStates()
  {
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
       this->printMCState(std::cout, i, this->StateEnergy->at(i), C, Cindx, NPrintC);  
    }
  }

  template <typename MatsT, typename IntsT>
  void MCWaveFunction<MatsT,IntsT>::printMCState(std::ostream & out, 
                                                 size_t i,
                                                 double energy,
                                                 MatsT * C,
                                                 std::vector<size_t> & sorted_CAddr,
                                                 size_t N,
                                                 const size_t n_item_per_row)
  {
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
 
  } // MCWaveFunction::printMCWaveFunctionState



}; // namespace ChronusQ
