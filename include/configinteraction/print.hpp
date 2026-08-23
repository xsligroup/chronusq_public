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

#include <configinteraction.hpp>
#include <util/matout.hpp>
#include <util/print.hpp>
#include <cxxapi/output.hpp>
#include <mointstransformer/moranges.hpp>

namespace ChronusQ {

void CISettings::print(bool fourComp) {
  
  std::cout << std::endl;
  FormattedLine(std::cout,"* Computation Parameters:");
  if(ciAlg == CIDiagonalizationAlgorithm::CI_FULL_MATRIX) {
    FormattedLine(std::cout,"  CI Algorithm:",  "Full Matrix");
  } else if (ciAlg == CIDiagonalizationAlgorithm::CI_DAVIDSON) {
    FormattedLine(std::cout,"  CI Algorithm:",  "Davidson");
    
    std::string ciSigma2eContAlgPrintStr;
    if (ciSigma2eContAlg == "NAIVE" or ciSigma2eContAlg == "NAIVELOOP"
        or ciSigma2eContAlg == "NL") {
      ciSigma2eContAlgPrintStr = "Naive Loop (Default)"; 
    } else if (ciSigma2eContAlg == "KNOWLESHANDY" or ciSigma2eContAlg == "KH" ) {
      ciSigma2eContAlgPrintStr = "Knowles-Handy"; 
    } else if (ciSigma2eContAlg == "OLSENROOS" or ciSigma2eContAlg == "OR") {
      ciSigma2eContAlgPrintStr = "Olsen-Roos"; 
    } else if (ciSigma2eContAlg == "FRISCHLI" or ciSigma2eContAlg == "FL") {
      ciSigma2eContAlgPrintStr = "Frisch-Li"; 
    } else if (ciSigma2eContAlg == "SMALLBLOCKH" or ciSigma2eContAlg == "SBH" 
        or ciSigma2eContAlg == "SMALLBLOCKHAMILTONIAN") {
      ciSigma2eContAlgPrintStr = "Small Block Hamiltonian"; 
    } else if (ciSigma2eContAlg == "LARGEBLOCKH" or ciSigma2eContAlg == "LBH" 
        or ciSigma2eContAlg == "LARGEBLOCKHAMILTONIAN") {
      ciSigma2eContAlgPrintStr = "Large Block Hamiltonian"; 
    } else {
      CErr("Unknown Sigma 2e Contraction Algorithm in DASCI");
    }
    FormattedLine(std::cout,"  CI Sigma 2e Contraction Algorithm:",  ciSigma2eContAlgPrintStr);
    
    FormattedLine(std::cout,"  CI Maximum Number of Iteration:",  maxCIIter);
    FormattedLine(std::cout,"  CI Vector Convergence Threshold:", ciVectorConv);
    FormattedLine(std::cout,"  Max Len of Davidson Subspace (x NRoots):", maxDavidsonSpace);
    FormattedLine(std::cout,"  Number of Davidson Guess(x NRoots):", nDavidsonGuess);
    FormattedLine(std::cout,"  Sparse Implementation of STP-DAS:", SparseDavidson ? "True" : "False");
    if(SparseDavidson)
      FormattedLine(std::cout,"  Sparse Davidson Screening Threshold:", SparseDavidsonEps);

  } else CErr("NYI CI Algorithm");
  
  if(this->doSCF) {
    std::cout << std::endl;
    ORSettings.print(fourComp);
    
    std::cout << std::endl;
    FormattedLine(std::cout,"  SCF Maximum Number of Iteration:",    maxSCFIter);
    FormattedLine(std::cout,"  SCF Energy Convergence Threshold:",   scfEnergyConv);
    FormattedLine(std::cout,"  SCF Gradient Convergence Threshold:", scfGradientConv);
  }

} // CISettings::print

template <typename MatsT, typename IntsT>
void ConfigurationInteraction<MatsT, IntsT>::printCIHeader() {
  
  auto & corrS = this->corrSpace;
  auto & ref = *this->reference();

  std::cout << std::endl;
  std::cout << bannerTop << std::endl;
  std::cout << "Configuration Interaction Settings:" << std::boolalpha << std::endl << std::endl;
   
  std::cout << std::left << std::setprecision(3) << std::scientific;
  
  std::string job_title = "* Job:  " + std::to_string(ref.nC) + "C-";
  
  // if (corrS.scheme == CAS) {
  //   job_title += "CAS(" + std::to_string(corrS.nCorrE) + "," + std::to_string(corrS.nCorrO)
  //               +")-";
  // } else if (corrS.scheme == RAS) {
  //   job_title += "RAS(" + std::to_string(corrS.nCorrE) + "," + std::to_string(corrS.nCorrO)+")-";
  // } else {
  //   CErr("Other than CAS has not been implemented yet "); 
  // }
  
  if (this->ciSettings.doSCF) {
    job_title += "SCF";
  } else {
    job_title += "CI";
  }
  
  FormattedLine(std::cout,job_title);
  
  std::cout << std::endl;
  if (ref.nC == 4)  {
    FormattedLine(std::cout, "* No-pair Approxmiation:", this->FourCompNoPair);
  }
  
  //std::cout << std::endl;
  //this->printCorrMOSpace();
  //FormattedLine(std::cout, "* Active Spaces:");
  //std::cout << ciSettings.activeSpaces << std::endl;
  
  if (ciSettings.refOcc.size() > 0) {
    //std::cout << "  * Reference Category Occupancy (From input) in Each Space:" << std::endl;
    //size_t counter = 1ul;
    //for (auto const & occs: ciSettings.refOcc) {
    //  std::cout << "    - Reference " << std::setw(3) << counter << ": ";
    //  for (auto const & occ: occs) std::cout << std::setw(3) << occ << " ";
    //  std::cout << std::endl;
    //  counter++;
    //}
    std::cout << "  * Maximum Excitation from Reference: " 
      << (ciSettings.maxInterSpaceEx < 0 ? "All Excitation" : std::to_string(ciSettings.maxInterSpaceEx) )
      << std::endl;
  }
  
  detFactory->ketCategoricalSpace()->output(std::cout, "Categories Generated in Configuration Interaction"); 
  
  ciSettings.print(ref.nC == 4);
   if (!this->ciSettings.energyRefs.empty()) {
    FormattedLine(std::cout,"  Energy specific settings:");
    size_t nLowRoots = this->NStates;
    for (auto & pair: this->ciSettings.energyRefs) {
      FormattedLine(std::cout, "  Energy threshold:", pair.first, " #Roots:", pair.second);        
      nLowRoots -= pair.second;
    }
    FormattedLine(std::cout,"  Number of low energy roots:", nLowRoots);
  }
 
  std::cout << std::endl;
  if(this->ciSettings.doSCF and this->StateAverage) {
    FormattedLine(std::cout,"  State Average is ON, with weights:");
    auto & weights = this->SAWeight;
    for (auto i = 0ul; i < this->NStates; i ++)
      std::cout << "        State " << std::setw(5) << std::right << i << ":" 
                << std::setw(16) <<  weights[i] << std::endl;
    
    std::cout << std::left << std::endl;
  }
  
  std::cout << std::endl << bannerTop << std::endl << std::endl;

} //ConfigInteraction::printCIHeader

template <typename MatsT, typename IntsT>
void ConfigurationInteraction<MatsT, IntsT>::printStateEnergy() {
  
  std::cout << std::left << std::endl;
  
  FormattedLine(std::cout, "Energy at this Cycle:");
  std::cout << std::right << std::setprecision(10) << std::scientific;
  
  for (auto i = 0ul; i < this->NStates; i++)
    std::cout << "      State " << std::setw(5) << i + 1 << ":" 
              << std::setw(20) << this->StateEnergy[i] << std::endl;;

  std::cout << std::left << std::endl;
} // printStateEnergy 

namespace {
template <typename MatsT>
void printCIState(std::shared_ptr<DeterminantFactory> detFactory, std::ostream &out, size_t i, double energy, 
  const std::vector<size_t>& kLargestCAddr, 
  const std::vector<MatsT>& kLargestC,
  const size_t n_item_per_row = 5,
  const bool expandAddressToOcc = false) {
  
  out << std::fixed << std::right<< std::setprecision(10);
  out.fill(' ');
  
  out << std::endl <<  "State:" << std::setw(4) << i + 1 << "  Energy (Hartree):" 
  << std::setw(16) << energy <<  std::endl;
  
  out << std::fixed << std::right<< std::setprecision(7);
  
  size_t C_length = 10;
  size_t nPrint = kLargestCAddr.size();

  for (auto i = 0ul; i < (nPrint - 1) / n_item_per_row + 1; ++i) {
    size_t jBegin = i * n_item_per_row;
    size_t jEnd = std::min((i + 1) * n_item_per_row, nPrint);
    
    for (auto j = jBegin; j < jEnd; ++j) {
      out << "(";
      out << std::setw(5) << kLargestCAddr[j];
      if (expandAddressToOcc){
        out << " occ: ";
        const auto ketSpace = detFactory->ketCategoricalSpace();
        const auto ketCatIdx = ketSpace->getCategoryIdx(kLargestCAddr[j]);
        //const auto& ketCat = ketSpace->getCategory(ketCatIdx);
        const auto& ketCat = dynamic_cast<const FullDeterminantCategory&>(*ketSpace->getCategory(ketCatIdx));
        const auto addresser = ketCat.template addresser<uint64_t>();
        auto occInfo = addresser.addressToOccInfo(kLargestCAddr[j]);
        //std::cout << "CAT " << ketCatIdx << " " << ketCat.offset()  << "," << ketCat.offset() + ketCat.nDeterminants() << "   " << j << std::endl;
        auto count = 0;
        for (auto occOrb : occInfo){
          if (count != 0)
            out << ", ";
          out << occOrb + 1;
          count++;
        }
        //const auto detGen = ketCat->generator();
        //const auto addresser = ketCat->addresser();
      }
      out << ") " << std::setw(C_length) << std::real(kLargestC[j]); 
      if(std::is_same<MatsT, dcomplex>::value)  
        out << " " << std::setw(C_length) << std::imag(kLargestC[j]);
      out << "  ";  
    }
    out << std::endl;
  }

} // printCIState
} // namespace

template <typename MatsT, typename IntsT>
void ConfigurationInteraction<MatsT, IntsT>::printCIFooter( ) {
  
  std::cout << std::endl << "Configuration Interaction Results:" << std::endl;
  std::cout << BannerTop << std::endl;

  this->printMOInfo(std::cout);
  
  std::cout << " *---------------------------------------------*" << std::endl;      
  std::cout << " * Configuration Interaction (CI) Eigen States *" << std::endl;      
  std::cout << " *---------------------------------------------*" << std::endl;      
  
  std::cout << std::endl << BannerTop << std::endl;
  
  size_t nDet  = this->nDeterminants();
  size_t nS    = this->NStates;
  size_t minnDetPerNode = 25;
  const auto ketSpace = detFactory->ketCategoricalSpace();
  for (auto& catLen : ketSpace->distributedCategoryLengths()) {
    if (minnDetPerNode > catLen) {
      minnDetPerNode = catLen;
    }
  }
  size_t nPrintC = std::min(nDet, size_t(minnDetPerNode));
  
  std::vector<size_t> kLargestCAddr; 
  std::vector<MatsT> kLargestC;


  if(!ciSettings.SparseDavidson) {
    std::shared_ptr<DistributedVectors<MatsT>> CIVectorsCast = std::dynamic_pointer_cast<DistributedVectors<MatsT>>(CIVectors);
    for (auto i = 0ul; i < nS; i++) {

     kLargestCAddr.clear();
     kLargestC.clear();
     std::vector<std::pair<double, size_t>> CWindows;
     CWindows.emplace_back(0.0, nPrintC);
     
     CIVectorsCast->getKIndicesAndValues(nPrintC, i, kLargestCAddr, kLargestC, CWindows,
         [](const MatsT& a, const MatsT& b) { return std::norm(a) > std::norm(b); }
     );

     auto largestVal = CIVectorsCast->get(kLargestCAddr[0], i);
     //Rotate each state such that the largest coefficient is real and positive
     for(auto& coeff : kLargestC) { 
       coeff *= std::abs(largestVal) / largestVal;
     }

     // only print k Largest coefficient
     size_t nDetPerLine = this->printDetailedCICoeffs ? 1 : 5;
     printCIState(this->detFactory, std::cout, i, this->StateEnergy[i], kLargestCAddr, kLargestC, nDetPerLine, this->printDetailedCICoeffs);  

    }
  }
  else {
#ifdef CQ_ENABLE_SPARSE
    std::shared_ptr<DistributedSparseVectors<MatsT>> CIVectorsCast = std::dynamic_pointer_cast<DistributedSparseVectors<MatsT>>(CIVectors);
    for (auto i = 0ul; i < nS; i++) {

     kLargestCAddr.clear();
     kLargestC.clear();
     CIVectorsCast->getKIndicesAndValues(nPrintC, i, kLargestCAddr, kLargestC,
         [](const MatsT& a, const MatsT& b) { return std::norm(a) > std::norm(b); }
     );

     auto largestVal = CIVectorsCast->get(kLargestCAddr[0], i);
     //Rotate each state such that the largest coefficient is real and positive
     for(auto& coeff : kLargestC) {
       coeff *= std::abs(largestVal) / largestVal;
     }

     // only print k Largest coefficient
     size_t nDetPerLine = this->printDetailedCICoeffs ? 1 : 5;
     printCIState(this->detFactory, std::cout, i, this->StateEnergy[i], kLargestCAddr, kLargestC, nDetPerLine, this->printDetailedCICoeffs);  

    }
#else
    CErr("Please ENABLE_SPARSE in compilation");
#endif
  }


  this->print1RDMs();
  
  std::cout << BannerTop << std::endl;

} //ConfigInteraction::printCIFooter
 
} // namespace ChronusQ
