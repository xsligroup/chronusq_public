/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2018 Li Research Group (University of Washington)
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

#ifndef DETFACTORY_HPP
#error This file may only be included from detfactory.hpp
#endif

#define DETFACTORY_USE_2E_PERMUTATIONAL_SYMMETRY

namespace ChronusQ {

namespace {
inline std::string generateFullCD1eExListId(
  const DeterminantGroup& tGroup, size_t tSpace,
  const DeterminantGroup& uGroup, size_t uSpace) {

  std::string id = "1e-("
                   + std::to_string(tGroup.nElectrons()) + "e,"
                   + std::to_string(tGroup.nOrbitals()) + "o)";
  if (tSpace != uSpace) {
    id += "-(" + std::to_string(uGroup.nElectrons()) + "e,"
          + std::to_string(uGroup.nOrbitals()) + "o)";
    if (tSpace > uSpace) id += "-R";
  } 

  return id; 
} // generateFullCD1eExListId
} // namespace

/*
 * Figure out whether it is intra- or inter-space full configuration-driven excitation list
 * Call the correct function to generate the excitation list
 */
inline std::shared_ptr<NewExcitationList>
DeterminantFactory::constructFullCD1eExList(
  const DeterminantGroup& tGroup, size_t tSpace,
  const DeterminantGroup& uGroup, size_t uSpace) {

  // Generate a unique ID for the excitation list
  // The ID of the excitation list is only dependent on the number of electrons and orbitals
  std::string id = generateFullCD1eExListId(tGroup, tSpace, uGroup, uSpace);
  std::shared_ptr<NewExcitationList> exList;
  
  if (exLists_.count(id) != 0) {
    // If an excitation list already exits for the same spaces
    exList = exLists_.at(id);
  } else {
    // Create a new excitation list
    try {
      if (tSpace == uSpace)  {
        if (tGroup.nElectrons() == 0) return nullptr;
        exList = constructIntraSpaceFullCD1eExList(tGroup);
      } else {
        exList = constructInterSpaceFullCD1eExList(tGroup, uGroup, tSpace > uSpace);
      }
      exLists_.emplace(id, exList);
    } catch (const std::invalid_argument& e) {
      exList = nullptr;
    }
  }
  
  return exList;
} // DetFactory::constructFullCD1eExList

/*
 * Create an excitation list for <K|a^\dagger_t a_u|L>, defined in terms of excitations between spaces
 * for each Bra and Ket category pair.
 * For one-electron excitation, there can only be one space from each category involved.
 */
inline void DeterminantFactory::constructOneEExcitation(
  const std::shared_ptr<const DeterminantCategory>& braDetsCat,
  size_t braCategoricalIndex, size_t tSpace, size_t uSpace, 
  size_t ketCategoricalIndex, std::string twoE_RITerms = "") {

  // Generate a unique Id based on the space Id for the one-electron excitation term.
  // Note the oneETerm Id is different from that for the excitation list which is only dependent
  // on the number of electrons and orbitals.
  std::string oneETerm = "h1e(" + std::to_string(tSpace) + ","
      + std::to_string(uSpace) + ")" + twoE_RITerms;
  
  // std::cout << "start to build 1e terms: " << oneETerm << std::endl;
  const auto& tGroup = braDetsCat->detGroups()[tSpace];
  const auto& uGroup = braDetsCat->detGroups()[uSpace];

  // Currently, we can only construct the full configuration-driven (FullCD) excitation list.
  auto exList = constructFullCD1eExList(tGroup, tSpace, uGroup, uSpace);
  if (exList == nullptr) {
    // std::cout << "skip 1e term: " << oneETerm << std::endl;
    return;
  }
  
  oneEExTerms_.emplace(oneETerm);
  
  DetsCat1eExcitation excitation({oneETerm});
  excitation.exLists[0] = exList;
  excitation.categoricalIndices = {braCategoricalIndex, ketCategoricalIndex};
  excitation.exSpaces = {tSpace, uSpace};
  excitation.symmetryFactor = braDetsCat->signOffsetBetween(tSpace, uSpace);

  oneEExcitations_.emplace_back(std::move(excitation));
  // std::cout << "build 1e terms Done" << std::endl;
  return;
} // constructOneEExcitation 

inline std::string DeterminantFactory::constructTwoEExcitationWithExRI(
  const std::shared_ptr<const DeterminantCategory>& braDetsCat,
  size_t braCategoricalIndex,
  size_t tSpace, size_t uSpace, size_t wSpace, size_t vSpace, 
  const std::shared_ptr<const DeterminantCategory>& ketDetsCat,
  size_t ketCategoricalIndex) {
  
  std::string twoETerm = "g2e(" + std::to_string(tSpace) 
      + "," + std::to_string(uSpace) + "," + std::to_string(wSpace) 
      + "," + std::to_string(vSpace) + ")";  
  
  // find out aux groups, t u from bra side, w v from ket side
  // These are local categories that include active spaces of interest with respect
  // to the excitation list.
  DeterminantGroup tAuxGroup(braDetsCat->detGroups()[tSpace]);
  DeterminantGroup uAuxGroup(braDetsCat->detGroups()[uSpace]);
  DeterminantGroup wAuxGroup(ketDetsCat->detGroups()[wSpace]);
  DeterminantGroup vAuxGroup(ketDetsCat->detGroups()[vSpace]);
  
  try {
    if (tSpace != uSpace) {
      // We need to make sure that w v excitation is already part of the Bra.
      // Therefore, we start from the Bra side of the category.
      // Since u and t are reversed for easy access to resolution of idenity.
      tAuxGroup.removeNElectrons(1);
      uAuxGroup.addNElectrons(1);
    }
    if (wSpace != vSpace) {
      // We need to make sure that u t are not mixed. Therefore, we need to
      // prepare the first excitation from the Ket side.
      wAuxGroup.addNElectrons(1);
      vAuxGroup.removeNElectrons(1);
    }
  } catch (const std::invalid_argument& e) {
    return ""; 
  }

  double symmFac = 1.;
#ifdef DETFACTORY_USE_2E_PERMUTATIONAL_SYMMETRY
  if (tSpace != wSpace or uSpace != vSpace) {
    symmFac = 2.;
  }
  if (tSpace != wSpace and uSpace != vSpace) {
    twoETerm += "-X"; // add exchange terms
  }
#endif 

  auto exList_ut = constructFullCD1eExList(uAuxGroup, uSpace, tAuxGroup, tSpace);
  auto exList_wv = constructFullCD1eExList(wAuxGroup, wSpace, vAuxGroup, vSpace);
  
  if (exList_ut == nullptr or exList_wv == nullptr) {
    return "";
  }
  
  twoEExTerms_.emplace(twoETerm);

  DetsCat2eRIExcitation excitation({twoETerm});
  excitation.exLists[0] = exList_ut;
  excitation.exLists[1] = exList_wv;
  excitation.categoricalIndices = {braCategoricalIndex, ketCategoricalIndex};
  excitation.exSpaces = {tSpace, uSpace, wSpace, vSpace};
  excitation.symmetryFactor = symmFac 
      * braDetsCat->signOffsetBetween(tSpace, uSpace) 
      * ketDetsCat->signOffsetBetween(wSpace, vSpace); 
  
  twoEExcitations_.emplace_back(std::move(excitation));
   
  return twoETerm;
} // DetFactory::constructTwoEExcitationWithExRI

inline void DeterminantFactory::generateComputingGraph(bool usingExRI) {

  if(not usingExRI) CErr("Two electron interaction not usingExRI is NYI");
  
  std::vector<size_t> ketSpaceOcc, braSpaceOcc;
  std::vector<size_t> iBraSpace, iKetSpace; // space indices that are associated with non-zero excitation

  // main loop
  for (auto i = 0ul; i < braCategoricalSpace_->nCategories(); ++i) {
    const auto& braCat = braCategoricalSpace_->getCategory(i);
      braSpaceOcc = braCat->SpaceOccupations();
    
    // MPI Parallelism 
    if (not braCategoricalSpace_->containsLocalCategory(i)) continue;

    for (auto j = 0ul; j < ketCategoricalSpace_->nCategories(); ++j) {
      const auto& ketCat = ketCategoricalSpace_->getCategory(j);
        ketSpaceOcc = ketCat->SpaceOccupations();
      
      // find excitations spaces 
      iBraSpace.clear();
      iKetSpace.clear();
      size_t nExcitations = 0ul;

      // compare occupation numbers between Bra and Ket spaces in different categories
      // use the total number of occupation difference to determine if there are non-zero
      // excitations between two categories. If so, figure out if they are one-electron
      // or two electron excitations.
      for (auto k = 0ul; k < ketSpaceOcc.size(); ++k) {
        size_t kSpaceExcitation = std::abs(int(braSpaceOcc[k]) - int(ketSpaceOcc[k]));

        // total number of excitation includes both exciting into and out of the space
        nExcitations += kSpaceExcitation;
        if (nExcitations > 4ul) break;
        
        // Add space indices that are associated with non-zero excitations
        // The space is saved once for every one-electron excitation, i.e., if there
        //  are two electrons excited, the space is saved twice.
        // Figure out creation space (Bra) and annihilation space
        if (braSpaceOcc[k] > ketSpaceOcc[k]) {
          // this is in the creation Bra space
          for (auto l = 0ul; l < kSpaceExcitation; ++l) {
            iBraSpace.push_back(k);
          }
        } else if (braSpaceOcc[k] < ketSpaceOcc[k]){
          // this is in the annihilation Ket space
          for (auto l = 0ul; l < kSpaceExcitation; ++l) {
            iKetSpace.push_back(k);
          }
        }
      }

      // total number of excitation includes both exciting into and out of the space
      if (nExcitations > 4ul) continue;

      // since the total number of excitation includes both exciting into and
      // out of the space, the effective excitation level should be divided by 2.
      nExcitations /= 2;
      
      if (nExcitations == 2ul) {
        // two-electron excitations <t(1)w(2)|u(1)v(2)> or (t(1)u(1)|w(2)v(2))
        // there could be several different arrangements among Bra and Ket spaces, respectively
        size_t tSpace = iBraSpace[0]; // creation
        size_t uSpace = iKetSpace[0]; // annihilation
        size_t wSpace = iBraSpace[1]; // creation
        size_t vSpace = iKetSpace[1]; // annihilation

        //XSLI, I do not know why the following statement has anything to do with permutation symmetry
        // I think they are just zero terms
#ifdef DETFACTORY_USE_2E_PERMUTATIONAL_SYMMETRY
        if ((tSpace == wSpace and braSpaceOcc[tSpace] < 2) or
            (braSpaceOcc[tSpace] < 1 or braSpaceOcc[wSpace] < 1) or
            (uSpace == vSpace and ketSpaceOcc[uSpace] < 2) or
            (ketSpaceOcc[uSpace] < 1 or ketSpaceOcc[vSpace] < 1)) continue;
#endif
         constructTwoEExcitationWithExRI(braCat, i, tSpace, uSpace, wSpace, vSpace, ketCat, j);
#ifndef DETFACTORY_USE_2E_PERMUTATIONAL_SYMMETRY
         if (tSpace != wSpace)  {
           constructTwoEExcitationWithExRI(braCat, i, wSpace, uSpace, tSpace, vSpace, ketCat, j);
         }
         if (uSpace != vSpace) {
           constructTwoEExcitationWithExRI(braCat, i, tSpace, vSpace, wSpace, uSpace, ketCat, j);
         }
         if (tSpace != wSpace and uSpace != vSpace) {
           constructTwoEExcitationWithExRI(braCat, i, wSpace, vSpace, tSpace, uSpace, ketCat, j);
         }
#endif  
      } else if (nExcitations == 1ul) {
        size_t tSpace = iBraSpace[0]; // creation
        size_t uSpace = iKetSpace[0]; // annihilation
#ifndef DETFACTORY_USE_2E_PERMUTATIONAL_SYMMETRY
        std::string twoE_RITerms = "";
        for (auto wSpace = 0ul; wSpace < nEs_ket.size(); ++wSpace) {
          constructTwoEExcitationWithExRI(braCat, i, tSpace, wSpace, wSpace, uSpace, ketCat, j);
          if (tSpace != wSpace) {
            constructTwoEExcitationWithExRI(braCat, i, wSpace, wSpace, tSpace, uSpace, ketCat, j);
          }
          if (uSpace != wSpace) {
            constructTwoEExcitationWithExRI(braCat, i, tSpace, uSpace, wSpace, wSpace, ketCat, j);
          }
          if (tSpace != wSpace and uSpace != wSpace) {
            constructTwoEExcitationWithExRI(braCat, i, wSpace, uSpace, tSpace, wSpace, ketCat, j);
          }
          twoE_RITerms += "+RI[g2e(" + std::to_string(tSpace) + "," 
              + std::to_string(wSpace) + "," + std::to_string(wSpace) + "," 
              + std::to_string(uSpace) + ")]";  
        }
        constructOneEExcitation(braCat, i, tSpace, uSpace, ketCat, j, twoE_RITerms);
#else          
        if (braSpaceOcc[tSpace] < 1 or ketSpaceOcc[uSpace] < 1) continue;

        for (auto wSpace = 0ul; wSpace < ketSpaceOcc.size(); ++wSpace) {
          // <K|a^\dagger_w a^\dagger_t w u|L>, defined in terms of excitations between spaces
          constructTwoEExcitationWithExRI(braCat, i, wSpace, uSpace, tSpace, wSpace, ketCat, j);
        }
        // <K|a^\dagger_t u|L>, defined in terms of excitations between spaces
        constructOneEExcitation(braCat, i, tSpace, uSpace, j);
#endif        
      } else { // == 0ul same space
        for (auto tSpace = 0ul; tSpace < ketSpaceOcc.size(); ++tSpace) {
#ifndef DETFACTORY_USE_2E_PERMUTATIONAL_SYMMETRY
          std::string twoE_RITerms = "";
          for (auto uSpace = 0ul; uSpace < nEs_ket.size(); ++uSpace) {
            if ( tSpace <= uSpace) {
              constructTwoEExcitationWithExRI(braCat, i, tSpace, uSpace, uSpace, tSpace, ketCat, j);
              if ( tSpace != uSpace) {
                constructTwoEExcitationWithExRI(braCat, i, uSpace, tSpace, tSpace, uSpace, ketCat, j);
                constructTwoEExcitationWithExRI(braCat, i, tSpace, tSpace, uSpace, uSpace, ketCat, j);
                constructTwoEExcitationWithExRI(braCat, i, uSpace, uSpace, tSpace, tSpace, ketCat, j);
              }
            }
            twoE_RITerms += "+RI[g2e(" + std::to_string(tSpace) + "," 
                + std::to_string(uSpace) + "," + std::to_string(uSpace) + "," 
                + std::to_string(tSpace) + ")]";  
          }
          constructOneEExcitation(braCat, i, tSpace, tSpace, j, twoE_RITerms);
#else
          if (braSpaceOcc[tSpace] < 1) continue;
          // XSLI: do we need to check if we have >1 electrons for two-e excitation?
          for (auto uSpace = tSpace; uSpace < ketSpaceOcc.size(); ++uSpace) {
            // <K|a^\dagger_t a^\dagger_u a_t a_u|L>, defined in terms of excitations between spaces
            constructTwoEExcitationWithExRI(braCat, i, tSpace, tSpace, uSpace, uSpace, ketCat, j);
          }
          std::string twoE_RITerms = "+RI[g2e(" + std::to_string(tSpace)
              + "," + std::to_string(tSpace) + "," + std::to_string(tSpace) 
              + "," + std::to_string(tSpace) + ")]";
          // <K|a^\dagger_t a_t|L>, defined in terms of excitations between spaces
          constructOneEExcitation(braCat, i, tSpace, tSpace, j, twoE_RITerms);
#endif
        }
      }
    } // braCat
  } // cat_ket
  
#ifdef CQ_ENABLE_MPI
  //taskScheduler_.generateTwoEExTaskMap(twoEExcitations_, activeSpaces_);
#endif

} // DetFactory::generateComputingGraph

inline void DeterminantFactory::output(std::ostream & out, const std::string & s) const {
   
  std::string outputStr;
    
  if (s == "")
    outputStr = "Determinant Factory";
  else
    outputStr = "Determinant Factory Of " + s;

  outputStr += " on node " + std::to_string(MPIRank(this->comm_));
  
  out << "---------------------" << std::endl;
  out << outputStr << std::endl;
  out << std::endl;
  
  out << "* Excitation Lists Needed: " << std::endl;
  
  double totalExListStorage; 
  size_t count = 0;
  for (auto & l: exLists_) {
    out << " # " << std::setw(5) << count << ", "; 
    double lstorage = l.second->storageSize() / 1e9;
    totalExListStorage += lstorage;
    out << std::setw(30) << l.first << ", needs storage: " << lstorage 
        << " GB" << std::endl;
    count++;
  }
  out << "----Total Storage Needed for Excitation Lists: " 
      << totalExListStorage << " GB" << std::endl;
  out << std::endl;
  
  auto printDetsCatExcitation = 
      [&] (const std::string& term, 
           const std::pair<size_t, size_t>& categoricalIndices, 
           double symmetryFactor) {
        out << "    Term: " << term << std::endl;
        out << "    Excitation Categories: " << std::setw(5) << categoricalIndices.first
            << " <- " << std::setw(5) << categoricalIndices.second << std::endl;
        out << "    Excitation Occupation Info: (";
        for (const auto& i : braCategoricalSpace()->getCategory(categoricalIndices.first)->SpaceOccupations()) {
          out << std::setw(3) << i;
        }
        out << ") <- (";
        for (const auto& i : ketCategoricalSpace()->getCategory(categoricalIndices.second)->SpaceOccupations()) {
          out << std::setw(3) << i;
        }
        out << ")" << std::endl;
        out << "    Symmetry Factor: " << symmetryFactor << std::endl;
        out << std::endl;
      };
  
  out << "* One-Body Interactions:" << std::endl;

#if 0
  count = 0ul;
  for (const auto& oneEEx : oneEExcitations_) {
    out << " # " << count << ":" << std::endl; 
    printDetsCatExcitation(oneEEx.term, oneEEx.categoricalIndices, oneEEx.symmetryFactor);
    count++;
  }
  
  out << std::endl;
  out << "* Two-Body Interactions:" << std::endl;
  count = 0;
  for (const auto& twoEEx : twoEExcitations_) {
    out << " # " << count << ":" << std::endl; 
    printDetsCatExcitation(twoEEx.term, twoEEx.categoricalIndices, twoEEx.symmetryFactor);
    count++;
  }
#endif

#if 0
  out << std::endl;
  out << "* Two-Body Interactions Summary By Terms: " << std::endl;
  
  std::unordered_map<std::string, std::pair<size_t, double>> twoETermSummary;
  for (const auto& twoEEx : twoEExcitations_) {
    if (twoETermSummary.count(twoEEx.term) == 0) {
      twoETermSummary.emplace(twoEEx.term, std::pair<size_t, double>({0ul, 0.}));
    }
    twoETermSummary[twoEEx.term].first++;
    twoETermSummary[twoEEx.term].second += twoEEx.estimatedComputationalCost();
  }
  
  const auto& activeSpaces = ketCategoricalSpace()->activeSpaces();
  std::vector<size_t> span;
  std::vector<std::string> termVec;
  double totalTwoETermStorage = 0.;
  for (const auto& [term, summary] : twoETermSummary) {
    parseTermSpan(term, termVec, span);
    double nt = activeSpaces[span[0]].nOrbitals; 
    double nu = activeSpaces[span[1]].nOrbitals; 
    double nw = activeSpaces[span[2]].nOrbitals; 
    double nv = activeSpaces[span[3]].nOrbitals; 
    double storageNeeded = nt * nw * nu * nv / 1e9; 

    out << "  - term " << std::left << std::setw(20) << term 
        << ", storage: sizeof(MatsT) x " << std::right << std::setw(10) << storageNeeded << " GB"
        << ", # of tasks:" << std::setw(4) << summary.first
        << ", estimated load: " << summary.second 
        << std::endl;
    totalTwoETermStorage += storageNeeded;
  }
  size_t nTOrb = 0ul;
  for (const auto& s: activeSpaces) {
    nTOrb += s.nOrbitals;
  }

  out << "---- Total number of terms: " <<  twoETermSummary.size() 
      << ", total storage: sizeof(MatsT) x" << totalTwoETermStorage << " GB" 
      << " (compared to " << std::pow(nTOrb, 4) / 1e9 
      << " GB in full space)" << std::endl;
  
#ifdef CQ_ENABLE_MPI
  out << std::endl;
  out << "* MPI Two Body Interactions Schedule Summary " << std::endl; 
  std::vector<double> estimatedLoad(MPISize(comm_), 0.);
  std::vector<std::vector<size_t>> taskPools(MPISize(comm_));
  
  for (auto i = 0ul; i < MPISize(comm_); ++i) {
    double estimatedLoad = 0;
    std::unordered_set<std::string> terms;
    size_t numberOfTask = 0ul;
    for (auto j = 0ul; j < twoEExcitations_.size(); ++j) {
      if (twoEExTaskMap[j] != i) continue;
      ++numberOfTask;
      estimatedLoad += twoEExcitations_[j].estimatedComputationalCost();
      terms.emplace(twoEExcitations_[j].term);
    }
    out << " # On Node " << i << ":" << std::endl;
    out << "   - Number of Task: " << numberOfTask << std::endl;
    out << "   - Total load: " << estimatedLoad << std::endl;
    out << "   - Terms needed: " << std::endl;
    double localTwoETermStorage = 0.;
    for (const auto& t : terms) {
      out << "     $ " << t << std::endl;
      parseTermSpan(t, termVec, span);
      double nt = activeSpaces[span[0]].nOrbitals; 
      double nu = activeSpaces[span[1]].nOrbitals; 
      double nw = activeSpaces[span[2]].nOrbitals; 
      double nv = activeSpaces[span[3]].nOrbitals; 
      localTwoETermStorage += nt * nw * nu * nv / 1e9; 
    }
    out << "   - Local two E term storage: sizeof(MatsT) x" << localTwoETermStorage << " GB" << std::endl; 
    out << std::endl;
  }

  out << "-------------------------" << std::endl;
#endif
#endif
} // DetFactory::output

} // namespace ChronusQ
