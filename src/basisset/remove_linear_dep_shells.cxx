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
#include <basisset/remove_linear_dep_shells.hpp>
#include <particleintegrals/onepints.hpp>
#include <physcon.hpp>
#include <optional>

//#define DEBUG_BASIS_LIN_DEP

namespace ChronusQ {

  void read_option_and_remove_linear_dependency(
      const CQInputFile &input, BasisSet &basis, Molecule &mol, std::ostream &output) {

    bool rm_lin_dep = false;
    if (input.containsData("SCF.REMOVELINEARDEP")) {
      rm_lin_dep = input.getData<bool>("SCF.REMOVELINEARDEP");
    }
    if ( rm_lin_dep ) {
      // First figure out if this is 4-component calculation
      std::string reference;
      try {
        reference = input.getData<std::string>("QM.REFERENCE");
      } catch(...) {
        CErr("QM.REFERENCE Keyword not found!", output);
      }
      // Digest reference string
      // Trim Spaces
      trim(reference);

      // Split into tokens
      std::vector<std::string> tokens;
      split(tokens,reference);
      for(auto &X : tokens) trim(X);
      bool is4C = tokens.back() == "4CHF";

      bool libcint = false;
      if (input.containsData("INTS.LIBCINT")) {
        libcint = input.getData<bool>("INTS.LIBCINT");
      }

      // Check for linear dependency threshold
      double linearDepTol = 1e-12;
      if ( input.containsData("SCF.LINEARDEPTOL")) {
        linearDepTol = input.getData<double>("SCF.LINEARDEPTOL");
      }
      // For NEO, it's possible to set different thresholds for electronic
      // and nuclear bases
      if ( basis.nucBasis && input.containsData("SCF.PROT_LINEARDEPTOL")){
        linearDepTol = input.getData<double>("SCF.PROT_LINEARDEPTOL");
      }

      bool atomicOnly = false;
      if ( input.containsData("SCF.REMOVEATOMICLINEARDEPONLY")) {
        atomicOnly = input.getData<bool>("SCF.REMOVEATOMICLINEARDEPONLY");
      }

      // Check for linear dependencies
      HamiltonianOptions hamiltonianOptions;
      hamiltonianOptions.Libcint = libcint;
      remove_linear_dependency<double>(basis, mol, hamiltonianOptions, linearDepTol, is4C, atomicOnly);
    }

  }

  template <typename IntsT>
  void remove_linear_dependency(
      BasisSet &basis, Molecule &mol, const HamiltonianOptions &hamiltonianOptions,
      double linearDependencyThreshold, bool is4C, bool atomicOnly) {

    size_t nBasis = basis.nBasis;
    std::shared_ptr<cqmatrix::Matrix<IntsT>> overlapMatrix;
    std::shared_ptr<cqmatrix::Matrix<IntsT>> kineticMatrix;

    // Compute overlap and kinetic matrices
    OnePInts<IntsT> overlapInts(nBasis);
    EMPerturbation emPert;
    overlapInts.computeAOInts(basis, mol, emPert, OPERATOR::OVERLAP, hamiltonianOptions);
    if (is4C) {
      overlapMatrix = std::make_shared<cqmatrix::Matrix<IntsT>>(overlapInts.matrix());
      overlapInts.computeAOInts(basis, mol, emPert, OPERATOR::KINETIC, hamiltonianOptions);
      kineticMatrix = std::make_shared<cqmatrix::Matrix<IntsT>>(std::move(overlapInts.matrix()));
    } else {
      overlapMatrix = std::make_shared<cqmatrix::Matrix<IntsT>>(std::move(overlapInts.matrix()));
    }

    // Determine shells to keep (center by center)
    std::vector<size_t> keptShells = two_step_remove_linearly_dependent_shells(
        basis, mol, overlapMatrix, kineticMatrix, linearDependencyThreshold, is4C, atomicOnly);

    // Find out which shells are removed in keptShells
    std::vector<libint2::Shell> removedShells;
    size_t nextIndex = 0;
    for (size_t curIndex : keptShells) {
      while (nextIndex < curIndex) {
        removedShells.push_back(basis.shells[nextIndex]);
        nextIndex++;
      }
      nextIndex++;
    }

    if (removedShells.empty()) { // No shell removed
      return;

    } else {
      std::cout << "Warning: Linear dependency in basisset detected, we removed some shells to avoid linear dependency." << std::endl;
      std::cout << "Removed " << removedShells.size() << " Shells:" << std::endl << bannerTop << std::endl;

      std::cout << "  " << "  " << std::left;
      std::cout << std::setw(5) << "#" ;
      std::cout << std::setw(5) << "L" ;
      std::cout << std::setw(15) << std::right << "Exponents";
      std::cout << std::setw(30) << std::right << "Basis Center";
      std::cout << std::endl << std::endl;

      size_t iShell = 0;
      for(auto& removedShell : removedShells){
        std::cout << "  " << "  " << std::left << std::setprecision(4)
                  << std::scientific;
        std::cout << std::setw(5) << iShell++;
        std::cout << std::setw(5) << removedShell.contr[0].l << std::right;
        std::cout << std::setw(15) << removedShell.alpha[0];
        std::cout << "          ( " << removedShell.O[0] << ", " << removedShell.O[1] << ", " << removedShell.O[2] << " )";

        std::cout << std::endl;
      }
      std::cout << std::endl << bannerEnd << std::endl << std::endl;

      // Update the basis set
      std::vector<libint2::Shell> newShells;
      for (size_t shellIndex : keptShells) {
        newShells.push_back(std::move(basis.shells[shellIndex]));
      }
      basis.shells = std::move(newShells);
      basis.update();

      std::cout << "Shells used in the calculation:" << std::endl;
      std::cout << basis << std::endl;
    }

  }

  template void remove_linear_dependency<double>(
      BasisSet &basis, Molecule &mol, const HamiltonianOptions &hamiltonianOptions,
      double linearDependencyThreshold, bool is4C, bool atomicOnly);
  template void remove_linear_dependency<dcomplex>(
      BasisSet &basis, Molecule &mol, const HamiltonianOptions &hamiltonianOptions,
      double linearDependencyThreshold, bool is4C, bool atomicOnly);

  template<typename IntsT>
  std::vector<size_t> two_step_remove_linearly_dependent_shells(
      BasisSet &originalBasis, Molecule &mol,
      std::shared_ptr<cqmatrix::Matrix<IntsT>> overlapMatrix,
      std::shared_ptr<cqmatrix::Matrix<IntsT>> kineticMatrix,
      double linearDependencyThreshold,
      bool is4C, bool atomicOnly) {

    // Record atoms finished linear dependency removal, and the corresponding shells
    std::vector<size_t> finishedAtomNumbers;
    std::vector<std::pair<size_t, size_t>> finishedAtomOriginalShellRanges;
    std::vector<std::vector<size_t>> finishedAtomKeptShells;
    // A function to check if a finished atom matches the current shell center
    auto match_finished_atom = [&finishedAtomNumbers, &finishedAtomOriginalShellRanges, &originalBasis, &mol](
        const std::pair<size_t, size_t> &curShellRange) -> std::optional<size_t> {
      auto shellCenter = originalBasis.shells[curShellRange.first].O;
      size_t atomIndex = std::distance(mol.atoms.begin(),
          std::find_if( mol.atoms.begin(),
                        mol.atoms.end(),
                        [&shellCenter](const Atom &a){
                          return a.coord == shellCenter;
                        }));
      size_t atomNumber = mol.atoms[atomIndex].atomicNumber;
      for (size_t i = 0; i < finishedAtomNumbers.size(); i++) {
        if (atomNumber != finishedAtomNumbers[i])
          continue; // Different element, skip
        std::pair<size_t, size_t> &finishedShellRangeI = finishedAtomOriginalShellRanges[i];
        if (finishedShellRangeI.second - finishedShellRangeI.first != curShellRange.second - curShellRange.first)
          continue; // Different number of shells, skip
        bool allMatch = true;
        size_t shellOffset = curShellRange.first - finishedShellRangeI.first;
        for (size_t j = finishedShellRangeI.first; j < finishedShellRangeI.second; j++) {
          libint2::Shell &curShell = originalBasis.shells[j + shellOffset];
          libint2::Shell &finishedShell = originalBasis.shells[j];
          if (curShell.alpha != finishedShell.alpha or curShell.contr != finishedShell.contr) {
            allMatch = false; // Shell not found in finished atom, mismatch
            break;
          }
        }
        if (allMatch) { // Same atom and same shells, match found
          return i; // Return the index of the finished atom
        }
      }
      finishedAtomNumbers.push_back(atomNumber);
      return std::nullopt; // No match found
    };

    // Function to remove linearly dependent shells for a given center, and record the kept shells for potential future matching
    auto remove_atomic_linearly_dependent_shells = [&](
        const std::pair<size_t, size_t> &curShellRange, std::vector<size_t> &keptShells) -> void {
      // Check if we have processed the same atom before, if so reuse the result
      std::optional<size_t> matchedFinishedAtomIndex = match_finished_atom({curShellRange.first, curShellRange.second});
      if (matchedFinishedAtomIndex.has_value()) { // Match found, reuse the kept shells from the finished atom
        size_t matchedIndex = matchedFinishedAtomIndex.value();
        size_t shellOffset = curShellRange.first - finishedAtomOriginalShellRanges[matchedIndex].first;
        for (size_t shellIndex : finishedAtomKeptShells[matchedIndex]) {
          keptShells.push_back(shellIndex + shellOffset); // Append the matched shells with the appropriate offset
        }

      } else { // No match found, process the current center shells and record the result for potential future reuse
        std::vector<size_t> centerShellsKept(curShellRange.second - curShellRange.first);
        std::iota(centerShellsKept.begin(), centerShellsKept.end(), curShellRange.first);
        remove_linearly_dependent_shells(
            originalBasis, *overlapMatrix, centerShellsKept, linearDependencyThreshold);
        if (is4C) {
          remove_linearly_dependent_shells(
              originalBasis, *kineticMatrix, centerShellsKept,
              2.0 * SpeedOfLight * SpeedOfLight * linearDependencyThreshold);
        }
        // Append to keptShells
        keptShells.insert(keptShells.end(), centerShellsKept.begin(), centerShellsKept.end());

        // Record finished atom shells for potential matching with future centers
        finishedAtomOriginalShellRanges.emplace_back(curShellRange.first, curShellRange.second);
        finishedAtomKeptShells.push_back(std::move(centerShellsKept));
      }
    };

    // Determine shells to keep (center by center)
    std::vector<size_t> keptShells;
    auto preShellCenter = originalBasis.shells[0].O;
    size_t curentCenterShellStart = 0;
    for (size_t i = 0; i < originalBasis.nShell; i++) {
      if (originalBasis.shells[i].O != preShellCenter) {
        // New center encountered, process the previous center shells
        remove_atomic_linearly_dependent_shells({curentCenterShellStart, i}, keptShells);

        // Update for new center
        preShellCenter = originalBasis.shells[i].O;
        curentCenterShellStart = i;
      }
      if (i == originalBasis.nShell - 1) {
        // Last shell, process the last center shells
        remove_atomic_linearly_dependent_shells({curentCenterShellStart, i+1}, keptShells);
      }
    }

    // Return if only atomic level removal is desired
    if (atomicOnly) {
      return keptShells;
    }

    // Final check on the full basis set
    remove_linearly_dependent_shells(
        originalBasis, *overlapMatrix, keptShells, linearDependencyThreshold);

    if (is4C) {
      remove_linearly_dependent_shells(
          originalBasis, *kineticMatrix, keptShells,
          2.0 * SpeedOfLight * SpeedOfLight * linearDependencyThreshold);
    }

    return keptShells;

  }

  template std::vector<size_t> two_step_remove_linearly_dependent_shells(
      BasisSet &originalBasis, Molecule &mol,
      std::shared_ptr<cqmatrix::Matrix<double>> overlapMatrix,
      std::shared_ptr<cqmatrix::Matrix<double>> kineticMatrix,
      double linearDependencyThreshold,
      bool is4C, bool atomicOnly);

  template std::vector<size_t> two_step_remove_linearly_dependent_shells(
      BasisSet &originalBasis, Molecule &mol,
      std::shared_ptr<cqmatrix::Matrix<dcomplex>> overlapMatrix,
      std::shared_ptr<cqmatrix::Matrix<dcomplex>> kineticMatrix,
      double linearDependencyThreshold,
      bool is4C, bool atomicOnly);

  template<typename IntsT>
  void remove_linearly_dependent_shells(
      const BasisSet &originalBasis,
      cqmatrix::Matrix<IntsT>& overlapMatrix,
      std::vector<size_t>& keptShells,
      double linearDependencyThreshold) {

    if (not overlapMatrix.isSquareMatrix())
      CErr("Error: Overlap matrix is not square in remove_linearly_dependent_shells()", std::cout);
    if (overlapMatrix.nRows() != originalBasis.nBasis)
      CErr("Error: Overlap matrix size does not match basis set in remove_linearly_dependent_shells()", std::cout);
    
    // Index of current to full basis functions
    std::vector<size_t> cur2full;
    double *eVal = nullptr;
    size_t NB = 0;
    double lowestOverlapEigenvalue = 0.0;

    auto subOverlapMatrix =
        [&overlapMatrix](const std::vector<size_t>& cur2full) -> cqmatrix::Matrix<IntsT> {
          size_t NB = cur2full.size();

          // Compute the overlap matrix
          cqmatrix::Matrix<IntsT> subOverlap(NB, NB);
          for (size_t j = 0; j < NB; j++) {
            for (size_t i = 0; i < NB; i++) {
              subOverlap(i, j) = overlapMatrix(cur2full[i], cur2full[j]);
            }
          }
          return subOverlap;
    };

    while (lowestOverlapEigenvalue < linearDependencyThreshold) {
      // Update the current to full primitive index
      cur2full.clear();
      for (const size_t shellIndex : keptShells) {
        size_t bfStart = originalBasis.mapSh2Bf[shellIndex];
        size_t bfEnd = shellIndex + 1 < originalBasis.mapSh2Bf.size() ? originalBasis.mapSh2Bf[shellIndex + 1] : originalBasis.nBasis;
        for (size_t bfIndex = bfStart; bfIndex < bfEnd; bfIndex++) {
          cur2full.push_back(bfIndex);
        }
      }

      // Compute the overlap matrix
      cqmatrix::Matrix<IntsT> currentOverlap = subOverlapMatrix(cur2full);
      cqmatrix::Matrix<IntsT> removeIOverlap(currentOverlap); // Make a copy for removing shells

      // Compute the eigenvalues (singular values) of the overlap matrix
      NB = cur2full.size();
      if (eVal == nullptr) eVal = CQMemManager::get().template malloc<double>(NB);
      IntsT *eVec = currentOverlap.pointer();
      int info = lapack::gesvd(lapack::Job::NoVec, lapack::Job::OverwriteVec, NB, NB, currentOverlap.pointer(), NB,
                               eVal, nullptr, NB, nullptr, NB);
      lowestOverlapEigenvalue = eVal[NB-1];

      // Print debug info
#ifdef DEBUG_BASIS_LIN_DEP
      std::cout << std::scientific << std::setprecision(12);
      std::cout << "Lowest Eigenvalue: " << lowestOverlapEigenvalue << std::endl;
#endif

      if (lowestOverlapEigenvalue < linearDependencyThreshold) {

        // Determine which shell to remove
        std::vector<double> smallestEigAfterRemoveShell(keptShells.size(), 0.0);

        // Loop over each shell and try removing it to compute the lowest eigenvalue after removing the shell
        size_t nBasisAhead = 0;
        currentOverlap = removeIOverlap; // Reset the overlap matrix
        for (size_t i = 0; i < keptShells.size(); i++) {

          // Determine number of basis functions to remove
          size_t removedNB = originalBasis.shells[keptShells[i]].size();
          size_t removeINB = NB - removedNB;
          size_t nBasisTailBegin = nBasisAhead + removedNB;
          size_t nBasisTail = NB - nBasisTailBegin;

          // Gather the overlap matrix after removing the shell
          SetMat('N', nBasisAhead, nBasisAhead, IntsT(1.0),
                 currentOverlap.pointer(), NB,
                 removeIOverlap.pointer(), removeINB);
          SetMat('N', nBasisTail, nBasisAhead, IntsT(1.0),
                 currentOverlap.pointer() + nBasisTailBegin, NB,
                 removeIOverlap.pointer() + nBasisAhead, removeINB);
          SetMat('N', nBasisAhead, nBasisTail, IntsT(1.0),
                 currentOverlap.pointer() + nBasisTailBegin * NB, NB,
                 removeIOverlap.pointer() + nBasisAhead * removeINB, removeINB);
          SetMat('N', nBasisTail, nBasisTail, IntsT(1.0),
                 currentOverlap.pointer() + nBasisTailBegin + nBasisTailBegin * NB, NB,
                 removeIOverlap.pointer() + nBasisAhead + nBasisAhead * removeINB, removeINB);

          // Compute the eigenvalues of the overlap matrix
          int info = lapack::gesvd(lapack::Job::NoVec, lapack::Job::OverwriteVec, removeINB, removeINB,
                                   removeIOverlap.pointer(), removeINB,
                                   eVal, nullptr, removeINB, nullptr, removeINB);
          smallestEigAfterRemoveShell[i] = eVal[removeINB - 1];

          nBasisAhead = nBasisTailBegin;
        }
#ifdef DEBUG_BASIS_LIN_DEP
        std::vector<size_t> shellIndices(keptShells.size());
        std::iota(shellIndices.begin(), shellIndices.end(), 0);
        std::sort(shellIndices.begin(), shellIndices.end(),
                  [&originalBasis, &keptShells](size_t kept1, size_t kept2) {
                    size_t i1 = keptShells[kept1], i2 = keptShells[kept2];
                    if (originalBasis.shells[i1].O == originalBasis.shells[i2].O) {
                      if (originalBasis.shells[i1].contr[0].l != originalBasis.shells[i2].contr[0].l) {
                        return originalBasis.shells[i1].contr[0].l < originalBasis.shells[i2].contr[0].l;
                      } else {
                        return originalBasis.shells[i1].alpha[0] > originalBasis.shells[i2].alpha[0];
                      }
                    } else
                      return i1 < i2;
                  });
        for (size_t argI = 0; argI < shellIndices.size(); argI++) {
          size_t i = keptShells[shellIndices[argI]];
          std::cout << "  Remove Shell " << i << " (L=" << originalBasis.shells[i].contr[0].l
                    << ", Exponent=" << originalBasis.shells[i].alpha[0] << ")"
                    << " -> Smallest Eigenvalue: " << std::scientific << std::setprecision(12)
                    << smallestEigAfterRemoveShell[shellIndices[argI]]
                    << std::endl;
        }
#endif

        // Find the shell with the largest eigenvalue
        size_t linearDepShellIndex = std::distance(smallestEigAfterRemoveShell.begin(),
                                                   std::max_element(smallestEigAfterRemoveShell.begin(), smallestEigAfterRemoveShell.end()));

#ifdef DEBUG_BASIS_LIN_DEP
        std::cout << "Removing Shell " << keptShells[linearDepShellIndex]
                  << " (L=" << originalBasis.shells[keptShells[linearDepShellIndex]].contr[0].l
                  << ", Exponent=" << originalBasis.shells[keptShells[linearDepShellIndex]].alpha[0] << ")"
                  << " with Smallest Eigenvalue: " << std::scientific << std::setprecision(12)
                  << smallestEigAfterRemoveShell[linearDepShellIndex]
                  << std::endl;
#endif

        // Remove the basis shell
        auto it = keptShells.begin();
        std::advance(it, linearDepShellIndex);
        keptShells.erase(it);

      }

    }

    // Free Scratch Space
    if (eVal) CQMemManager::get().free(eVal);
  }

  template void remove_linearly_dependent_shells(
      const BasisSet &originalBasis,
      cqmatrix::Matrix<double>& overlapMatrix,
      std::vector<size_t>& keptShells,
      double linearDependencyThreshold);
  template void remove_linearly_dependent_shells(
      const BasisSet &originalBasis,
      cqmatrix::Matrix<dcomplex>& overlapMatrix,
      std::vector<size_t>& keptShells,
      double linearDependencyThreshold);

}; // namespace ChronusQ