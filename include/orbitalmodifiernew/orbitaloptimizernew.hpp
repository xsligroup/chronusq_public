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

#include <orbitalmodifiernew.hpp>

namespace ChronusQ {

struct SCFConvergenceNew {

    double deltaEnergy;   ///< Convergence of Energy
    double currentEnergy;
    double previousEnergy;
    double rmsdP;         ///< RMS change in density
    double maxdP;         ///< Max change in density
    size_t nSCFIter = 0;  ///< Number of SCF Iterations
};   // SCFConvergence struct

/*
 *     Brief: Object to holds the overlapping functions for orbital optimization
 *            such as printing the header and iteration progress. Thus, getNewOrbitals
 *            is modified depending on which step is used. ie. SCF or Newton-Raphson
 */
template <template <typename, typename> class singleSlaterT, typename MatsT, typename IntsT>
class OrbitalOptimizerNew : public OrbitalModifierNew<singleSlaterT, MatsT, IntsT> {
  private:
    std::vector<cqmatrix::Matrix<MatsT>> prevOnePDM;   ///< Previous density used to test convergence

  public:

    SCFControls scfControls;
    SCFConvergenceNew scfConv;
    bool doingDamp;                                ///< Whether damping is currently on or off (only used for printing)

    // Constructor
    OrbitalOptimizerNew(SCFControls sC, singleSlaterT<MatsT,IntsT> &referenceSS, MPI_Comm comm):
      scfControls(sC), OrbitalModifierNew<singleSlaterT,MatsT,IntsT>(referenceSS, comm) {

        // Allocate prevOnePDM
        vecShrdPtrMat<MatsT> onePDM = this->singleSlaterSystem.getOnePDM();
        for( size_t a = 0; a < onePDM.size(); a++ ) {
          prevOnePDM.emplace_back(onePDM[a]->dimension());
          prevOnePDM[a] = *onePDM[a];
        }
    };

    // Destructor
    ~OrbitalOptimizerNew() {};

    // Perform an SCF procedure (see include/singleslater/scf.hpp for docs)
    void run(EMPerturbation&) override;

    // Evaluate convergence
    bool evaluateProgress(EMPerturbation&);
    double computeDensityConv();
    virtual double computeFDCConv() { return 0.; };

    //   Print SCF header, footer and progress
    void printRunHeader(EMPerturbation&) override;
    void printHeaderFinal() const;
    void printIteration(bool printDiff = true) override;

    // Common SCF Functions
    void computeEigenvalues(EMPerturbation& pert);

    virtual void initialize(size_t maxPoints = 0) override {};
};
};   // namespace ChronusQ
