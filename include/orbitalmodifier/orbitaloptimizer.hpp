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

#include <orbitalmodifier.hpp>

namespace ChronusQ {

struct SCFConvergence {

    double deltaEnergy;   ///< Convergence of Energy
    double currentEnergy;
    double previousEnergy;
    double rmsdP;         ///< RMS change in density
    double maxdP;         ///< Max change in density
    double nrmFDC;        ///< 2-Norm of [F,D]
    double maxFDC;        ///< Maximum element of [F,D]
    size_t nSCFIter = 0;  ///< Number of SCF Iterations

};   // SCFConvergence struct

/*
 *     Brief: Object to holds the overlapping functions for orbital optimization
 *            such as printing the header and iteration progress. Thus, getNewOrbitals
 *            is modified depending on which step is used. ie. SCF or Newton-Raphson
 */
template<typename MatsT>
class OrbitalOptimizer : public OrbitalModifier<MatsT> {
  private:
    std::vector<cqmatrix::Matrix<MatsT>> prevOnePDM;   ///< Previous density used to test convergence
    double prevEnergy;                             ///< Previous Energy to test convergence

  public:

    SCFControls scfControls;
    SCFConvergence scfConv;
    bool doingDamp;                                ///< Whether damping is currently on or off (only used for printing)

    // Constructor
    OrbitalOptimizer(SCFControls sC, MPI_Comm comm, OrbitalModifierDrivers<MatsT> modOpt):
      scfControls(sC), OrbitalModifier<MatsT>(comm, modOpt) {

        // Allocate prevOnePDM
        vecShrdPtrMat<MatsT> onePDM = this->orbitalModifierDrivers.getOnePDM();
        for( size_t a = 0; a < onePDM.size(); a++ ) {
          prevOnePDM.emplace_back(onePDM[a]->dimension());
          prevOnePDM[a] = *onePDM[a];
        }
    };

    // Destructor
    ~OrbitalOptimizer() {};

    // Perform an SCF procedure (see include/singleslater/scf.hpp for docs)
    void runOrbitalModifier(EMPerturbation&, vecMORef<MatsT>&, vecEPtr&) override;

    // Evaluate convergence
    bool evaluateProgress(EMPerturbation&);
    double computeDensityConv();
    virtual double computeFDCConv() { return 0.; };

    //   Print SCF header, footer and progress
    void printRunHeader(std::ostream& out, EMPerturbation&) const override;
    void printHeaderFinal(std::ostream& out) const;
    void printIteration(std::ostream& out = std::cout, bool printDiff = true) const override;

    // Common SCF Functions
    void computeEigenvalues(EMPerturbation& pert, vecMORef<MatsT>&, vecEPtr&);
};
};   // namespace ChronusQ
