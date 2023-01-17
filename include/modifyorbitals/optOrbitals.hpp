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

#include <modifyorbitals.hpp>

namespace ChronusQ {

struct SCFConvergence {

    double deltaEnergy;   ///< Convergence of Energy
    double RMSDen;        ///< RMS change in Scalar density
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
class OptimizeOrbitals : public ModifyOrbitals<MatsT> {
  private:
    std::vector<SquareMatrix<MatsT>> prevOnePDM;   ///< Previous density used to test convergence
    double prevEnergy;                             ///< Previous Energy to test convergence

  public:
//xslis
    SCFControls scfControls;
//xslie
    SCFConvergence scfConv;
    bool doingDamp;                                ///< Whether damping is currently on or off (only used for printing)

    // Constructor
    OptimizeOrbitals(SCFControls sC, MPI_Comm comm, ModifyOrbitalsOptions<MatsT> modOpt, CQMemManager& mem):
        scfControls(sC), ModifyOrbitals<MatsT>(comm, modOpt, mem) {

        // Allocate prevOnePDM
        VecShrdPtrMat<MatsT> den = this->modOrbOpt.getOnePDM();
        for( size_t a = 0; a < den.size(); a++ ) {
            prevOnePDM.emplace_back(den[a]->memManager(), den[a]->dimension());
            prevOnePDM[a] = *den[a];
        }
    };

    // Destructor
    ~OptimizeOrbitals() {};

    // Perform an SCF procedure (see include/singleslater/scf.hpp for docs)
    void runModifyOrbitals(EMPerturbation&, VecMORef<MatsT>&, VecEPtr&);

    // Evaluate convergence
    bool evalProgress(EMPerturbation&);
    double computeDensityConv();
    virtual double computeFDCConv() { return 0.; };

    //   Print SCF header, footer and progress
    void printRunHeader(std::ostream& out, EMPerturbation&) const;
    void printHeaderFinal(std::ostream& out) const;
    void printIteration(std::ostream& out = std::cout, bool printDiff = true) const;

    // Common SCF Functions
    void computeEigenvalues(VecMORef<MatsT>&, VecEPtr&);
};
};   // namespace ChronusQ
