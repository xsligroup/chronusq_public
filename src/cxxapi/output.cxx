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

#include <cxxapi/output.hpp>
#include <util/mpi.hpp>
#include <util/timer.hpp>

namespace ChronusQ {

  template <typename Duration>
  void printTimerInternals(std::ostream& out, std::string unit) {

    ROOT_ONLY(MPI_COMM_WORLD);

    Timer& timer = *ProgramTimer::instance();

    // Magic formatting numbers
    size_t sectionWidth =  36;
    size_t durationWidth = 20;
    size_t precision = 4;
    size_t limit = std::pow(10, -(precision+1));

    // Printing helpers
    auto printBase = [&](std::string message, double dur, double avg) {

      // Don't print if it would be all 0
      if (dur < 1e-5)
        return;

      // Print the base
      out << std::left << std::setw(sectionWidth) << message << ":";
      out << std::fixed << std::setprecision(precision) << std::right;
      // Print the total
      out << std::setw(durationWidth) <<  dur;
      // Optionally print the average if it differs from the total
      if ( std::abs(avg - dur) > limit )
        out << "  " << std::setw(durationWidth) << avg;
      out << std::endl;
    };

    auto printReg = [&](std::string message, std::string label, size_t id) {
      auto summary = timer.getDurationSummary<Duration>(label, id);
      auto dur = summary.first.count();
      auto avg = summary.second.count();

      printBase(message, dur, avg);
    };

    auto printRegWithAvgScale = [&](std::string message, std::string label, size_t id,
        double avgScale) {
      auto summary = timer.getDurationSummary<Duration>(label, id);
      auto dur = summary.first.count();
      auto avg = summary.second.count();

      printBase(message, dur, avg * avgScale);
    };
    
    auto printSub = [&](std::string message, std::string label,
                        std::string filter, size_t id) {
      auto summary = timer.getDurationSummary<Duration>(label, id, filter);
      printBase(message, summary.first.count(), summary.second.count());
    };

    auto printAgg = [&](std::string message, std::string label,
                        std::string parent, size_t id) {
      auto summary = timer.getAggDurSummary<Duration>(label, parent, id);
      auto dur = summary.first.count();
      auto avg = summary.second.count();

      printBase(message, dur, avg);
    };


    // Print header
    out << std::setfill(' ');
    out << "\nTiming Summary:\n" << BannerTop << std::endl;

    std::string titleTotal = "Total Duration (" + unit + ")";
    std::string titleAvg = "Average Duration (" + unit + ")";
    out << std::left << std::setw(sectionWidth) << "  Section" << " ";
    out << std::right << std::setw(durationWidth) << titleTotal << "  ";
    out << std::right << std::setw(durationWidth) << titleAvg << "\n";
    out << bannerTop << std::endl;


    //
    // Print actual information
    //

    // Top level / general
    auto cqId = timer.getLabelId("Chronus Quantum");
    printReg("  Total program time", "Chronus Quantum", cqId);
    printReg("  - Memory allocation", "Memory Allocation", cqId);
    printAgg("  - Core Hamiltonian", "Form Core H", "Chronus Quantum", cqId);
    printAgg("  - Two Body Integrals", "Incore 2e Ints", "Chronus Quantum", cqId);
    printAgg("  - Guess Formation", "Form Guess", "Chronus Quantum", cqId);

    // SCF specific
    auto scfId = timer.getLabelId("SCF Total");
    if ( scfId != 0 ) {
      printReg("  - SCF", "SCF Total", scfId);
      printReg("    - SCF Iteration", "SCF Iter", scfId);
      printReg("      - Fock Formation", "Form Fock", scfId);
      printReg("        - J Contraction", "J Contract", scfId);
      printReg("        - K Contraction", "K Contract", scfId);
      printReg("        - ERI Contraction", "Contract Total", scfId);
      printReg("        - Vxc Formation", "Form VXC", scfId);
      // printAgg("          - Integrals", "Direct Int Form",
      //          "Contract Total", scfId);
      // printAgg("          - Contraction", "Direct Den Contract",
      //          "Contract Total", scfId);
      printReg("      - Property Evaluation", "Compute Properties", scfId);
    }

    // RT specific
    auto rtId = timer.getLabelId("Real Time Total");
    if ( rtId != 0 ) {
      printReg("  - RT Propagation", "Real Time Total", rtId);
      printReg("    - RT Iteration", "Real Time Iter", rtId);
      printReg("      - Fock Formation", "Form Fock", rtId);
      printReg("      - Propagator Formation", "Propagator Formation", rtId);
      printReg("      - Propagation", "Propagate Density", rtId);
    }

    // RESP specific
    auto respId = timer.getLabelId("Response Total");
    if ( respId != 0 ) {
      printReg("  - Response", "Response Total", respId);
      printReg("    - Full Hessian Formation", "Full Hessian Form", respId);
      //Types
      printReg("    - Full Diagonalization", "Full Diagonalize", respId);
      printReg("    - Partial Diagonalization", "Iter Diagonalize", respId);
      printReg("    - Full Frequency Response", "Full FDR", respId);
      printReg("    - Iterative Freq Response", "Iter FDR", respId);
      // (Iterative) Residue only
      printReg("      - Diagonalize Iteration", "Diagonalize Iter", respId);
      printSub("        - Hessian Contraction", "Direct Hessian Contract", "Diagonalize Iter", respId);
      // FDR only
      printReg("      - Frequencies", "Omega", respId);
      printReg("        - Linear Solution", "Solve Linear System", respId);
      printReg("        - Linear Solve Iteration", "Lin Solve Iter", respId);
      printSub("          - Hessian Contraction", "Direct Hessian Contract", "Lin Solve Iter", respId);
      // All
      printReg("    - Property Evaluation", "Property Eval", respId);
    }

    // MCSCF specific
    auto mcscfId = timer.getLabelId("MCSCF Total");
    if ( mcscfId != 0 ) {
      printReg("  - MCSCF", "MCSCF Total", mcscfId);
      // CI
      printReg("    - Solving CI", "Solve CI", mcscfId);
      printSub("      - Integral Transformation", "Integral Trans", "Solve CI",mcscfId);
      printSub("      - Diagonalization", "Diagonalization", "Solve CI",mcscfId);
      printReg("        - Full Matrix Formation", "Full Matrix", mcscfId);
      printReg("        - Sigma Formation", "Sigma", mcscfId);
      // Orbital Rotation
      printReg("    - Orbital Rotation", "Orbital Rotation", mcscfId);
      printReg("      - Gradient Formation", "Form Gradient", mcscfId);
      printReg("      - Hessian Formation", "Form Hessian", mcscfId);
      printReg("      - MO Rotation", "Rotate MO", mcscfId);
      printReg("      - IVO Generation", "Gen IVOs", mcscfId);
      printReg("    - Property Evaluation", "Property Eval", mcscfId);
    }

    // Visualization
    auto cubeEvalId = timer.getLabelId("Cube Eval");
    if  ( cubeEvalId != 0 ) {
      printReg("  - Cubegen Evaluation", "Cube Eval", cubeEvalId);
    }
    
    // PERTURB specific
    auto perturbId = timer.getLabelId("PERTURB Total");
    if ( perturbId != 0 ) {
      printReg("  - PERTURB", "PERTURB Total", perturbId);
      printReg("    - Integral Transformation", "Integral Trans", perturbId);
      printReg("    - XMS rotation", "Extend MS", perturbId);
      printReg("    - Fock build", "Fock build", perturbId);
      printReg("    - Full LHS build", "Full LHS build", perturbId);
      printReg("    - LHS diagonal build", "LHS diagonal build", perturbId);
      printReg("    - RHS build with Sigma", "RHS build", perturbId);
      printReg("    - Solve Linear System", "Solve Linear System", perturbId);
      printReg("    - LHS sigma build", "LHS sigma build", perturbId);
      printReg("    - compute HV with Sigma", "compute HV", perturbId);
      printReg("    - MS H effective", "MS Heff", perturbId);
      printSub("      -compute HV with Sigma", "compute HV", "MS Heff", perturbId); 
      printSub("      -diagonalization", "Heff diag", "MS Heff", perturbId);
    }


    // Print footer
    out << BannerEnd << std::endl;

  }; // printTimerInternals

  void printTimerSummary(std::ostream& out) {

    Timer& timer = *ProgramTimer::instance();

    if ( !timer.options.doSummary )
      return;

    if ( timer.options.unit == MILLISECONDS )
      printTimerInternals<CQMillisecond>( out, "ms" );
    else if ( timer.options.unit == SECONDS )
      printTimerInternals<CQSecond>( out, "s" );
    else if ( timer.options.unit == MINUTES )
      printTimerInternals<CQMinute>( out, "min" );
  };

}
