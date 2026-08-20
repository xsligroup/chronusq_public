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

#include <chronusq_sys.hpp>
#include <cerr.hpp>
#include <util/files.hpp>
#include <configinteraction.hpp>
#include <posthartreefock.hpp>
#include <singleslater.hpp>

//#define _DEBUG_MRPT

namespace ChronusQ {

  struct MRPTSettings {

    bool STATEAVERAGE = false;  // < use state-averaged density matrix to build Fock.
    bool ENPT = true; // < Classic Epstein-Nesbet Style PT Order 2
    bool GVVPT = false;  // < build H_eff within the whole primary space
    double LEVELSHIFT = 0.0; // < Level Shift operator
    double EPS = 0.0; // < PT2 threshold contribution to one-electron and two-electron terms

    //GVVPT2 settings:
    size_t SECONDARYROOTS = 0;

    // frozen orbital settings
    size_t FROZENCORE = 0;
    size_t FROZENVIRTUAL = 0;
    std::string SELECTVIRTUAL;

    // Use only Scalar two-pints for PT2
    bool SCALAR = false;

  };


  template <typename MatsT, typename IntsT>
  class DasPerturb : public PostHartreeFock<MatsT,IntsT> {

  protected:

    std::vector<double> E0_; // < General Object to Zero-Order Energies
    std::vector<double> E2_; // < General Object to PT2 Correlation Energies 
    std::shared_ptr<ConfigurationInteraction<MatsT,IntsT>> RefMCWfn_ = nullptr; // < MCSCF reference
    std::vector<size_t> Target_States_; // < State of Interest
    std::shared_ptr<DeterminantFactory> PTFactory_ = nullptr; // < PT determinant factory
    std::shared_ptr<DASCIBuilder<MatsT>> ptBuilder_ = nullptr; // < Ptr to CI builder
    std::shared_ptr<cqmatrix::Matrix<MatsT>> fockDiag_ = nullptr; // < Diagonal Fock for GVVPT2
    std::shared_ptr<cqmatrix::Matrix<MatsT>> oneRDMSA_ = nullptr; // State-Averaged one RDM

    // Caching the one- and two-electron integrals:
    std::vector<double> hDiagCache_;
    std::vector<double> ttuu2eCache_;

  public:

    MRPTSettings PTopts;

    // Disable default, copy and move constructors
    DasPerturb()                   = delete;
    DasPerturb(const DasPerturb &) = delete;
    DasPerturb(DasPerturb &&)      = delete;

    // Constructors

    /**
     *  \brief DasPerturb Constructor
     *  
     *  Obtain references from a "reference" CASCI/CASSCF object
     *
     */
    DasPerturb(std::shared_ptr<ConfigurationInteraction<MatsT,IntsT>> RefMCWfn, std::vector<size_t> Target_States) : 
        PostHartreeFock<MatsT,IntsT>(RefMCWfn->reference(), Target_States.size()), 
        Target_States_(Target_States), RefMCWfn_(RefMCWfn) { 

      // MRPT with 1C reference not implemented
      if (RefMCWfn->reference()->nC==1)
        CErr("1c + MRPT not yet implemented");
    }

    ~DasPerturb() { dealloc(); }

    // PERTURB procedural functions
    void run(EMPerturbation &);
    void PTInitialize();

    // EN2 specific functions:
    void CIEnergy();
    void buildDiagIntCache();
    void computeAmplitudes(size_t &, DistributedVectors<MatsT> &, double &);
    void computeEN2(size_t &);
#ifdef CQ_ENABLE_SPARSE 
    void computeSparseAmplitudes(size_t &, DistributedSparseVectors<MatsT> &, double &);
    void computeEN2Sparse(size_t &);
#endif

    // GVVPT specific functions:
    void semiCanonicalize();
    void computeSAOneRDM();
    void formFockDiag(); 
    void formFockDiagStateAvg(); 
    void computeZeroEnergy(); 
    void computeSAZeroEnergy();
    void buildHX(std::pair<MatsT,MatsT> &, 
      const std::shared_ptr<DistributedVectors<MatsT>> &, size_t);
    void buildHX(std::pair<MatsT,MatsT> &, 
      const std::shared_ptr<DistributedVectors<MatsT>> &,
      const std::shared_ptr<DistributedVectors<MatsT>> &, size_t, size_t);
#ifdef CQ_ENABLE_SPARSE
    void buildHX(std::pair<MatsT,MatsT> &,
      DistributedSparseVectors<MatsT> &, size_t);
    void buildHX(std::pair<MatsT,MatsT> &,
      DistributedSparseVectors<MatsT> &, size_t, size_t);
#endif
    void formEffectiveHSparse(std::shared_ptr<cqmatrix::Matrix<MatsT>> &);
    void diagEffHSparse();
    void formEffectiveH(std::shared_ptr<cqmatrix::Matrix<MatsT>> &);
    void diagEffH();

   
    //Printing functions
    void printMRPTHeader();
    void printMRPTFooter();
    void buildX();

    // Save Data to HDF5
    void saveCurrentStates();

    // Implement pure virtual functions from PostHartreeFock
    void computeTDM(size_t, size_t, std::shared_ptr<cqmatrix::Matrix<MatsT>>) 
      override {};
    void compute2TDM(size_t, size_t, std::shared_ptr<InCore4indexTPI<MatsT>>) 
      override {};
    void compute2RDM(size_t, size_t, std::shared_ptr<InCore4indexTPI<MatsT>>) 
      override {};

    // Memory functions
    void alloc();
    void dealloc();

  }; // class DasPerturb

}; // namespace ChronusQ
