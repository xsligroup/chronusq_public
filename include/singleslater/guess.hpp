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

#include <singleslater.hpp>
#include <singleslater/base.hpp>
#include <singleslater/guessdiff.hpp>
#include <cqlinalg.hpp>
#include <util/matout.hpp>
#include <fockbuilder/matrixfock.hpp>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <wavefunction/base.hpp>

namespace ChronusQ {

  /**
   *  \brief List of default atomic multiplicities for
   *  atomic SCF.
   *
   *  Atomic multiplicities obtained from physics.nist.gov
   */
  static std::map<size_t,size_t> defaultMult(
    {
      { 1  , 2 }, // H
      { 2  , 1 }, // He
      { 3  , 2 }, // Li
      { 4  , 1 }, // Be
      { 5  , 2 }, // B
      { 6  , 3 }, // C
      { 7  , 4 }, // N
      { 8  , 3 }, // O
      { 9  , 2 }, // F
      { 10 , 1 }, // Ne
      { 11 , 2 }, // Na
      { 12 , 1 }, // Mg
      { 13 , 2 }, // Al
      { 14 , 3 }, // Si
      { 15 , 4 }, // P
      { 16 , 3 }, // S
      { 17 , 2 }, // Cl
      { 18 , 1 }, // Ar
      { 19 , 2 }, // K
      { 20 , 1 }, // Ca
      { 21 , 2 }, // Sc
      { 22 , 3 }, // Ti
      { 23 , 4 }, // V
      { 24 , 7 }, // Cr
      { 25 , 6 }, // Mn
      { 26 , 5 }, // Fe
      { 27 , 4 }, // Co
      { 28 , 3 }, // Ni
      { 29 , 2 }, // Cu
      { 30 , 1 }, // Zn
      { 31 , 2 }, // Ga
      { 32 , 3 }, // Ge
      { 33 , 4 }, // As
      { 34 , 3 }, // Se
      { 35 , 2 }, // Br
      { 36 , 1 }, // Kr
      { 37 , 2 }, // Rb
      { 38 , 1 }, // Sr
      { 39 , 2 }, // Y
      { 40 , 3 }, // Zr
      { 41 , 6 }, // Nb
      { 42 , 7 }, // Mo
      { 43 , 6 }, // Tc
      { 44 , 5 }, // Ru
      { 45 , 4 }, // Rh
      { 46 , 1 }, // Pd
      { 47 , 2 }, // Ag
      { 48 , 1 }, // Cd
      { 49 , 2 }, // In
      { 50 , 3 }, // Sn
      { 51 , 4 }, // Sb
      { 52 , 3 }, // Te
      { 53 , 2 }, // I
      { 54 , 1 }, // Xe
      { 55 , 2 }, // Cs
      { 56 , 1 }, // Ba
      { 57 , 2 }, // La
      { 58 , 1 }, // Ce
      { 59 , 4 }, // Pr
      { 60 , 5 }, // Nd
      { 61 , 6 }, // Pm
      { 62 , 7 }, // Sm
      { 63 , 8 }, // Eu
      { 64 , 9 }, // Gd
      { 65 , 6 }, // Tb
      { 66 , 5 }, // Dy
      { 67 , 4 }, // Ho
      { 68 , 3 }, // Er
      { 69 , 2 }, // Tm
      { 70 , 1 }, // Yb
      { 71 , 2 }, // Lu
      { 72 , 3 }, // Hf
      { 73 , 4 }, // Ta
      { 74 , 5 }, // W
      { 75 , 6 }, // Re
      { 76 , 5 }, // Os
      { 77 , 4 }, // Ir
      { 78 , 3 }, // Pt
      { 79 , 2 }, // Au
      { 80 , 1 }, // Hg
      { 81 , 2 }, // Tl
      { 82 , 3 }, // Pb
      { 83 , 4 }, // Bi
      { 84 , 3 }, // Po
      { 85 , 2 }, // At
      { 86 , 1 }, // Rn
      { 87 , 2 }, // Fr
      { 88 , 1 }, // Ra
      { 89 , 2 }, // Ac
      { 90 , 3 }, // Th
      { 91 , 4 }, // Pa
      { 92 , 5 }, // U
      { 93 , 6 }, // Np
      { 94 , 7 }, // Pu
      { 95 , 8 }, // Am
      { 96 , 9 }, // Cm
      { 97 , 6 }, // Bk
      { 98 , 5 }, // Cf
      { 99 , 4 }, // Es
    }
  );

  // Printing out reference types
  std::unordered_map<int,std::string> refMap(
    {
      { 1    ,  "1-Component Restricted"            },
      { 2    ,  "1-Component Unrestricted"          },
      { 3    ,  "1-Component Restricted Open Shell" },
      { 4    ,  "2-Component Unrestricted"          },
      { 5    ,  "4-Component Unrestricted"          },
    }
  );


  /**
   *  \brief Forms a set of guess orbitals for a single slater
   *  determininant SCF in various ways
   */
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::formGuess(const SingleSlaterOptions &ssOptions) {

    ProgramTimer::tick("Form Guess");
     
    // populate AO-fock Matrix
    if( printLevel > 0 )
      std::cout << "  *** Forming Initial Guess Density for SCF Procedure ***"
                << std::endl << std::endl;

    // Form initial density
    // NOTE: The only postcondition that is necessarily satisfied by these
    //   methods is that onePDM has been populated by a guess density. It does
    //   NOT necessarily form an initial Fock matrix (sp. ReadMO/Read1PDM)
    //
    //   MO's are also initialized because they are needed before the first
    //   fock build for RI/Cholesky.
    if (std::dynamic_pointer_cast<MatrixFock<MatsT,IntsT>>(fockBuilder)) {

      EMPerturbation emPert;
      fockBuilder->formFock(*this, emPert, false, 0.0);
      getNewOrbitals();

    } else if( scfControls.guess == RANDOM ) RandomGuess();
      else if( scfControls.guess == READMO ) ReadGuessMO();
      else if( scfControls.guess == READDEN ) ReadGuess1PDM();
      else if( scfControls.guess == FCHKMO ) FchkGuessMO();
      else if( scfControls.guess == NEOTightProton ) NEOTightProtonGuess();
      else if( scfControls.guess == NEOConvergeClassical ) NEOConvergeClassicalGuess(ssOptions);
      else if( this->molecule().nAtoms == 1  and scfControls.guess == SAD ) CoreGuess();
      else if( scfControls.guess == CORE ) CoreGuess();
      else if( scfControls.guess == TIGHT ) TightGuess();
      else if( scfControls.guess == SAD ) SADGuess(ssOptions);
      else CErr("Unknown choice for SCF.GUESS",std::cout);

    // If RANDOM guess, scale the densites appropriately
    // *** Replicates on all MPI processes ***
    if( scfControls.guess == RANDOM ) {

      if (nC == 4)
        CErr("Random guess not implemented for 4-component",std::cout);

      size_t NB = this->basisSet().nBasis;

      double TS =
        this->template computeOBProperty<SCALAR>(this->aoints_->overlap->pointer());

      double TZ =
        this->template computeOBProperty<MZ>(this->aoints_->overlap->pointer());
      double TY =
        this->template computeOBProperty<MY>(this->aoints_->overlap->pointer());
      double TX =
        this->template computeOBProperty<MX>(this->aoints_->overlap->pointer());

      double magNorm = std::sqrt(TZ*TZ + TY*TY + TX*TX);

      this->onePDM->S() *= MatsT(this->nO)/MatsT(TS);

      double factorXYZ = magNorm > 1e-10 ? (this->nOA - this->nOB)/magNorm : 0.0;
      if (this->onePDM->hasZ())
        this->onePDM->Z() *= factorXYZ;
      if (this->onePDM->hasXY()) {
        this->onePDM->Y() *= factorXYZ;
        this->onePDM->X() *= factorXYZ;
      }

    }

    ProgramTimer::tock("Form Guess");

  }; // SingleSlater<T>::formGuess

  /**
   * \brief Populates the initial density as single occupied AOs on the
   *        tightest basis functions for each atom
   */
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::TightGuess() {
    std::vector<size_t> occ;
    BasisSet& basis = this->basisSet();

    this->onePDM->clear();

    // Get the shells on this atom
    std::vector<size_t> atShells;
    for( auto iSh = 0; iSh < basis.shells.size(); iSh++ ) {
      atShells.push_back(iSh);
    }

    // Shells in descending order of tightest primative
    // XXX: This is not necessarily the "tightest" basis function;
    //      We should also account for the coefficients for generalized contractions
    std::stable_sort(atShells.begin(), atShells.end(), 
      [&basis](size_t sh1, size_t sh2) {
        libint2::Shell& shell1 = basis.shells[sh1];
        libint2::Shell& shell2 = basis.shells[sh2];
        auto sh1m = std::max_element(shell1.alpha.begin(), shell1.alpha.end());
        auto sh2m = std::max_element(shell2.alpha.begin(), shell2.alpha.end());
        return *sh1m > *sh2m;
      }
    );

    // Occupy tightest basis functions
    size_t npart = this->nOA;
    size_t ishell = 0;
    while( npart > 0 ) {
      std::cout << "npart: " << npart << std::endl;
      size_t loc = basis.mapSh2Bf[atShells[ishell]];
      size_t nbas = basis.shells[atShells[ishell]].size();
      size_t nadd = std::min(nbas, npart);

      for(size_t iadd = 0; iadd < nadd; iadd++, npart--) {
          this->onePDM->S()(loc+iadd, loc+iadd) = 1.;
      }
      ishell += 1;
    }

    // Make charge and multiplicity correct
    if( this->onePDM->hasZ() and MPIRank(comm) == 0 ) {
      this->onePDM->Z() = (MatsT(this->nOA - this->nOB) / MatsT(this->nOA)) * this->onePDM->S();
    }

    this->onePDM->output(std::cout, "Guess PDM", true);
  };

  /**
   *  \brief Populates the initial Fock matrix with the core
   *  hamiltonian.
   *
   *  I.e. Neglecting electron-electron interaction. While this
   *  is the simplest choice for an initial guess, it is often
   *  a very poor guess.
   */
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::CoreGuess() {

    if( printLevel > 0 )
      std::cout << "    * Forming the Core Hamiltonian Guess (F = H)\n\n";

    // Form the core guess Fock matrix on the root MPI process
    if( MPIRank(comm) == 0 )
      // Copy over the Core Hamiltonian
      *fockMatrix = *coreH;

#ifdef CQ_ENABLE_MPI
    // BCast random fock to all MPI processes
    if( MPISize(comm) > 1 ) {
      std::cerr  << "  *** Scattering the CORE GUESS Fock ***\n";
      size_t NB = this->fockMatrix->dimension();
      for(auto mat : this->fockMatrix->SZYXPointers())
        MPIBCast(mat,NB*NB,0,comm);
    }
#endif

    // Forming Orbitals from core guess
    getNewOrbitals();

  }; // SingleSlater::CoreGuess

  /**
   *
   */
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::SADGuess(const SingleSlaterOptions &ssOptions) {



    // ROOT Communicator
    int color = (MPIRank(comm) == 0) ? 1 : MPI_UNDEFINED;
    MPI_Comm rcomm = MPICommSplit(comm,color,0);

    //std::cerr << "RCOMM " << rcomm << std::endl;

    if( printLevel > 0 )
      std::cout << "    * Forming the Superposition of Atomic Densities Guess"
                << " (SAD)\n\n";

    size_t NB = this->basisSet().nBasis;
    EMPerturbation pert;

    // Zero out the densities (For all MPI processes)
    this->onePDM->clear();


    // Determine the unique atoms
    std::vector<Atom> uniqueElements(this->molecule().atoms);

    // Sort the Atoms by atomic number for std::unique
    std::sort(uniqueElements.begin(),uniqueElements.end(),
      [](Atom &a, Atom &b) {
        return a.atomicNumber < b.atomicNumber;
      });

    // Obtain unique elements on sorted list
    auto it = std::unique(uniqueElements.begin(),uniqueElements.end(),
                [](Atom &a, Atom &b){
                  return a.atomicNumber == b.atomicNumber;
                });

    // Remove excess elements
    uniqueElements.resize(std::distance(uniqueElements.begin(),it));

    if( printLevel > 0 )
      std::cout << "  *** Found " << uniqueElements.size()
                << " unique Atoms in Molecule Specification ***\n";


    // Set all atomic centers to origin and force classical
    for(auto &u : uniqueElements) {
      u.coord = {0.,0.,0.};
      u.quantum = false;
    }

    // Create a map from atoms in molecule to unique atom index
    std::map<size_t,size_t> mapAtom2Uniq;
    for(auto iAtm = 0ul; iAtm < this->molecule().nAtoms; iAtm++) {
      size_t curAtomicNumber = this->molecule().atoms[iAtm].atomicNumber;

      auto el = std::find_if(uniqueElements.begin(),uniqueElements.end(),
                  [&](Atom &a){return a.atomicNumber == curAtomicNumber;});

      if( el == uniqueElements.end() )
        CErr("Whoops! Can't find atom in unique list! Contact a developer!");

      mapAtom2Uniq[iAtm] = std::distance(uniqueElements.begin(),el);
    }

    if( printLevel > 0 )
      std::cout << "  *** Running " << uniqueElements.size()
                << " Atomic SCF calculations for SAD Guess ***\n\n";


    if( MPIRank(comm) == 0 )
    for(auto iUn = 0; iUn < uniqueElements.size(); iUn++) {

      // Get default multiplicity for the Atom
      size_t defaultMultip;
      try {
        defaultMultip = defaultMult[uniqueElements[iUn].atomicNumber];
      } catch(...) {
        CErr("AtomZ = " + std::to_string(uniqueElements[iUn].atomicNumber)
             + " not supported for SAD Guess");
      }


      std::string multipName;
      switch( defaultMultip ) {

        case 1: multipName = "Singlet"; break;
        case 2: multipName = "Doublet"; break;
        case 3: multipName = "Triplet"; break;
        case 4: multipName = "Quadruplet"; break;
        case 5: multipName = "Quintuplet"; break;
        case 6: multipName = "Sextuplet"; break;
        case 7: multipName = "Septuplet"; break;
        case 8: multipName = "Octuplet"; break;

        default: multipName = "UNKNOWN"; break;

      }


      if( printLevel > 0 )
        std::cout << "    * Running AtomZ = "
                  << uniqueElements[iUn].atomicNumber << " as a "
                  << multipName << std::endl;

      BASIS_FUNCTION_TYPE basisType = this->basisSet().basisType;

      Molecule atom(0,defaultMultip,{ uniqueElements[iUn] });
      BasisSet basis(this->basisSet().basisName,
        this->basisSet().basisDef, this->basisSet().inputDef,
        atom, basisType, this->basisSet().forceCart, false);

      std::shared_ptr<Integrals<IntsT>> aointsAtom =
          std::make_shared<Integrals<IntsT>>();
      aointsAtom->TPI = std::make_shared<InCore4indexTPI<IntsT>>(
          basis.nBasis);

      SingleSlaterOptions guessSSOptions(ssOptions);
      guessSSOptions.refOptions.iCS = defaultMultip == 1;

      std::shared_ptr<SingleSlater<MatsT,IntsT>> ss =
          std::dynamic_pointer_cast<SingleSlater<MatsT,IntsT>>(
              guessSSOptions.buildSingleSlater(std::cout, 
                  atom, basis, aointsAtom));

      ss->comm = rcomm;
      ss->printLevel = 0;
      ss->scfControls.scfAlg = _CONVENTIONAL_SCF;
      ss->scfControls.doIncFock = false;
      ss->scfControls.dampError = 1e-4;
      ss->scfControls.nKeep     = 8;
      ss->buildOrbitalModifierOptions();

      ss->formCoreH(pert, false);
      aointsAtom->TPI->computeAOInts(basis, atom, pert,
          ELECTRON_REPULSION, guessSSOptions.hamiltonianOptions);


      ss->formGuess(ssOptions);
      ss->runSCF(pert);

      size_t NBbasis = basis.nBasis;

      // Place the atomic densities into the guess density
      // *** This only places it into root MPI process ***

      for(auto iAtm = 0; iAtm < mapAtom2Uniq.size(); iAtm++)
        if( mapAtom2Uniq[iAtm] == iUn ) {
          SetMat('N',NBbasis,NBbasis,MatsT(1.),
              ss->onePDM->S().pointer(),NBbasis,
              this->onePDM->S().pointer() + this->basisSet().mapCen2BfSt[iAtm]*(1+NB), NB);
        }

    }

    // Spin-Average the SAD density (on MPI root)
    if( this->onePDM->hasZ() and MPIRank(comm) == 0 ) {
      this->onePDM->Z() = (MatsT(this->nOA - this->nOB) / MatsT(this->nO)) * this->onePDM->S();
    }


#ifdef CQ_ENABLE_MPI
    // Broadcast the 1PDM to all MPI processes
    if( MPISize(comm) > 1 ) {
      std::cerr  << "  *** Scattering the SAD Density ***\n";
      for(MatsT *mat : this->onePDM->SZYXPointers())
        MPIBCast(mat,NB*NB,0,comm);
    }

    // Free the ROOT communicator
    if( rcomm != MPI_COMM_NULL) MPICommFree(rcomm);
#endif


    if( printLevel > 0 )
      std::cout << std::endl
                << "  *** Forming Initial Fock Matrix from SAD Density ***\n\n";

    ao2orthoDen();
    //computeNaturalOrbitals();

  }; // SingleSlater<T>::SADGuess



  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::RandomGuess() {

    size_t NB = this->fockMatrix->dimension();

    // Set up random number generator
    std::random_device rd;
    std::mt19937 gen(rd());
//  std::uniform_real_distribution<> dis(-5,5);
    std::normal_distribution<> dis(0,5);

    // Form the random Fock matrix on the root MPI process
    if( MPIRank(comm) == 0 ) {

      // Copy over the Core Hamiltonian
      *fockMatrix = *coreH;

      // Randomize the Fock matricies
      for(auto F : this->fockMatrix->SZYXPointers()) {
        for(auto k = 0ul; k < NB*NB; k++) F[k] += dis(gen);
        HerMat('L',NB,F,NB);
      }

    }

#ifdef CQ_ENABLE_MPI
    // BCast random fock to all MPI processes
    if( MPISize(comm) > 1 ) {
      std::cerr  << "  *** Scattering the RANDOM Fock ***\n";
      for(auto mat : this->fockMatrix->SZYXPointers())
        MPIBCast(mat,NB*NB,0,comm);
    }
#endif

    // Forming Orbitals from random guess
    EMPerturbation pert; // Dummy EM perturbation
    getNewOrbitals();

  }

  /**
   *  \brief Reads in 1PDM from bin file.
   *
   **/
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::ReadGuess1PDM() {

    //Check if 1PDM comes from save file or scratch file
    if( MPIRank(comm) == 0 ) {

      if( scrBinFileName.empty() ){

        if( printLevel > 0 )
          std::cout << "    * Reading in guess density (restart file) from file "
            << savFile.fName() << "\n";

        readSameTypeDenBin();

      } else {

        if( printLevel > 0 )
          std::cout << "    * Reading in guess density (scratch file) from file "
            << scrBinFileName << "\n";

        readDiffTypeDenBin(scrBinFileName);

        if( printLevel > 0 )
          std::cout << "    * Saving prepared 1-PDMs to file "
            << savFile.fName() << "\n";

        // Saving post-transformed 1-PDMs to restart file
        if( savFile.exists() ) {

          std::string prefix = "SCF/";
          if( this->particle.charge == 1.0 ) prefix = "PROT_" + prefix;

          savFile.safeWriteData(prefix + "1PDM", *this->onePDM);

        }

      }

    }

    // dimension of 1PDM
    auto NB = basisSet().nBasis;
    if( this->nC == 4 ) NB=2*NB;
    auto NB2 = NB*NB;

#ifdef CQ_ENABLE_MPI
    // BCast onePDM to all MPI processes
    if( MPISize(comm) > 1 ) {
      std::cerr  << "  *** Scattering the onePDM ***\n";
      for(auto mat : this->onePDM->SZYXPointers())
        MPIBCast(mat,NB2,0,comm);
    }
#endif

    ao2orthoDen();
    //computeNaturalOrbitals();

  } // SingleSlater<T>::ReadGuess1PDM()

  /**
   *  \brief Reads in 1PDM from bin file
   *  of the same type as calculation.
   *
   **/
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::readSameTypeDenBin() {

    if( MPIRank(comm) == 0 ) {

      size_t t_hash = std::is_same<MatsT,double>::value ? 1 : 2;
      size_t d_hash = 1;
      size_t c_hash = 2;

      size_t savHash;
      std::string prefix = "/SCF/";
      if (this->particle.charge == 1.0)
        prefix = "/PROT_SCF/";

      try{
        savFile.readData(prefix + "FIELD_TYPE", &savHash);
      } catch (...) {
        CErr("Cannot find " + prefix + "FIELD_TYPE on rstFile!",std::cout);
      }


      if( t_hash != savHash ) {

        bool t_is_double  = t_hash == d_hash;
        bool t_is_complex = t_hash == c_hash;

        bool s_is_double  = savHash == d_hash;
        bool s_is_complex = savHash == c_hash;

        std::string t_field = t_is_double ? "REAL" : "COMPLEX";
        std::string s_field = s_is_double ? "REAL" : "COMPLEX";

        std::string message = prefix + "FIELD_TYPE on disk (" + s_field +
          ") is incompatible with current FIELD_TYPE (" + t_field + ")";

        CErr(message,std::cout);
      }


      // dimension of 1PDM
      auto NB = basisSet().nBasis;
      if( this->nC == 4 ) NB=2*NB;
      auto NB2 = NB*NB;

      auto DSdims = savFile.getDims( prefix + "1PDM_SCALAR" );
      auto DZdims = savFile.getDims( prefix + "1PDM_MZ" );
      auto DYdims = savFile.getDims( prefix + "1PDM_MY" );
      auto DXdims = savFile.getDims( prefix + "1PDM_MX" );
  
      bool hasDS = DSdims.size() != 0;
      bool hasDZ = DZdims.size() != 0;
      bool hasDY = DYdims.size() != 0;
      bool hasDX = DXdims.size() != 0;

      bool r2DS = DSdims.size() == 2;
      bool r2DZ = DZdims.size() == 2;
      bool r2DY = DYdims.size() == 2;
      bool r2DX = DXdims.size() == 2;


      // Errors in 1PDM SCALAR
      if( not hasDS )
        CErr(prefix+"1PDM_SCALAR does not exist in " + savFile.fName(), std::cout);

      else if( not r2DS )
        CErr(prefix + "1PDM_SCALAR not saved as a rank-2 tensor in " +
            savFile.fName(), std::cout);

      else if( DSdims[0] != NB or DSdims[1] != NB ) {

        std::cout << "    * Incompatible " + prefix + "1PDM_SCALAR:";
        std::cout << "  Recieved (" << DSdims[0] << "," << DSdims[1] << ")"
          << " :";
        std::cout << "  Expected (" << NB << "," << NB << ")";
        CErr("Wrong dimension of 1PDM SCALAR!",std::cout);

      }

      // Read in 1PDM SCALAR
      std::cout << "    * Found " + prefix + "1PDM_SCALAR !" << std::endl;
      savFile.readData(prefix + "1PDM_SCALAR",this->onePDM->S().pointer());


      // Oddities in Restricted
      if( this->nC == 1 and this->iCS ) {

        if( hasDZ )
          std::cout << "    * WARNING: Reading in " + prefix + "1PDM_SCALAR as "
            << "restricted guess but " << savFile.fName()
            << " contains " + prefix + "1PDM_MZ" << std::endl;

        if( hasDY )
          std::cout << "    * WARNING: Reading in " + prefix + "1PDM_SCALAR as "
            << "restricted guess but " << savFile.fName()
            << " contains SCF/1PDM_MY" << std::endl;

        if( hasDX )
          std::cout << "    * WARNING: Reading in " + prefix + "1PDM_SCALAR as "
            << "restricted guess but " << savFile.fName()
            << " contains SCF/1PDM_MX" << std::endl;

      }


      // MZ
      if( this->nC == 2 or not this->iCS ) {

        if( not hasDZ ) {

          std::cout <<  "    * WARNING: " + prefix + "1PDM_MZ does not exist in "
            << savFile.fName() << " -- Zeroing out " + prefix + "1PDM_MZ" << std::endl;

          this->onePDM->Z().clear();


        } else if( not r2DZ )
          CErr(prefix + "1PDM_MZ not saved as a rank-2 tensor in " +
              savFile.fName(), std::cout);

        else if( DZdims[0] != NB or DZdims[1] != NB ) {

          std::cout << "    * Incompatible " + prefix + "1PDM_MZ:";
          std::cout << "  Recieved (" << DZdims[0] << "," << DZdims[1] << ")"
            << " :";
          std::cout << "  Expected (" << NB << "," << NB << ")";
          CErr("Wrong dimension of 1PDM MZ!",std::cout);

        } else {

          std::cout << "    * Found " + prefix + "1PDM_MZ !" << std::endl;
          savFile.readData(prefix + "1PDM_MZ",this->onePDM->Z().pointer());

        }

        // Oddities in Unrestricted
        if( this->nC == 2 ) {

          if( hasDY )
            std::cout << "    * WARNING: Reading in " + prefix + "1PDM_MZ as "
              << "unrestricted guess but " << savFile.fName()
              << " contains " + prefix + "1PDM_MY" << std::endl;

          if( hasDX )
            std::cout << "    * WARNING: Reading in " + prefix + "1PDM_MZ as "
              << "unrestricted guess but " << savFile.fName()
              << " contains " + prefix + "1PDM_MX" << std::endl;

        }

      }


      if( this->nC == 2 or this->nC == 4 ) {

        if( not hasDY ) {

          std::cout <<  "    * WARNING: " + prefix + "1PDM_MY does not exist in "
            << savFile.fName() << " -- Zeroing out " + prefix + "1PDM_MY" << std::endl;

          this->onePDM->Y().clear();


        } else if( not r2DY )
          CErr(prefix + "1PDM_MY not saved as a rank-2 tensor in " +
              savFile.fName(), std::cout);

        else if( DYdims[0] != NB or DYdims[1] != NB ) {

          std::cout << "    * Incompatible " + prefix + "1PDM_MY:";
          std::cout << "  Recieved (" << DYdims[0] << "," << DYdims[1] << ")"
            << " :";
          std::cout << "  Expected (" << NB << "," << NB << ")";
          CErr("Wrong dimension of 1PDM MY!",std::cout);

        } else {

          std::cout << "    * Found " + prefix + "1PDM_MY !" << std::endl;
          savFile.readData(prefix + "1PDM_MY",this->onePDM->Y().pointer());

        }


        if( not hasDX ) {

          std::cout <<  "    * WARNING: " + prefix + "1PDM_MX does not exist in "
            << savFile.fName() << " -- Zeroing out " + prefix + "1PDM_MX" << std::endl;

          this->onePDM->X().clear();


        } else if( not r2DX )
          CErr(prefix + "1PDM_MX not saved as a rank-2 tensor in " +
              savFile.fName(), std::cout);

        else if( DXdims[0] != NB or DXdims[1] != NB ) {

          std::cout << "    * Incompatible " + prefix + "1PDM_MX:";
          std::cout << "  Recieved (" << DXdims[0] << "," << DXdims[1] << ")"
            << " :";
          std::cout << "  Expected (" << NB << "," << NB << ")";
          CErr("Wrong dimension of 1PDM MX!",std::cout);

        } else {

          std::cout << "    * Found " + prefix + "1PDM_MX !" << std::endl;
          savFile.readData(prefix + "1PDM_MX",this->onePDM->X().pointer());

        }

      }

    }

  } // SingleSlater<T>::readSameTypeDenBin()

  /**
   *  \brief Reads in 1PDM from bin file
   *  of different type as calculation
   *  and uses it as initial guess.
   *
   **/
  template <typename MatsT, typename IntsT>
  template <typename ScrMatsT>
  void SingleSlater<MatsT,IntsT>::getScr1PDM(SafeFile& scrBin) {

    if( MPIRank(comm) == 0 ) {

      // dimension of 1PDM
      auto NB = basisSet().nBasis;
      if( this->nC == 4 ) NB=2*NB;
      auto NB2 = NB*NB;

      auto DSdims = scrBin.getDims( "SCF/1PDM_SCALAR" );
      auto DZdims = scrBin.getDims( "SCF/1PDM_MZ" );
      auto DYdims = scrBin.getDims( "SCF/1PDM_MY" );
      auto DXdims = scrBin.getDims( "SCF/1PDM_MX" );

      bool hasDS = DSdims.size() != 0;
      bool hasDZ = DZdims.size() != 0;
      bool hasDY = DYdims.size() != 0;
      bool hasDX = DXdims.size() != 0;

      bool r2DS = DSdims.size() == 2;
      bool r2DZ = DZdims.size() == 2;
      bool r2DY = DYdims.size() == 2;
      bool r2DX = DXdims.size() == 2;

      int scrRefType, binRefType;
      scrBin.readData("REF/REFTYPE",&scrRefType);
      savFile.readData("REF/REFTYPE",&binRefType);

      std::cout << "    * Converting from " << refMap[scrRefType] << " to "
        << refMap[binRefType] << std::endl;

      // onePDM on scr bin file
      // assume square and same dimension between S,X,Y,Z
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<ScrMatsT>> onePDMtmp;
      onePDMtmp = std::make_shared<cqmatrix::PauliSpinorMatrices<ScrMatsT>>(DSdims[0],hasDY,hasDZ);

      // Errors in 1PDM SCALAR
      if( not hasDS )
        CErr("SCF/1PDM_SCALAR does not exist in " + scrBin.fName(), std::cout);

      else if( not r2DS )
        CErr("SCF/1PDM_SCALAR not saved as a rank-2 tensor in " +
            scrBin.fName(), std::cout);

      // Read in 1PDM SCALAR
      std::cout << "    * Looking for SCF/1PDM_SCALAR !" << std::endl;
      scrBin.readData("/SCF/1PDM_SCALAR",onePDMtmp->S().pointer());

      // MZ
      if( onePDMtmp->hasZ() ){

        std::cout << "    * Looking for SCF/1PDM_MZ !" << std::endl;
        if( not r2DZ )
          CErr("SCF/1PDM_MZ not saved as a rank-2 tensor in " +
            scrBin.fName(), std::cout);
        scrBin.readData("SCF/1PDM_MZ",onePDMtmp->Z().pointer());

      }

      // MY
      if( onePDMtmp->hasXY() ){

        std::cout << "    * Looking for SCF/1PDM_MX !" << std::endl;
        if( not r2DX )
          CErr("SCF/1PDM_MX not saved as a rank-2 tensor in " +
            scrBin.fName(), std::cout);
        scrBin.readData("SCF/1PDM_MX",onePDMtmp->X().pointer());

        std::cout << "    * Looking for SCF/1PDM_MY !" << std::endl;
        if( not r2DY )
          CErr("SCF/1PDM_MY not saved as a rank-2 tensor in " +
            scrBin.fName(), std::cout);
        scrBin.readData("SCF/1PDM_MY",onePDMtmp->Y().pointer());

      }

      // Initialize onePDM
      auto scr1PDMSize = onePDMtmp->dimension();
      // Guess 1PDM same size as calculation 1PDM
      if( scr1PDMSize == NB ) *this->onePDM = *onePDMtmp;
      // Guess 1PDM smaller than 1PDM
      else if( scr1PDMSize < NB ){
        auto p1Comps = this->onePDM->SZYXPointers();
        auto p2Comps = onePDMtmp->SZYXPointers();
        auto nComp = p1Comps.size();
        auto n2Comp = p2Comps.size();
        for( auto iComp=0; iComp<nComp; iComp++ ){
          if( iComp < n2Comp )
            SetMat('N',scr1PDMSize,scr1PDMSize,MatsT(1.),
               p2Comps[iComp],scr1PDMSize,p1Comps[iComp],NB);
        }
      } else CErr("Cannot use a guess of larger size.");

      std::cout << "\n" << std::endl;
      onePDMtmp = nullptr;

    }

  } // SingleSlater<T>::getScr1PDM()

  template <>
  template <>
  void SingleSlater<double,double>::getScr1PDM<dcomplex>(SafeFile& scrBin) {

    CErr("Cannot do complex guess density for real calculation.");

  }

  template <>
  template <>
  void SingleSlater<double,dcomplex>::getScr1PDM<dcomplex>(SafeFile& scrBin) {

    CErr("Cannot do complex guess density for real calculation.");

  }

  /**
   *  \brief Reads in 1PDM from bin file
   *  of different type as calculation.
   *
   **/
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::readDiffTypeDenBin(std::string binName) {

    if( MPIRank(comm) == 0 ) {

      bool scrBinExists;
      SafeFile binFile(binName, scrBinExists);

      size_t t_hash = std::is_same<MatsT,double>::value ? 1 : 2;
      size_t d_hash = 1;
      size_t c_hash = 2;
      size_t savHash;

      std::string prefix = "/SCF/";
      if (this->particle.charge == 1.0)
        prefix = "/PROT_SCF/";

      try{
        binFile.readData("/SCF/FIELD_TYPE", &savHash);
      } catch (...) {
        CErr("Cannot find /SCF/FIELD_TYPE on rstFile!",std::cout);
      }

      // type of 1PDM
      bool t_is_double  = t_hash == d_hash;
      bool t_is_complex = t_hash == c_hash;

      bool s_is_double  = savHash == d_hash;
      bool s_is_complex = savHash == c_hash;

      std::string t_field = t_is_double ? "REAL" : "COMPLEX";
      std::string s_field = s_is_double ? "REAL" : "COMPLEX";

      std::string message = "    * Going from /SCF/FIELD_TYPE on disk (" + s_field +
        ") to current FIELD_TYPE (" + t_field + ")";

      std::cout << message << std::endl;

      // Determine storage of 1PDM on scr bin file
      // Assumes scalar 1PDM is same size as MX, MY, MZ
      // Assumes square 1PDM
      if( s_is_double ){

        getScr1PDM<double>(binFile);

      } else if( s_is_complex ){

        getScr1PDM<dcomplex>(binFile);

      } else CErr("Could not determine type of scratch bin file");

    }

  } // SingleSlater<T>::readDiffTypeDenBin()


  /**
   *  \brief Reads in MOs from bin file.
   *
   **/
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::ReadGuessMO() {

    //Check if MOs come from save file or scratch file
    if( MPIRank(comm) == 0 ) {

      if( scrBinFileName.empty() ){

        if( printLevel > 0 )
          std::cout << "    * Reading in guess MOs (restart file) from file "
            << savFile.fName() << "\n";

        readSameTypeMOBin();

      } else {

        if( printLevel > 0 )
          std::cout << "    * Reading in guess MOs (scratch file) from file "
            << scrBinFileName << "\n";

        readDiffTypeMOBin(scrBinFileName);

        if( printLevel > 0 )
          std::cout << "    * Saving prepared MOs to file "
            << savFile.fName() << "\n";

        // Saving post-transformed MOs to restart file
        if( savFile.exists() ) {

          size_t NB  = this->nAlphaOrbital();
          size_t NBC = this->nC * NB;
          std::string prefix = "SCF/";
          if( this->particle.charge == 1.0 ) prefix = "PROT_" + prefix;

          savFile.safeWriteData(prefix + "MO1", this->mo[0].pointer(), {NBC, NBC});
          if (this->nC == 1 and not this->iCS) savFile.safeWriteData(prefix + "MO2", this->mo[1].pointer(), {NBC, NBC});

        }

      }

    }

    // MO coefficients from AO to othonormalized basis
    orthoAOMO();

    // MO swapping if requested
    if( this->moPairs[0].size() != 0 ){

      this->swapMOs(this->moPairs,SpinType::isAlpha);

      if( printLevel > 0 )
        std::cout << "    * Saving swapped MOs to file "
          << savFile.fName() << "\n";

      // Saving post-transformed MOs to restart file
      if( savFile.exists() ) {

        size_t NB  = this->nAlphaOrbital();
        size_t NBC = this->nC * NB;
        std::string prefix = "SCF/";
        if( this->particle.charge == 1.0 ) prefix = "PROT_" + prefix;

        savFile.safeWriteData(prefix + "MO1", this->mo[0].pointer(), {NBC, NBC});

      }

    }
    if( this->moPairs[1].size() != 0 ){

      this->swapMOs(this->moPairs,SpinType::isBeta);

      if( printLevel > 0 )
        std::cout << "    * Saving swapped beta MOs to file "
          << savFile.fName() << "\n";

      // Saving post-transformed MOs to restart file
      if( savFile.exists() ) {

        size_t NB  = this->nAlphaOrbital();
        size_t NBC = this->nC * NB;
        std::string prefix = "SCF/";
        if( this->particle.charge == 1.0 ) prefix = "PROT_" + prefix;

        savFile.safeWriteData(prefix + "MO2", this->mo[1].pointer(), {NBC, NBC});

      }

    }

    // Form density from MOs
    formDensity();

  }

  /**
   *  \brief Reads in MOs from bin file.
   *
   **/
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::readSameTypeMOBin() {

    if( MPIRank(comm) == 0 ) {

      size_t t_hash = std::is_same<MatsT,double>::value ? 1 : 2;
      size_t d_hash = 1;
      size_t c_hash = 2;

      size_t savHash;

      std::string prefix = "/SCF/";
      if (this->particle.charge == 1.0)
        prefix = "/PROT_SCF/";

      try{
        savFile.readData(prefix + "FIELD_TYPE", &savHash);
      } catch (...) {
        CErr("Cannot find /SCF/FIELD_TYPE on rstFile!",std::cout);
      }


      if( t_hash != savHash ) {

        bool t_is_double  = t_hash == d_hash;
        bool t_is_complex = t_hash == c_hash;

        bool s_is_double  = savHash == d_hash;
        bool s_is_complex = savHash == c_hash;

        std::string t_field = t_is_double ? "REAL" : "COMPLEX";
        std::string s_field = s_is_double ? "REAL" : "COMPLEX";

        std::string message = prefix + "FIELD_TYPE on disk (" + s_field +
          ") is incompatible with current FIELD_TYPE (" + t_field + ")";

        CErr(message,std::cout);
      }

      // dimension of mo1 and mo2
      auto NB = this->nC * this->nAlphaOrbital();
      auto NB2 = NB*NB;

      auto MO1dims = savFile.getDims( prefix + "MO1" );
      auto MO2dims = savFile.getDims( prefix + "MO2" );


      // Find errors in MO1
      if( MO1dims.size() == 0 )
        CErr(prefix + "MO1 does not exist in " + savFile.fName(), std::cout);

      if( MO1dims.size() != 2 )
        CErr(prefix + "MO1 not saved as a rank-2 tensor in " + savFile.fName(),
            std::cout);

      if( MO1dims[0] != NB or MO1dims[1] != NB ) {

        std::cout << "    * Incompatible SCF/MO1:";
        std::cout << "  Recieved (" << MO1dims[0] << "," << MO1dims[1] << ")"
          << " :";
        std::cout << "  Expected (" << NB << "," << NB << ")";
        CErr("Wrong number of MO coefficients!",std::cout);

      }



      // MO2 + RHF is odd, print warning
      if( MO2dims.size() != 0 and this->nC == 1 and this->iCS )
        std::cout << "    * WARNING: Reading in SCF/MO1 as restricted guess "
                  << "but " << savFile.fName() << " contains SCF/MO2"
                  << std::endl;


      // Read in MO1
      std::cout << "    * Found SCF/MO1 !" << std::endl;
      savFile.readData(prefix + "MO1",this->mo[0].pointer());

      // Unrestricted calculations
      if( this->nC == 1 and not this->iCS ) {

        if( MO2dims.size() == 0 )
          std::cout << "    * WARNING: SCF/MO2 does not exist in "
            << savFile.fName() << " -- Copying SCF/MO1 -> SCF/MO2 " << std::endl;

        if( MO2dims.size() > 2  )

          CErr("SCF/MO2 not saved as a rank-2 tensor in " + savFile.fName(),
              std::cout);

        else if( MO2dims[0] != NB or MO2dims[1] != NB ) {

          std::cout << "    * Incompatible SCF/MO2:";
          std::cout << "  Recieved (" << MO2dims[0] << "," << MO2dims[1] << ")"
            << " :";
          std::cout << "  Expected (" << NB << "," << NB << ")";
          CErr("Wrong number of MO coefficients!",std::cout);

        }


        // Read in MO2
        if( MO2dims.size() == 0 )
          this->mo[1] = this->mo[0];
        else {
          std::cout << "    * Found SCF/MO2 !" << std::endl;
          savFile.readData(prefix + "MO2",this->mo[1].pointer());
        }

      }

    }

  } // SingleSlater<T>::readSameTypeMOBin()

  /**
   *  \brief Reads in MOs from bin file.
   *
   **/
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::readDiffTypeMOBin(std::string binName) {

    if( MPIRank(comm) == 0 ) {

      bool scrBinExists;
      SafeFile binFile(binName, scrBinExists);

      size_t t_hash = std::is_same<MatsT,double>::value ? 1 : 2;
      size_t d_hash = 1;
      size_t c_hash = 2;

      size_t savHash;

      std::string prefix = "/SCF/";
      if (this->particle.charge == 1.0)
        prefix = "/PROT_SCF/";

      try{
        binFile.readData(prefix + "FIELD_TYPE", &savHash);
      } catch (...) {
        CErr("Cannot find /SCF/FIELD_TYPE on rstFile!",std::cout);
      }

      // type of MO
      bool t_is_double  = t_hash == d_hash;
      bool t_is_complex = t_hash == c_hash;

      bool s_is_double  = savHash == d_hash;
      bool s_is_complex = savHash == c_hash;

      std::string t_field = t_is_double ? "REAL" : "COMPLEX";
      std::string s_field = s_is_double ? "REAL" : "COMPLEX";

      std::string message = "    * Going from /SCF/FIELD_TYPE on disk (" + s_field +
        ") to current FIELD_TYPE (" + t_field + ")";

      std::cout << message << std::endl;

      // Determine storage of MOs on scr bin file
      if( s_is_double ){

        getScrMO<double>(binFile);

      } else if( s_is_complex ){

        getScrMO<dcomplex>(binFile);

      } else CErr("Could not determine type of scratch bin file");

    }

  }

  /**
   *  \brief Reads in MOs from fchk file.
   *
   **/
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::FchkGuessMO() {

    if( printLevel > 0 ){
      std::cout << "    * Reading in guess orbitals from file " << fchkFileName << "\n";
      std::cout << "      Please check that IOp(3/60=-1) was included in your Gaussian calculation." << "\n";
      std::cout << "      This functionality works best when the number of primitives is the same." << "\n";
    }

    if( this->nC == 4 ) CErr("FCHKMO NYI for 4c",std::cout);

    std::vector<int> shellList;

    // Parsing fchk file and populating mo1(mo2)
    // Outputs shell listing needed for angular momentum reorganization
    shellList=fchkToCQMO();

    // Dimension of mo1 and mo2
    auto NB = this->nC * this->basisSet().nBasis;
    auto NB2 = NB*NB;

    // Reorders shells of mo1 to Chronus ordering
    MatsT* mo1tmp = CQMemManager::get().malloc<MatsT>(NB2);
    SetMat('N',NB,NB,MatsT(1.),this->mo[0].pointer(),NB,mo1tmp,NB);
    reorderAngMO(shellList,mo1tmp,0);
    CQMemManager::get().free(mo1tmp);

    // Reorders shells of mo2 to Chronus ordering
    if( this->nC == 1 and not this->iCS){
      MatsT* mo2tmp = CQMemManager::get().malloc<MatsT>(NB2);
      SetMat('N',NB,NB,MatsT(1.),this->mo[1].pointer(),NB,mo2tmp,NB);
      reorderAngMO(shellList,mo2tmp,1);
      CQMemManager::get().free(mo2tmp);
    }

    // Reorder spin components
    if( this->nC == 2 ) reorderSpinMO();

    // Orthogonalize MOs
    orthoAOMO();

    // Form density from MOs
    formDensity();

  } // SingleSlater<T>::FchkGuessMO()

  /**
   *  \brief For NEO Protons Only: 
   *         Occupy the tightest orbital for each quantum pron
   **/
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::NEOTightProtonGuess(){
    
    if(this->particle.charge<0)
      CErr("NEOTightProtonGuess Doesn't Apply to Electronic Wavefunction");

    size_t numProt = this->nOA, NB = this->basisSet().nBasis, NB_per_prot = NB / numProt;
    
    this->onePDM->S().clear();
    this->mo[0].clear();
    this->mo[1].clear();
    for(int i = 0; i < numProt; i++) {
      this->onePDM->S()(i*NB_per_prot, i*NB_per_prot) = 1;
      this->mo[0](i*NB_per_prot, i) = 1;
    }
    this->onePDM->Z() = this->onePDM->S();

    std::cout << "      Each quantum proton occupies the tightest orbital. " << std::endl;
    std::cout << std::endl;

  } // SingleSlater<MatsT,IntsT>::NEOTightProtonGuess()

  /**
   *  \brief For NEO Electronic Wavefunction Only: 
   *         Converge a classical SCF, then use converged density as guess for NEO electronic subsystem
   **/
  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT>::NEOConvergeClassicalGuess(const SingleSlaterOptions &ssOptions){


    std::cout << "    * Converging a classical SCF calculation " << std::endl;
    std::cout << "      The converged density will be used as the guess the electronic subsystem " << "\n" << std::endl;

    size_t NB = this->basisSet().nBasis;
    EMPerturbation pert;

    // Zero out the densities (For all MPI processes)
    this->onePDM->clear();

    SingleSlaterOptions classicalSSOptions(ssOptions);

    // Create a temporary molecule and set quantum atoms to be classical
    Molecule tempMol(this->molecule());
    for ( Atom& atom : tempMol.atoms ) if (atom.quantum) atom.quantum = false;
    tempMol.update();

    // Grab the pointer of the pre-built and computed aoints 
    std::shared_ptr<Integrals<IntsT>> tempaoints = this->aoints_;

    // Create temporary singleslater object to converge a classical calculation
    std::shared_ptr<SingleSlater<MatsT,IntsT>> classicalSS =
        std::dynamic_pointer_cast<SingleSlater<MatsT,IntsT>>(
            classicalSSOptions.buildSingleSlater(std::cout,  tempMol, 
                this->basisSet(), tempaoints));

    classicalSS->printLevel = 1;
    classicalSS->scfControls.scfAlg = _CONVENTIONAL_SCF;
    classicalSS->scfControls.guess =   SAD;
    classicalSS->scfControls.diisAlg =   CDIIS;
    classicalSS->buildOrbitalModifierOptions();

    classicalSS->formCoreH(pert, false);
    classicalSS->formGuess(classicalSSOptions);
    classicalSS->formFock(pert, false);
    classicalSS->runSCF(pert);

    std::cout << "\nClassical SCF Calculation Converged!" << std::endl;
    std::cout << "\nWill use this converged density as the guess density\n"  
        << "for the electronic subsystem in subsequent NEO calculation\n" << std::endl;

    // Place the converged SCF densities into the guess density
    // *** This only places it into root MPI process ***
    SetMat('N',NB,NB,MatsT(1.),classicalSS->onePDM->S().pointer(),NB,
      this->onePDM->S().pointer(),NB);

    if( this->onePDM->hasZ() ) 
      this->onePDM->Z() = (MatsT(this->nOA - this->nOB) / MatsT(this->nO)) * this->onePDM->S();

    ao2orthoDen();
    computeNaturalOrbitals();

  } // SingleSlater<MatsT,IntsT>::NEOConvergeClassicalGuess

  /*
   * Brief: Computes the Natural orbitals from the orthogonal
   *        density
   */
  template<typename MatsT,typename IntsT>
  void SingleSlater<MatsT,IntsT> :: computeNaturalOrbitals() {

    if( MPIRank(comm) == 0) {

    //if( printLevel > 0 ) std::cout << "  *** Computing Natural Orbitals from Guess Density ***" << std::endl << std::endl;

    size_t NBC = this->nC*basisSet().nBasis;
    bool iRO = (std::dynamic_pointer_cast<ROFock<MatsT, IntsT>>(this->fockBuilder) != nullptr);

    if( this->nC == 1 ){
      // Allocate Local Matrices
      std::vector<cqmatrix::Matrix<MatsT>> SCR = this->onePDMOrtho->template spinGatherToBlocks<MatsT>(false);
      double* eVals = CQMemManager::get().malloc<double>(NBC);


      // Diagonalize Density
      int INFO  = HermetianEigen('V', 'L', NBC, SCR[0].pointer(), NBC, eVals);
      if( INFO != 0 )
        CErr("HermetianEigen failed in computing Natural Orbitals", std::cout);

      // Copy in reverse order to MO's
      // Because the highest occupation numbers are last
      for( size_t i=0; i<NBC; i++ ){
         size_t disp = ((NBC - 1) - i)*NBC;
         std::copy_n(SCR[0].pointer()+disp,NBC,this->mo[0].pointer()+i*NBC);
      }
      if( iRO ) std::copy_n(this->mo[0].pointer(),NBC*NBC,this->mo[1].pointer());

#ifdef _SINGLESLATER_NATURAL_ORBITALS
      std::cout << "Alpha Natural Orbital Occupations:" << std::endl;
      for( size_t i=0; i<NBC; i++ )
        std::cout << "  " << eVals[i] << std::endl;
      prettyPrintSmart(std::cout, "Natural Orbitals (Alpha)", this->mo[0].pointer(),NBC,NBC,NBC);
#endif

      // Compute Beta Orbitals for Unrestricted
      if( not (this->iCS) ){
        INFO  = HermetianEigen('V', 'L', NBC, SCR[1].pointer(), NBC, eVals);
        if( INFO != 0 )
          CErr("HermetianEigen failed in computing Natural Orbitals", std::cout);

        // Copy in reverse order
        for( size_t i=0; i<NBC; i++ ){
           size_t disp = ((NBC - 1) - i)*NBC;
           std::copy_n(SCR[1].pointer()+disp,NBC,this->mo[1].pointer()+i*NBC);
        }

#ifdef _SINGLESLATER_NATURAL_ORBITALS
      std::cout << "Beta Natural Orbital Occupations:" << std::endl;
      for( size_t i=0; i<NBC; i++ )
        std::cout << "  " << eVals[i] << std::endl;
      prettyPrintSmart(std::cout, "Natural Orbitals (Beta)", this->mo[1].pointer(),NBC,NBC,NBC);
#endif
      }
      CQMemManager::get().free(eVals);

    // 2C and 4C
    } else {

      cqmatrix::Matrix<MatsT> SCR = this->onePDMOrtho->template spinGather<MatsT>();
      double* eVals = CQMemManager::get().malloc<double>(NBC);

      // Diagonalize Density
      int INFO  = HermetianEigen('V', 'L', NBC, SCR.pointer(), NBC, eVals);
      if( INFO != 0 )
        CErr("HermetianEigen failed in computing Natural Orbitals", std::cout);

      // Copy in reverse order to MOs
      if( this->nC == 4 ){
        size_t nPos = NBC/2;
        // Copy Positive energy orbitals
        for( size_t i=0; i<nPos; ++i ){
            size_t disp = ((NBC-1) - i)*NBC;
           std::copy_n(SCR.pointer()+disp,NBC,this->mo[0].pointer()+(i+nPos)*NBC);
        }
        // Copy Negative energy orbitals
        for( size_t i=nPos; i<NBC; ++i ){
          size_t disp = ((NBC-1) - i)*NBC;
          std::copy_n(SCR.pointer()+disp,NBC,this->mo[0].pointer()+(i-nPos)*NBC);
        }
      } else {
        // 2C Copy natural orbitals
        for( size_t i=0; i<NBC; ++i ){
           size_t disp = ((NBC - 1) - i)*NBC;
           std::copy_n(SCR.pointer()+disp,NBC,this->mo[0].pointer()+i*NBC);
        }
      }

#ifdef _SINGLESLATER_NATURAL_ORBITALS
      std::cout << "Natural Orbital Occupations:" << std::endl;
      for( size_t i=0; i<NBC; i++ )
        std::cout << "  " << eVals[i] << std::endl;
      prettyPrintSmart(std::cout, "Natural Orbitals", this->mo[0].pointer(),NBC,NBC,NBC);
#endif
      CQMemManager::get().free(eVals);
    }
    } // MPIRank == 0
    ortho2aoMOs();
  }

  template <typename MatsT, typename IntsT>
  void SingleSlater<MatsT,IntsT> :: getNewOrbitals() {
      ao2orthoFock();
      diagOrthoFock();
      ortho2aoMOs();
      formDensity();
      saveCurrentState();
  }

}; // namespace ChronusQ

