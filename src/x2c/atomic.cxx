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
#include <corehbuilder/x2c/atomic.hpp>
#include <corehbuilder/nonrel.hpp>
#include <particleintegrals/onepints/relativisticints.hpp>
#include <matrix.hpp>
#include <cqlinalg.hpp>
#include <physcon.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  AtomicX2C<MatsT,IntsT>::AtomicX2C(const AtomicX2C<MatsT,IntsT> &other) :
    AtomicX2C(other,0) {}

  template <typename MatsT, typename IntsT>
  AtomicX2C<MatsT,IntsT>::AtomicX2C(AtomicX2C<MatsT,IntsT> &&other) :
    AtomicX2C(std::move(other),0) {}

  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  AtomicX2C<MatsT,IntsT>::AtomicX2C(const AtomicX2C<MatsU,IntsT> &other, int dummy) :
    X2C<MatsT,IntsT>(other),
    type_(other.type_), atomIdx_(other.atomIdx_) {
    for (const auto& atom : other.atoms_)
      atoms_.emplace_back(atom);
  }

  template <typename MatsT, typename IntsT>
  template <typename MatsU>
  AtomicX2C<MatsT,IntsT>::AtomicX2C(AtomicX2C<MatsU,IntsT> &&other, int dummy) :
    X2C<MatsT,IntsT>(other),
    type_(other.type_), atomIdx_(other.atomIdx_) {
    for (const auto& atom : other.atoms_)
      atoms_.emplace_back(atom);
  }

  template <typename MatsT, typename IntsT>
  void AtomicX2C<MatsT,IntsT>::dealloc() {
    X2C<MatsT,IntsT>::dealloc();
  }

  /**
   *  \brief Compute the AtomicX2C Core Hamiltonian
   */
  template <typename MatsT, typename IntsT>
  void AtomicX2C<MatsT, IntsT>::computeOneEX2C(EMPerturbation &emPert,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> coreH) {
    size_t NP = this->uncontractedBasis_.nPrimitive;
    size_t NB = this->basisSet_.nBasis;

    size_t Natom = this->molecule_.atoms.size();

    std::map<size_t, std::vector<size_t>> uniqueElements;
    std::vector<size_t> cumeNBs;

    if (type_.diagonalOnly)
      NRCoreH<MatsT, IntsT>(this->aoints_, this->ssOptions_.hamiltonianOptions)
          .computeNRCH(emPert, coreH);

    std::vector<MatsT*> CH(coreH->SZYXPointers());

    atoms_.clear();
    atoms_.reserve(Natom);

    size_t maxAtomNB = 0, maxAtomNP = 0;

    size_t cumeNB = 0;
    for (size_t i = 0; i < Natom; i++) {

      Atom atom(this->molecule_.atoms[i]);
      if (type_.isolateAtom) {
        size_t aN = atom.atomicNumber;
        uniqueElements[aN].push_back(i);
        if (uniqueElements[aN][0] == i) {
          atomIdx_.push_back(atoms_.size());
          atom.coord = {0., 0., 0.};
        } else {
          atomIdx_.push_back(atomIdx_[uniqueElements[aN][0]]);
          cumeNBs.push_back(cumeNB);
          cumeNB += atoms_[atomIdx_.back()].basisSet_.nBasis;
          continue;
        }
      } else {
        uniqueElements[i].push_back(i);
        atomIdx_.push_back(i);
      }

      Molecule atomMol(0, atom.atomicNumber % 2 + 1, { atom });
      BasisSet basis(this->basisSet_.basisName,
        this->basisSet_.basisDef, this->basisSet_.inputDef,
        atomMol, this->basisSet_.basisType,
        this->basisSet_.forceCart, false);

      cumeNBs.push_back(cumeNB);
      cumeNB += basis.nBasis;
      maxAtomNB = std::max(basis.nBasis, maxAtomNB);
      maxAtomNP = std::max(basis.nPrimitive, maxAtomNP);

      std::shared_ptr<Integrals<IntsT>> aointsAtom;
      if (type_.isolateAtom) {
        aointsAtom = std::make_shared<Integrals<IntsT>>();
        atoms_.emplace_back(*aointsAtom, atomMol, basis, this->ssOptions_);
      } else {
        aointsAtom = std::make_shared<Integrals<IntsT>>();
        atoms_.emplace_back(*aointsAtom, this->molecule_, basis, this->ssOptions_);
      }

    }

    for (size_t k = 0; k < atoms_.size(); k++) {
      size_t atomNB = atoms_[k].basisSet_.nBasis;
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<MatsT>> atomCoreH;
      if (coreH->hasXY())
        atomCoreH = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(atomNB, true);
      else if (coreH->hasZ())
        atomCoreH = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(atomNB, false);
      else
        atomCoreH = std::make_shared<cqmatrix::PauliSpinorMatrices<MatsT>>(atomNB, false, false);
      atomCoreH->clear();
      if (type_.diagonalOnly) {
        atoms_[k].computeOneEX2C_corr(emPert, atomCoreH);
        size_t aN = k;
        if (type_.isolateAtom)
          aN = atoms_[k].molecule_.atoms[0].atomicNumber;
        std::vector<MatsT*> atomCH(atomCoreH->SZYXPointers());
        for (size_t i : uniqueElements[aN])
          for (size_t j=0; j<atomCH.size(); j++)
            MatAdd('N','N', atomNB, atomNB,
              MatsT(1.), atomCH[j], atomNB,
              MatsT(1.), CH[j] + NB * cumeNBs[i] + cumeNBs[i], NB,
              CH[j] + NB * cumeNBs[i] + cumeNBs[i], NB);
      } else {
        atoms_[k].computeOneEX2C(emPert, atomCoreH);
#ifndef REAL_SPACE_X2C_ALGORITHM
        atoms_[k].computeOneEX2C_Umatrix();
#endif
      }
    }

    assert(cumeNB == NB);

    if (type_.diagonalOnly) return;

    /// DLU/ALU O(N^3) Algorithm (commented out):
    ///   1> Combile diagonal matrix U = diag(U_A, U_B, U_C, ...)
    ///   2> Compute Hx2c = U * D * U
#ifdef UDU_ATOMIC_X2C_ALGORITHM
    computeOneEX2C_Umatrix();
    this->uncontractedInts_.computeAOOneP(
        this->molecule_, this->uncontractedBasis_, emPert,
        {{OVERLAP,0}, {KINETIC,0}, {NUCLEAR_POTENTIAL,0}},
        this->ssOptions_.hamiltonianOptions);
    this->computeOneEX2C_UDU(emPert, coreH);
    return;
#endif


    /// DLU/ALU O(N^2) Algorithm:
    ///   Hx2c_AB = U_A * D_AB * U_B

    // Allocate memory
    IntsT *T2c = CQMemManager::get().malloc<IntsT>(4*maxAtomNP*maxAtomNP);
    IntsT *V2c = CQMemManager::get().malloc<IntsT>(4*maxAtomNP*maxAtomNP);
    MatsT *W2c = CQMemManager::get().malloc<MatsT>(4*maxAtomNP*maxAtomNP);
    MatsT *SCR = CQMemManager::get().malloc<MatsT>(4*maxAtomNP*maxAtomNB);
    std::fill_n(SCR, 4*maxAtomNP*maxAtomNB, MatsT(0.));
    MatsT *Hx2c = CQMemManager::get().malloc<MatsT>(4*maxAtomNB*maxAtomNB);
    std::fill_n(Hx2c, 4*maxAtomNB*maxAtomNB, MatsT(0.));

    // Allocate memory for the uncontracted spin components
    // of the 2C CH
    MatsT *HUnS = CQMemManager::get().malloc<MatsT>(maxAtomNB*maxAtomNB);
    MatsT *HUnZ = CQMemManager::get().malloc<MatsT>(maxAtomNB*maxAtomNB);
    MatsT *HUnX = CQMemManager::get().malloc<MatsT>(maxAtomNB*maxAtomNB);
    MatsT *HUnY = CQMemManager::get().malloc<MatsT>(maxAtomNB*maxAtomNB);

    this->uncontractedInts_.computeAOOneP(
        this->molecule_, this->uncontractedBasis_, emPert,
        {{KINETIC,0}, {NUCLEAR_POTENTIAL,0}},
        this->ssOptions_.hamiltonianOptions);

    this->W = std::make_shared<cqmatrix::Matrix<MatsT>>(
        std::dynamic_pointer_cast<OnePRelInts<IntsT>>(
            this->uncontractedInts_.potential)->template formW<MatsT>());

    size_t cumeINP = 0;
    size_t cumeINB = 0;
    for (size_t i = 0; i < Natom; i++) {
      size_t I = atomIdx_[i];

      size_t atomINP = atoms_[I].uncontractedBasis_.nPrimitive;
      size_t atomINB = atoms_[I].basisSet_.nBasis;

      size_t cumeJNP = 0;
      size_t cumeJNB = 0;
      for (size_t j = 0; j <= i; j++) {
        size_t J = atomIdx_[j];

        size_t atomJNP = atoms_[J].uncontractedBasis_.nPrimitive;
        size_t atomJNB = atoms_[J].basisSet_.nBasis;

        size_t TVshift = cumeINP + NP*cumeJNP;
        // T2c_AB = [ T_AB   0   ]
        //          [  0    T_AB ]
        SetMatDiag(atomINP,atomJNP,
            this->uncontractedInts_.kinetic->pointer() + TVshift,
            NP,T2c,2*atomINP);

        // V2c_AB = [ V_AB   0   ]
        //          [  0    V_AB ]
        SetMatDiag(atomINP,atomJNP,
            this->uncontractedInts_.potential->pointer() + TVshift,
            NP,V2c,2*atomINP);

        // W2c_AB = [ W11_AB  W12_AB]
        //          [ W21_AB  W22_AB]
        size_t WsubShift = cumeINP + 2*NP*cumeJNP;
        std::fill_n(W2c,4*atomINP*atomJNP,MatsT(0.));
        SetMat('N',atomINP,atomJNP,MatsT(1.),
          this->W->pointer() + WsubShift,2*NP,W2c,2*atomINP);
        SetMat('N',atomINP,atomJNP,MatsT(1.),
          this->W->pointer() + 2*NP*NP + WsubShift,2*NP,
          W2c + 2*atomINP*atomJNP,2*atomINP);
        SetMat('N',atomINP,atomJNP,MatsT(1.),
          this->W->pointer() + NP + WsubShift,2*NP,
          W2c + atomINP,2*atomINP);
        SetMat('N',atomINP,atomJNP,MatsT(1.),
          this->W->pointer() + 2*NP*NP + NP + WsubShift,2*NP,
          W2c + 2*atomINP*atomJNP + atomINP,2*atomINP);

        // Hx2c_AB = U_AA D_AB U_BB
        // Hx2c = UL^H * T2c * US
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*atomINP,2*atomJNB,2*atomJNP,MatsT(1.),
          T2c,2*atomINP,atoms_[J].US,2*atomJNP,MatsT(0.),SCR,2*atomINP);
        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*atomINB,2*atomJNB,2*atomINP,MatsT(1.),
          atoms_[I].UL,2*atomINP,SCR,2*atomINP,MatsT(0.),Hx2c,2*atomINB);
        // Hx2c += US^H * T2c * UL
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*atomINP,2*atomJNB,2*atomJNP,MatsT(1.),
          T2c,2*atomINP,atoms_[J].UL,2*atomJNP,MatsT(0.),SCR,2*atomINP);
        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*atomINB,2*atomJNB,2*atomINP,MatsT(1.),
          atoms_[I].US,2*atomINP,SCR,2*atomINP,MatsT(1.),Hx2c,2*atomINB);
        // Hx2c -= US^H * T2c * US
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*atomINP,2*atomJNB,2*atomJNP,MatsT(1.),
          T2c,2*atomINP,atoms_[J].US,2*atomJNP,MatsT(0.),SCR,2*atomINP);
        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*atomINB,2*atomJNB,2*atomINP,MatsT(-1.),
          atoms_[I].US,2*atomINP,SCR,2*atomINP,MatsT(1.),Hx2c,2*atomINB);
        // Hx2c += UL^H * V2c * UL
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*atomINP,2*atomJNB,2*atomJNP,MatsT(1.),
          V2c,2*atomINP,atoms_[J].UL,2*atomJNP,MatsT(0.),SCR,2*atomINP);
        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*atomINB,2*atomJNB,2*atomINP,MatsT(1.),
          atoms_[I].UL,2*atomINP,SCR,2*atomINP,MatsT(1.),Hx2c,2*atomINB);
        // Hx2c += 1/(4*C**2) US^H * W * US
        blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,2*atomINP,2*atomJNB,2*atomJNP,
          MatsT(0.25/SpeedOfLight/SpeedOfLight),
          W2c,2*atomINP,atoms_[J].US,2*atomJNP,MatsT(0.),SCR,2*atomINP);
        blas::gemm(blas::Layout::ColMajor,blas::Op::ConjTrans,blas::Op::NoTrans,2*atomINB,2*atomJNB,2*atomINP,MatsT(1.),
          atoms_[I].US,2*atomINP,SCR,2*atomINP,MatsT(1.),Hx2c,2*atomINB);

        if (this->ssOptions_.hamiltonianOptions.OneESpinOrbit)
          SpinScatter(atomINB,atomJNB,Hx2c,2*atomINB,HUnS,atomINB,
            HUnZ,atomINB,HUnY,atomINB,HUnX,atomINB);
        else {
          MatAdd('N','N',atomINB,atomJNB,MatsT(1.),Hx2c,2*atomINB,MatsT(1.),
            Hx2c+atomINB+2*atomINB*atomJNB,2*atomINB,HUnS,atomINB);
          if (CH.size() > 1) {
            std::fill_n(HUnZ,atomINB*atomJNB,MatsT(0.));
            std::fill_n(HUnY,atomINB*atomJNB,MatsT(0.));
            std::fill_n(HUnX,atomINB*atomJNB,MatsT(0.));
          }
        }

        size_t CHshift = cumeINB + NB*cumeJNB;
        SetMat('N',atomINB,atomJNB,MatsT(1.),HUnS,atomINB,CH[0]+CHshift,NB);
        if (CH.size() > 1) {
          SetMat('N',atomINB,atomJNB,MatsT(1.),HUnZ,atomINB,CH[1]+CHshift,NB);
          SetMat('N',atomINB,atomJNB,MatsT(1.),HUnY,atomINB,CH[2]+CHshift,NB);
          SetMat('N',atomINB,atomJNB,MatsT(1.),HUnX,atomINB,CH[3]+CHshift,NB);
        }

        if (i != j) {
          CHshift = cumeJNB + NB*cumeINB;
          SetMat('C',atomINB,atomJNB,MatsT(1.),HUnS,atomINB,CH[0]+CHshift,NB);
          if (CH.size() > 1) {
            SetMat('C',atomINB,atomJNB,MatsT(1.),HUnZ,atomINB,CH[1]+CHshift,NB);
            SetMat('C',atomINB,atomJNB,MatsT(1.),HUnY,atomINB,CH[2]+CHshift,NB);
            SetMat('C',atomINB,atomJNB,MatsT(1.),HUnX,atomINB,CH[3]+CHshift,NB);
          }
        }

        cumeJNP += atomJNP;
        cumeJNB += atomJNB;
      }

      cumeINP += atomINP;
      cumeINB += atomINB;
    }

    CQMemManager::get().free(T2c, V2c, W2c, SCR, Hx2c, HUnS, HUnZ, HUnY, HUnX);

  }

  template<> void AtomicX2C<dcomplex,dcomplex>::computeOneEX2C(EMPerturbation&,
      std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>>) {
    CErr("X2C + Complex Ints NYI",std::cout);
  }

  /**
   *  \brief Compute the picture change matrices UL, US
   */
  template <typename MatsT, typename IntsT>
  void AtomicX2C<MatsT, IntsT>::computeOneEX2C_Umatrix() {

    size_t NP = this->uncontractedBasis_.nPrimitive;
    size_t NB = this->basisSet_.nBasis;

    this->UL = CQMemManager::get().malloc<MatsT>(4*NP*NB);
    std::fill_n(this->UL,4*NP*NB,MatsT(0.));
    this->US = CQMemManager::get().malloc<MatsT>(4*NP*NB);
    std::fill_n(this->US,4*NP*NB,MatsT(0.));

    size_t cumeNP = 0, cumeNB = 0;

    for (size_t k : atomIdx_) {
      const auto &at = atoms_[k];

      size_t atomNP = at.uncontractedBasis_.nPrimitive;
      size_t atomNB = at.basisSet_.nBasis;

      SetMat('N', atomNP, atomNB, MatsT(1.),
        at.UL, 2*atomNP,
        this->UL + 2*NP*cumeNB + cumeNP, 2*NP);
      SetMat('N', atomNP, atomNB, MatsT(1.),
        at.UL + atomNP, 2*atomNP,
        this->UL + 2*NP*cumeNB + NP + cumeNP, 2*NP);
      SetMat('N', atomNP, atomNB, MatsT(1.),
        at.UL + 2*atomNP*atomNB, 2*atomNP,
        this->UL + 2*NP*NB + 2*NP*cumeNB + cumeNP, 2*NP);
      SetMat('N', atomNP, atomNB, MatsT(1.),
        at.UL + 2*atomNP*atomNB + atomNP, 2*atomNP,
        this->UL + 2*NP*NB + 2*NP*cumeNB + NP + cumeNP, 2*NP);

      SetMat('N', atomNP, atomNB, MatsT(1.),
        at.US, 2*atomNP,
        this->US + 2*NP*cumeNB + cumeNP, 2*NP);
      SetMat('N', atomNP, atomNB, MatsT(1.),
        at.US + atomNP, 2*atomNP,
        this->US + 2*NP*cumeNB + NP + cumeNP, 2*NP);
      SetMat('N', atomNP, atomNB, MatsT(1.),
        at.US + 2*atomNP*atomNB, 2*atomNP,
        this->US + 2*NP*NB + 2*NP*cumeNB + cumeNP, 2*NP);
      SetMat('N', atomNP, atomNB, MatsT(1.),
        at.US + 2*atomNP*atomNB + atomNP, 2*atomNP,
        this->US + 2*NP*NB + 2*NP*cumeNB + NP + cumeNP, 2*NP);

      cumeNP += atomNP;
      cumeNB += atomNB;

    }

  }

  template<> void AtomicX2C<dcomplex,dcomplex>::computeOneEX2C_Umatrix() {
    CErr("X2C + Complex Ints NYI",std::cout);
  }

  template class AtomicX2C<double,double>;
  template class AtomicX2C<dcomplex,double>;
  template class AtomicX2C<dcomplex,dcomplex>;

  // Instantiate copy constructors
  template AtomicX2C<dcomplex,double>::AtomicX2C(const AtomicX2C<double,double> &, int);
  template AtomicX2C<dcomplex,dcomplex>::AtomicX2C(const AtomicX2C<dcomplex,dcomplex> &, int);

}; // namespace ChronusQ
