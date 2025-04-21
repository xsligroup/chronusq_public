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

#include <coupledcluster/CCSD.hpp>
#include <coupledcluster/DFCCSD.hpp>
#include <coupledcluster/CCSDT.hpp>
#include <coupledcluster/EOMCC.hpp>
#include <coupledcluster/MBExpansionImpl.hpp>
#include <coupledcluster/TAERI.hpp>
#include <singleslater.hpp>
#include <cxxapi/options.hpp>
#include <cxxapi/boilerplate.hpp>
#include <intermediates.hpp>

namespace ChronusQ {

  template class CCSD<dcomplex>;
  template class DFCCSD<dcomplex>;
  template class CCSDT<dcomplex>;
  template class EOMCCBase<dcomplex>;
  template class MBExpansion<dcomplex>;
  template class MBExpansionSet<dcomplex>;
  template class MBExpansionSetDebug<dcomplex>;

  std::string orbitalSelectionToString(std::vector<size_t> orbitals);

  std::shared_ptr<EOMCCBase<dcomplex>> build_CVSEOMCCSD(
                                              const SafeFile &savFile,
                                              CCIntermediates<dcomplex> &intermediates, const EOMSettings &eomSettings,
                                              const CoupledClusterSettings &ccSettings) ;
  std::shared_ptr<EOMCCBase<dcomplex>> build_EOMCCSD(
                                              const SafeFile &savFile,
                                              CCIntermediates<dcomplex> &intermediates, const EOMSettings &eomSettings,
                                              const CoupledClusterSettings &ccSettings) ;
  std::shared_ptr<EOMCCBase<dcomplex>> build_EOMDIP_3h1p(
                                              const SafeFile &savFile,
                                              CCIntermediates<dcomplex> &intermediates, const EOMSettings &eomSettings,
                                              const CoupledClusterSettings &ccSettings) ;
  std::shared_ptr<EOMCCBase<dcomplex>> build_EOMDIP_4h2pCCSDT(
                                              const SafeFile &savFile,
                                              CCIntermediates<dcomplex> &intermediates, const EOMSettings &eomSettings,
                                              const CoupledClusterSettings &ccSettings) ;
  std::shared_ptr<EOMCCBase<dcomplex>> build_EOMEA(
                                              const SafeFile &savFile,
                                              CCIntermediates<dcomplex> &intermediates, const EOMSettings &eomSettings,
                                              const CoupledClusterSettings &ccSettings) ;
  std::shared_ptr<EOMCCBase<dcomplex>> build_EOMIP_2h1p(
                                              const SafeFile &savFile,
                                              CCIntermediates<dcomplex> &intermediates, const EOMSettings &eomSettings,
                                              const CoupledClusterSettings &ccSettings) ;
  std::shared_ptr<EOMCCBase<dcomplex>> build_EOMIP_3h2p(
                                              const SafeFile &savFile,
                                              CCIntermediates<dcomplex> &intermediates, const EOMSettings &eomSettings,
                                              const CoupledClusterSettings &ccSettings) ;

  void EOMSettings::printEOMCCSettings(std::ostream &out) {
    out << "Equation of Motion Coupled Cluster (EOMCC) Settings:" << std::endl << std::endl;

    out << std::setw(45) << std::left << "  Type:"
        << "Singles and Doubles" << std::endl;

    out << std::setw(45) << std::left << "  EOM type:";
    switch (eom_type) {
      case EOM_TYPE::EE:
        out << "EE";
        break;
      case EOM_TYPE::IP:
        out << "IP";
        break;
      case EOM_TYPE::EA:
        out << "EA";
        break;
      case EOM_TYPE::DIP:
        switch (ip_level) {
          case 4:
            out << "DIP (up to 4 hole 2 particle)";
            break;
          case 3:
            out << "DIP (up to 3 hole 1 particle)";
            break;
        }
        break;
    }
    out << std::endl;

    out << std::setw(45) << std::left << "  Hbar matrix type:";
    switch (hbar_type) {
      case EOM_HBAR_TYPE::EXPLICIT:
        out << "Explicit full matrix";
        break;
      case EOM_HBAR_TYPE::IMPLICIT:
        out << "Implicit many-body tensor contraction";
        break;
      case EOM_HBAR_TYPE::DEBUG:
        out << "Debug: Both explicit full matrix and implicit many-body tensor contraction";
        break;
    }
    out << std::endl;

    out << std::setw(45) << std::left << "  Diagonalization method:";
    switch (diag_method) {
      case EOM_DIAG_METHOD::FULL:
        out << "Full diagonalization";
        break;
      case EOM_DIAG_METHOD::DAVIDSON:
        out << "Davidson-Liu";
        break;
      case EOM_DIAG_METHOD::GPLHR:
        out << "GPLHR";
        break;
    }
    out << std::endl;

    out << std::setw(45) << std::left << "  Compute Oscillator Strength:"
        << (oscillator_strength ? "ON" : "OFF") << std::endl;

    out << std::setw(45) << std::left << "  Save Hbar Matrix:"
        << (save_hamiltonian ? "ON" : "OFF") << std::endl;

    out << std::setw(45) << std::left << "  Contain Active Space:"
        << (containActive() ? "ON" : "OFF") << std::endl;

    if (containActive()) {
      if (not cvs_core.empty())
        out << "    * External Occupied / CVS Core orbitals         : " << orbitalSelectionToString(cvs_core) << std::endl;
      if (not external_virtual.empty())
        out << "    * External Virtual / CVS Continuum orbitals         : " << orbitalSelectionToString(external_virtual) << std::endl;
    }

    std::pair<double, char> mem_postfix = memSize(estimate_mem_peak());
    std::cout << std::setw(45) << std::left << "  Estimated TiledArray memory requirement: " << std::fixed << std::setprecision(1)
    << mem_postfix.first << mem_postfix.second << "B" << std::endl;

  }

  size_t EOMSettings::estimate_mem_peak() const {

    size_t count = 0;

    count += intermediate_mem();

    size_t nVec = n_MBExpansion();
    size_t size_per_MBExpansion = MBExpansionSize();
    count += nVec * size_per_MBExpansion;
    
    return count * sizeof(dcomplex);
  }
  size_t EOMSettings::n_MBExpansion() const {
    size_t nVec = 0;
    if (hbar_type == EOM_HBAR_TYPE::EXPLICIT) {
      nVec += nroots * (oscillator_strength ? 2 : 0);
    } else {
      nVec += nroots * ((davidson_whenSc > 1 ? davidson_guess_multiplier : 1) *2
          + davidson_subspace_multiplier * 2 + davidson_guess_multiplier
          + (oscillator_strength ? 2 : 1));
    }
    return nVec;
  }
  EOM_IMPLEMENTATION EOMSettings::find_eom_implementation() {
    if (not containActive()) {
      if (eom_type == EOM_TYPE::EE){
          eom_implementation = EOM_IMPLEMENTATION::EOMCCSD;
          return EOM_IMPLEMENTATION::EOMCCSD;
      } else if (eom_type == EOM_TYPE::EA){
          eom_implementation = EOM_IMPLEMENTATION::EOMEA;
          return EOM_IMPLEMENTATION::EOMEA;
      } else if (eom_type == EOM_TYPE::IP){
          if (ip_level == 2){
            eom_implementation = EOM_IMPLEMENTATION::EOMIP_2h1p;
            return EOM_IMPLEMENTATION::EOMIP_2h1p;
          } else if (ip_level == 3) {
            eom_implementation = EOM_IMPLEMENTATION::EOMIP_3h2p;
            return EOM_IMPLEMENTATION::EOMIP_3h2p;
          }
      } else if (eom_type == EOM_TYPE::DIP){
          if (ip_level == 3) {
            eom_implementation = EOM_IMPLEMENTATION::EOMDIP_3h1p;
            return EOM_IMPLEMENTATION::EOMDIP_3h1p;
          } else if (ip_level == 4) {
            eom_implementation = EOM_IMPLEMENTATION::EOMDIP_4h2p;
            return EOM_IMPLEMENTATION::EOMDIP_4h2p;
          }
      }
    } else {
      if (eom_type == EOM_TYPE::EE){
        eom_implementation = EOM_IMPLEMENTATION::CVSEOMCCSD;
        return EOM_IMPLEMENTATION::CVSEOMCCSD;
      }
    }
    return EOM_IMPLEMENTATION::SOMETHING_WRONG;
  }

  size_t EOMSettings::MBExpansionSize() const {
    TAManager &TAmanager = TAManager::get();
    size_t count = 0;
    switch (static_cast<int>(eom_implementation)) {
      case static_cast<int>(EOM_IMPLEMENTATION::CVSEOMCCSD):
        count += TAmanager.elem_per_TA("rc");
        count += TAmanager.elem_per_TA("rrcc");
        count += TAmanager.elem_per_TA("rrch");
        break;
      case static_cast<int>(EOM_IMPLEMENTATION::EOMCCSD):
        count += TAmanager.elem_per_TA("ov");
        count += TAmanager.elem_per_TA("oovv");
        break;
      case static_cast<int>(EOM_IMPLEMENTATION::EOMDIP_3h1p):
        count += TAmanager.elem_per_TA("oo");
        count += TAmanager.elem_per_TA("ooov");
        break;
      case static_cast<int>(EOM_IMPLEMENTATION::EOMDIP_4h2p):
        count += TAmanager.elem_per_TA("oo");
        count += TAmanager.elem_per_TA("ooov");
        count += TAmanager.elem_per_TA("oooovv");
        break;
      case static_cast<int>(EOM_IMPLEMENTATION::EOMEA):
        count += TAmanager.elem_per_TA("v");
        count += TAmanager.elem_per_TA("ovv");
        break;
      case static_cast<int>(EOM_IMPLEMENTATION::EOMIP_2h1p):
        count += TAmanager.elem_per_TA("o");
        count += TAmanager.elem_per_TA("oov");
        break;
      case static_cast<int>(EOM_IMPLEMENTATION::EOMIP_3h2p):
        count += TAmanager.elem_per_TA("o");
        count += TAmanager.elem_per_TA("oov");
        count += TAmanager.elem_per_TA("ooovv");
        break;
      default:
        CErr("Unknown EOM implementation");
    }
    return count;
  }

  size_t EOMSettings::intermediate_mem() const {
    TAManager &TAmanager = TAManager::get();
    size_t count = 0;
    switch (static_cast<int>(eom_implementation)) {
      case static_cast<int>(EOM_IMPLEMENTATION::CVSEOMCCSD):
        // W
        count += TAmanager.elem_per_TA("ooov");
        count += TAmanager.elem_per_TA("vovv");
        count += TAmanager.elem_per_TA("ovoo");
        count += TAmanager.elem_per_TA("vvvo");
        //tmp
        count += TAmanager.elem_per_TA("ovvo");
        count += TAmanager.elem_per_TA("rrcc");
        count += TAmanager.elem_per_TA("co");
        count += TAmanager.elem_per_TA("rr");
        count += TAmanager.elem_per_TA("rrch");
        count += TAmanager.elem_per_TA("co");
        count += TAmanager.elem_per_TA("ho");
        count += TAmanager.elem_per_TA("rr");
        count += TAmanager.elem_per_TA("rrcc");
        count += TAmanager.elem_per_TA("rrch");
        //G
        count += TAmanager.elem_per_TA("rr");
        count += TAmanager.elem_per_TA("oo");
        // MOints
        count += TAmanager.elem_per_TA("vooo");
        count += TAmanager.elem_per_TA("vvoo");
        count += TAmanager.elem_per_TA("vovo");
        count += TAmanager.elem_per_TA("vvvo");
        if (oscillator_strength) {
          // tmp
          count += 1 * TAmanager.elem_per_TA("oo");
          count += 2 * TAmanager.elem_per_TA("ov");
          count += 1 * TAmanager.elem_per_TA("vv");
          // Density
          count += 1 * TAmanager.elem_per_TA("oo");
          count += 2 * TAmanager.elem_per_TA("ov");
          count += 1 * TAmanager.elem_per_TA("vv");
          // Right hand ground state
          count += MBExpansionSize();
        }
        break;
      case static_cast<int>(EOM_IMPLEMENTATION::EOMCCSD):
        //W
        count += TAmanager.elem_per_TA("ooov");
        count += TAmanager.elem_per_TA("vovv");
        count += TAmanager.elem_per_TA("ovoo");
        count += TAmanager.elem_per_TA("vvvo");
        //tmp
        count += TAmanager.elem_per_TA("ovvo");
        count += TAmanager.elem_per_TA("ovvo");
        count += TAmanager.elem_per_TA("vv");
        count += TAmanager.elem_per_TA("oo");
        //G
        count += TAmanager.elem_per_TA("vv");
        count += TAmanager.elem_per_TA("oo");
        // MOints
        count += TAmanager.elem_per_TA("vooo");
        count += TAmanager.elem_per_TA("vvoo");
        count += TAmanager.elem_per_TA("vovo");
        count += TAmanager.elem_per_TA("vvvo");
        if (oscillator_strength) {
          // tmp
          count += 1 * TAmanager.elem_per_TA("oo");
          count += 2 * TAmanager.elem_per_TA("ov");
          count += 1 * TAmanager.elem_per_TA("vv");
          // Density
          count += 1 * TAmanager.elem_per_TA("oo");
          count += 2 * TAmanager.elem_per_TA("ov");
          count += 1 * TAmanager.elem_per_TA("vv");
          // Right hand ground state
          count += MBExpansionSize();
        }
        break;
      case static_cast<int>(EOM_IMPLEMENTATION::EOMIP_3h2p):
        // W
        count += TAmanager.elem_per_TA("ooov");
        count += TAmanager.elem_per_TA("vovv");
        count += TAmanager.elem_per_TA("ovoo");
        count += TAmanager.elem_per_TA("vvvo");
        // tmp
        count += TAmanager.elem_per_TA("ovvo");
        count += TAmanager.elem_per_TA("ovvo");
        count += TAmanager.elem_per_TA("voo");
        count += TAmanager.elem_per_TA("v");
        count += TAmanager.elem_per_TA("vvooo");
        count += TAmanager.elem_per_TA("ooo");
        count += TAmanager.elem_per_TA("ooo");
        count += TAmanager.elem_per_TA("vov");
        count += TAmanager.elem_per_TA("v");
        count += TAmanager.elem_per_TA("vvv");
        count += TAmanager.elem_per_TA("oov");
        // MOints
        count += TAmanager.elem_per_TA("vooo");
        count += TAmanager.elem_per_TA("vvoo");
        count += TAmanager.elem_per_TA("vovo");
        count += TAmanager.elem_per_TA("vvvo");
        break;
      case static_cast<int>(EOM_IMPLEMENTATION::EOMIP_2h1p):
        // W
        count += TAmanager.elem_per_TA("ooov");
        count += TAmanager.elem_per_TA("vovv");
        count += TAmanager.elem_per_TA("ovoo");
        // tmp
        count += TAmanager.elem_per_TA("ovvo");
        count += TAmanager.elem_per_TA("voo");
        count += TAmanager.elem_per_TA("v");
        // MOints
        count += TAmanager.elem_per_TA("vooo");
        count += TAmanager.elem_per_TA("vvoo");
        count += TAmanager.elem_per_TA("vovo");
        count += TAmanager.elem_per_TA("vvvo");
        break;
      case static_cast<int>(EOM_IMPLEMENTATION::EOMDIP_3h1p):
        count += 14 * TAmanager.elem_per_TA("vo");
        count +=  6 * TAmanager.elem_per_TA("oo");
        count +=  4 * TAmanager.elem_per_TA("vooo");
        count +=  3 * TAmanager.elem_per_TA("oooo");
        count +=  2 * TAmanager.elem_per_TA("vv");
        count +=  1 * TAmanager.elem_per_TA("vvoo");
        // recycled temporary memory space from tmp TA objects
	    count += TAmanager.elem_per_TA("oo");
        count += TAmanager.elem_per_TA("vooo");
        count += TAmanager.elem_per_TA("vvvo");
        // MOints
        count += TAmanager.elem_per_TA("oooo");
        count += TAmanager.elem_per_TA("vooo");
        count += TAmanager.elem_per_TA("vovo");
        count += TAmanager.elem_per_TA("vvoo");
        count += TAmanager.elem_per_TA("vvvo");
		break;
      case static_cast<int>(EOM_IMPLEMENTATION::EOMDIP_4h2p):
        count += 17 * TAmanager.elem_per_TA("vvoo");
        count += 15 * TAmanager.elem_per_TA("vo");
        count += 11 * TAmanager.elem_per_TA("vooo");
        count +=  8 * TAmanager.elem_per_TA("voov");
        count +=  8 * TAmanager.elem_per_TA("ovvo");
        count +=  6 * TAmanager.elem_per_TA("vovo");
        count +=  6 * TAmanager.elem_per_TA("oo");
        count +=  5 * TAmanager.elem_per_TA("oovv");
        count +=  4 * TAmanager.elem_per_TA("vv");
        count +=  4 * TAmanager.elem_per_TA("voovoo");
        count +=  4 * TAmanager.elem_per_TA("ovov");
        count +=  4 * TAmanager.elem_per_TA("ov");
        count +=  4 * TAmanager.elem_per_TA("oooo");
        count +=  3 * TAmanager.elem_per_TA("vovv");
        count +=  3 * TAmanager.elem_per_TA("oovvoo");
        count +=  3 * TAmanager.elem_per_TA("oovo");
        count +=  2 * TAmanager.elem_per_TA("vvvo");
        count +=  2 * TAmanager.elem_per_TA("vvov");
        count +=  2 * TAmanager.elem_per_TA("vvoooo");
        count +=  2 * TAmanager.elem_per_TA("oovoov");
        count +=  2 * TAmanager.elem_per_TA("ooov");
        count +=  1 * TAmanager.elem_per_TA("vovooo");
        count +=  1 * TAmanager.elem_per_TA("vooooo");
        count +=  1 * TAmanager.elem_per_TA("ovvooo");
        count +=  1 * TAmanager.elem_per_TA("ovoo");
        // recycled temporary memory space from tmp TA objects
	    count += TAmanager.elem_per_TA("oooooo");
	    count += TAmanager.elem_per_TA("oooo");
	    count += TAmanager.elem_per_TA("oooovv");
	    count += TAmanager.elem_per_TA("ooovov");
	    count += TAmanager.elem_per_TA("ooovvo");
	    count += TAmanager.elem_per_TA("oo");
	    count += TAmanager.elem_per_TA("oovo");
	    count += TAmanager.elem_per_TA("oovovo");
	    count += TAmanager.elem_per_TA("ovovoo");
	    count += TAmanager.elem_per_TA("ovvooo");
	    count += TAmanager.elem_per_TA("vooo");
	    count += TAmanager.elem_per_TA("vooovo");
	    count += TAmanager.elem_per_TA("voovoo");
	    count += TAmanager.elem_per_TA("vovooo");
	    count += TAmanager.elem_per_TA("vvoooo");
        break;
      case static_cast<int>(EOM_IMPLEMENTATION::EOMEA):
        count += 14 * TAmanager.elem_per_TA("vo");
        count +=  5 * TAmanager.elem_per_TA("vvoo");
        count +=  5 * TAmanager.elem_per_TA("oo");
        count +=  3 * TAmanager.elem_per_TA("vvvo");
        count +=  3 * TAmanager.elem_per_TA("vv");
        count +=  2 * TAmanager.elem_per_TA("vooo");
        count +=  1 * TAmanager.elem_per_TA("vvvv");
        // recycled temporary memory space from tmp TA objects
    	count += TAmanager.elem_per_TA("ooo");
    	count += TAmanager.elem_per_TA("o");
    	count += TAmanager.elem_per_TA("voo");
    	count += TAmanager.elem_per_TA("vvo");
    	count += TAmanager.elem_per_TA("vvvo");
        // MOints
        count += TAmanager.elem_per_TA("vooo");
        count += TAmanager.elem_per_TA("vvoo");
        count += TAmanager.elem_per_TA("vovo");
        count += TAmanager.elem_per_TA("vvvo");
        count += TAmanager.elem_per_TA("vvvv");
		break;
      default:
        CErr("Unknown EOM implementation");
    }
    return count;
  }

  template <typename MatsT>
  void swapVectorToFirst(size_t groundIndex, MatsT* M, size_t ldm) {
    MatsT* tmpVec = CQMemManager::get().malloc<MatsT>(ldm);

    SetMat('N', ldm, 1, 1.0, M + groundIndex * ldm, ldm, tmpVec, ldm);
    while (groundIndex > 0) {
      SetMat('N', ldm, 1, 1.0, M + (groundIndex - 1) * ldm, ldm, M + groundIndex * ldm, ldm);
      groundIndex--;
    }
    SetMat('N', ldm, 1, 1.0, tmpVec, ldm, M, ldm);

    CQMemManager::get().free(tmpVec);
  }

  template <typename MatsT>
  void findTrueGroundStateEOMCCEigen(size_t Hbar_dim_w0,
                                     MatsT* theta_w0, MatsT* VL_w0, MatsT* VR_w0, double e_conv) {

    size_t groundIndex = 0;
    double absGroundL0 = std::abs(VL_w0[0]);
    double secondLargestL0 = std::abs(VL_w0[Hbar_dim_w0]);
    if (secondLargestL0 > absGroundL0) {
      groundIndex = 1;
      std::swap(secondLargestL0, absGroundL0);
    }
    for (size_t i = 2; i < Hbar_dim_w0; i++) {
      if (std::abs(VL_w0[i * Hbar_dim_w0]) > absGroundL0) {
        groundIndex = i;
        absGroundL0 = std::abs(VL_w0[i * Hbar_dim_w0]);
      }
    }

    std::cout << "Found true ground state at index " << groundIndex
              << " with abs(L0) = " << absGroundL0
              <<", the next largest abs(L0) = " << secondLargestL0 << std::endl;

    if (groundIndex != 0) {

      std::cout << "Swap ground state to index 0 ..." << std::endl;

      swapVectorToFirst(groundIndex, theta_w0, 1);
      swapVectorToFirst(groundIndex, VL_w0, Hbar_dim_w0);
      swapVectorToFirst(groundIndex, VR_w0, Hbar_dim_w0);

      std::cout << "Swap finished" << std::endl;

    }
  }

  template <typename _F>
  void biOrthoNormalize(size_t N, size_t nR, RawVectors<_F> &VL, RawVectors<_F> &VR) {

    // Biorthonomalize VR_ and VL_
    // VL_^\dagger*VR_ = P*L*U
    // VL_^\dagger = L^{-1}*P^{T}*VL_\dagger
    // VR_ = VR_ * U^{-1}
    std::vector<int64_t> IPIV(nR);
    cqmatrix::Matrix<_F> LUMat(nR);

//    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans,
//               nR,nR,N,_F(1.),VL,N,VR,N,_F(0.),LUMat,nR);
    VL.dot_product(0, VR, 0, nR, nR, LUMat.pointer(), nR, false);
//    prettyPrintSmart(std::cout,"LUMat:  ",LUMat,nR,nR,nR);

    if (MPIRank() == 0) {
      lapack::getrf(nR, nR, LUMat.pointer(), nR, IPIV.data());
//    prettyPrintSmart(std::cout,"LUMat after LU:  ",LUMat,nR,nR,nR);


      // Compute the inverse of lower and upper triangular matrices in-place
      lapack::trtri(lapack::Uplo::Upper, lapack::Diag::NonUnit, nR, LUMat.pointer(), nR);
      lapack::trtri(lapack::Uplo::Lower, lapack::Diag::Unit, nR, LUMat.pointer(), nR);

      // VR_ = VR_ * U^{-1}
      blas::trmm(blas::Layout::ColMajor, blas::Side::Right, blas::Uplo::Upper,
                 blas::Op::NoTrans, lapack::Diag::NonUnit, N, nR, _F(1.0), LUMat.pointer(), nR, VR.getPtr(), N);

      // Apply P^{T}
      //IPIV represents elementary permutation matrices, whose transpose are themselves.
      //P = P1 * P2 * ... * Pn
      // P^{T} = Pn^{T} * ... * P2^{T} * P1^{T}
      // = Pn * ... * P2 * P1
      Eigen::Map<
          Eigen::Matrix<_F,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>
      > VLMap(VL.getPtr(),N,nR);

      for (int i = 0; i < nR; i++){
        if (i != IPIV[i] - 1){
//          std::cout << i << "   "  << IPIV[i] - 1 << std::endl;
          VLMap.col(IPIV[i] - 1).swap(VLMap.col(i));
        }
      }


      // VL_^\dagger = L^{-1}*P^{T}*VL_\dagger
      blas::trmm(blas::Layout::ColMajor, blas::Side::Right, blas::Uplo::Lower,
                 blas::Op::Trans, lapack::Diag::Unit, N, nR, _F(1.0), LUMat.pointer(), nR, VL.getPtr(), N);

//    prettyPrintSmart(std::cout,"New VL:  ",VL,N,nR,N);
//    prettyPrintSmart(std::cout,"New VR:  ",VR,N,nR,N);
    } // END ROOT_ONLY Section

    for (size_t i = 0; i < nR; i++) {
//      double norm = blas::nrm2(N, VR + i * N, 1);
//      blas::scal(N, 1.0/norm, VR + i * N, 1);
//      blas::scal(N, norm, VL + i * N, 1);
      double norm = VR.norm2F(i, 1);
      VR.scale(1.0/norm, i, 1);
      VL.scale(norm, i, 1);
    }

//    prettyPrintSmart(std::cout,"New VL after scale:  ",VL,N,nR,N);
//    prettyPrintSmart(std::cout,"New VR after scale:  ",VR,N,nR,N);

  }

  template <typename _F>
  void biOrthoNormalize(size_t nR, MBExpansionSet<_F> &VL, MBExpansionSet<_F> &VR) {

    // Biorthonomalize VR_ and VL_
    // VL_^\dagger*VR_ = P*L*U
    // VL_^\dagger = L^{-1}*P^{T}*VL_\dagger
    // VR_ = VR_ * U^{-1}
    std::vector<int64_t> IPIV(nR);
    cqmatrix::Matrix<_F> LUMat(nR);

    //    blas::gemm(blas::Layout::ColMajor, blas::Op::ConjTrans, blas::Op::NoTrans,
    //               nR,nR,N,_F(1.),VL,N,VR,N,_F(0.),LUMat,nR);
    VL.dot_product(0, VR, 0, nR, nR, LUMat.pointer(), nR, false);
    //    prettyPrintSmart(std::cout,"LUMat:  ",LUMat,nR,nR,nR);

    if (MPIRank() == 0) {
      lapack::getrf(nR, nR, LUMat.pointer(), nR, IPIV.data());
      //    prettyPrintSmart(std::cout,"LUMat after LU:  ",LUMat,nR,nR,nR);


      // Compute the inverse of lower and upper triangular matrices in-place
      lapack::trtri(lapack::Uplo::Upper, lapack::Diag::NonUnit, nR, LUMat.pointer(), nR);
      lapack::trtri(lapack::Uplo::Lower, lapack::Diag::Unit, nR, LUMat.pointer(), nR);
    }

    TA::get_default_world().gop.fence();
    LUMat.broadcast();
    MPIBCast(IPIV.data(), nR, 0, MPI_COMM_WORLD);

    cqmatrix::Matrix<_F> UMat(LUMat);
    UMat.setTriangle(blas::Uplo::Lower, 0.0, false);
    cqmatrix::Matrix<_F> LMat(LUMat);
    LMat.setTriangle(blas::Uplo::Upper, 0.0, true, 1.0);

    // VR_ = VR_ * U^{-1}
//    blas::trmm(blas::Layout::ColMajor, blas::Side::Right, blas::Uplo::Upper,
//               blas::Op::NoTrans, lapack::Diag::NonUnit, N, nR, _F(1.0), LUMat, nR, VR.getPtr(), N);
    MBExpansionSet<_F> Vcopy(VR);
    Vcopy.multiply_matrix(0, blas::Op::NoTrans, nR, nR, _F(1.0), UMat.pointer(), nR, _F(0.0), VR, 0);

    // Apply P^{T}
    //IPIV represents elementary permutation matrices, whose transpose are themselves.
    //P = P1 * P2 * ... * Pn
    // P^{T} = Pn^{T} * ... * P2^{T} * P1^{T}
    // = Pn * ... * P2 * P1
    Eigen::Map<
        Eigen::Matrix<_F,Eigen::Dynamic,Eigen::Dynamic,Eigen::ColMajor>
    > LMap(LMat.pointer(),nR,nR);

    for (int i = nR; i > 0; i--){
      if (i != IPIV[i - 1]){
//        std::cout << i - 1 << "   "  << IPIV[i - 1] - 1 << std::endl;
        LMap.col(IPIV[i - 1] - 1).swap(LMap.col(i - 1));
      }
    }

    LMat = LMat.T();

    // VL_^\dagger = L^{-1}*P^{T}*VL_\dagger
    Vcopy.set_data(0, nR, VL, 0);
    Vcopy.multiply_matrix(0, blas::Op::NoTrans, nR, nR, _F(1.0), LMat.pointer(), nR, _F(0.0), VL, 0);

    //    prettyPrintSmart(std::cout,"New VL:  ",VL,N,nR,N);
    //    prettyPrintSmart(std::cout,"New VR:  ",VR,N,nR,N);

    for (size_t i = 0; i < nR; i++) {
      //      double norm = blas::nrm2(N, VR + i * N, 1);
      //      blas::scal(N, 1.0/norm, VR + i * N, 1);
      //      blas::scal(N, norm, VL + i * N, 1);
      double norm = VR.norm2F(i, 1);
      VR.scale(1.0/norm, i, 1);
      VL.scale(norm, i, 1);
    }

    //    prettyPrintSmart(std::cout,"New VL after scale:  ",VL,N,nR,N);
    //    prettyPrintSmart(std::cout,"New VR after scale:  ",VR,N,nR,N);

  }


  std::vector<size_t> getGuessIndices(size_t nGuess, size_t length, const EOMSettings& eomSettings,
                                      const dcomplex *eomDiag, MPI_Comm comm) { 

    std::vector<size_t> guessIndices;
    guessIndices.reserve(nGuess);
    if (MPIRank(comm) == 0) {
      std::vector<size_t> diagSort(length, 0);
      std::iota(diagSort.begin(), diagSort.end(), 0);

      std::stable_sort(diagSort.begin(), diagSort.end(),
          [&] (size_t i , size_t j) { return std::real(eomDiag[i]) < std::real(eomDiag[j]); });

      if (eomSettings.davidson_Eref.empty()) {
        std::copy_n(diagSort.begin(), nGuess, std::back_inserter(guessIndices));

      } else {
        std::vector<size_t>::iterator curIterBegin = diagSort.begin()
            + eomSettings.davidson_guess_multiplier * eomSettings.davidson_nLowRoots;
        std::copy(diagSort.begin(), curIterBegin, std::back_inserter(guessIndices));
        double Eoffset = eomSettings.davidson_ErefAbs? 0. : std::real(eomDiag[diagSort[0]]);

        for (auto & pair: eomSettings.davidson_Eref) {
          double curERef = pair.first + Eoffset;
          size_t curNGuess = eomSettings.davidson_guess_multiplier * pair.second;
          curIterBegin = std::lower_bound(curIterBegin, diagSort.end(), curERef,
                                          [&eomDiag](size_t i, double x){ return std::real(eomDiag[i]) < x; });
          if (curIterBegin <= diagSort.end() - curNGuess) {
            std::copy_n(curIterBegin, curNGuess, std::back_inserter(guessIndices));
            curIterBegin += curNGuess;
          } else {
            CErr("No enough element above the reference energy to select.");
          }
        }
      }
    }
    MPIBCast(guessIndices.data(), nGuess, 0, comm);

    return guessIndices;

  } // getGuessIndices

  std::vector<size_t> LinearRange(size_t size, size_t start_index, size_t blksize){
    size_t blocks = size % blksize == 0 ? size / blksize : size / blksize + 1;
    std::vector<size_t> blk;
    blk.reserve(blocks);

    for (auto i = 0 ; i < blocks; i++)
      blk.push_back(blksize * i + start_index);
    return blk;
  }


  template <typename MatsT>
  template <typename IntsT>
  void CCIntermediates<MatsT>::initializeIntegrals(const cqmatrix::PauliSpinorMatrices<MatsT> &aoCoreH,
                                                   const cqmatrix::PauliSpinorMatrices<MatsT> &aoFock,
                                                   const cqmatrix::PauliSpinorMatrices<MatsT> &aoTwoeH,
                                                   const TwoPInts<IntsT> &aoTPI,
                                                   const MultipoleInts<IntsT> &lenElectric,
                                                   CoupledClusterSettings& ccSettings,
                                                   EOMSettings& eomSettings,
                                                   MatsT *mo, size_t nO, size_t nV,
                                                   size_t blksize, double nucRepEnergy,
                                                   CC_TYPE cctype, double denomshift_, bool pertT3, bool rebuildFock) {

    auto initIntStart = tick();

    size_t nMO = nO + nV, nAO = nMO/2;

    bool isRI = false;

    TAManager &TAmanager = TAManager::get();

    // Initialize ranges
    nOcc = nO;
    nVir = nV;

    TAERI<IntsT> taERI(aoTPI, blksize);
    TAmanager.addRangeType(aoLabel, taERI.getAOrange());

    for (auto it = ccSettings.frozen_occupied.rbegin(); it != ccSettings.frozen_occupied.rend(); it++) {
      if (*it >= nO + nV)
        CErr("EOMCC: Orbital index in input EOMCC.FROZENOCCUPIED out of range");
      if (*it >= nO)
        CErr("EOMCC: Virtual orbital appears in input EOMCC.FROZENOCCUPIED");
    }

    for (size_t i : eomSettings.cvs_core) {
      if (i >= nO + nV)
        CErr("EOMCC: Orbital index in input EOMCC.CVSCORE out of range");
      if (i >= nO)
        CErr("EOMCC: Virtual orbital appears in input EOMCC.CVSCORE");
      if (std::find(ccSettings.frozen_occupied.begin(), ccSettings.frozen_occupied.end(), i) != ccSettings.frozen_occupied.end() )
        CErr("EOMCC: Orbital appear in both EOMCC.CVSCORE and EOMCC.FROZENOCCUPIED inputs");
    }

    for (auto it = ccSettings.frozen_virtual.rbegin(); it != ccSettings.frozen_virtual.rend(); it++) {
      if (*it >= nO + nV)
        CErr("EOMCC: Orbital index in input EOMCC.FROZENVRITUAL out of range");
      if (*it < nO)
        CErr("EOMCC: Occupied orbital appears in input EOMCC.FROZENVRITUAL");
    }

    for (size_t i : eomSettings.external_virtual) {
      if (i >= nO + nV)
        CErr("EOMCC: Orbital index in input EOMCC.CVSVIRTUAL out of range");
      if (i < nO)
        CErr("EOMCC: Occupied orbital appears in input EOMCC.CVSVIRTUAL");
      if (std::find(ccSettings.frozen_virtual.begin(), ccSettings.frozen_virtual.end(), i) != ccSettings.frozen_virtual.end() )
        CErr("EOMCC: Orbital appear in both EOMCC.CVSVIRTUAL and EOMCC.FROZENVRITUAL inputs");
    }


    //if (!eomSettings.contain_active_space() && ! eomSettings.containActive()) {
    //  std::vector<size_t> V_blk = LinearRange(nV, 0, blksize);
    //  std::vector<size_t> O_blk = LinearRange(nO, 0, blksize);
    //  TAmanager.addRangeType(VLabel, TA::TiledRange1(V_blk.begin(), V_blk.end()));
    //  TAmanager.addRangeType(OLabel, TA::TiledRange1(O_blk.begin(), O_blk.end()));
    //}
    //else {
    // reorder mo by space
      eomSettings.nO = nO;
      // full space calculation
      if (eomSettings.cvs_core.empty() and eomSettings.active_virtual.empty() and eomSettings.active_occupied.empty()){
        for (size_t i = 0; i < nO; i++){
          if (std::find(ccSettings.frozen_occupied.begin(), ccSettings.frozen_occupied.end(), i) == ccSettings.frozen_occupied.end())
            eomSettings.active_occupied.push_back(i);
        }
        for (size_t i = nO; i < nMO; i++){
          if (std::find(ccSettings.frozen_virtual.begin(), ccSettings.frozen_virtual.end(), i) == ccSettings.frozen_virtual.end())
            eomSettings.active_virtual.push_back(i);
        }
      }
      // CVS calculations
      if (not eomSettings.cvs_core.empty() and eomSettings.active_virtual.empty() and eomSettings.active_occupied.empty()){
        // CVS calculation, set all nonfrozen virtual to external virtual
        for (size_t i = nO; i < nMO; i++){
          if (std::find(ccSettings.frozen_virtual.begin(), ccSettings.frozen_virtual.end(), i) == ccSettings.frozen_virtual.end())
            eomSettings.external_virtual.push_back(i);
        }
        // set active occupied space to non CVS occupied space
        for (size_t i = 0; i < nO; i++){
          if (std::find(ccSettings.frozen_occupied.begin(), ccSettings.frozen_occupied.end(), i) == ccSettings.frozen_occupied.end() and
              std::find(eomSettings.cvs_core.begin(), eomSettings.cvs_core.end(), i) == eomSettings.cvs_core.end() )
            eomSettings.active_occupied.push_back(i);
        }
      }
      // active calculations
      if (eomSettings.cvs_core.empty() and not eomSettings.active_occupied.empty() and not eomSettings.active_virtual.empty()) {
        // set cvs core to non active occupied space
        for (size_t i = 0; i < nO; i++) {
          if (std::find(ccSettings.frozen_occupied.begin(), ccSettings.frozen_occupied.end(), i) == ccSettings.frozen_occupied.end() and
              std::find(eomSettings.active_occupied.begin(), eomSettings.active_occupied.end(), i) == eomSettings.active_occupied.end() )
            eomSettings.cvs_core.push_back(i);
        }
        for (size_t i = nO; i < nMO; i++){
          if (std::find(ccSettings.frozen_virtual.begin(), ccSettings.frozen_virtual.end(), i) == ccSettings.frozen_virtual.end() and
              std::find(eomSettings.active_virtual.begin(), eomSettings.active_virtual.end(), i) == eomSettings.active_virtual.end() )
            eomSettings.external_virtual.push_back(i);
        }
      }

      int nFZC = ccSettings.frozen_occupied.size();
      int nCVSCore = eomSettings.cvs_core.size();
      int nCVSOValence = eomSettings.active_occupied.size(); 
      int nCVSVValence = eomSettings.active_virtual.size();
      int nCVSVirtual = eomSettings.external_virtual.size();
      int nFZV = ccSettings.frozen_virtual.size();

      nOcc = nO - nFZC;
      nVir = nV - nFZV;

      if (nCVSCore != nO || nCVSVirtual != nV ) {
        MatsT * mo_by_space = CQMemManager::get().malloc<MatsT>(nMO * nAO * 2);
        for (size_t i = 0; i < nFZC; i++) {
          size_t mo_index = ccSettings.frozen_occupied[i];
          memcpy(mo_by_space + i * nMO, mo + mo_index * nMO, nMO * sizeof(MatsT));
        }
        for (size_t i = nFZC; i < nFZC+nCVSCore; i++) {
          size_t mo_index = eomSettings.cvs_core[i-nFZC];
          memcpy(mo_by_space + i * nMO, mo + mo_index * nMO, nMO * sizeof(MatsT));
        }
        for (size_t i = nFZC+nCVSCore; i < nO; i++) {
          size_t mo_index = eomSettings.active_occupied[i-nFZC-nCVSCore];
          memcpy(mo_by_space + i * nMO, mo + mo_index * nMO, nMO * sizeof(MatsT));
        }
        for (size_t i = nO; i < nO+nCVSVValence; i++) {
          size_t mo_index = eomSettings.active_virtual[i-nO];
          memcpy(mo_by_space + i * nMO, mo + mo_index * nMO, nMO * sizeof(MatsT));
        }
        for (size_t i = nO+nCVSVValence; i < nMO-nFZV; i++) {
          size_t mo_index = eomSettings.external_virtual[i-nO-nCVSVValence];
          memcpy(mo_by_space + i * nMO, mo + mo_index * nMO, nMO * sizeof(MatsT));
        }
        for (size_t i = nMO-nFZV; i < nMO; i++) {
          size_t mo_index = ccSettings.frozen_virtual[i+nFZV-nMO];
          memcpy(mo_by_space + i * nMO, mo + mo_index * nMO, nMO * sizeof(MatsT));
        }

        memcpy(mo, mo_by_space, nMO * nMO * sizeof(MatsT));
        CQMemManager::get().free(mo_by_space);
      }

      // reorder CVS orbital indicies

      ccSettings.frozen_occupied.clear();
      eomSettings.cvs_core.clear();
      eomSettings.external_virtual.clear();
      ccSettings.frozen_virtual.clear();
      for (size_t i = 0; i < nFZC; i++) ccSettings.frozen_occupied.push_back(i);
      for (size_t i = nFZC; i < nFZC+nCVSCore; i++) eomSettings.cvs_core.push_back(i);
      for (size_t i = nO+nCVSVValence; i < nMO-nFZV; i++) eomSettings.external_virtual.push_back(i);
      for (size_t i = nMO-nFZV; i < nMO; i++) ccSettings.frozen_virtual.push_back(i);

      // define occupied and virtual TA ranges by subspace
      std::vector<size_t> V_blk, O_blk, v_blk, o_blk;
      std::vector<size_t> deep_blk, core_blk, homo_blk, lumo_blk, rydb_blk, free_blk;
      if (nFZC        ) deep_blk = LinearRange(nFZC, 0, blksize);
      if (nFZV        ) free_blk = LinearRange(nFZV, nCVSVirtual+nCVSVValence, blksize);

      if (nCVSCore    ) core_blk = LinearRange(nCVSCore, 0, blksize);
      if (nCVSOValence) homo_blk = LinearRange(nCVSOValence, nCVSCore, blksize);
      if (nCVSVValence) lumo_blk = LinearRange(nCVSVValence, 0, blksize);
      if (nCVSVirtual ) rydb_blk = LinearRange(nCVSVirtual, nCVSVValence, blksize);

      O_blk.insert(O_blk.end(), deep_blk.begin(), deep_blk.end());
      O_blk.insert(O_blk.end(), core_blk.begin(), core_blk.end());
      O_blk.insert(O_blk.end(), homo_blk.begin(), homo_blk.end());
      O_blk.push_back(nO);
      V_blk.insert(V_blk.end(), lumo_blk.begin(), lumo_blk.end());
      V_blk.insert(V_blk.end(), rydb_blk.begin(), rydb_blk.end());
      V_blk.insert(V_blk.end(), free_blk.begin(), free_blk.end());
      V_blk.push_back(nV);
      o_blk.insert(o_blk.end(), core_blk.begin(), core_blk.end());
      o_blk.insert(o_blk.end(), homo_blk.begin(), homo_blk.end());
      o_blk.push_back(nCVSCore + nCVSOValence);
      v_blk.insert(v_blk.end(), lumo_blk.begin(), lumo_blk.end());
      v_blk.insert(v_blk.end(), rydb_blk.begin(), rydb_blk.end());
      v_blk.push_back(nCVSVValence + nCVSVirtual);

      if (nFZC) {
        deep_blk.push_back(nFZC);
        TAmanager.addRangeType(dLabel, TA::TiledRange1(deep_blk.begin(),deep_blk.end()));
      }

      //TAmanager.addRangeType(VLabel, TA::TiledRange1(V_blk.begin(),V_blk.end()));
      //TAmanager.addRangeType(OLabel, TA::TiledRange1(O_blk.begin(),O_blk.end()));
      TAmanager.addRangeType(vLabel, TA::TiledRange1(v_blk.begin(),v_blk.end()));
      TAmanager.addRangeType(oLabel, TA::TiledRange1(o_blk.begin(),o_blk.end()));

      //// tiles within the full occ/vir space
      //TAmanager.addBlockRangeType(dLabel, 0, deep_blk.size()); 
      //TAmanager.addBlockRangeType(oLabel, deep_blk.size(), O_blk.size()-1); 
      //TAmanager.addBlockRangeType(fLabel, homo_blk.size()+rydb_blk.size(), V_blk.size()-1); 
      //TAmanager.addBlockRangeType(vLabel, 0, V_blk.size()-free_blk.size()-1); 

      // tiles within active occ/vir space
      TAmanager.addBlockRangeType(cLabel, 0, core_blk.size()); 
      TAmanager.addBlockRangeType(hLabel, core_blk.size(), o_blk.size()-1); 
      TAmanager.addBlockRangeType(oLabel, 0, o_blk.size()-1); 
      TAmanager.addBlockRangeType(lLabel, 0, lumo_blk.size());
      TAmanager.addBlockRangeType(rLabel, lumo_blk.size(), v_blk.size()-1); 
      if (core_blk.size()) core_blk.push_back(nCVSCore); 
      if (homo_blk.size()) homo_blk.push_back(nCVSOValence + nCVSCore); 
      if (lumo_blk.size()) lumo_blk.push_back(nCVSVValence); 
      if (rydb_blk.size()) rydb_blk.push_back(nCVSVirtual + nCVSVValence);
      for (auto it = homo_blk.begin(); it < homo_blk.end(); it++) *it -= nCVSCore;
      for (auto it = rydb_blk.begin(); it < rydb_blk.end(); it++) *it -= nCVSVValence;
      if (core_blk.size()) TAmanager.addRangeType(cLabel, TA::TiledRange1(core_blk.begin(), core_blk.end())); 
      if (homo_blk.size()) TAmanager.addRangeType(hLabel, TA::TiledRange1(homo_blk.begin(), homo_blk.end())); 
      if (lumo_blk.size()) TAmanager.addRangeType(lLabel, TA::TiledRange1(lumo_blk.begin(), lumo_blk.end())); 
      if (rydb_blk.size()) TAmanager.addRangeType(rLabel, TA::TiledRange1(rydb_blk.begin(), rydb_blk.end())); 

    //}

    try {
      const InCoreRITPI<IntsT> &aoRITPI = dynamic_cast<const InCoreRITPI<IntsT>&>(aoTPI);
      isRI = true;
      TAmanager.addRangeType(auxLabel, TAERI<IntsT>::LinRange(aoRITPI.nRIBasis(), blksize));

      // So DFCCSD knows how many RI basis functions
      nRI = aoRITPI.nRIBasis();
    } catch (const std::bad_cast&) {}

    std::map<std::string,TArray> ao2mo;
    std::vector<std::string> ao2moTypes{"ad","bd","ao","bo","av","bv"};
    for(const auto& ao2moType : ao2moTypes){

      std::vector<size_t> offset(2, 0);
      offset[0] = ao2moType[0] == 'b' ? nAO : 0;
      switch ( ao2moType[1] ) {
          case 'd': offset[1] = 0; break;
          case 'o': offset[1] = ccSettings.frozen_occupied.size(); break;
          case 'v': offset[1] = nO; break;
      }

      std::string rangeStr(ao2moType);
      rangeStr[0] = 'a';
      if (nFZC || rangeStr[1] != 'd') {
        TArray tmp = TAmanager.malloc_fresh<dcomplex>(rangeStr);

        tmp.init_elements([mo, offset, nMO](const typename TArray::index &i){
          return mo[i[0] + offset[0] + (i[1] + offset[1]) * nMO];
        });

        ao2mo[ao2moType] = tmp;
      }
    }
#ifdef DEBUG_CCSD
    prettyPrintSmart(std::cout, "MO", mo, nMO, nMO, nMO);
    if (nFZC) std::cout << "ao2mo[ad]:" << ao2mo["ad"] << std::endl;
    if (nFZC) std::cout << "ao2mo[bd]:" << ao2mo["bd"] << std::endl;
    std::cout << "ao2mo[ao]:" << ao2mo["ao"] << std::endl;
    std::cout << "ao2mo[bo]:" << ao2mo["bo"] << std::endl;
    std::cout << "ao2mo[av]:" << ao2mo["av"] << std::endl;
    std::cout << "ao2mo[bv]:" << ao2mo["bv"] << std::endl;
#endif

    // Create MO TPI
    TArray aoTPIta = taERI.template generateAOERI<MatsT>();
#ifdef DEBUG_CCSD
    std::cout << "aoTPIta:" << std::endl << aoTPIta << std::endl;
#endif

    if (isRI) { // RI 3-index aoTPIta
      std::array<char, 3> occvir{'d', 'o', 'v'};
      for ( char p : occvir)
        for ( char q : occvir) {
          if (!nFZC && (p=='d' || q=='d')) continue;
          std::string k = std::string("b")+p+q;
          riMoInts[k] = TAmanager.malloc<dcomplex>(k);
          riMoInts[k]("L,p,q")  = aoTPIta("L,m,n") * conj(ao2mo[std::string("a")+p]("m,p")) * ao2mo[std::string("a")+q]("n,q");
          riMoInts[k]("L,p,q") += aoTPIta("L,m,n") * conj(ao2mo[std::string("b")+p]("m,p")) * ao2mo[std::string("b")+q]("n,q");

        }
      TAmanager.free("baa", std::move(aoTPIta), true);

      // Construct oooo, vooo, and vovo in case rebuildFock is required
      if (rebuildFock || ccSettings.cctype != CC_TYPE::DFCCSD) {
        // oooo
        antiSymMoInts["oooo"] = TAmanager.malloc<dcomplex>("oooo");
        antiSymMoInts["oooo"]("p,r,q,s")  = riMoInts["boo"]("L,p,q") * riMoInts["boo"]("L,r,s");
        antiSymMoInts["oooo"]("p,q,r,s") -= antiSymMoInts["oooo"]("p,q,s,r");

        // vooo
        antiSymMoInts["vooo"] = TAmanager.malloc<dcomplex>("vooo");
        antiSymMoInts["vooo"]("p,r,q,s")  = riMoInts["bvo"]("L,p,q") * riMoInts["boo"]("L,r,s");
        antiSymMoInts["vooo"]("p,q,r,s") -= antiSymMoInts["vooo"]("p,q,s,r");

        // vovo
        antiSymMoInts["vovo"] = TAmanager.malloc<dcomplex>("vovo");
        antiSymMoInts["vovo"]("p,r,q,s")  = riMoInts["bvv"]("L,p,q") * riMoInts["boo"]("L,r,s");
        TArray tmpvoov = TAmanager.malloc<dcomplex>("voov");
        tmpvoov("p,r,q,s")  = riMoInts["bvo"]("L,p,q") * riMoInts["bov"]("L,r,s");
        antiSymMoInts["vovo"]("p,q,r,s") -= tmpvoov("p,q,s,r");
        TAmanager.free("voov", std::move(tmpvoov), true);
      }

      // Skip building these slices for DFCCSD because they will be built on the fly
      if (ccSettings.cctype != CC_TYPE::DFCCSD) {
        // vvoo
        antiSymMoInts["vvoo"] = TAmanager.malloc<dcomplex>("vvoo");
        antiSymMoInts["vvoo"]("p,r,q,s") = riMoInts["bvo"]("L,p,q") * riMoInts["bvo"]("L,r,s");
        antiSymMoInts["vvoo"]("p,q,r,s") -= antiSymMoInts["vvoo"]("p,q,s,r");

        // vvvo
        antiSymMoInts["vvvo"] = TAmanager.malloc<dcomplex>("vvvo");
        antiSymMoInts["vvvo"]("p,r,q,s") = riMoInts["bvv"]("L,p,q") * riMoInts["bvo"]("L,r,s");
        antiSymMoInts["vvvo"]("p,q,r,s") -= antiSymMoInts["vvvo"]("q,p,r,s");

        // vvvv
        antiSymMoInts["vvvv"] = TAmanager.malloc<dcomplex>("vvvv");
        antiSymMoInts["vvvv"]("p,r,q,s") = riMoInts["bvv"]("L,p,q") * riMoInts["bvv"]("L,r,s");
        antiSymMoInts["vvvv"]("p,q,r,s") -= antiSymMoInts["vvvv"]("p,q,s,r");
      }

      if (nFZC) {
        // dddd
        antiSymMoInts["dddd"] = TAmanager.malloc<dcomplex>("dddd");
        antiSymMoInts["dddd"]("p,r,q,s")  = riMoInts["bdd"]("L,p,q") * riMoInts["bdd"]("L,r,s");
        antiSymMoInts["dddd"]("p,q,r,s") -= antiSymMoInts["dddd"]("p,q,s,r");

        // dodo
        antiSymMoInts["dodo"] = TAmanager.malloc<dcomplex>("dodo");
        antiSymMoInts["dodo"]("p,r,q,s")  = riMoInts["bdd"]("L,p,q") * riMoInts["boo"]("L,r,s");
        TArray tmpdood = TAmanager.malloc<dcomplex>("dood");
        tmpdood("p,r,q,s")  = riMoInts["bdo"]("L,p,q") * riMoInts["bod"]("L,r,s");
        antiSymMoInts["dodo"]("p,q,r,s") -= tmpdood("p,q,s,r");
        TAmanager.free("dood", std::move(tmpdood), true);

        // vdod
        antiSymMoInts["vdod"] = TAmanager.malloc<dcomplex>("vdod");
        antiSymMoInts["vdod"]("p,r,q,s")  = riMoInts["bvo"]("L,p,q") * riMoInts["bdd"]("L,r,s");
        TArray tmpvddo = TAmanager.malloc<dcomplex>("vddo");
        tmpvddo("p,r,q,s")  = riMoInts["bvd"]("L,p,q") * riMoInts["bdo"]("L,r,s");
        antiSymMoInts["vdod"]("p,q,r,s") -= tmpvddo("p,q,s,r");
        TAmanager.free("vddo", std::move(tmpvddo), true);

        // vdvd
        antiSymMoInts["vdvd"] = TAmanager.malloc<dcomplex>("vdvd");
        antiSymMoInts["vdvd"]("p,r,q,s")  = riMoInts["bvv"]("L,p,q") * riMoInts["bdd"]("L,r,s");
        TArray tmpvddv = TAmanager.malloc<dcomplex>("vddv");
        tmpvddv("p,r,q,s")  = riMoInts["bvd"]("L,p,q") * riMoInts["bdv"]("L,r,s");
        antiSymMoInts["vdvd"]("p,q,r,s") -= tmpvddv("p,q,s,r");
        TAmanager.free("vddv", std::move(tmpvddv), true);
      }

      for ( char p : occvir)
        for ( char q : occvir) {
          if (!nFZC && (p=='d' || q=='d')) continue;
          std::string k = std::string("b")+p+q;
          if (ccSettings.cctype != CC_TYPE::DFCCSD) TAmanager.free(k, std::move(riMoInts[k]), true);
        }


    } else { // 4-index aoTPIta

    if (nFZC) {
      // dddd
      antiSymMoInts["dddd"] = TAmanager.malloc<dcomplex>("dddd");
      antiSymMoInts["dddd"]("p,r,q,s")  = aoTPIta("m,n,l,g") * conj(ao2mo["ad"]("m,p")) * ao2mo["ad"]("n,q") * conj(ao2mo["ad"]("l,r")) * ao2mo["ad"]("g,s");
      antiSymMoInts["dddd"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bd"]("m,p")) * ao2mo["bd"]("n,q") * conj(ao2mo["bd"]("l,r")) * ao2mo["bd"]("g,s");
      antiSymMoInts["dddd"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["ad"]("m,p")) * ao2mo["ad"]("n,q") * conj(ao2mo["bd"]("l,r")) * ao2mo["bd"]("g,s");
      antiSymMoInts["dddd"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bd"]("m,p")) * ao2mo["bd"]("n,q") * conj(ao2mo["ad"]("l,r")) * ao2mo["ad"]("g,s");
      antiSymMoInts["dddd"]("p,q,r,s") -= antiSymMoInts["dddd"]("p,q,s,r");
      // dodo
      antiSymMoInts["dodo"] = TAmanager.malloc<dcomplex>("dodo");
      antiSymMoInts["dodo"]("p,r,q,s")  = aoTPIta("m,n,l,g") * conj(ao2mo["ad"]("m,p")) * ao2mo["ad"]("n,q") * conj(ao2mo["ao"]("l,r")) * ao2mo["ao"]("g,s");
      antiSymMoInts["dodo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bd"]("m,p")) * ao2mo["bd"]("n,q") * conj(ao2mo["bo"]("l,r")) * ao2mo["bo"]("g,s");
      antiSymMoInts["dodo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["ad"]("m,p")) * ao2mo["ad"]("n,q") * conj(ao2mo["bo"]("l,r")) * ao2mo["bo"]("g,s");
      antiSymMoInts["dodo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bd"]("m,p")) * ao2mo["bd"]("n,q") * conj(ao2mo["ao"]("l,r")) * ao2mo["ao"]("g,s");
      TArray tmpdood = TAmanager.malloc<dcomplex>("dood");
      tmpdood("p,r,q,s")  = aoTPIta("m,n,l,g") * conj(ao2mo["ad"]("m,p")) * ao2mo["ao"]("n,q") * conj(ao2mo["ao"]("l,r")) * ao2mo["ad"]("g,s");
      tmpdood("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bd"]("m,p")) * ao2mo["bo"]("n,q") * conj(ao2mo["bo"]("l,r")) * ao2mo["bd"]("g,s");
      tmpdood("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["ad"]("m,p")) * ao2mo["ao"]("n,q") * conj(ao2mo["bo"]("l,r")) * ao2mo["bd"]("g,s");
      tmpdood("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bd"]("m,p")) * ao2mo["bo"]("n,q") * conj(ao2mo["ao"]("l,r")) * ao2mo["ad"]("g,s");
      antiSymMoInts["dodo"]("p,q,r,s") -= tmpdood("p,q,s,r");
      TAmanager.free("dood", std::move(tmpdood), true);
      // vdod
      antiSymMoInts["vdod"] = TAmanager.malloc<dcomplex>("vdod");
      antiSymMoInts["vdod"]("p,r,q,s")  = aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["ao"]("n,q") * conj(ao2mo["ad"]("l,r")) * ao2mo["ad"]("g,s");
      antiSymMoInts["vdod"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bo"]("n,q") * conj(ao2mo["bd"]("l,r")) * ao2mo["bd"]("g,s");
      antiSymMoInts["vdod"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["ao"]("n,q") * conj(ao2mo["bd"]("l,r")) * ao2mo["bd"]("g,s");
      antiSymMoInts["vdod"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bo"]("n,q") * conj(ao2mo["ad"]("l,r")) * ao2mo["ad"]("g,s");
      TArray tmpvddo = TAmanager.malloc<dcomplex>("vddo");
      tmpvddo("p,r,q,s")  = aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["ad"]("n,q") * conj(ao2mo["ad"]("l,r")) * ao2mo["ao"]("g,s");
      tmpvddo("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bd"]("n,q") * conj(ao2mo["bd"]("l,r")) * ao2mo["bo"]("g,s");
      tmpvddo("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["ad"]("n,q") * conj(ao2mo["bd"]("l,r")) * ao2mo["bo"]("g,s");
      tmpvddo("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bd"]("n,q") * conj(ao2mo["ad"]("l,r")) * ao2mo["ao"]("g,s");
      antiSymMoInts["vdod"]("p,q,r,s") -= tmpvddo("p,q,s,r");
      TAmanager.free("vddo", std::move(tmpvddo), true);
      // vdvd
      antiSymMoInts["vdvd"] = TAmanager.malloc<dcomplex>("vdvd");
      antiSymMoInts["vdvd"]("p,r,q,s")  = aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["av"]("n,q") * conj(ao2mo["ad"]("l,r")) * ao2mo["ad"]("g,s");
      antiSymMoInts["vdvd"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bv"]("n,q") * conj(ao2mo["bd"]("l,r")) * ao2mo["bd"]("g,s");
      antiSymMoInts["vdvd"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["av"]("n,q") * conj(ao2mo["bd"]("l,r")) * ao2mo["bd"]("g,s");
      antiSymMoInts["vdvd"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bv"]("n,q") * conj(ao2mo["ad"]("l,r")) * ao2mo["ad"]("g,s");
      TArray tmpvddv = TAmanager.malloc<dcomplex>("vddv");
      tmpvddv("p,r,q,s")  = aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["ad"]("n,q") * conj(ao2mo["ad"]("l,r")) * ao2mo["av"]("g,s");
      tmpvddv("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bd"]("n,q") * conj(ao2mo["bd"]("l,r")) * ao2mo["bv"]("g,s");
      tmpvddv("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["ad"]("n,q") * conj(ao2mo["bd"]("l,r")) * ao2mo["bv"]("g,s");
      tmpvddv("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bd"]("n,q") * conj(ao2mo["ad"]("l,r")) * ao2mo["av"]("g,s");
      antiSymMoInts["vdvd"]("p,q,r,s") -= tmpvddv("p,q,s,r");
      TAmanager.free("vddv", std::move(tmpvddv), true);
    }

    // oooo
    antiSymMoInts["oooo"] = TAmanager.malloc<dcomplex>("oooo");
    antiSymMoInts["oooo"]("p,r,q,s")  = aoTPIta("m,n,l,g") * conj(ao2mo["ao"]("m,p")) * ao2mo["ao"]("n,q") * conj(ao2mo["ao"]("l,r")) * ao2mo["ao"]("g,s");
    antiSymMoInts["oooo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bo"]("m,p")) * ao2mo["bo"]("n,q") * conj(ao2mo["bo"]("l,r")) * ao2mo["bo"]("g,s");
    antiSymMoInts["oooo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["ao"]("m,p")) * ao2mo["ao"]("n,q") * conj(ao2mo["bo"]("l,r")) * ao2mo["bo"]("g,s");
    antiSymMoInts["oooo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bo"]("m,p")) * ao2mo["bo"]("n,q") * conj(ao2mo["ao"]("l,r")) * ao2mo["ao"]("g,s");
    antiSymMoInts["oooo"]("p,q,r,s") -= antiSymMoInts["oooo"]("p,q,s,r");


    // vooo
    antiSymMoInts["vooo"] = TAmanager.malloc<dcomplex>("vooo");
    antiSymMoInts["vooo"]("p,r,q,s")  = aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["ao"]("n,q") * conj(ao2mo["ao"]("l,r")) * ao2mo["ao"]("g,s");
    antiSymMoInts["vooo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bo"]("n,q") * conj(ao2mo["bo"]("l,r")) * ao2mo["bo"]("g,s");
    antiSymMoInts["vooo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["ao"]("n,q") * conj(ao2mo["bo"]("l,r")) * ao2mo["bo"]("g,s");
    antiSymMoInts["vooo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bo"]("n,q") * conj(ao2mo["ao"]("l,r")) * ao2mo["ao"]("g,s");
    antiSymMoInts["vooo"]("p,q,r,s") -= antiSymMoInts["vooo"]("p,q,s,r");

    // vovo
    antiSymMoInts["vovo"] = TAmanager.malloc<dcomplex>("vovo");
    antiSymMoInts["vovo"]("p,r,q,s")  = aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["av"]("n,q") * conj(ao2mo["ao"]("l,r")) * ao2mo["ao"]("g,s");
    antiSymMoInts["vovo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bv"]("n,q") * conj(ao2mo["bo"]("l,r")) * ao2mo["bo"]("g,s");
    antiSymMoInts["vovo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["av"]("n,q") * conj(ao2mo["bo"]("l,r")) * ao2mo["bo"]("g,s");
    antiSymMoInts["vovo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bv"]("n,q") * conj(ao2mo["ao"]("l,r")) * ao2mo["ao"]("g,s");
    TArray tmpvoov = TAmanager.malloc<dcomplex>("voov");
    tmpvoov("p,r,q,s")  = aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["ao"]("n,q") * conj(ao2mo["ao"]("l,r")) * ao2mo["av"]("g,s");
    tmpvoov("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bo"]("n,q") * conj(ao2mo["bo"]("l,r")) * ao2mo["bv"]("g,s");
    tmpvoov("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["ao"]("n,q") * conj(ao2mo["bo"]("l,r")) * ao2mo["bv"]("g,s");
    tmpvoov("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bo"]("n,q") * conj(ao2mo["ao"]("l,r")) * ao2mo["av"]("g,s");
    antiSymMoInts["vovo"]("p,q,r,s") -= tmpvoov("p,q,s,r");
    TAmanager.free("voov", std::move(tmpvoov), true);

    // vvoo
    antiSymMoInts["vvoo"] = TAmanager.malloc<dcomplex>("vvoo");
    antiSymMoInts["vvoo"]("p,r,q,s")  = aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["ao"]("n,q") * conj(ao2mo["av"]("l,r")) * ao2mo["ao"]("g,s");
    antiSymMoInts["vvoo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bo"]("n,q") * conj(ao2mo["bv"]("l,r")) * ao2mo["bo"]("g,s");
    antiSymMoInts["vvoo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["ao"]("n,q") * conj(ao2mo["bv"]("l,r")) * ao2mo["bo"]("g,s");
    antiSymMoInts["vvoo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bo"]("n,q") * conj(ao2mo["av"]("l,r")) * ao2mo["ao"]("g,s");
    antiSymMoInts["vvoo"]("p,q,r,s") -= antiSymMoInts["vvoo"]("p,q,s,r");

    // vvvo
    antiSymMoInts["vvvo"] = TAmanager.malloc<dcomplex>("vvvo");
    antiSymMoInts["vvvo"]("p,r,q,s")  = aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["av"]("n,q") * conj(ao2mo["av"]("l,r")) * ao2mo["ao"]("g,s");
    antiSymMoInts["vvvo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bv"]("n,q") * conj(ao2mo["bv"]("l,r")) * ao2mo["bo"]("g,s");
    antiSymMoInts["vvvo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["av"]("n,q") * conj(ao2mo["bv"]("l,r")) * ao2mo["bo"]("g,s");
    antiSymMoInts["vvvo"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bv"]("n,q") * conj(ao2mo["av"]("l,r")) * ao2mo["ao"]("g,s");
    antiSymMoInts["vvvo"]("p,q,r,s") -= antiSymMoInts["vvvo"]("q,p,r,s");

    // vvvv
    antiSymMoInts["vvvv"] = TAmanager.malloc<dcomplex>("vvvv");
    antiSymMoInts["vvvv"]("p,r,q,s")  = aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["av"]("n,q") * conj(ao2mo["av"]("l,r")) * ao2mo["av"]("g,s");
    antiSymMoInts["vvvv"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bv"]("n,q") * conj(ao2mo["bv"]("l,r")) * ao2mo["bv"]("g,s");
    antiSymMoInts["vvvv"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["av"]("m,p")) * ao2mo["av"]("n,q") * conj(ao2mo["bv"]("l,r")) * ao2mo["bv"]("g,s");
    antiSymMoInts["vvvv"]("p,r,q,s") += aoTPIta("m,n,l,g") * conj(ao2mo["bv"]("m,p")) * ao2mo["bv"]("n,q") * conj(ao2mo["av"]("l,r")) * ao2mo["av"]("g,s");
    antiSymMoInts["vvvv"]("p,q,r,s") -= antiSymMoInts["vvvv"]("p,q,s,r");

    TAmanager.free("aaaa", std::move(aoTPIta), true);
    }
    for (auto ta : ao2mo) {
      std::string rangeStr(ta.first);
      rangeStr[0] = 'a';
      TAmanager.free(rangeStr, std::move(ta.second), true);
    }

    // Create MO Density matrics
    TArray moDen = TAmanager.template malloc_fresh<dcomplex>("oo");
    TArray moDen_dd;
    if (nFZC) {
        moDen_dd = TAmanager.template malloc_fresh<dcomplex>("dd");
        moDen_dd.init_elements([](const typename TArray::index &i){
          return i[0] == i[1] ? 1.0 : 0.0;
        });
    }
    moDen.init_elements([](const typename TArray::index &i) {
      return i[0] == i[1] ? 1.0 : 0.0;
    });


    std::map<std::string, TArray> twoeHta;

    if (rebuildFock) {
      /*
       * Rebuild the Fock matrix from coreH and ERI after ao2mo transformation
       * This block **will** be problematic for mmfX2C because the new Fock matrix
       * will not have all the 2-electron relativistic effects captured in the
       * original Fock matrix coming out of the 4c->2c transformation. Thus,
       * one should only use it with caution.
       */
      // Create MO H
      cqmatrix::Matrix<MatsT> moCoreH = aoCoreH.template spinGather<MatsT>().transform('N', mo, nMO, nMO);
      std::map<std::string, TArray> coreHta;

      // moCoreH.output(std::cout, "moCoreH", true);

      // Build Fock from coreH and TPI to TA blocks
      std::vector<std::string> onePTypes{"oo", "vo", "vv", "ov"};
      for (const auto &onePType: onePTypes) {

        std::vector<size_t> offset;

        for (const auto &otype: onePType) {
          if (otype == 'o') {
            offset.push_back(nFZC);
          } else {
            offset.push_back(nO);
          }
        }

        coreHta[onePType] = TAmanager.template malloc_fresh<dcomplex>(onePType);
        coreHta[onePType].init_elements([&moCoreH, offset](const typename TArray::index &i) {
          return moCoreH(i[0] + offset[0], i[1] + offset[1]);
        });

        fockMatrix[onePType] = TAmanager.template malloc<dcomplex>(onePType);
      }

#ifdef DEBUG_CCSD
      std::cout << "Hvv:" << coreHta["vv"] << std::endl;
      std::cout << "Hov:" << coreHta["ov"] << std::endl;
      std::cout << "Hvo:" << coreHta["vo"] << std::endl;
      std::cout << "Hoo:" << coreHta["oo"] << std::endl;
#endif

      fockMatrix["oo"]("p,q") = coreHta["oo"]("p,q") + antiSymMoInts["oooo"]("p,i,q,j") * moDen("i,j");
      fockMatrix["vo"]("p,q") = coreHta["vo"]("p,q") + antiSymMoInts["vooo"]("p,i,q,j") * moDen("i,j");
      fockMatrix["vv"]("p,q") = coreHta["vv"]("p,q") + antiSymMoInts["vovo"]("p,i,q,j") * moDen("i,j");
      fockMatrix["ov"]("p,q") = coreHta["ov"]("p,q") + conj(antiSymMoInts["vooo"]("q,j,p,i")) * moDen("i,j");

      // Free these slices as soon as fockMatrix is built
      if (ccSettings.cctype == CC_TYPE::DFCCSD) {
        // oooo still needed further down for rebuildFock with frozen core
        TAmanager.free("vooo", std::move(antiSymMoInts["vooo"]), true);
        TAmanager.free("vovo", std::move(antiSymMoInts["vovo"]), true);
      }

      if (nFZC) {
        coreHta["dd"] = TAmanager.template malloc_fresh<dcomplex>("dd");
        coreHta["dd"].init_elements([&moCoreH](const typename TArray::index &i){
          return moCoreH(i[0], i[1]);
        });
#ifdef DEBUG_CCSD
      if (nFZC) std::cout << "Hdd:" << coreHta["dd"] << std::endl;
#endif
        fockMatrix["dd"] = TAmanager.template malloc<dcomplex>("dd");
        fockMatrix["dd"]("p,q") = coreHta["dd"]("p,q") 
            + antiSymMoInts["dodo"]("p,i,q,j") * moDen("i,j") 
            + antiSymMoInts["dddd"]("p,i,q,j") * moDen_dd("i,j");
        //EG += 0.5 * (antiSymMoInts["dddd"]("i,k,j,l") * moDen_dd("i,j")).dot(moDen_dd("k,l")).get();
        //EG += 0.5 * (antiSymMoInts["dodo"]("i,k,j,l") * moDen_dd("i,j")).dot(moDen("k,l")).get();
        //EG += 0.5 * (antiSymMoInts["dodo"]("k,i,l,j") * moDen("i,j")).dot(moDen_dd("k,l")).get();
        fockMatrix["oo"]("p,q") += antiSymMoInts["dodo"]("i,p,j,q") * moDen_dd("i,j");
        fockMatrix["vo"]("p,q") += antiSymMoInts["vdod"]("p,i,q,j") * moDen_dd("i,j");
        fockMatrix["vv"]("p,q") += antiSymMoInts["vdvd"]("p,i,q,j") * moDen_dd("i,j");
        fockMatrix["ov"]("p,q") += conj(antiSymMoInts["vdod"]("q,j,p,i")) * moDen_dd("i,j");
      }

      for (auto ta: coreHta)
        TAmanager.free(ta.first, std::move(ta.second), true);
    } else {
      /*
       * Instead of rebuilding Fock matrix when employing frozen core approximation,
       * we should instead grab the Fock matrix from SingleSlater and slice it afterward
       * to obtain the appropriate spaces. Recomputing E_ref is unnecesssary because it
       * lives in SingleSlater too, but it can be useful to leave as is for checking.
       */

      cqmatrix::Matrix<MatsT> moFock = aoFock.template spinGather<MatsT>().transform('N', mo, nMO, nMO);

      std::vector<std::string> onePTypes{"oo", "vo", "vv", "ov"};
      // temporary fix just in case both rebuildFock and frozen core are in use
      if (nFZC) {
        onePTypes.push_back("dd");
      }
      for (const auto &onePType: onePTypes) {

        std::vector<size_t> offset;

        for (const auto &otype: onePType) {
          if (otype == 'd') {
            offset.push_back(0);
          } else if (otype == 'o') {
            offset.push_back(nFZC);
          } else {
            offset.push_back(nO);
          }
        }

        fockMatrix[onePType] = TAmanager.template malloc_fresh<dcomplex>(onePType);
        fockMatrix[onePType].init_elements([&moFock, offset](const typename TArray::index &i) {
          return moFock(i[0] + offset[0], i[1] + offset[1]);
        });
        
      }
      TA::get_default_world().gop.fence();
    }
#ifdef DEBUG_CCSD
      std::cout << "Fvv:" << fockMatrix["vv"] << std::endl;
      std::cout << "Fov:" << fockMatrix["ov"] << std::endl;
      std::cout << "Fvo:" << fockMatrix["vo"] << std::endl;
      std::cout << "Foo:" << fockMatrix["oo"] << std::endl;
      if(nFZC) std::cout << "Fdd:" << fockMatrix["dd"] << std::endl;
#endif


    // Create MO lenElectric multipoles
    MultipoleInts<MatsT> moMU = lenElectric.template spatialToSpinBlock<IntsT>().transform('N', mo, nMO, nMO);

    // Build Fock from coreH and TPI to TA blocks
    std::vector<std::string> onePTypes{"oo","vo","vv","ov"};
    for(const auto& onePType:onePTypes){

      std::vector<size_t> offset;

      for(const auto& otype:onePType){
        if (otype == 'o') {
          offset.push_back(nFZC);
        } else {
          offset.push_back(nO);
        }
      }

      for (size_t j = 0; j < 3; j++) {
        muMatrix[static_cast<char>('X' + j) + onePType] = TAmanager.template malloc_fresh<dcomplex>(onePType);

        muMatrix[static_cast<char>('X' + j) + onePType].init_elements([&moMU, offset, j](const typename TArray::index &i){
          return (*moMU[std::string()+static_cast<char>('X' + j)])(i[0] + offset[0], i[1] + offset[1]);
        });
      }
    }


    // Compute diagonal Fock (orbital energies)
    eps.clear();
    eps.resize(nMO-nFZC-nFZV, 0.0);
    std::vector<double> eps_d(nFZC, 0.0);
    if (nFZC) {
      foreach_inplace(fockMatrix["dd"],[&](TA::Tensor<MatsT> &tile) {

        const auto& lobound = tile.range().lobound();
        if (lobound[0] == lobound[1]) {
          const auto& upbound = tile.range().upbound();

          std::size_t x[] = {0, 0};
          for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
            x[1] = x[0];
            eps_d[x[0]] = std::real(tile[x]);
          }
        }
      });
    }
    foreach_inplace(fockMatrix["oo"],[&](TA::Tensor<MatsT> &tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] == lobound[1]) {
        const auto& upbound = tile.range().upbound();

        std::size_t x[] = {0, 0};
        for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
          x[1] = x[0];
          eps[x[0]] = std::real(tile[x]);
        }
      }
    });
    foreach_inplace(fockMatrix["vv"], [&](TA::Tensor<MatsT> &tile){

      const auto& lobound = tile.range().lobound();
      if (lobound[0] == lobound[1]) {
        const auto& upbound = tile.range().upbound();

        std::size_t x[] = {0, 0};
        for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
          x[1] = x[0];
          eps[nOcc + x[0]] = std::real(tile[x]);
        }
      }
    });
    TA::get_default_world().gop.fence();

    double *eps_copy = CQMemManager::get().malloc<double>(nMO-nFZV);
    std::copy_n(eps.data(), nMO-nFZC-nFZV, eps_copy);
    std::copy_n(eps_d.data(), nFZC, eps_copy+nMO-nFZC-nFZV);
    std::fill_n(eps.data(), nMO-nFZC-nFZV, double(0.0));
    std::fill_n(eps_d.data(), nFZC, double(0.0));
    MPIAllReduce(eps_copy, nMO-nFZC-nFZV, eps.data(), MPI_COMM_WORLD);
    MPIAllReduce(eps_copy+nMO-nFZC-nFZV, nFZC, eps_d.data(), MPI_COMM_WORLD);
    CQMemManager::get().free(eps_copy);

    TA::get_default_world().gop.fence();    

#ifdef DEBUG_CCSD
    for (size_t i = 0; i < nMO-nFZC-nFZV; i++) {
      std::cout << "Orbital " << i << " : " << eps[i] << std::endl;
    }
#endif
#ifdef DEBUG_CCSD
    if (nFZC) std::cout << "Fdd:" << fockMatrix["dd"] << std::endl;
    std::cout << "Fvv:" << fockMatrix["vv"] << std::endl;
    std::cout << "Fov:" << fockMatrix["ov"] << std::endl;
    std::cout << "Fvo:" << fockMatrix["vo"] << std::endl;
    std::cout << "Foo:" << fockMatrix["oo"] << std::endl;
#endif


    double EF = 0.0;
    double EF_fzc = 0.0;
    MatsT EG = 0.0;
    MatsT EG_fzc = 0.0;

    for (size_t i = 0; i < nFZC; i++)
      EF_fzc += eps_d[i];
    for (size_t i = 0; i < nOcc; i++)
      EF += eps[i];

    if (rebuildFock) {
      EG += 0.5 * (antiSymMoInts["oooo"]("i,k,j,l") * moDen   ("i,j")).dot(moDen   ("k,l")).get();
      TA::get_default_world().gop.fence();    
      if (nFZC) {
        EG_fzc += 0.5 * (antiSymMoInts["dodo"]("i,k,j,l") * moDen_dd("i,j")).dot(moDen   ("k,l")).get();
        TA::get_default_world().gop.fence();    
        EG_fzc += 0.5 * (antiSymMoInts["dodo"]("k,i,l,j") * moDen   ("i,j")).dot(moDen_dd("k,l")).get();
        TA::get_default_world().gop.fence();    
        EG_fzc += 0.5 * (antiSymMoInts["dddd"]("i,k,j,l") * moDen_dd("i,j")).dot(moDen_dd("k,l")).get();
        TA::get_default_world().gop.fence();

        // In DF-CCSD, free this slice as soon as fockMatrix is built
        if (ccSettings.cctype == CC_TYPE::DFCCSD) {
          TAmanager.free("oooo", std::move(antiSymMoInts["oooo"]), true);
        }
      }
    }
    else {
      cqmatrix::Matrix<MatsT> moTwoeH = aoTwoeH.template spinGather<MatsT>().transform('N', mo, nMO, nMO);
      //coreHta["oo"]("p,q") = fockMatrix["oo"]("p,q") - moTwoeH_TA(p,q); //coreH with relativistic folded in
                                                                        //build HF energy based on 1/2(coreH+F)
      for (size_t i = nFZC; i < nO; i++)
        EG += 0.5 * moTwoeH(i, i);
      for (size_t i = 0; i < nFZC; i++)
        EG_fzc += 0.5 * moTwoeH(i, i);
    }

    TAmanager.free("oo", std::move(moDen), true);
    if (nFZC) {
      TAmanager.free("dddd", std::move(antiSymMoInts["dddd"]), true);
      TAmanager.free("dodo", std::move(antiSymMoInts["dodo"]), true);
      TAmanager.free("vdod", std::move(antiSymMoInts["vdod"]), true);
      TAmanager.free("vdvd", std::move(antiSymMoInts["vdvd"]), true);
      antiSymMoInts.erase("dddd");
      antiSymMoInts.erase("dodo");
      antiSymMoInts.erase("vdod");
      antiSymMoInts.erase("vdvd");
      TAmanager.free("dd", std::move(moDen_dd), true);
      TAmanager.free("dd", std::move(fockMatrix["dd"]), true);
      fockMatrix.erase("dd");
    }

    //E_fzc = coreHta + EG_fzc
    E_fzc = EF_fzc - std::real(EG_fzc) + nucRepEnergy;
    E_ref = EF - std::real(EG) + E_fzc;

    // Build diagonal elements of moFock
    fockMatrix["oo_diag"] = TAmanager.template malloc_fresh<dcomplex>("oo");
    fockMatrix["oo_diag"].init_elements([this](const typename TArray::index &i){
      if (i[0] == i[1])
        return eps[i[0]];
      return 0.0;
    });
    fockMatrix["vv_diag"] = TAmanager.template malloc_fresh<dcomplex>("vv");
    fockMatrix["vv_diag"].init_elements([this](const typename TArray::index &i){
      if (i[0] == i[1])
        return eps[nOcc + i[0]];
      return 0.0;
    });


    D_abij = TAmanager.template malloc_fresh<dcomplex>("vvoo");
    D_abij.init_elements([this,&denomshift_](const typename TArray::index &i){
      return 1.0/(eps[i[2]] + eps[i[3]] - eps[i[0] + nOcc] - eps[i[1] + nOcc] - denomshift_);
    });

    D_ai = TAmanager.template malloc_fresh<dcomplex>("vo");
    D_ai.init_elements([this,&denomshift_](const typename TArray::index &i){
      return 1.0 / (eps[i[1]] - eps[i[0] + nOcc] - denomshift_);
    });

    if (cctype == CC_TYPE::CCSDT){
      std::vector<std::string> tmp;
      tmp.push_back(std::string({vLabel}));
      tmp.push_back(std::string({oLabel}));
      tmp.push_back(std::string("OneBody"));
      tmp.push_back(std::string({vLabel,vLabel}));
      tmp.push_back(std::string({oLabel,oLabel}));
      tmp.push_back(std::string("TwoBody"));
      tmp.push_back(std::string({vLabel,vLabel,vLabel}));
      tmp.push_back(std::string({oLabel,oLabel,oLabel}));
      tmp.push_back(std::string("ThreeBody"));
      T = std::make_shared<MBExpansion<dcomplex>>(tmp);
    } else if (cctype == CC_TYPE::CCSD or cctype == CC_TYPE::DFCCSD) {

      std::vector<std::string> tmp;
      tmp.push_back(std::string({vLabel}));
      tmp.push_back(std::string({oLabel}));
      tmp.push_back(std::string("OneBody"));
      tmp.push_back(std::string({vLabel,vLabel}));
      tmp.push_back(std::string({oLabel,oLabel}));
      tmp.push_back(std::string("TwoBody"));
      T = std::make_shared<MBExpansion<dcomplex>>(tmp);
    }

    TA::get_default_world().gop.fence();


    std::cout << "    * Initialize MO integrals for coupled cluster took "
              << std::setw(10) << std::right << std::setprecision(6) << std::fixed
              << tock(initIntStart) << " s." << std::endl;
  } // CCIntermediates::initializeIntegrals


  void runCoupledCluster(JobType jobType, Molecule &mol, std::shared_ptr<SingleSlaterBase> ss,
                         std::shared_ptr<IntegralsBase> aoints,
                         SafeFile &rstFile, CQInputFile &input, std::ostream &output) {
    int rank = MPIRank();

    // Convert ss dynamic type
    std::shared_ptr<SingleSlater<dcomplex,double>> ccref =
        std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss);

    if (not ccref) {
      CErr("CC only support complex-matrix real-integral reference.", output);
    }

    if (ccref->nC != 2) {
      CErr("Only GHF/X2C-CCSD is supported!", output);
    }

    if (not TA::initialized()) {
      // Initialize TA
      int argc = 1;
      char **argv = NULL;
      auto& world = madness::initialize(argc, argv, MPI_COMM_WORLD, GetNumThreads(), /* quiet = */ true);
      TA::initialize(argc, argv, world.mpi.comm());
      initialized_tiledarray(true);
    }
    SetLAThreads(1);
    TAManager::get().reset();


    // Read CC options
    CoupledClusterSettings ccSettings = CQCCOptions(output, input);

    // check for things that cannot be done
//    if ((ccSettings.cctype != CC_TYPE::CCSD and ccSettings.cctype != CC_TYPE::DFCCSD and ccSettings.save) or
//        (ccSettings.cctype != CC_TYPE::CCSD and ccSettings.cctype != CC_TYPE::DFCCSD and ccSettings.restart)) {
//      CErr("Only CCSD have save and restart implemented so far", output);
//    }


    // Initialize integrals
    CCIntermediates<dcomplex> intermediates;

    std::shared_ptr<MultipoleInts<double>> &aoMU =
        std::dynamic_pointer_cast<Integrals<double>>(aoints)->lenElectric;
    if (aoMU == nullptr) {
      aoMU = std::make_shared<MultipoleInts<double>>(ccref->nAlphaOrbital(), 3, true);
    }
    aoMU->broadcast();

    // Read EOMCC options
    EOMSettings eomSettings;
    if (jobType == JobType::EOMCC) {
      eomSettings = CQEOMCCOptions(output, input);
    }

//    if(ccSettings.frozen_occupied.size() && !ccSettings.rebuildFock) 
//        CErr("RebuildFock option must be used when frozen core orbitals exist.", output);
    //if(jobType == JobType::CC && eomSettings.contain_active_space()) 
    //    CErr("You have requested CC calculation, please remove the active space related keywords in the EOMCC section.", output);
    //if(jobType == JobType::EOMCC && ccSettings.frozen_occupied != eomSettings.frozen_occupied)
    //    CErr("Frozen core space in CC disagree with EOMCC section.", output);
    //if(jobType == JobType::EOMCC && ccSettings.frozen_virtual != eomSettings.frozen_virtual)
    //    CErr("Frozen virtual space in CC disagree with EOMCC section.", output);

    std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>> coreHAO, fockMatrixAO, twoeHAO;
    std::shared_ptr<cqmatrix::Matrix<dcomplex>> mo1;
    size_t nAO = ccref->nAlphaOrbital();
    size_t nMO = ss->nOrbital();
    double Eref_ = 0.0; // dummy holder to read data
    coreHAO = std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(nAO);
    fockMatrixAO = std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(nAO);
    twoeHAO = std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(nAO);
    mo1 = std::make_shared<cqmatrix::Matrix<dcomplex>>(nMO);
    if (MPIRank() == 0) {
      if (ccSettings.skipSCF) {
        // If SCF is skipped, read from file
        std::cout << "  * Reading quantities from Bin File." << std::endl;
        rstFile.readData("INTS/CORE_HAMILTONIAN", *coreHAO);
        rstFile.readData("SCF/FOCK", *fockMatrixAO);
        rstFile.readData("SCF/TWOEH", *twoeHAO);
        rstFile.readData("SCF/MO1", mo1->pointer());
        rstFile.readData("SCF/TOTAL_ENERGY", &Eref_); // read reference energy from file
      } else {
        // If SCF is not skipped, read from intermediates
        coreHAO = CQIntermediates::getInstance().getData<cqmatrix::PauliSpinorMatrices<dcomplex>>(
            "INTS/CORE_HAMILTONIAN");
        fockMatrixAO = CQIntermediates::getInstance().getData<cqmatrix::PauliSpinorMatrices<dcomplex>>("SCF/FOCK");
        twoeHAO = CQIntermediates::getInstance().getData<cqmatrix::PauliSpinorMatrices<dcomplex>>("SCF/TWOEH");
        mo1 = CQIntermediates::getInstance().getData<cqmatrix::Matrix<dcomplex>>("SCF/MO1");
        Eref_ = ss->totalEnergy; // read reference energy from SingleSlater object
      }
    }
    coreHAO->broadcast();
    fockMatrixAO->broadcast();
    twoeHAO->broadcast();
    mo1->broadcast();
    MPIBCast(Eref_, 0, MPI_COMM_WORLD);
//    std::shared_ptr<cqmatrix::PauliSpinorMatrices<dcomplex>> coreHAO, fockMatrixAO, twoeHAO;
//    std::shared_ptr<cqmatrix::Matrix<dcomplex>> mo1;
//    if (ccSettings.skipSCF) {
//      size_t nAO = ccref->nAlphaOrbital();
//      size_t nMO = ss->nOrbital();
//      coreHAO = std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(nAO);
//      fockMatrixAO = std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(nAO);
//      twoeHAO = std::make_shared<cqmatrix::PauliSpinorMatrices<dcomplex>>(nAO);
//      mo1 = std::make_shared<cqmatrix::Matrix<dcomplex>>(nMO);
//      if (MPIRank() == 0) {
//        std::cout << "  * Reading quantities from Bin File." << std::endl;
//        rstFile.readData("INTS/CORE_HAMILTONIAN", *coreHAO);
//        rstFile.readData("SCF/FOCK", *fockMatrixAO);
//        rstFile.readData("SCF/TWOEH", *twoeHAO);
//        rstFile.readData("SCF/MO1", mo1->pointer());
//      }
////      mo1->output(std::cout, "mo1 read 0", true);
//
//      coreHAO->broadcast();
//      fockMatrixAO->broadcast();
//      twoeHAO->broadcast();
//      mo1->broadcast();
//      CQIntermediates::getInstance().addData("INTS/CORE_HAMILTONIAN", coreHAO);
//      CQIntermediates::getInstance().addData("SCF/FOCK", fockMatrixAO);
//      CQIntermediates::getInstance().addData("SCF/TWOEH", twoeHAO);
//      CQIntermediates::getInstance().addData("SCF/MO1", mo1);
//    }
//
//    coreHAO = CQIntermediates::getInstance().getData<cqmatrix::PauliSpinorMatrices<dcomplex>>("INTS/CORE_HAMILTONIAN");
//    fockMatrixAO = CQIntermediates::getInstance().getData<cqmatrix::PauliSpinorMatrices<dcomplex>>("SCF/FOCK");
//    twoeHAO = CQIntermediates::getInstance().getData<cqmatrix::PauliSpinorMatrices<dcomplex>>("SCF/TWOEH");
//    mo1 = CQIntermediates::getInstance().getData<cqmatrix::Matrix<dcomplex>>("SCF/MO1");
//
    intermediates.initializeIntegrals(*coreHAO,
                                      *fockMatrixAO,
                                      *twoeHAO,
                                      *std::dynamic_pointer_cast<Integrals<double>>(aoints)->TPI,
                                      *aoMU,
                                      ccSettings,
                                      eomSettings, 
                                      mo1->pointer(),
                                      ccref->nO + ccSettings.nEvariation,
                                      ccref->nV - ccSettings.nEvariation,
                                      ccSettings.blksize, mol.nucRepEnergy, ccSettings.cctype,
                                      ccSettings.denomshift, ccSettings.pertT3, ccSettings.rebuildFock);


    // Create CC object
    std::shared_ptr<CCBase<dcomplex>> cc = nullptr;
    cc = intermediates.build_cc(rank == 0 ? rstFile : SafeFile(),
                               ccSettings);

    // Stupid redundant flag to override E_ref if using mmfX2C and frozen core
//#define DEBUG_EREF
#ifdef DEBUG_EREF
    std::cout << "  Reference energy from initializeIntegrals is " << intermediates.E_ref << std::endl;
    std::cout << "  Reference energy from SingleSlater object is " << Eref_ << std::endl;
#endif
    if ( not ccSettings.rebuildFock ) {
      intermediates.E_ref = Eref_;
    }

    if (MPIRank() == 0) {
      cc->printBanner(intermediates.E_ref);
    }
    cc->run();

    // CRCC job
    if(ccSettings.crcc == true){
      // Exit out if asking EOMCCSDT
      if (ccSettings.cctype != CC_TYPE::CCSD){
        CErr("Only CR-CC(2,3) can be run for now. Other options NYI.", output);
      }

      std::cout << BannerTop << std::endl;
      eomSettings.printEOMCCSettings(std::cout);
      std::cout << BannerMid << std::endl << std::endl;

      auto beginIntermediates = tick();

      // Build CC intermediates
      cc->buildIntermediates();
      
      // Create EOMCCSD object
      EOMCCSD<dcomplex> eomcc(rank == 0 ? rstFile : SafeFile(),
                                      intermediates, eomSettings, ccSettings);

      eomcc.initializeGroundStateLambda();

      // Build EOMCC intermediates
      eomcc.prepEOMCC();
      std::cout << "  * Form EOMCC Intermediates spent "
                << std::setw(10) << std::right << std::setprecision(6) << std::fixed
                << tock(beginIntermediates) << " s." << std::endl;

      std::cout << BannerMid << std::endl << std::endl;

      // Run CC lambda equations
      eomcc.runLambda();

      // Run CR-CC procedure
      eomcc.runCR(cc->CorrE);

    }

    // EOMCC job
    if(jobType == JobType::EOMCC){

      // Exit out if asking EOMCCSDT
//      if ((ccSettings.cctype != CC_TYPE::CCSD) && (ccSettings.cctype != CC_TYPE::DFCCSD)){
//        CErr("EOMCC can only be run at CCSD level of theory.", output);
//      }
      if ((ccSettings.cctype == CC_TYPE::CCSDT) && 
             (!((eomSettings.eom_type == EOM_TYPE::IP) && (eomSettings.ip_level == 3))) && 
             //(!((eomSettings.eom_type == EOM_TYPE::DIP) && (eomSettings.ip_level == 4)))) {
             (!(eomSettings.eom_type == EOM_TYPE::DIP) )) {
        CErr("CCSDT can be run with IP-EOMCCSDT(3h2p) and DIP-EOMCCSDT(4h2p) type only.", output);
      }
      
      //// Read EOMCC options
      //EOMSettings eomSettings = CQEOMCCOptions(output, input);

      std::cout << BannerTop << std::endl;
      eomSettings.printEOMCCSettings(std::cout);
      std::cout << BannerMid << std::endl << std::endl;

      if (eomSettings.eom_type != EOM_TYPE::EE && eomSettings.eom_type != EOM_TYPE::DIP && eomSettings.eom_type != EOM_TYPE::EA && eomSettings.eom_type != EOM_TYPE::IP){
        CErr("Only EE, EA, IP and DIP types of EOMCC are implemented.", output);
      }

      auto beginIntermediates = tick();

      // Build CC intermediates for EOM-EE or IP calcs (not using pq-generated intermediates)
      if (eomSettings.eom_type == EOM_TYPE::EE or eomSettings.eom_type == EOM_TYPE::IP) {
        cc->buildIntermediates();
      }
      cc->cleanMemory();

      // Create EOMCC object
      std::shared_ptr<EOMCCBase<dcomplex>> eomcc = nullptr;
      switch (static_cast<int>(eomSettings.eom_implementation)) {
        case static_cast<int>(EOM_IMPLEMENTATION::CVSEOMCCSD):
          eomcc = build_CVSEOMCCSD(
                  rank == 0 ? rstFile : SafeFile(),
                  intermediates, eomSettings, ccSettings);
          break;
        case static_cast<int>(EOM_IMPLEMENTATION::EOMCCSD):
          eomcc = build_EOMCCSD(
                  rank == 0 ? rstFile : SafeFile(),
                  intermediates, eomSettings, ccSettings);
          break;
        case static_cast<int>(EOM_IMPLEMENTATION::EOMDIP_3h1p):
          eomcc = build_EOMDIP_3h1p(
                  rank == 0 ? rstFile : SafeFile(),
                  intermediates, eomSettings, ccSettings);
          break;
        case static_cast<int>(EOM_IMPLEMENTATION::EOMDIP_4h2p):
          eomcc = build_EOMDIP_4h2pCCSDT(
                  rank == 0 ? rstFile : SafeFile(),
                  intermediates, eomSettings, ccSettings);
          break;
        case static_cast<int>(EOM_IMPLEMENTATION::EOMEA):
          eomcc = build_EOMEA(
                  rank == 0 ? rstFile : SafeFile(),
                  intermediates, eomSettings, ccSettings);
          break;
        case static_cast<int>(EOM_IMPLEMENTATION::EOMIP_2h1p):
          eomcc = build_EOMIP_2h1p(
                  rank == 0 ? rstFile : SafeFile(),
                  intermediates, eomSettings, ccSettings);
          break;
        case static_cast<int>(EOM_IMPLEMENTATION::EOMIP_3h2p):
          eomcc = build_EOMIP_3h2p(
                  rank == 0 ? rstFile : SafeFile(),
                  intermediates, eomSettings, ccSettings);
          break;
        default:
          CErr("Unknown EOM implementation");
      }

      if (eomSettings.oscillator_strength or eomSettings.diag_method == EOM_DIAG_METHOD::FULL) {
        eomcc->initializeGroundStateLambda();
      }
      // Build EOMCC intermediates
      eomcc->prepEOMCC();
      std::cout << "  * Form EOMCC Intermediates spent "
                << std::setw(10) << std::right << std::setprecision(6) << std::fixed
                << tock(beginIntermediates) << " s." << std::endl;

      std::cout << BannerMid << std::endl << std::endl;

      // Full diagonalization algorithm case
      if (eomSettings.diag_method == EOM_DIAG_METHOD::FULL) {
        eomcc->full_diagonalization();
      } else {
      // Davidson diagonalization case

        // Run CC lambda equations
        if (eomSettings.oscillator_strength) {
          eomcc->runLambda(); // has CVS version
        }

        eomcc->davidsonSolve();
      }


      std::cout << BannerEnd << std::endl;
    }


    // Reset LAThreads
    SetLAThreads(GetNumThreads());

  }

  template <typename MatsT>
  CCIntermediates<MatsT>::~CCIntermediates() {

    TAManager &TAmanager = TAManager::get();

    for (auto ta : fockMatrix)
      if (ta.second) TAmanager.free(ta.first.substr(0,2), std::move(ta.second), true);
    fockMatrix.clear();

    // Clean antiSymMoInts if it is not empty for accurate accounting
    if (!antiSymMoInts.empty()) {
      for (auto ta: antiSymMoInts)
        if (ta.second) TAmanager.free(ta.first, std::move(ta.second), true);
      antiSymMoInts.clear();
    }

    // Clean riMoInts if it is not empty for accurate accounting
    if (!riMoInts.empty()) {
      for (auto ta: riMoInts)
        if (ta.second) TAmanager.free(ta.first.substr(0, 3), std::move(ta.second), true);
      riMoInts.clear();
    }

    for (auto ta : muMatrix)
      if (ta.second) TAmanager.free(ta.first.substr(1), std::move(ta.second), true);
    muMatrix.clear();

    if (D_ai) TAmanager.free("vo", std::move(D_ai), true);
    if (D_abij) TAmanager.free("vvoo", std::move(D_abij), true);
    if (tau) TAmanager.free("vvoo", std::move(tau), true);
    if (tilde_tau) TAmanager.free("vvoo", std::move(tilde_tau), true);
    if (F_ae) TAmanager.free("vv", std::move(F_ae), true);
    if (F_mi) TAmanager.free("oo", std::move(F_mi), true);
    if (F_me) TAmanager.free("ov", std::move(F_me), true);
    if (W_mnij) TAmanager.free("oooo", std::move(W_mnij), true);
    if (W_abef) TAmanager.free("vvvv", std::move(W_abef), true);
    if (W_mbej) TAmanager.free("ovvo", std::move(W_mbej), true);
    if (W_mnie) TAmanager.free("ooov", std::move(W_mnie), true);
    if (W_amef) TAmanager.free("vovv", std::move(W_amef), true);
    if (W_mbij) TAmanager.free("ovoo", std::move(W_mbij), true);
    if (W_abei) TAmanager.free("vvvo", std::move(W_abei), true);
    if (G_ae) TAmanager.free("vv", std::move(G_ae), true);
    if (G_mi) TAmanager.free("oo", std::move(G_mi), true);

  }
}; // namespace ChronusQ
