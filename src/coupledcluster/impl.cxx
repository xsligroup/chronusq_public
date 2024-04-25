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
#include <coupledcluster/EOMCCSD.hpp>
#include <coupledcluster/CVSEOMCCSD.hpp>
#include <coupledcluster/EOMCCSDLambda.hpp>
#include <coupledcluster/EOMCCSDDensity.hpp>
#include <coupledcluster/EOMCCSDVectorImpl.hpp>
#include <coupledcluster/TAERI.hpp>
#include <singleslater.hpp>
#include <cxxapi/options.hpp>
#include <cxxapi/boilerplate.hpp>

namespace ChronusQ {

  template class CCSD<dcomplex,double>;
  template class EOMCCSD<dcomplex,double>;
  template class EOMCCSDVector<dcomplex>;
  template class EOMCCSDVectorSet<dcomplex>;
  template class EOMCCSDVectorSetDebug<dcomplex>;

  std::string orbitalSelectionToString(std::vector<size_t> orbitals);

  void EOMSettings::printEOMCCSettings(std::ostream &out) {
    out << "Equation of Motion Coupled Cluster (EOMCC) Settings:" << std::endl << std::endl;

    out << std::setw(45) << std::left << "  Type:"
        << "Singles and Doubles" << std::endl;

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

    out << std::setw(45) << std::left << "  Save Hbar matrix:"
        << (save_hamiltonian ? "ON" : "OFF") << std::endl;

    out << std::setw(45) << std::left << "  Core-Valence Separation approximation:"
        << (doCVS() ? "ON" : "OFF") << std::endl;

    if (doCVS()) {
      if (not frozen_occupied.empty())
        out << "    * Frozen Occupied orbitals  : " << orbitalSelectionToString(frozen_occupied) << std::endl;
      if (not cvs_core.empty())
        out << "    * CVS Core orbitals         : " << orbitalSelectionToString(cvs_core) << std::endl;
      if (not cvs_virtual.empty())
        out << "    * CVS Continuum orbitals    : " << orbitalSelectionToString(cvs_virtual) << std::endl;
      if (not frozen_virtual.empty())
        out << "    * Frozen Unoccupied orbitals: " << orbitalSelectionToString(frozen_virtual) << std::endl;
    }

    std::pair<double, char> mem_postfix = memSize(estimate_mem_peak());
    std::cout << std::setw(45) << std::left << "  Estimated TiledArray memory requirement: " << std::fixed << std::setprecision(1)
    << mem_postfix.first << mem_postfix.second << "B" << std::endl;

  }

  size_t EOMSettings::estimate_mem_peak() const {
    TAManager &TAmanager = TAManager::get();

    size_t count = 0;
    count += 8 * TAmanager.elem_per_TA("oo");
    count += 2 * TAmanager.elem_per_TA("oooo");
    count += 1 * TAmanager.elem_per_TA("ooov");
    count += 5 * TAmanager.elem_per_TA("ov");
    count += 1 * TAmanager.elem_per_TA("ovoo");
    count += 1 * TAmanager.elem_per_TA("ovvo");
    count += 6 * TAmanager.elem_per_TA("vo");
    count += 1 * TAmanager.elem_per_TA("vooo");
    count += 1 * TAmanager.elem_per_TA("vovo");
    count += 2 * TAmanager.elem_per_TA("vovv");
    count += 8 * TAmanager.elem_per_TA("vv");
    count += 5 * TAmanager.elem_per_TA("vvoo");
    count += 2 * TAmanager.elem_per_TA("vvvo");
    count += 2 * TAmanager.elem_per_TA("vvvv");

    size_t nVec = 0;
    if (hbar_type == EOM_HBAR_TYPE::EXPLICIT) {
      nVec += nroots * (oscillator_strength ? 2 : 0);
    } else {
      nVec += nroots * ((davidson_whenSc > 1 ? davidson_guess_multiplier : 1) *2
          + davidson_subspace_multiplier * 2 + davidson_guess_multiplier
          + (oscillator_strength ? 2 : 1));
    }
    count += nVec * TAmanager.elem_per_TA("vo");
    count += nVec * TAmanager.elem_per_TA("vvoo");

    if (oscillator_strength) {
      count += 1 * TAmanager.elem_per_TA("oo");
      count += 2 * TAmanager.elem_per_TA("ov");
      count += 1 * TAmanager.elem_per_TA("vv");
    }

    return count * sizeof(dcomplex);

  }

  template <typename MatsT>
  size_t CCIntermediates<MatsT>::estimate_mem_peak(size_t nDIIS) const {
    TAManager &TAmanager = TAManager::get();

    size_t count = 0;
    count += 8 * TAmanager.elem_per_TA("oo");
    count += 2 * TAmanager.elem_per_TA("oooo");
    count += 5 * TAmanager.elem_per_TA("ov");
    count += 1 * TAmanager.elem_per_TA("ovoo");
    count += 1 * TAmanager.elem_per_TA("ovvo");
    count += 7 * TAmanager.elem_per_TA("vo");
    count += 1 * TAmanager.elem_per_TA("vooo");
    count += 1 * TAmanager.elem_per_TA("vovo");
    count += 8 * TAmanager.elem_per_TA("vv");
    count += 7 * TAmanager.elem_per_TA("vvoo");
    count += 1 * TAmanager.elem_per_TA("vvvo");
    count += 2 * TAmanager.elem_per_TA("vvvv");

    if (nDIIS) {
      count += (nDIIS + 1) * 2 * TAmanager.elem_per_TA("vo");
      count += (nDIIS + 1) * 2 * TAmanager.elem_per_TA("vvoo");
    }
    return count * sizeof(MatsT);
  }

  template <typename MatsT>
  template <typename IntsT>
  void CCIntermediates<MatsT>::initializeIntegrals(const cqmatrix::PauliSpinorMatrices<MatsT> &aoCoreH,
                                                   const cqmatrix::PauliSpinorMatrices<MatsT> &aoFock,
                                                   const cqmatrix::PauliSpinorMatrices<MatsT> &aoTwoeH,
                                                   const TwoPInts<IntsT> &aoTPI,
                                                   const MultipoleInts<IntsT> &lenElectric,
                                                   MatsT *mo, size_t nO, size_t nV,
                                                   size_t blksize, double nucRepEnergy,
                                                   bool rebuildFock) {

    auto initIntStart = tick();

    nOcc = nO;
    nVir = nV;

    size_t nMO = nO + nV, nAO = nMO/2;

    TAManager &TAmanager = TAManager::get();

    // Initialize ranges
    TAERI<IntsT> taERI(aoTPI, blksize);
    TAmanager.addRangeType(aoLabel, taERI.getAOrange());
    TAmanager.addRangeType(vLabel, TAERI<IntsT>::LinRange(nV, blksize));
    TAmanager.addRangeType(oLabel, TAERI<IntsT>::LinRange(nO, blksize));

    std::map<std::string,TArray> ao2mo;
    std::vector<std::string> ao2moTypes{"ao","bo","av","bv"};
    for(const auto& ao2moType : ao2moTypes){

      std::vector<size_t> offset(2, 0);
      offset[0] = ao2moType[0] == 'b' ? nAO : 0;
      offset[1] = ao2moType[1] == 'v' ? nO : 0;

      std::string rangeStr(ao2moType);
      rangeStr[0] = 'a';
      TArray tmp = TAmanager.malloc_fresh<dcomplex>(rangeStr);

      tmp.init_elements([mo, offset, nMO](const typename TArray::index &i){
        return mo[i[0] + offset[0] + (i[1] + offset[1]) * nMO];
      });

      ao2mo[ao2moType] = tmp;
    }
#ifdef DEBUG_CCSD
    prettyPrintSmart(std::cout, "MO", mo, nMO, nMO, nMO);
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

    for (auto ta : ao2mo) {
      std::string rangeStr(ta.first);
      rangeStr[0] = 'a';
      TAmanager.free(rangeStr, std::move(ta.second), true);
    }
    TAmanager.free("aaaa", std::move(aoTPIta), true);

    // Create MO Density matrics
    TArray moDen = TAmanager.template malloc_fresh<dcomplex>("oo");
    moDen.init_elements([](const typename TArray::index &i) {
      return i[0] == i[1] ? 1.0 : 0.0;
    });


    std::map<std::string, TArray> twoeHta;

    if (rebuildFock) {
      // Create MO H
      cqmatrix::Matrix<MatsT> moCoreH = aoCoreH.template spinGather<MatsT>().transform('N', mo, nMO, nMO);
      std::map<std::string, TArray> coreHta;

//    moCoreH.output(std::cout, "moCoreH", true);

      // Build Fock from coreH and TPI to TA blocks
      std::vector<std::string> onePTypes{"oo", "vo", "vv", "ov"};
      for (const auto &onePType: onePTypes) {

        std::vector<size_t> offset;

        for (const auto &otype: onePType) {
          if (otype == 'o') {
            offset.push_back(0);
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

      for (auto ta: coreHta)
        TAmanager.free(ta.first, std::move(ta.second), true);
    } else {

      cqmatrix::Matrix<MatsT> moFock = aoFock.template spinGather<MatsT>().transform('N', mo, nMO, nMO);

      std::vector<std::string> onePTypes{"oo", "vo", "vv", "ov"};
      for (const auto &onePType: onePTypes) {

        std::vector<size_t> offset;

        for (const auto &otype: onePType) {
          if (otype == 'o') {
            offset.push_back(0);
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
#ifdef DEBUG_CCSD
      std::cout << "Fvv:" << fockMatrix["vv"] << std::endl;
      std::cout << "Fov:" << fockMatrix["ov"] << std::endl;
      std::cout << "Fvo:" << fockMatrix["vo"] << std::endl;
      std::cout << "Foo:" << fockMatrix["oo"] << std::endl;
#endif
    }


    // Create MO lenElectric multipoles
    MultipoleInts<MatsT> moMU = lenElectric.template spatialToSpinBlock<IntsT>().transform('N', mo, nMO, nMO);

//    moMU.getByOrder(1).output(std::cout, "moMU", true);

    // Build Fock from coreH and TPI to TA blocks
    std::vector<std::string> onePTypes{"oo","vo","vv","ov"};
    for(const auto& onePType:onePTypes){

      std::vector<size_t> offset;

      for(const auto& otype:onePType){
        if(otype == 'o'){
          offset.push_back(0);
        }
        else{
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
    eps.resize(nMO, 0.0);
    foreach_inplace(fockMatrix["oo"], [&](TA::Tensor<MatsT> &tile){

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
          eps[nO + x[0]] = std::real(tile[x]);
        }
      }
    });
    TA::get_default_world().gop.fence();
    TA::get_default_world().gop.template reduce(eps.data(), nMO, std::plus<double>());

#ifdef DEBUG_CCSD
    for (size_t i = 0; i < nO + nV; i++) {
      std::cout << "Orbital " << i << " : " << eps[i] << std::endl;
    }
#endif

    double EF = 0.0;
    MatsT EG = 0.0;

    for (size_t i = 0; i < nO; i++)
      EF += eps[i];

    if (rebuildFock) {
      EG = 0.5 * (antiSymMoInts["oooo"]("i,k,j,l") * moDen("i,j")).dot(moDen("k,l")).get();
      TA::get_default_world().gop.fence();
    } else {

      cqmatrix::Matrix<MatsT> moTwoeH = aoTwoeH.template spinGather<MatsT>().transform('N', mo, nMO, nMO);
      for (size_t i = 0; i < nO; i++)
        EG += 0.5 * moTwoeH(i, i);
    }
    TAmanager.free("oo", std::move(moDen), true);

    E_ref = EF - std::real(EG) + nucRepEnergy;

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

#ifdef DEBUG_CCSD
    std::cout << "Fvv:" << fockMatrix["vv"] << std::endl;
    std::cout << "Fov:" << fockMatrix["ov"] << std::endl;
    std::cout << "Fvo:" << fockMatrix["vo"] << std::endl;
    std::cout << "Foo:" << fockMatrix["oo"] << std::endl;
#endif

    // Compute denominators
    D_abij = TAmanager.template malloc_fresh<dcomplex>("vvoo");
    D_abij.init_elements([this](const typename TArray::index &i){
      return 1.0/(eps[i[2]] + eps[i[3]] - eps[i[0] + nOcc] - eps[i[1] + nOcc]);
    });

    D_ai = TAmanager.template malloc_fresh<dcomplex>("vo");
    D_ai.init_elements([this](const typename TArray::index &i){
      return 1.0 / (eps[i[1]] - eps[i[0] + nOcc]);
    });

    TA::get_default_world().gop.fence();

    std::cout << "    * Initialize MO integrals for coupled cluster took "
              << std::setw(10) << std::right << std::setprecision(6) << std::fixed
              << tock(initIntStart) << " s." << std::endl;
  } // CCIntermediates::initializeIntegrals


  std::vector<size_t> getGuessIndices(size_t nGuess, size_t length, const EOMSettings& eomSettings,
                                      const dcomplex *eomDiag, MPI_Comm comm) {//, bool sortByDistance = false) {

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

  void runCoupledCluster(JobType jobType, Molecule &mol, std::shared_ptr<SingleSlaterBase> ss,
                         std::shared_ptr<IntegralsBase> aoints,
                         SafeFile &rstFile, CQInputFile &input, std::ostream &output) {
    int rank = MPIRank();

    // Convert ss dynamic type
    std::shared_ptr<SingleSlater<dcomplex,double>> ccref =
        std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ss);

    if (not ccref) {
      CErr("CCSD only support complex-matrix real-integral reference.", output);
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

    // Initialize integrals
    CCIntermediates<dcomplex> intermediates;
    std::shared_ptr<MultipoleInts<double>> &aoMU =
        std::dynamic_pointer_cast<Integrals<double>>(aoints)->lenElectric;
    if (aoMU == nullptr) {
      aoMU = std::make_shared<MultipoleInts<double>>(ccref->nAlphaOrbital(), 3, true);
    }
    aoMU->broadcast();
    intermediates.initializeIntegrals(*ccref->coreH,
                                      *ccref->fockMatrix,
                                      *ccref->twoeH,
                                      *std::dynamic_pointer_cast<Integrals<double>>(aoints)->TPI,
                                      *aoMU,
                                      ccref->mo[0].pointer(),
                                      ccref->nO + ccSettings.nEvariation,
                                      ccref->nV - ccSettings.nEvariation,
                                      ccSettings.blksize, mol.nucRepEnergy,
                                      ccSettings.rebuildFock);
    intermediates.T = std::make_shared<EOMCCSDVector<dcomplex>>(intermediates.vLabel, intermediates.oLabel);

    // Create CCSD object
    CCSD<dcomplex, double> cc(rank == 0 ? rstFile : SafeFile(),
                              intermediates, ccSettings);
    if (MPIRank() == 0) {
      cc.printBanner(intermediates.E_ref);
    }
    cc.run();

    // EOMCC job
    if(jobType == JobType::EOMCC){

      // Read EOMCC options
      EOMSettings eomSettings = CQEOMCCOptions(output, input);

      std::cout << BannerTop << std::endl;
      eomSettings.printEOMCCSettings(std::cout);
      std::cout << BannerMid << std::endl << std::endl;

      auto beginIntermediates = tick();

      // Build CC intermediates
      cc.buildIntermediates();
      if (eomSettings.oscillator_strength or eomSettings.diag_method == EOM_DIAG_METHOD::FULL)
        intermediates.Lg = std::make_shared<EOMCCSDVector<dcomplex>>(intermediates.vLabel, intermediates.oLabel);

      // Create CCSD object
      EOMCCSD<dcomplex, double> eomcc(rank == 0 ? rstFile : SafeFile(),
                                      intermediates, eomSettings, ccSettings);

      // Build EOMCC intermediates
      eomcc.run();
      std::cout << "  * Form EOMCC Intermediates spent "
                << std::setw(10) << std::right << std::setprecision(6) << std::fixed
                << tock(beginIntermediates) << " s." << std::endl;

      std::cout << BannerMid << std::endl << std::endl;

      // Run CC lambda equations
      if (eomSettings.oscillator_strength and not eomSettings.doCVS())
        eomcc.runLambda();

      // Full diagonalization algorithm case
      if (eomSettings.diag_method == EOM_DIAG_METHOD::FULL) {
        eomcc.full_diagonalization();
        return;
      }

      size_t Hbar_dim = eomcc.getHbarDim();
      dcomplex* eomDiag = CQMemManager::get().malloc<dcomplex>(Hbar_dim);

      typename Davidson<dcomplex>::VecsGen_t vecsGenerator = Davidson<dcomplex>::VecsGen_t(); // Generator for new vector sets
      typename Davidson<dcomplex>::LinearTrans_t sigmaBuilder; // Sigma vector builder
      typename Davidson<dcomplex>::LinearTrans_t preConditioner; // Preconditioner

      /// Functions for Davidson
      EOMCCEigenVecType eigenVecType = EOMCCEigenVecType::RIGHT;

      size_t nGuess = eomSettings.nroots * eomSettings.davidson_guess_multiplier;
      nGuess = std::min(nGuess, eomcc.getHbarDim());
      dcomplex * curEig = CQMemManager::get().malloc<dcomplex>(nGuess);
      double PCsmall = eomSettings.davidson_preCond_small;

      // Algorithm with implicit Hbar matrix
      typename Davidson<dcomplex>::VecsGen_t vecsGenEOM;
      typename Davidson<dcomplex>::LinearTrans_t funcEOM;
      typename Davidson<dcomplex>::LinearTrans_t PCEOM;
      if (eomSettings.hbar_type == EOM_HBAR_TYPE::IMPLICIT
          or eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
        vecsGenEOM = [&intermediates](size_t nVec)->std::shared_ptr<SolverVectors<dcomplex>> {
          return std::make_shared<EOMCCSDVectorSet<dcomplex>>(intermediates.vLabel, intermediates.oLabel, nVec);
        }; // implicit vecsGenerator

        funcEOM = [&eomcc, &eigenVecType]( size_t nVec, SolverVectors<dcomplex> &V,
            SolverVectors<dcomplex> &AV) {

          EOMCCSDVectorSet<dcomplex> *V_ptr = nullptr, *AV_ptr = nullptr;
          size_t Vshift = 0, AVshift = 0;
          try {
            V_ptr = &dynamic_cast<EOMCCSDVectorSet<dcomplex>&>(V);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<dcomplex>& V_view = dynamic_cast<SolverVectorsView<dcomplex>&>(V);
            V_ptr = &dynamic_cast<EOMCCSDVectorSet<dcomplex>&>(V_view.getVecs());
            Vshift = V_view.shift();
          }

          try {
            AV_ptr = &dynamic_cast<EOMCCSDVectorSet<dcomplex>&>(AV);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<dcomplex>& AV_view = dynamic_cast<SolverVectorsView<dcomplex>&>(AV);
            AV_ptr = &dynamic_cast<EOMCCSDVectorSet<dcomplex>&>(AV_view.getVecs());
            AVshift = AV_view.shift();
          }

          for (size_t i = 0; i < nVec; i++) {
            const EOMCCSDVector<dcomplex> &Vi = V_ptr->get(i + Vshift);
            EOMCCSDVector<dcomplex> &AVi = AV_ptr->get(i + AVshift);
            eomcc.buildSigma(Vi.oneBody(), Vi.twoBody(), AVi.oneBody(), AVi.twoBody(), eigenVecType);
            AVi.enforceTwoBodySymmetry();
          }

        }; // implicit sigmaBuilder

        PCEOM = [&intermediates, &eomcc, eomDiag, curEig, PCsmall]( size_t nVec, SolverVectors<dcomplex> &V,
            SolverVectors<dcomplex> &AV) {

          AV.set_data(0, nVec, V, 0);

          EOMCCSDVectorSet<dcomplex> *AV_ptr = nullptr;
          size_t AVshift = 0;

          try {
            AV_ptr = &dynamic_cast<EOMCCSDVectorSet<dcomplex>&>(AV);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<dcomplex>& AV_view = dynamic_cast<SolverVectorsView<dcomplex>&>(AV);
            AV_ptr = &dynamic_cast<EOMCCSDVectorSet<dcomplex>&>(AV_view.getVecs());
            AVshift = AV_view.shift();
          }

          for (size_t iVec = 0; iVec < nVec; iVec++) {

            EOMCCSDVector<dcomplex> &curB = AV_ptr->get(iVec + AVshift);

            TA::foreach_inplace(curB.oneBody(), [iVec, curEig, eomDiag, &eomcc, PCsmall](TA::TensorZ &tile){
              const auto& lobound = tile.range().lobound();
              const auto& upbound = tile.range().upbound();

              dcomplex denom = 0.0;
              std::vector<std::size_t> x{0, 0};
              for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
                for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
                  denom = curEig[iVec] - eomDiag[eomcc.CVStoCompoundS(x[0], x[1])];
                  if (std::abs(denom) >= PCsmall) tile[x] /= denom;
                }
            });

            dcomplex *diagD = eomDiag + intermediates.nVir * intermediates.nOcc;
            TA::foreach_inplace(curB.twoBody(), [iVec, curEig, diagD, &eomcc, PCsmall](TA::TensorZ &tile){
              const auto& lobound = tile.range().lobound();
              const auto& upbound = tile.range().upbound();

              dcomplex denom = 0.0;
              std::vector<std::size_t> x{0, 0, 0, 0};
              for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0])
                for(x[1] = lobound[1]; x[1] < upbound[1]; ++x[1]) {
                  if (x[0] == x[1])
                    continue;
                  size_t a = x[0], b = x[1];
                  for(x[2] = lobound[2]; x[2] < upbound[2]; ++x[2])
                    for(x[3] = lobound[3]; x[3] < upbound[3]; ++x[3]) {
                      if (x[2] == x[3])
                        continue;
                      size_t i = x[2], j = x[3];
                      eomcc.signD(a,b,i,j);
                      denom = curEig[iVec] - diagD[eomcc.CVStoCompoundD(a, b, i, j)];
                      if (std::abs(denom) >= PCsmall) tile[x] /= denom;
                    }
                }
            });
            TA::get_default_world().gop.fence();

            curB.enforceTwoBodySymmetry();
          }
        }; // implicit preConditioner

        vecsGenerator = vecsGenEOM;
        sigmaBuilder = funcEOM;
        preConditioner = PCEOM;
      }

      // Algorithm with explicit Hbar matrix
      std::shared_ptr<cqmatrix::Matrix<dcomplex>> fullMat = nullptr;
      typename Davidson<dcomplex>::LinearTrans_t funcRaw;
      typename Davidson<dcomplex>::LinearTrans_t PCRaw;
      if (eomSettings.hbar_type == EOM_HBAR_TYPE::EXPLICIT
          or eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {

        std::cout << "  *** Start building the full matrix for explicit diagonalization ***" << std::endl;

        auto beginBuildHbar = tick();
        fullMat = std::make_shared<cqmatrix::Matrix<dcomplex>>(eomcc.buildHbarCVS(false));
        std::cout << "    * Build Hbar spent "
                  << std::setw(10) << std::right << std::setprecision(6) << std::fixed
                  << tock(beginBuildHbar) << " s." << std::endl;

        funcRaw = [&fullMat, &Hbar_dim, &eigenVecType]( size_t nVec, SolverVectors<dcomplex> &V,
            SolverVectors<dcomplex> &AV) {
          ROOT_ONLY(MPI_COMM_WORLD);

          size_t N = Hbar_dim;
          auto V_ptr = tryGetRawVectorsPointer(V);
          auto AV_ptr = tryGetRawVectorsPointer(AV);

          switch(eigenVecType) {
            case EOMCCEigenVecType::RIGHT:
              blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
                         N,nVec,N,dcomplex(1.),fullMat->pointer(),N,V_ptr,N,dcomplex(0.),AV_ptr,N);
              break;
            case EOMCCEigenVecType::LEFT:
              blas::gemm(blas::Layout::ColMajor,blas::Op::Trans,blas::Op::NoTrans,
                         nVec,N,N,dcomplex(1.),V_ptr,N,fullMat->pointer(),N,dcomplex(0.),AV_ptr,nVec);
              IMatCopy('T', nVec, N, 1.0, AV_ptr,nVec, N);
              break;
          }

        }; // explicit sigmaBuilder

        PCRaw = [Hbar_dim, eomDiag, curEig, PCsmall]( size_t nVec, SolverVectors<dcomplex> &V,
            SolverVectors<dcomplex> &AV) {

          ROOT_ONLY(MPI_COMM_WORLD);

          //            prettyPrintSmart(std::cout, "eomDiag", eomDiag, Hbar_dim, 1, Hbar_dim);

          dcomplex nom = 0.0, denom = 0.0;
          const dcomplex *Vptr = tryGetRawVectorsPointer(V);
          dcomplex *AVptr = tryGetRawVectorsPointer(AV);

          // Scale by inverse diagonals
          for (size_t i = 0; i < nVec; i++) {
            for(auto k = 0ul; k < Hbar_dim; k++ ) {
              nom = Vptr[k];
              denom = curEig[i] - eomDiag[k];
              if (std::abs(denom) >= PCsmall)
                AVptr[k] = nom / denom;
              else
                AVptr[k] = nom;
            }
            Vptr += Hbar_dim;
            AVptr += Hbar_dim;
          }

        }; // explicit preConditioner

        sigmaBuilder = funcRaw;
        preConditioner = PCRaw;

      }

      // Algorithm for debug, comparing implicit and explicit
      if (eomSettings.hbar_type == EOM_HBAR_TYPE::DEBUG) {
        typename Davidson<dcomplex>::VecsGen_t vecsGenRaw = [&Hbar_dim](size_t nVec)->std::shared_ptr<SolverVectors<dcomplex>> {
          return std::make_shared<RawVectors<dcomplex>>(
              MPI_COMM_WORLD, Hbar_dim, nVec
              );
        };

        typename Davidson<dcomplex>::VecsGen_t vecsGenDebug =
            [&intermediates, &eomcc](size_t nVec)->std::shared_ptr<SolverVectors<dcomplex>> {
          return std::make_shared<EOMCCSDVectorSetDebug<dcomplex>>(
              intermediates.vLabel, intermediates.oLabel, nVec, MPI_COMM_WORLD
              );
        };

        typename Davidson<dcomplex>::LinearTrans_t funcDebug = [&funcRaw, &funcEOM]( size_t nVec, SolverVectors<dcomplex> &V,
            SolverVectors<dcomplex> &AV) {

          EOMCCSDVectorSetDebug<dcomplex> *V_ptr = nullptr, *AV_ptr = nullptr;
          size_t Vshift = 0, AVshift = 0;

          try {
            V_ptr = &dynamic_cast<EOMCCSDVectorSetDebug<dcomplex>&>(V);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<dcomplex>& V_view = dynamic_cast<SolverVectorsView<dcomplex>&>(V);
            V_ptr = &dynamic_cast<EOMCCSDVectorSetDebug<dcomplex>&>(V_view.getVecs());
            Vshift = V_view.shift();
          }

          try {
            AV_ptr = &dynamic_cast<EOMCCSDVectorSetDebug<dcomplex>&>(AV);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<dcomplex>& AV_view = dynamic_cast<SolverVectorsView<dcomplex>&>(AV);
            AV_ptr = &dynamic_cast<EOMCCSDVectorSetDebug<dcomplex>&>(AV_view.getVecs());
            AVshift = AV_view.shift();
          }

          SolverVectorsView<dcomplex> V_EOM(V_ptr->getEOMCCSet(), Vshift);
          SolverVectorsView<dcomplex> V_Raw(V_ptr->getRawSet(), Vshift);
          SolverVectorsView<dcomplex> AV_EOM(AV_ptr->getEOMCCSet(), AVshift);
          SolverVectorsView<dcomplex> AV_Raw(AV_ptr->getRawSet(), AVshift);

          std::cout << "procedural.cxx::funcDebug before error = "
          << V_ptr->compareDebug(Vshift, nVec) << std::endl;

          funcEOM(nVec, V_EOM, AV_EOM);
          funcRaw(nVec, V_Raw, AV_Raw);

          std::cout << "procedural.cxx::funcDebug error = "
          << AV_ptr->compareDebug(AVshift, nVec) << std::endl;

        };

        typename Davidson<dcomplex>::LinearTrans_t PCDebug = [&PCRaw, &PCEOM]( size_t nVec, SolverVectors<dcomplex> &V,
            SolverVectors<dcomplex> &AV) {

          EOMCCSDVectorSetDebug<dcomplex> *V_ptr = nullptr, *AV_ptr = nullptr;
          size_t Vshift = 0, AVshift = 0;

          try {
            V_ptr = &dynamic_cast<EOMCCSDVectorSetDebug<dcomplex>&>(V);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<dcomplex>& V_view = dynamic_cast<SolverVectorsView<dcomplex>&>(V);
            V_ptr = &dynamic_cast<EOMCCSDVectorSetDebug<dcomplex>&>(V_view.getVecs());
            Vshift = V_view.shift();
          }

          try {
            AV_ptr = &dynamic_cast<EOMCCSDVectorSetDebug<dcomplex>&>(AV);
          } catch(const std::bad_cast& e) {
            SolverVectorsView<dcomplex>& AV_view = dynamic_cast<SolverVectorsView<dcomplex>&>(AV);
            AV_ptr = &dynamic_cast<EOMCCSDVectorSetDebug<dcomplex>&>(AV_view.getVecs());
            AVshift = AV_view.shift();
          }

          SolverVectorsView<dcomplex> V_EOM(V_ptr->getEOMCCSet(), Vshift);
          SolverVectorsView<dcomplex> V_Raw(V_ptr->getRawSet(), Vshift);
          SolverVectorsView<dcomplex> AV_EOM(AV_ptr->getEOMCCSet(), AVshift);
          SolverVectorsView<dcomplex> AV_Raw(AV_ptr->getRawSet(), AVshift);

          std::cout << "procedural.cxx::PCDebug before error = "
          << V_ptr->compareDebug(Vshift, nVec) << std::endl;

          PCEOM(nVec, V_EOM, AV_EOM);
          PCRaw(nVec, V_Raw, AV_Raw);

          std::cout << "procedural.cxx::PCDebug error = "
          << AV_ptr->compareDebug(AVshift, nVec) << std::endl;

        };

        vecsGenerator = vecsGenDebug;
        sigmaBuilder = funcDebug;
        preConditioner = PCDebug;
      }

      // Clear cached TA objects
      TAManager::get().discard_cache();

      std::cout << std::endl << "Davidson-Liu algorithm for EOMCC:" << std::endl;
      std::cout << BannerMid << std::endl << std::endl;

      // Compute diagonal elements
      auto beginBuildDiag = tick();
      eomcc.buildDiag(eomDiag);
      std::cout << "  * Build diagonal elements for iterative solver spent "
                << std::setw(10) << std::right << std::setprecision(6) << std::fixed
                << tock(beginBuildDiag) << " s." << std::endl << std::endl;

      std::cout << "Right eigensolver iterations:" << std::endl << std::endl;
      auto beginRightEig = tick();

      Davidson<dcomplex> davidson(MPI_COMM_WORLD, eomcc.getHbarDim(),
                                  eomSettings.davidson_max_macro_iter,
                                  eomSettings.davidson_max_micro_iter,
                                  eomSettings.davidson_residual_conv, eomSettings.nroots,
                                  sigmaBuilder, preConditioner, vecsGenerator);

      davidson.setWhenSc(eomSettings.davidson_whenSc);
      davidson.setM(eomSettings.davidson_subspace_multiplier);
      davidson.setkG(eomSettings.davidson_guess_multiplier);
      davidson.setGramSchmidtRepeat(eomSettings.GramSchmidt_NRe);
      davidson.setGramSchmidtEps(eomSettings.GramSchmidt_eps);
      davidson.setResidueConvCheck(eomSettings.davidson_check_residual);
      davidson.setEigenVectorConvCheck(eomSettings.davidson_check_eigen_vector);
      davidson.setEigenValueConvCheck(eomSettings.davidson_check_eigen_value);
      davidson.setConvOnGramSchmidt(eomSettings.davidson_conv_on_GramSchmidt);
      davidson.setEigenVectorConvCriteria(eomSettings.davidson_eigen_vector_conv);
      davidson.setEigenValueConvCriteria(eomSettings.davidson_eigen_value_conv);
      davidson.setEigForT(curEig);

      if (!eomSettings.davidson_Eref.empty())
        davidson.setEnergySpecific(eomSettings.davidson_Eref, eomSettings.davidson_ErefAbs);

//      if (eomSettings.davidson_Eref > 0.0) {
//        davidson.useEnergySpecific(eomSettings.nroots, 0.0, eomSettings.davidson_Eref);
//        if (eomSettings.davidson_sort_by_distance)
//          davidson.setSortByDistance();
//      }

//      bool sortByDistance = eomSettings.davidson_sort_by_distance;
      std::vector<size_t> guessIndices = getGuessIndices(nGuess, eomcc.getHbarDim(), eomSettings, eomDiag, MPI_COMM_WORLD);//, sortByDistance);

      davidson.setGuess(nGuess, [&guessIndices] (size_t nGuess, SolverVectors<dcomplex> &guessVec, size_t length) {
        guessVec.clear();
        for(size_t i = 0; i < nGuess; i++) guessVec.set(guessIndices[i], i, 1.0);
      });

      davidson.run();

      if (not davidson.hasConverged()) {
        CErr("EOMCC right Davidson iteration failed to converge!");
      }

      eomcc.setR(davidson.VR());

      std::cout << "  * EOMCC right eigensolver spent "
                << std::setw(10) << std::right << std::setprecision(6) << std::fixed
                << tock(beginRightEig) << " s." << std::endl;

      if (eomSettings.oscillator_strength) {
        eigenVecType = EOMCCEigenVecType::LEFT;

        std::cout << BannerMid << std::endl << std::endl;
        std::cout << "Left eigensolver iterations:" << std::endl << std::endl;
        auto beginLeftEig = tick();

        Davidson<dcomplex> davidsonLeft(MPI_COMM_WORLD, eomcc.getHbarDim(),
                                    eomSettings.davidson_max_macro_iter,
                                    eomSettings.davidson_max_micro_iter,
                                    eomSettings.davidson_residual_conv,
                                    eomSettings.nroots,
                                    sigmaBuilder, preConditioner, vecsGenerator);

        davidsonLeft.setWhenSc(eomSettings.davidson_whenSc);
        davidsonLeft.setM(eomSettings.davidson_subspace_multiplier);
        davidsonLeft.setkG(1);
        davidsonLeft.setGramSchmidtRepeat(eomSettings.GramSchmidt_NRe);
        davidsonLeft.setGramSchmidtEps(eomSettings.GramSchmidt_eps);
        davidsonLeft.setResidueConvCheck(eomSettings.davidson_check_residual);
        davidsonLeft.setEigenVectorConvCheck(eomSettings.davidson_check_eigen_vector);
        davidsonLeft.setEigenValueConvCheck(eomSettings.davidson_check_eigen_value);
        davidsonLeft.setConvOnGramSchmidt(eomSettings.davidson_conv_on_GramSchmidt);
        davidsonLeft.setEigenVectorConvCriteria(eomSettings.davidson_eigen_vector_conv);
        davidsonLeft.setEigenValueConvCriteria(eomSettings.davidson_eigen_value_conv);
        davidsonLeft.setEigForT(curEig);

        // Reuse scratch spaces for right Davidson
        davidsonLeft.setGuessScratch(davidson.getGuessScratch());
        davidsonLeft.setSubspaceScratch(davidson.getSubspaceScratch());
        davidsonLeft.setSigmaVecScratch(davidson.getSigmaVecScratch());
        davidsonLeft.setScratchR(davidson.getScratchR());
        davidsonLeft.setScratchS(davidson.getScratchS());

        davidson.clear_scratch();

        if (!eomSettings.davidson_Eref.empty())
          davidsonLeft.setEnergySpecific(eomSettings.davidson_Eref, eomSettings.davidson_ErefAbs);
//        if (eomSettings.davidson_Eref > 0.0) {
//          davidsonLeft.useEnergySpecific(eomSettings.nroots, 0.0, eomSettings.davidson_Eref);
//          if (eomSettings.davidson_sort_by_distance)
//            davidsonLeft.setSortByDistance();
//        }

        davidsonLeft.setGuess(eomSettings.nroots, [&davidson] (size_t nGuess, SolverVectors<dcomplex> &guessVec, size_t length) {
          guessVec.clear();
          guessVec.set_data(0, nGuess, *davidson.VR(), 0);
          guessVec.conjugate(0, nGuess);
        });

        davidsonLeft.run();
        if (not davidsonLeft.hasConverged()) {
          CErr("EOMCC left Davidson iteration failed to converge!");
        }
        davidsonLeft.clear_scratch();

        std::cout << "  * EOMCC left eigensolver spent "
                  << std::setw(10) << std::right << std::setprecision(6) << std::fixed
                  << tock(beginLeftEig) << " s." << std::endl;
        std::cout << BannerMid << std::endl;

        std::shared_ptr<EOMCCSDVectorSet<dcomplex>> VL, VR;

        switch (eomSettings.hbar_type) {
          case EOM_HBAR_TYPE::IMPLICIT:
            VR = std::dynamic_pointer_cast<EOMCCSDVectorSet<dcomplex>>(davidson.VR());
            VL = std::dynamic_pointer_cast<EOMCCSDVectorSet<dcomplex>>(davidsonLeft.VR());
            break;
          case EOM_HBAR_TYPE::DEBUG:
            VR = std::make_shared<EOMCCSDVectorSet<dcomplex>>(
                std::dynamic_pointer_cast<EOMCCSDVectorSetDebug<dcomplex>>(davidson.VR())->getEOMCCSet());
            VL = std::make_shared<EOMCCSDVectorSet<dcomplex>>(
                std::dynamic_pointer_cast<EOMCCSDVectorSetDebug<dcomplex>>(davidsonLeft.VR())->getEOMCCSet());
            break;
          case EOM_HBAR_TYPE::EXPLICIT:
            VR = std::make_shared<EOMCCSDVectorSet<dcomplex>>(intermediates.vLabel, intermediates.oLabel, eomSettings.nroots);
            VR->fromRaw(MPI_COMM_WORLD,
                        *std::dynamic_pointer_cast<RawVectors<dcomplex>>(davidson.VR()),
                        eomcc, false, 0, 0, eomSettings.nroots);
            VL = std::make_shared<EOMCCSDVectorSet<dcomplex>>(intermediates.vLabel, intermediates.oLabel, eomSettings.nroots);
            VL->fromRaw(MPI_COMM_WORLD,
                        *std::dynamic_pointer_cast<RawVectors<dcomplex>>(davidsonLeft.VR()),
                        eomcc, false, 0, 0, eomSettings.nroots);
            break;

        }

        eomcc.setR(VR);
        eomcc.setL(VL);
        eomcc.setTheta(davidson.eigVal(), eomSettings.nroots);

        eomcc.buildRightZeroBody(eomSettings.nroots);

        if (eomSettings.davidson_biortho) {
          // BiOrthonormalization
          std::cout << "  *** BiOrthonormalize left and right eigenvectors ***" << std::endl;

          biOrthoNormalize(eomSettings.nroots, *VL, *VR);
//          RawVectors<dcomplex> rawLex(VL->toRaw(MPI_COMM_WORLD, eomcc, true, 0, eomSettings.nroots));
//          RawVectors<dcomplex> rawRex(VR->toRaw(MPI_COMM_WORLD, eomcc, true, 0, eomSettings.nroots));
//
//          biOrthoNormalize(eomcc.getHbarDim(true), eomSettings.nroots, rawLex, rawRex);
//
//          VL->fromRaw(MPI_COMM_WORLD, rawLex, eomcc, true, 0, 0, eomSettings.nroots);
//          VR->fromRaw(MPI_COMM_WORLD, rawRex, eomcc, true, 0, 0, eomSettings.nroots);

        }

        std::cout << BannerMid << std::endl;

        std::cout << "EOMCC results:" << std::endl << std::endl;

        std::cout << std::setw(18) << std::left <<  "  Excited states";
        std::cout << std::setw(34) << std::left << "Excitation Energy (Eh)";
        std::cout << std::setw(19) << std::left << "Oscillator Strength";
        std::cout << std::endl;
        std::cout << std::setw(18) << std::left <<  "  -------------";
        std::cout << std::setw(34) << std::left << "-----------------";
        std::cout << std::setw(18) << std::left << "-----------------";
        std::cout << std::endl;

        auto beginOsc = tick();

        eomcc.initilizeDensity();
        std::cout << "----------------------------------------------" << std::endl;
        std::vector<double> oscStrength;
        std::vector<double> excitationE;
        for (size_t i = 0; i < eomSettings.nroots ; i++){
          dcomplex f = eomcc.calcOscillatorStrength(i);

          std::cout << std::setprecision(12) << std::fixed;
          std::cout << "      State "  << std::setw(6) << std::left << i+1;
          std::cout << std::setw(34) << std::left;
          if (std::abs(davidsonLeft.eigVal()[i]) > 1e-6)
            std::cout << std::fixed << std::setprecision(12);
          else
            std::cout << std::scientific << std::setprecision(6);
          std::cout << davidsonLeft.eigVal()[i];
          std::cout << std::setw(34) << std::left;
          if (std::abs(f) > 1e-6)
            std::cout << std::fixed << std::setprecision(12);
          else
            std::cout << std::scientific << std::setprecision(6);
          std::cout << f;
          std::cout << std::endl;
          
          oscStrength.push_back(std::real(f));
          excitationE.push_back(std::real(davidson.eigVal()[i]));
        }

        TA::get_default_world().gop.fence();
        // Write data to bin file
        if (rstFile.exists()) {
          rstFile.safeWriteData("/CC/EXCITATION_ENERGIES",excitationE.data(), {eomSettings.nroots});
          rstFile.safeWriteData("/CC/OSCILLATOR_STRENGTHS", oscStrength.data(), {eomSettings.nroots});
        }

        std::cout << std::endl << "  * Compute oscillator strength spent "
                  << std::setw(10) << std::right << std::setprecision(6) << std::fixed
                  << tock(beginOsc) << " s." << std::endl;

      }

      davidson.clear_scratch();

      CQMemManager::get().free(eomDiag);
      if (curEig) CQMemManager::get().free(curEig);

      std::cout << BannerEnd << std::endl;
    }


    // Reset LAThreads
    SetLAThreads(GetNumThreads());

  }

  template <typename MatsT>
  CCIntermediates<MatsT>::~CCIntermediates() {

    TAManager &TAmanager = TAManager::get();

    for (auto ta : fockMatrix)
      TAmanager.free(ta.first.substr(0,2), std::move(ta.second), true);
    fockMatrix.clear();

    for (auto ta : antiSymMoInts)
      TAmanager.free(ta.first, std::move(ta.second), true);
    antiSymMoInts.clear();

    for (auto ta : muMatrix)
      TAmanager.free(ta.first.substr(1), std::move(ta.second), true);
    fockMatrix.clear();

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
