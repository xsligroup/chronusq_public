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
#include <util/threads.hpp>
#include <coupledcluster.hpp>
#include <coupledcluster/TAERI.hpp>

namespace ChronusQ {

  std::vector<size_t> LinearRange(size_t size, size_t start_index, size_t blksize);


  template <>
  template <>
  void CCIntermediates<dcomplex>::initializeIntegrals_rhfRef(const cqmatrix::PauliSpinorMatrices<dcomplex> &aoCoreH,
                                                   const cqmatrix::PauliSpinorMatrices<dcomplex> &aoFock,
                                                   const cqmatrix::PauliSpinorMatrices<dcomplex> &aoTwoeH,
                                                   const TwoPInts<double> &aoTPI,
                                                   const MultipoleInts<double> &lenElectric,
                                                   CoupledClusterSettings& ccSettings,
                                                   EOMSettings& eomSettings,
                                                   dcomplex *mo, size_t nO, size_t nV,
                                                   size_t blksize, double nucRepEnergy,
                                                   CC_TYPE cctype, double denomshift_, bool pertT3, bool rebuildFock) {
    CErr("RCCSD is not compatible with complex-valued matrices. Please use CCSD instead.");
  }

  template <>
  template <>
  void CCIntermediates<double>::initializeIntegrals_rhfRef(const cqmatrix::PauliSpinorMatrices<double> &aoCoreH,
                                                   const cqmatrix::PauliSpinorMatrices<double> &aoFock,
                                                   const cqmatrix::PauliSpinorMatrices<double> &aoTwoeH,
                                                   const TwoPInts<double> &aoTPI,
                                                   const MultipoleInts<double> &lenElectric,
                                                   CoupledClusterSettings& ccSettings,
                                                   EOMSettings& eomSettings,
                                                   double *mo, size_t nO, size_t nV,
                                                   size_t blksize, double nucRepEnergy,
                                                   CC_TYPE cctype, double denomshift_, bool pertT3, bool rebuildFock) {

    auto initIntStart = tick();

    oLabel = 'o';
    vLabel = 'v';

    nO /= 2;
    nV /= 2; // In RHF reference, nO and nV read from SingleSlater are actually 2*nO and 2*nV,
             // so we need to divide by 2 here to get the correct number of spatial orbitals.
             // This is a bit hacky but it avoids changing the interface of initializeIntegrals().

    size_t nAO = nO + nV, nMO = nAO; // In RHF reference, nMO is equal to nAO, and both are equal to the number of spatial orbitals.

    TAManager &TAmanager = TAManager::get();

    // Initialize ranges
    nOcc = nO;
    nVir = nV;

    TAERI<double> taERI(aoTPI, blksize);
    TAmanager.addRangeType(aoLabel, taERI.getAOrange());

    // Validate the frozen and CVS orbital indices in input
    eomSettings.validateOrbitalSpaces(nO, nV, ccSettings);

    // Fill in unspecified occupied orbitals into empty categories
    eomSettings.assignUnspecifiedOrbitalToSpaces(nO, nV, ccSettings);

    eomSettings.nO = nO;

    int nFZC = ccSettings.frozen_occupied.size();
    int nCVSCore = eomSettings.cvs_core.size();
    int nCVSVirtual = eomSettings.external_virtual.size();
    int nFZV = ccSettings.frozen_virtual.size();

    if (nCVSCore > 0 or nCVSVirtual > 0) {
      CErr("CVS-RCCSD not yet implemented. Please use GCCSD instead.");
    }

    reorderMOs(mo, nO, nV, ccSettings, eomSettings);

    nOcc = nO - nFZC;
    nVir = nV - nFZV;

    std::vector<size_t> deep_blk;
    if (nFZC) {
      deep_blk = LinearRange(nFZC, 0, blksize);
      deep_blk.push_back(nFZC);
      TAmanager.addRangeType(dLabel, TA::TiledRange1(deep_blk.begin(),deep_blk.end()));
    }
    std::vector<size_t> v_blk = LinearRange(nVir, 0, blksize);
    v_blk.push_back(nVir); // add the end point for convenience in slicing
    TAmanager.addRangeType(vLabel, TA::TiledRange1(v_blk.begin(), v_blk.end()));
    std::vector<size_t> o_blk = LinearRange(nOcc, 0, blksize);
    o_blk.push_back(nOcc); // add the end point for convenience in slicing
    TAmanager.addRangeType(oLabel, TA::TiledRange1(o_blk.begin(), o_blk.end()));

    eomSettings.nO = nO;

    std::map<std::string,TArray> ao2mo;
    std::vector<std::string> ao2moTypes{"ao","av"};
    if (nFZC and rebuildFock) ao2moTypes.push_back("ad");
    for(const auto& ao2moType : ao2moTypes){

      std::vector<size_t> offset(2, 0);
      offset[0] = 0;
      switch ( ao2moType[1] ) {
        case 'd': offset[1] = 0; break;
        case 'o': offset[1] = nFZC; break;
        case 'v': offset[1] = nO; break;
      }

      TArray tmp = TAmanager.malloc_fresh<double>(ao2moType);

      tmp.init_elements([mo, offset, nAO](const typename TArray::index &i){
        return mo[i[0] + offset[0] + (i[1] + offset[1]) * nAO];
      });

      ao2mo[ao2moType] = tmp;
    }
#ifdef DEBUG_CCSD
    prettyPrintSmart(std::cout, "MO", mo, nMO, nMO, nMO);
    std::cout << "ao2mo[ao]:" << ao2mo["ao"] << std::endl;
    std::cout << "ao2mo[av]:" << ao2mo["av"] << std::endl;
#endif

    // Create MO TPI
    TArray aoTPIta = taERI.template generateAOERI<double>();
#ifdef DEBUG_CCSD
    std::cout << "aoTPIta:" << std::endl << aoTPIta << std::endl;
#endif

    if (nFZC and rebuildFock) {
      // dddd
      moInts["dddd"] = TAmanager.malloc<double>("dddd");
      moInts["dddd"]("p,r,q,s")  = aoTPIta("m,n,l,g") * ao2mo["ad"]("m,p") * ao2mo["ad"]("n,q") * ao2mo["ad"]("l,r") * ao2mo["ad"]("g,s");

      // dodo
      moInts["dodo"] = TAmanager.malloc<double>("dodo");
      moInts["dodo"]("p,r,q,s")  = aoTPIta("m,n,l,g") * ao2mo["ad"]("m,p") * ao2mo["ad"]("n,q") * ao2mo["ao"]("l,r") * ao2mo["ao"]("g,s");

      // dood
      moInts["dood"] = TAmanager.malloc<double>("dood");
      moInts["dood"]("p,r,q,s")  = aoTPIta("m,n,l,g") * ao2mo["ad"]("m,p") * ao2mo["ao"]("n,q") * ao2mo["ao"]("l,r") * ao2mo["ad"]("g,s");

      // vdod
      moInts["vdod"] = TAmanager.malloc<double>("vdod");
      moInts["vdod"]("p,r,q,s")  = aoTPIta("m,n,l,g") * ao2mo["av"]("m,p") * ao2mo["ao"]("n,q") * ao2mo["ad"]("l,r") * ao2mo["ad"]("g,s");

      // vddo
      moInts["vddo"] = TAmanager.malloc<double>("vddo");
      moInts["vddo"]("p,r,q,s")  = aoTPIta("m,n,l,g") * ao2mo["av"]("m,p") * ao2mo["ad"]("n,q") * ao2mo["ad"]("l,r") * ao2mo["ao"]("g,s");

      // vdvd
      moInts["vdvd"] = TAmanager.malloc<double>("vdvd");
      moInts["vdvd"]("p,r,q,s")  = aoTPIta("m,n,l,g") * ao2mo["av"]("m,p") * ao2mo["av"]("n,q") * ao2mo["ad"]("l,r") * ao2mo["ad"]("g,s");

      // vddv
      moInts["vddv"] = TAmanager.malloc<double>("vddv");
      moInts["vddv"]("p,r,q,s")  = aoTPIta("m,n,l,g") * ao2mo["av"]("m,p") * ao2mo["ad"]("n,q") * ao2mo["ad"]("l,r") * ao2mo["av"]("g,s");
    }

    // oooo
    moInts["oooo"] = TAmanager.malloc<double>("oooo");
    moInts["oooo"]("p,r,q,s")  = aoTPIta("m,n,l,g") * ao2mo["ao"]("m,p") * ao2mo["ao"]("n,q") * ao2mo["ao"]("l,r") * ao2mo["ao"]("g,s");

    // vooo
    moInts["vooo"] = TAmanager.malloc<double>("vooo");
    moInts["vooo"]("p,r,q,s")  = aoTPIta("m,n,l,g") * ao2mo["av"]("m,p") * ao2mo["ao"]("n,q") * ao2mo["ao"]("l,r") * ao2mo["ao"]("g,s");

    // voov
    moInts["voov"] = TAmanager.malloc<double>("voov");
    moInts["voov"]("p,r,q,s")  = aoTPIta("m,n,l,g") * ao2mo["av"]("m,p") * ao2mo["ao"]("n,q") * ao2mo["ao"]("l,r") * ao2mo["av"]("g,s");

    // vovo
    moInts["vovo"] = TAmanager.malloc<double>("vovo");
    moInts["vovo"]("p,r,q,s")  = aoTPIta("m,n,l,g") * ao2mo["av"]("m,p") * ao2mo["av"]("n,q") * ao2mo["ao"]("l,r") * ao2mo["ao"]("g,s");

    // vvoo
    moInts["vvoo"] = TAmanager.malloc<double>("vvoo");
    moInts["vvoo"]("p,r,q,s")  = aoTPIta("m,n,l,g") * ao2mo["av"]("m,p") * ao2mo["ao"]("n,q") * ao2mo["av"]("l,r") * ao2mo["ao"]("g,s");

    // vvvo
    moInts["vvvo"] = TAmanager.malloc<double>("vvvo");
    moInts["vvvo"]("p,r,q,s")  = aoTPIta("m,n,l,g") * ao2mo["av"]("m,p") * ao2mo["av"]("n,q") * ao2mo["av"]("l,r") * ao2mo["ao"]("g,s");

    // vvvv
    moInts["vvvv"] = TAmanager.malloc<double>("vvvv");
    moInts["vvvv"]("p,r,q,s")  = aoTPIta("m,n,l,g") * ao2mo["av"]("m,p") * ao2mo["av"]("n,q") * ao2mo["av"]("l,r") * ao2mo["av"]("g,s");

    TAmanager.free("aaaa", std::move(aoTPIta), true);

    for (auto ta : ao2mo) {
      TAmanager.free(ta.first, std::move(ta.second), true);
    }

    // Create MO Density matrics
    TArray moDen, moDen_dd;
    if (rebuildFock) {
      if (nFZC) {
        moDen_dd = TAmanager.template malloc_fresh<double>("dd");
        moDen_dd.init_elements([](const typename TArray::index &i){
          return i[0] == i[1] ? 1.0 : 0.0;
        });
      }
      moDen = TAmanager.template malloc_fresh<double>("oo");
      moDen.init_elements([](const typename TArray::index &i) {
        return i[0] == i[1] ? 1.0 : 0.0;
      });

      /*
       * Rebuild the Fock matrix from coreH and ERI after ao2mo transformation
       */
      // Create MO H
      cqmatrix::Matrix<double> moCoreH = 0.5 * aoCoreH.S().transform('N', mo, nAO, nAO);
      std::map<std::string, TArray> coreHta;

      // moCoreH.output(std::cout, "moCoreH", true);

      // Build Fock from coreH and TPI to TA blocks
      std::vector<std::string> onePTypes{"oo", "vo", "vv", "ov"};
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

        coreHta[onePType] = TAmanager.template malloc_fresh<double>(onePType);
        coreHta[onePType].init_elements([&moCoreH, offset](const typename TArray::index &i) {
          return moCoreH(i[0] + offset[0], i[1] + offset[1]);
        });

        fockMatrix[onePType] = TAmanager.template malloc<double>(onePType);
      }

#ifdef DEBUG_CCSD
      if (nFZC) std::cout << "Hdd:" << coreHta["dd"] << std::endl;
      std::cout << "Hvv:" << coreHta["vv"] << std::endl;
      std::cout << "Hov:" << coreHta["ov"] << std::endl;
      std::cout << "Hvo:" << coreHta["vo"] << std::endl;
      std::cout << "Hoo:" << coreHta["oo"] << std::endl;
#endif

      fockMatrix["oo"]("p,q") = coreHta["oo"]("p,q") + 2.0 * moInts["oooo"]("p,i,q,j") * moDen("i,j");
      fockMatrix["oo"]("p,q") -= moInts["oooo"]("p,i,j,q") * moDen("i,j");
      fockMatrix["vo"]("p,q") = coreHta["vo"]("p,q") + 2.0 * moInts["vooo"]("p,i,q,j") * moDen("i,j");
      fockMatrix["vo"]("p,q") -= moInts["vooo"]("p,i,j,q") * moDen("i,j");
      fockMatrix["vv"]("p,q") = coreHta["vv"]("p,q") + 2.0 * moInts["vovo"]("p,i,q,j") * moDen("i,j");
      fockMatrix["vv"]("p,q") -= moInts["voov"]("p,i,j,q") * moDen("i,j");
      fockMatrix["ov"]("p,q") = coreHta["ov"]("p,q") + 2.0 * moInts["vooo"]("q,j,p,i") * moDen("i,j");
      fockMatrix["ov"]("p,q") -= moInts["vooo"]("q,j,i,p") * moDen("i,j");

      if (nFZC) {
        fockMatrix["dd"] = TAmanager.template malloc<double>("dd");
        fockMatrix["dd"]("p,q")  = coreHta["dd"]("p,q") + 2.0 * moInts["dddd"]("p,i,q,j") * moDen_dd("i,j");
        fockMatrix["dd"]("p,q") -= moInts["dddd"]("p,i,j,q") * moDen_dd("i,j");
        fockMatrix["dd"]("p,q") += 2.0 * moInts["dodo"]("p,i,q,j") * moDen("i,j");
        fockMatrix["dd"]("p,q") -= moInts["dood"]("p,i,j,q") * moDen("i,j");

        fockMatrix["oo"]("p,q") += 2.0 * moInts["dodo"]("i,p,j,q") * moDen_dd("i,j");
        fockMatrix["oo"]("p,q") -= moInts["dood"]("i,p,q,j") * moDen_dd("i,j");
        fockMatrix["vo"]("p,q") += 2.0 * moInts["vdod"]("p,i,q,j") * moDen_dd("i,j");
        fockMatrix["vo"]("p,q") -= moInts["vddo"]("p,i,j,q") * moDen_dd("i,j");
        fockMatrix["vv"]("p,q") += 2.0 * moInts["vdvd"]("p,i,q,j") * moDen_dd("i,j");
        fockMatrix["vv"]("p,q") -= moInts["vddv"]("p,i,j,q") * moDen_dd("i,j");
        fockMatrix["ov"]("p,q") += 2.0 * moInts["vdod"]("q,j,p,i") * moDen_dd("i,j");
        fockMatrix["ov"]("p,q") -= moInts["vddo"]("q,j,i,p") * moDen_dd("i,j");
      }

      for (auto ta: coreHta)
        TAmanager.free(ta.first, std::move(ta.second), true);

    } else {
      /*
       * We grab the Fock matrix from SingleSlater and slice it afterward
       * to obtain the appropriate spaces. Recomputing E_ref is unnecesssary because it
       * lives in SingleSlater too, but it can be useful to leave as is for checking.
       */

      cqmatrix::Matrix<double> moFock = 0.5 * aoFock.S().transform('N', mo, nAO, nAO);

      std::vector<std::string> onePTypes{"oo", "vo", "vv", "ov"};
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

        fockMatrix[onePType] = TAmanager.template malloc_fresh<double>(onePType);
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
    MultipoleInts<double> moMU = lenElectric.transform('N', mo, nAO, nAO);

    // Build Fock from coreH and TPI to TA blocks
    std::vector<std::string> onePTypes{"oo", "vo", "vv", "ov"};
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
        muMatrix[static_cast<char>('X' + j) + onePType] = TAmanager.template malloc_fresh<double>(onePType);

        muMatrix[static_cast<char>('X' + j) + onePType].init_elements([&moMU, offset, j](const typename TArray::index &i){
          return (*moMU[std::string()+static_cast<char>('X' + j)])(i[0] + offset[0], i[1] + offset[1]);
        });
      }
    }
    TA::get_default_world().gop.fence();

    // std::cout << "Frozen core dipole moments:";
    for (size_t j = 0; j < 3; j++) {
      Mu_fzc[j] = 0.0;
      for (size_t i = 0; i < nFZC; i++)
        Mu_fzc[j] += 2.0 * (*moMU[std::string()+static_cast<char>('X' + j)])(i, i);
      // std::cout << " " << Mu_fzc[j];
    }
    // std::cout << std::endl;
    MPIBCast(Mu_fzc.data(), 3, 0, MPI_COMM_WORLD);


    // Compute diagonal Fock (orbital energies)
    eps.clear();
    eps.resize(nMO-nFZC-nFZV, 0.0);
    std::vector<double> eps_d(nFZC, 0.0);
    if (nFZC) {
      foreach_inplace(fockMatrix["dd"],[&](TA::Tensor<double> &tile) {

        const auto& lobound = tile.range().lobound();
        if (lobound[0] == lobound[1]) {
          const auto& upbound = tile.range().upbound();

          std::size_t x[] = {0, 0};
          for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
            x[1] = x[0];
            eps_d[x[0]] = tile[x];
          }
        }
      });
    }
    foreach_inplace(fockMatrix["oo"],[&](TA::Tensor<double> &tile) {

      const auto& lobound = tile.range().lobound();
      if (lobound[0] == lobound[1]) {
        const auto& upbound = tile.range().upbound();

        std::size_t x[] = {0, 0};
        for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
          x[1] = x[0];
          eps[x[0]] = tile[x];
        }
      }
    });
    foreach_inplace(fockMatrix["vv"], [&](TA::Tensor<double> &tile){

      const auto& lobound = tile.range().lobound();
      if (lobound[0] == lobound[1]) {
        const auto& upbound = tile.range().upbound();

        std::size_t x[] = {0, 0};
        for(x[0] = lobound[0]; x[0] < upbound[0]; ++x[0]) {
          x[1] = x[0];
          eps[nOcc + x[0]] = tile[x];
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
    for (size_t i = 0; i < nAO; i++) {
      std::cout << "Orbital " << i << " : " << eps[i] << std::endl;
    }
#endif
#ifdef DEBUG_CCSD
    std::cout << "Fvv:" << fockMatrix["vv"] << std::endl;
    std::cout << "Fov:" << fockMatrix["ov"] << std::endl;
    std::cout << "Fvo:" << fockMatrix["vo"] << std::endl;
    std::cout << "Foo:" << fockMatrix["oo"] << std::endl;
#endif


    double EF = 0.0;
    double EF_fzc = 0.0;
    double EG = 0.0;
    double EG_fzc = 0.0;

    for (size_t i = 0; i < nFZC; i++)
      EF_fzc += 2.0 * eps_d[i];
    for (size_t i = 0; i < nOcc; i++)
      EF += 2.0 * eps[i];

    if (rebuildFock) {
      EG += 2.0 * (moInts["oooo"]("i,k,j,l") * moDen   ("i,j")).dot(moDen   ("k,l")).get();
      EG -=       (moInts["oooo"]("i,k,l,j") * moDen   ("i,j")).dot(moDen   ("k,l")).get();
      TA::get_default_world().gop.fence();

      if (nFZC) {
        EG_fzc += 2.0 * (moInts["dodo"]("i,k,j,l") * moDen_dd("i,j")).dot(moDen   ("k,l")).get();
        EG_fzc -=       (moInts["dood"]("i,k,l,j") * moDen_dd("i,j")).dot(moDen   ("k,l")).get();
        EG_fzc += 2.0 * (moInts["dodo"]("k,i,l,j") * moDen   ("i,j")).dot(moDen_dd("k,l")).get();
        EG_fzc -=       (moInts["dood"]("k,i,j,l") * moDen   ("i,j")).dot(moDen_dd("k,l")).get();
        EG_fzc += 2.0 * (moInts["dddd"]("i,k,j,l") * moDen_dd("i,j")).dot(moDen_dd("k,l")).get();
        EG_fzc -=       (moInts["dddd"]("i,k,l,j") * moDen_dd("i,j")).dot(moDen_dd("k,l")).get();
        TA::get_default_world().gop.fence();

        // Free and erase all moInts with frozen core ('d') indices
        for (auto it = moInts.begin(); it != moInts.end(); ) {
          if (it->first.find('d') != std::string::npos) {
            TAmanager.free(it->first, std::move(it->second), true);
            it = moInts.erase(it);
          } else {
            ++it;
          }
        }
      }
    } else {
      cqmatrix::Matrix<double> moTwoeH = 0.5 * aoTwoeH.S().transform('N', mo, nAO, nAO);
      //coreHta["oo"]("p,q") = fockMatrix["oo"]("p,q") - moTwoeH_TA(p,q); //coreH with relativistic folded in
      //build HF energy based on 1/2(coreH+F)
      for (size_t i = nFZC; i < nO; i++)
        EG += moTwoeH(i, i);
      for (size_t i = 0; i < nFZC; i++)
        EG_fzc += moTwoeH(i, i);
    }

    if (rebuildFock) TAmanager.free("oo", std::move(moDen), true);
    if (nFZC) {
      TAmanager.free("dd", std::move(fockMatrix["dd"]), true);
      fockMatrix.erase("dd");
    }

    //E_fzc = coreHta + EG_fzc
    E_fzc = EF_fzc - EG_fzc + nucRepEnergy;
    E_ref = EF - EG + E_fzc;

    D_abij = TAmanager.template malloc_fresh<double>("vvoo");
    D_abij.init_elements([this,&denomshift_](const typename TArray::index &i){
      return 1.0/(eps[i[2]] + eps[i[3]] - eps[i[0] + nOcc] - eps[i[1] + nOcc] - denomshift_);
    });

    D_ai = TAmanager.template malloc_fresh<double>("vo");
    D_ai.init_elements([this,&denomshift_](const typename TArray::index &i){
      return 1.0 / (eps[i[1]] - eps[i[0] + nOcc] - denomshift_);
    });

    std::vector<std::string> tmp;
    tmp.push_back(std::string({vLabel}));
    tmp.push_back(std::string({oLabel}));
    tmp.push_back(std::string("OneBody"));
    tmp.push_back(std::string({vLabel,vLabel}));
    tmp.push_back(std::string({oLabel,oLabel}));
    tmp.push_back(std::string("TwoBody"));
    T = std::make_shared<MBExpansion<double>>(tmp, true, MBTensorSymmetry::RCCSD);


    TA::get_default_world().gop.fence();


    std::cout << "    * Initialize MO integrals for coupled cluster took "
              << std::setw(10) << std::right << std::setprecision(6) << std::fixed
              << tock(initIntStart) << " s." << std::endl;
  } // CCIntermediates::initializeIntegrals


  template <typename MatsT>
  void RCCSD<MatsT>::convertTto2C(const MBExpansion<MatsT> &T_RHF, MBExpansion<MatsT> &T_2C) {
    // Convert T_RHF to raw
    size_t size_rhf = T_RHF.length();
    MatsT * t_amp_rhf = CQMemManager::get().malloc<MatsT>(size_rhf);
    TA::get_default_world().gop.fence();
    T_RHF.toRaw(t_amp_rhf, false);
    TA::get_default_world().gop.fence();

    size_t size_2c = T_2C.length();
    MatsT * t_amp_2c = CQMemManager::get().malloc<MatsT>(size_2c);

    if (MPIRank() == 0) {
      const MBTensor<MatsT>& rhfT1 = T_RHF.V_[0], &twoCT1 = T_2C.V_[0];
      const std::vector<size_t>& rhfT1dims = rhfT1.dim();
      size_t nV = rhfT1dims[0], nO = rhfT1dims[1];
      std::fill_n(t_amp_2c, size_2c, MatsT(0.0));

      // One-body
      std::vector<size_t> ai(2, 0);
      for (size_t a = 0; a < nV; a++) {
        for (size_t i = 0; i < nO; i++) {
          size_t idx_rhf = rhfT1.toCompoundIdx({a,i}, true);
          MatsT elem = t_amp_rhf[idx_rhf];
          size_t idx_2c = twoCT1.toCompoundIdx({a,i}, true);
          t_amp_2c[idx_2c] = elem;
          idx_2c = twoCT1.toCompoundIdx({a + nV,i + nO}, true);
          t_amp_2c[idx_2c] = elem;
        }
      }

      // Two-body
      const MBTensor<MatsT>& rhfT2 = T_RHF.V_[1], &twoCT2 = T_2C.V_[1];
      MatsT * t_amp_rhf_2body = t_amp_rhf + rhfT1.size(); // pointer to the start of T_RHF two-body amplitudes
      MatsT * t_amp_2c_2body = t_amp_2c + twoCT1.size(); // pointer to the start of T_2C two-body amplitudes
      size_t idx_2c = 0;
      for (size_t b = 0; b < nV; b++) {
        for (size_t j = 0; j < nO; j++) {
          for (size_t a = 0; a <= b; a++) {
            for (size_t i = 0, iMax = (a < b? nO : j + 1); i < iMax; i++) {
              MatsT abij = t_amp_rhf_2body[rhfT2.toCompoundIdx({a,b,i,j}, true)];
              MatsT abji = t_amp_rhf_2body[rhfT2.toCompoundIdx({a,b,j,i}, true)];

              // <αα||αα> and <ββ||ββ>: need a<b and i<j
              if (a < b and i < j) {
                // <αα||αα>
                idx_2c = twoCT2.toCompoundIdx({a,b,i,j}, true);
                t_amp_2c_2body[idx_2c] = abij - abji;
                // <ββ||ββ>
                idx_2c = twoCT2.toCompoundIdx({a+nV,b+nV,i+nO,j+nO}, true);
                t_amp_2c_2body[idx_2c] = abij - abji;
              }

              // <αβ||αβ>: need every combination of a,b and i,j
              // <βα||βα>: not stored because a>b and i>j
              idx_2c = twoCT2.toCompoundIdx({a,b+nV,i,j+nO}, true);
              t_amp_2c_2body[idx_2c] = abij;
              idx_2c = twoCT2.toCompoundIdx({b,a+nV,j,i+nO}, true);
              t_amp_2c_2body[idx_2c] = abij;

              // <αβ||βα> and <βα||αβ>: not stored because a>b or i>j
            }
          }
        }
      }
    }

    // Broadcast raw amplitudes to all ranks
    MPIBCast(t_amp_2c, size_2c, 0, MPI_COMM_WORLD);

    // Convert raw to T_2C
    TA::get_default_world().gop.fence();
    T_2C.fromRaw(t_amp_2c, false);
    TA::get_default_world().gop.fence();
    T_2C.zeroBody() = T_RHF.zeroBody();
    CQMemManager::get().free(t_amp_rhf, t_amp_2c);
  }

  template <typename MatsT>
  void RCCSD<MatsT>::convertTtoRHF(const MBExpansion<MatsT> &T_2C, MBExpansion<MatsT> &T_RHF) {
    // Convert T_2C to raw
    size_t size_2c = T_2C.length();
    MatsT * t_amp_2c = CQMemManager::get().malloc<MatsT>(size_2c);
    TA::get_default_world().gop.fence();
    T_2C.toRaw(t_amp_2c, false);
    TA::get_default_world().gop.fence();

    size_t size_rhf = T_RHF.length();
    MatsT * t_amp_rhf = CQMemManager::get().malloc<MatsT>(size_rhf);

    if (MPIRank() == 0) {
      const MBTensor<MatsT>& twoCT1 = T_2C.V_[0], &rhfT1 = T_RHF.V_[0];
      const std::vector<size_t>& rhfT1dims = rhfT1.dim();
      size_t nV = rhfT1dims[0], nO = rhfT1dims[1];
      std::fill_n(t_amp_rhf, size_rhf, MatsT(0.0));

      // One-body
      std::vector<size_t> ai(2, 0);
      for (size_t a = 0; a < nV; a++) {
        for (size_t i = 0; i < nO; i++) {
          size_t idx_2c = twoCT1.toCompoundIdx({a,i}, true);
          MatsT elem = t_amp_2c[idx_2c];
          size_t idx_rhf = rhfT1.toCompoundIdx({a,i}, true);
          t_amp_rhf[idx_rhf] = elem;
        }
      }

      // Two-body
      const MBTensor<MatsT>& twoCT2 = T_2C.V_[1], &rhfT2 = T_RHF.V_[1];
      MatsT * t_amp_2c_2body = t_amp_2c + twoCT1.size(); // pointer to the start of T_2C two-body amplitudes
      MatsT * t_amp_rhf_2body = t_amp_rhf + rhfT1.size(); // pointer to the start of T_RHF two-body amplitudes
      for (size_t b = 0; b < nV; b++) {
        for (size_t j = 0; j < nO; j++) {
          for (size_t a = 0; a <= b; a++) {
            for (size_t i = 0, iMax = (a < b? nO : j + 1); i < iMax; i++) {
              // Get data from <αβ||αβ>
              size_t idx_2c = twoCT2.toCompoundIdx({a,b+nV,i,j+nO}, true);
              MatsT abij = t_amp_2c_2body[idx_2c];
              size_t idx_rhf = rhfT2.toCompoundIdx({a,b,i,j}, true);
              t_amp_rhf_2body[idx_rhf] = abij;
              // t_amp_rhf_2body[rhfT2.toCompoundIdx({a,b,i,j}, true)]
              //     = t_amp_2c_2body[twoCT2.toCompoundIdx({a,b+nV,i,j+nO}, true)];
            }
          }
        }
      }
    }

    // Broadcast raw amplitudes to all ranks
    MPIBCast(t_amp_rhf, size_rhf, 0, MPI_COMM_WORLD);

    // Convert raw to T_RHF
    TA::get_default_world().gop.fence();
    T_RHF.fromRaw(t_amp_rhf, false);
    TA::get_default_world().gop.fence();
    T_RHF.zeroBody() = T_2C.zeroBody();
    CQMemManager::get().free(t_amp_2c, t_amp_rhf);
  }


  template <typename MatsT>
  RCCSD<MatsT>::RCCSD(const SafeFile &savFile,
                          CCIntermediates<MatsT> &intermediatesRref,
                          const CoupledClusterSettings &ccSettings):
      CCBase<MatsT>(savFile, intermediatesRref, ccSettings),
      moInts(intermediatesRref.moInts),
      tau_RHF(intermediatesRref.tau),
      tilde_tau_RHF(intermediatesRref.tilde_tau),
      Fae_RHF(intermediatesRref.F_ae),
      Fmi_RHF(intermediatesRref.F_mi),
      Fme_RHF(intermediatesRref.F_me),
      Wmnij_RHF(intermediatesRref.W_mnij),
      Wabef_RHF(intermediatesRref.W_abef),
      Wmbej_RHF_baab(intermediatesRref.W_mbej_baab),
      Wmbej_RHF_baba(intermediatesRref.W_mbej_baba),
      fockMatrix_ta_RHF(intermediatesRref.fockMatrix),
      T_RHF(*intermediatesRref.T),
      T1_RHF(intermediatesRref.T->get_tensor("OneBody")),
      T2_RHF(intermediatesRref.T->get_tensor("TwoBody")),
      Dai_RHF(intermediatesRref.D_ai),
      Dabij_RHF(intermediatesRref.D_abij),
      T_RHF_old(T_RHF) {}


  template <typename MatsT>
  RCCSD<MatsT>::~RCCSD() {

    cleanMemory();

  }

  template <typename MatsT>
  void RCCSD<MatsT>::cleanMemory(){
    TAManager &TAmanager = TAManager::get();
    if(tilde_tau_RHF) TAmanager.free("vvoo", std::move(tilde_tau_RHF), true);
  }

  template <typename MatsT>
  void RCCSD<MatsT>::initIntermediates() {
    TAManager &TAmanager = TAManager::get();
    if (not tau_RHF.is_initialized()){
      tau_RHF = TAmanager.malloc<MatsT>("vvoo");
    }
    if (not tilde_tau_RHF.is_initialized()){
      tilde_tau_RHF = TAmanager.malloc<MatsT>("vvoo");
    }
    if (not Fae_RHF.is_initialized()){
      Fae_RHF = TAmanager.malloc<MatsT>("vv");
    }
    if (not Fmi_RHF.is_initialized()){
      Fmi_RHF = TAmanager.malloc<MatsT>("oo");
    }
    if (not Fme_RHF.is_initialized()){
      Fme_RHF = TAmanager.malloc<MatsT>("ov");
    }
    if (not Wmnij_RHF.is_initialized()){
      Wmnij_RHF = TAmanager.malloc<MatsT>("oooo");
    }
    if (not Wabef_RHF.is_initialized()){
      Wabef_RHF = TAmanager.malloc<MatsT>("vvvv");
    }
    if (not Wmbej_RHF_baab.is_initialized()){
      Wmbej_RHF_baab = TAmanager.malloc<MatsT>("ovvo");
    }
    if (not Wmbej_RHF_baba.is_initialized()){
      Wmbej_RHF_baba = TAmanager.malloc<MatsT>("ovvo");
    }

  }

  template <typename MatsT>
  void RCCSD<MatsT>::initAmplitudes() {
    if(this->ccSettings_.restart){
      size_t size = this->T_RHF.length();
      MatsT * t_amp = CQMemManager::get().malloc<MatsT>(size);
      TA::get_default_world().gop.fence();
      if (MPIRank() == 0) this->savFile_.readData("/CC/T_AMPLITUDE", t_amp);
      if (MPIRank() == 0) this->savFile_.readData("/CC/REFERENCE_ENERGY",   &this->intermediates_.E_ref);
      if (MPIRank() == 0) this->savFile_.readData("/CC/CORRELATION_ENERGY", &this->CorrE);
      MPIBCast(t_amp, size, 0, MPI_COMM_WORLD);
      MPIBCast(&this->intermediates_.E_ref, 1, 0, MPI_COMM_WORLD);
      MPIBCast(&this->CorrE               , 1, 0, MPI_COMM_WORLD);
      TA::get_default_world().gop.fence();
      this->T_RHF.fromRaw(t_amp, false);
      TA::get_default_world().gop.fence();
      CQMemManager::get().free(t_amp);
    } else {
      this->T_RHF.scale(0.0);
    }

  }

  template <typename MatsT>
  void RCCSD<MatsT>::doDIIS(MBExpansion<MatsT> &T_old, std::shared_ptr<DIISTA<MatsT>> diis ){

    // give solution vector to diis
    diis->WriteVector(this->T_RHF);

    // Compute difference of old amplitudes from new amplitudes, write difference into old amplitudes
    T_old.scale(-1.0);
    T_old.axpy(1.0, this->T_RHF);

    //set error vector in DIIS
    diis->WriteErrorVector(T_old);

    // extrapolate new amplitudes from previous amplitudes and their errors. Overwrites solution vector
    diis->Extrapolate(this->T_RHF);

  }

  template <typename MatsT>
  void RCCSD<MatsT>::doDIIS(MBExpansion<MatsT> &T_old, MBExpansion<MatsT> &T_new, std::shared_ptr<DIISTA<MatsT>> diis ){

    // give solution vector to diis
    diis->WriteVector(T_new);

    // Compute difference of old amplitudes from new amplitudes, write difference into old amplitudes
    T_old.scale(-1.0);
    T_old.axpy(1.0, T_new);

    //set error vector in DIIS
    diis->WriteErrorVector(T_old);

    // extrapolate new amplitudes from previous amplitudes and their errors. Overwrites solution vector
    diis->Extrapolate(T_new);

  }

  template <typename MatsT>
  void RCCSD<MatsT>::run() {
    runConventional();
  }

  template <typename MatsT>
  void RCCSD<MatsT>::runConventional(){

    auto cc_start = tick();
    TAManager &TAmanager = TAManager::get();

    std::shared_ptr<DIISTA<MatsT> > diis_RHF = nullptr;
    if(this->ccSettings_.useDIIS){
      diis_RHF = std::make_shared<DIISTA<MatsT>>(this->ccSettings_.nDIIS);
    }

    this->CorrE = 0.0; // make public and not initialized in coupledcluster.hpp

    initIntermediates();
    initAmplitudes();

    if (this->ccSettings_.skipCC && this->ccSettings_.restart) return;

    MBExpansion<MatsT> T_RHF_old(this->T_RHF);

    std::cout << std::setw(18) << std::left <<  "  CC Iterations";
    std::cout << std::setw(34) << std::left << "Corr. Energy (Eh)";
    std::cout << std::setw(19) << std::right << "\u0394Ec (Eh)";
    std::cout << std::setw(19) << std::right << "|\u0394T|";
    std::cout << std::endl;
    std::cout << std::setw(18) << std::left <<  "  -------------";
    std::cout << std::setw(34) << std::left << "-----------------";
    std::cout << std::setw(18) << std::right << "--------";
    std::cout << std::setw(18) << std::right << "----";
    std::cout << std::endl << std::endl;

    for (auto iter = 0; iter < this->ccSettings_.maxiter; iter++){

      T_RHF_old = this->T_RHF;

      buildIntermediates();
#ifdef DEBUG_CCSD
      std::cout << "tau_RHF:" << tau_RHF << std::endl;
      std::cout << "tilde_tau_RHF:" << tilde_tau_RHF << std::endl;
      std::cout << "Fae_RHF:" << Fae_RHF << std::endl;
      std::cout << "Fmi_RHF:" << Fmi_RHF << std::endl;
      std::cout << "Fme_RHF:" << Fme_RHF << std::endl;
      std::cout << "Wmnij_RHF:" << Wmnij_RHF << std::endl;
      std::cout << "Wabef_RHF:" << Wabef_RHF << std::endl;
      std::cout << "Wmbej_RHF:" << Wmbej_RHF << std::endl;
#endif

      updateT1(T_RHF_old.get_tensor("OneBody"), T_RHF_old.get_tensor("TwoBody"));
      updateT2(T_RHF_old.get_tensor("OneBody"), T_RHF_old.get_tensor("TwoBody"));
#ifdef DEBUG_CCSD
      std::cout << "T1_RHF:" << this->T1_RHF << std::endl;
      std::cout << "T2_RHF:" << this->T2_RHF << std::endl;
#endif

      if (this->ccSettings_.useDIIS){
      // diis
        this->doDIIS(T_RHF_old, diis_RHF);
      } else {
        T_RHF_old.axpy(-1, this->T_RHF);
      }

      MatsT Eold = this->CorrE;
      this->getCorrEnergy();
      this->intermediates_.E_cc = this->intermediates_.E_ref + std::real(this->CorrE);
      double dE = std::abs(this->CorrE - Eold);

      double dT = T_RHF_old.norm();

      std::cout << std::setprecision(12) << std::fixed;
      std::cout << "  Iteration "  << std::setw(6) << std::left << iter;
      std::cout << std::setw(34) << std::left << std::fixed << this->CorrE;
      std::cout << std::setw(18) << std::right << std::fixed << dE;
      std::cout << std::setw(18) << std::right << std::fixed << std::abs(dT);
      std::cout << std::endl;

      if (dE < this->ccSettings_.eConv and dT < this->ccSettings_.tConv) {

        std::cout << std::endl << "  CC Completed: Corr. E is "<< std::setw(18) << std::right
                  << std::setprecision(12) << this->CorrE << " Eh" << std::endl;
        std::cout << std::endl << "  CC Completed: Total E is "<< std::setw(18) << std::right
                  << std::setprecision(12) << this->intermediates_.E_ref + this->CorrE << " Eh" << std::endl;
        std::cout << std::endl << "  CC Completed: Iteration total time "<< std::setw(10) << std::right
                  << std::setprecision(6) << tock(cc_start) << " s" << std::endl;


        if(this->ccSettings_.save){
          size_t size = this->T_RHF.length();
          MatsT * t_amp = CQMemManager::get().malloc<MatsT>(size);
          TA::get_default_world().gop.fence();
          this->T_RHF.toRaw(t_amp, false);
          TA::get_default_world().gop.fence();
          if(this->savFile_.exists()){
            this->savFile_.safeWriteData("/CC/T_AMPLITUDE", t_amp, {size});
          }
          CQMemManager::get().free(t_amp);
        }

        if (this->savFile_.exists()) {
          this->savFile_.safeWriteData("/CC/REFERENCE_ENERGY",&this->intermediates_.E_ref, {1});
          this->savFile_.safeWriteData("/CC/CORRELATION_ENERGY",&this->CorrE, {1});
        }

        std::cout << BannerEnd << std::endl;

        break;
      }

      if(iter == this->ccSettings_.maxiter - 1){
        CErr(std::string("CC iterations didn't converge in ") + std::to_string(this->ccSettings_.maxiter) + " steps." );
      }
    }

  }

  template <typename MatsT>
  void RCCSD<MatsT>::getCorrEnergy() {
    // Piecuch, et al. Comput Phys. Commun. 149, 71 (2002) Eq.(50)
    MatsT CorrEOneBody_RHF = 2.0 * fockMatrix_ta_RHF["ov"]("i,a").dot(T1_RHF("a,i")).get();
    MatsT CorrETwoBodyT2_RHF = 2.0 * conj(moInts["vvoo"]("c,d,k,l")).dot(T2_RHF("c,d,k,l")).get();
    CorrETwoBodyT2_RHF -= conj(moInts["vvoo"]("d,c,k,l")).dot(T2_RHF("c,d,k,l")).get();
    MatsT CorrETwoBodyT1_RHF = 2.0 * conj(moInts["vvoo"]("c,d,k,l")).dot(T1_RHF("c,k") * T1_RHF("d,l")).get();
    CorrETwoBodyT1_RHF -= conj(moInts["vvoo"]("d,c,k,l")).dot(T1_RHF("c,k") * T1_RHF("d,l")).get();
    this->CorrE = CorrEOneBody_RHF + CorrETwoBodyT2_RHF + CorrETwoBodyT1_RHF;
  }

  /**
   * Build Eq. III(d.1) and Eq. III(d.2)
   */
  template <typename MatsT>
  void RCCSD<MatsT>::build_tau_and_tilde_tau() {
    tau_RHF("a,b,i,j") = T1_RHF("a,i") * T1_RHF("b,j");
    tilde_tau_RHF("a,b,i,j") = 0.5 * tau_RHF("a,b,i,j") + T2_RHF("a,b,i,j");
    tau_RHF("a,b,i,j") += T2_RHF("a,b,i,j");
  }

  /**
   * Build Eq. III(a.1)
   */
  template <typename MatsT>
  void RCCSD<MatsT>::build_tilde_Fae() {
    Fae_RHF("a,e") = this->fockMatrix_ta_RHF["vv"]("a,e");
    Fae_RHF("a,e") -= 0.5 * this->fockMatrix_ta_RHF["ov"]("m,e") * this->T1_RHF("a,m");

    Fae_RHF("a,e") += 2.0 * this->T1_RHF("f,m") * conj(moInts["vvvo"]("e,f,a,m"));
    Fae_RHF("a,e") -= this->T1_RHF("f,m") * conj(moInts["vvvo"]("f,e,a,m"));

    Fae_RHF("a,e") -= 2.0 * tilde_tau_RHF("a,f,m,n") * conj(moInts["vvoo"]("e,f,m,n"));
    Fae_RHF("a,e") += tilde_tau_RHF("a,f,m,n") * conj(moInts["vvoo"]("e,f,n,m"));
  }

  /**
   * Build Eq. III(a.2)
   */
  template <typename MatsT>
  void RCCSD<MatsT>::build_tilde_Fmi() {
    Fmi_RHF("m,i") = this->fockMatrix_ta_RHF["oo"]("m,i");
    Fmi_RHF("m,i") += 0.5 * this->fockMatrix_ta_RHF["ov"]("m,e") * this->T1_RHF("e,i");

    Fmi_RHF("m,i") += 2.0 * this->T1_RHF("e,n") * conj(moInts["vooo"]("e,i,n,m"));
    Fmi_RHF("m,i") -= this->T1_RHF("e,n") * conj(moInts["vooo"]("e,i,m,n"));

    Fmi_RHF("m,i") += 2.0 * tilde_tau_RHF("e,f,i,n") * conj(moInts["vvoo"]("e,f,m,n"));
    Fmi_RHF("m,i") -= tilde_tau_RHF("e,f,i,n") * conj(moInts["vvoo"]("f,e,m,n"));
  }

  /**
   * Build Eq. III(a.3)
   */
  template <typename MatsT>
  void RCCSD<MatsT>::build_tilde_Fme() {
    Fme_RHF("m,e") = this->fockMatrix_ta_RHF["ov"]("m,e");
    Fme_RHF("m,e") += 2.0 * this->T1_RHF("f,n") * conj(moInts["vvoo"]("e,f,m,n"));
    Fme_RHF("m,e") -= this->T1_RHF("f,n") * conj(moInts["vvoo"]("e,f,n,m"));
  }

  /**
   * Build Eq. III(a.4)
   */
  template <typename MatsT>
  void RCCSD<MatsT>::build_tilde_Wmnij() {
    Wmnij_RHF("m,n,i,j") = this->T1_RHF("e,i") * conj(moInts["vooo"]("e,j,m,n"));
    Wmnij_RHF("m,n,i,j") += this->T1_RHF("e,j") * conj(moInts["vooo"]("e,i,n,m"));
    Wmnij_RHF("m,n,i,j") += moInts["oooo"]("m,n,i,j");
    Wmnij_RHF("m,n,i,j") += 0.5 * tau_RHF("e,f,i,j") * conj(moInts["vvoo"]("e,f,m,n"));
  }

  /**
   * Build Eq. III(a.5)
   */
  template <typename MatsT>
  void RCCSD<MatsT>::build_tilde_Wabef() {
    Wabef_RHF("a,b,e,f") = -this->T1_RHF("b,m") * conj(moInts["vvvo"]("e,f,a,m"));
    Wabef_RHF("a,b,e,f") -= this->T1_RHF("a,m") * conj(moInts["vvvo"]("f,e,b,m"));
    Wabef_RHF("a,b,e,f") += moInts["vvvv"]("a,b,e,f");
    Wabef_RHF("a,b,e,f") += 0.5 * tau_RHF("a,b,m,n") * conj(moInts["vvoo"]("e,f,m,n"));
  }

  /**
   * Build Eq. III(a.6)
   */
  template <typename MatsT>
  void RCCSD<MatsT>::build_tilde_Wmbej() {
    Wmbej_RHF_baab("m,b,e,j") = -moInts["vovo"]("b,m,e,j");
    Wmbej_RHF_baba("m,b,e,j") = moInts["voov"]("b,m,j,e");

    Wmbej_RHF_baab("m,b,e,j") -= this->T1_RHF("f,j") * conj(moInts["vvvo"]("e,f,b,m"));
    Wmbej_RHF_baba("m,b,e,j") += this->T1_RHF("f,j") * conj(moInts["vvvo"]("f,e,b,m"));
    Wmbej_RHF_baab("m,b,e,j") += this->T1_RHF("b,n") * conj(moInts["vooo"]("e,j,n,m"));
    Wmbej_RHF_baba("m,b,e,j") -= this->T1_RHF("b,n") * conj(moInts["vooo"]("e,j,m,n"));

    TArray tmpRHF = TAManager::get().malloc<MatsT>("vvoo");
    tmpRHF("f,b,j,n") = 0.5 * this->T2_RHF("f,b,j,n");
    tmpRHF("f,b,j,n") += this->T1_RHF("f,j") * this->T1_RHF("b,n");
    Wmbej_RHF_baab("m,b,e,j") += tmpRHF("f,b,j,n") * conj(moInts["vvoo"]("f,e,m,n"));

    tmpRHF("f,b,j,n") -= 0.5 * this->T2_RHF("b,f,j,n");
    Wmbej_RHF_baba("m,b,e,j") -= tmpRHF("f,b,j,n") * conj(moInts["vvoo"]("e,f,m,n"));
    Wmbej_RHF_baba("m,b,e,j") += 0.5 * this->T2_RHF("b,f,j,n") * conj(moInts["vvoo"]("e,f,m,n"));
    Wmbej_RHF_baba("m,b,e,j") -= 0.5 * this->T2_RHF("b,f,j,n") * conj(moInts["vvoo"]("f,e,m,n"));
    TAManager::get().free("vvoo", std::move(tmpRHF));
  }

  /**
   * Build Eq. III(a.1-6)
   */
  template <typename MatsT>
  void RCCSD<MatsT>::buildIntermediates() {
    build_tau_and_tilde_tau();
    build_tilde_Fae();
    build_tilde_Fmi();
    build_tilde_Fme();
    build_tilde_Wmnij();
    build_tilde_Wabef();
    build_tilde_Wmbej();
  }

  /**
   * Build Eq. I(a)
   */
  template <typename MatsT>
  void RCCSD<MatsT>::updateT1(const TArray T1_old, const TArray T2_old) {

    this->T1_RHF("a,i") = this->fockMatrix_ta_RHF["vo"]("a,i");
    this->T1_RHF("a,i") += Fae_RHF("a,e") * T1_old("e,i");
    this->T1_RHF("a,i") -= T1_old("a,m") * Fmi_RHF("m,i");
    this->T1_RHF("a,i") += 2.0 * Fme_RHF("m,e") * T2_old("a,e,i,m");
    this->T1_RHF("a,i") -= Fme_RHF("m,e") * T2_old("e,a,i,m");
    this->T1_RHF("a,i") += 2.0 * T1_old("e,m") * moInts["voov"]("a,m,i,e");
    this->T1_RHF("a,i") -= T1_old("e,m") * moInts["vovo"]("a,m,e,i");
    this->T1_RHF("a,i") += 2.0 * T2_old("e,f,i,m") * conj(moInts["vvvo"]("e,f,a,m"));
    this->T1_RHF("a,i") -= T2_old("e,f,i,m") * conj(moInts["vvvo"]("f,e,a,m"));
    this->T1_RHF("a,i") -= 2.0 * T2_old("a,e,m,n") * conj(moInts["vooo"]("e,i,n,m"));
    this->T1_RHF("a,i") += T2_old("a,e,m,n") * conj(moInts["vooo"]("e,i,m,n"));

    this->T1_RHF("a,i") = T1_old("a,i") + this->T1_RHF("a,i") * this->Dai_RHF("a,i");
  }

  /**
   * Build Eq. I(b)
   */
  template <typename MatsT>
  void RCCSD<MatsT>::updateT2(const TArray T1_old, const TArray T2_old) {

    TAManager &TAmanager = TAManager::get();

    this->T2_RHF("a,b,i,j") = moInts["vvoo"]("a,b,i,j");

    TArray TMPbe_RHF = TAmanager.malloc<MatsT>("vv");
    TMPbe_RHF("b,e") = Fae_RHF("b,e");
    TMPbe_RHF("b,e") -= 0.5 * T1_old("b,m") * Fme_RHF("m,e");
    this->T2_RHF("a,b,i,j") += T2_old("a,e,i,j") * TMPbe_RHF("b,e");
    this->T2_RHF("a,b,i,j") += T2_old("b,e,j,i") * TMPbe_RHF("a,e");
    TAmanager.free("vv", std::move(TMPbe_RHF));

    TArray TMPmj_RHF = TAmanager.malloc<MatsT>("oo");
    TMPmj_RHF("m,j") = Fmi_RHF("m,j");
    TMPmj_RHF("m,j") += 0.5 * T1_old("e,j") * Fme_RHF("m,e");
    this->T2_RHF("a,b,i,j") -= T2_old("a,b,i,m") * TMPmj_RHF("m,j");
    this->T2_RHF("a,b,i,j") -= T2_old("b,a,j,m") * TMPmj_RHF("m,i");
    TAmanager.free("oo", std::move(TMPmj_RHF));

    this->T2_RHF("a,b,i,j") += tau_RHF("a,b,m,n") * Wmnij_RHF("m,n,i,j");
    this->T2_RHF("a,b,i,j") += tau_RHF("e,f,i,j") * Wabef_RHF("a,b,e,f");

    TArray Pabij_RHF = TAmanager.malloc<MatsT>("vvoo");
    Pabij_RHF("a,b,i,j") = - moInts["voov"]("b,m,j,e") * T1_old("e,i") * T1_old("a,m");
    Pabij_RHF("a,b,i,j") -= moInts["vovo"]("a,m,e,j") * T1_old("e,i") * T1_old("b,m");
    Pabij_RHF("a,b,i,j") += 2.0 * T2_old("a,e,i,m") * Wmbej_RHF_baba("m,b,e,j");
    Pabij_RHF("a,b,i,j") -= T2_old("a,e,m,i") * Wmbej_RHF_baba("m,b,e,j");
    Pabij_RHF("a,b,i,j") += T2_old("a,e,i,m") * Wmbej_RHF_baab("m,b,e,j");
    Pabij_RHF("a,b,i,j") += T2_old("b,e,m,i") * Wmbej_RHF_baab("m,a,e,j");
    this->T2_RHF("a,b,i,j") += Pabij_RHF("a,b,i,j");
    this->T2_RHF("a,b,i,j") += Pabij_RHF("b,a,j,i");
    TAmanager.free("vvoo", std::move(Pabij_RHF));

    this->T2_RHF("a,b,i,j") += T1_old("e,i") * moInts["vvvo"]("a,b,e,j");
    this->T2_RHF("a,b,i,j") += T1_old("e,j") * moInts["vvvo"]("b,a,e,i");

    this->T2_RHF("a,b,i,j") -= T1_old("a,m") * moInts["vooo"]("b,m,j,i");
    this->T2_RHF("a,b,i,j") -= T1_old("b,m") * moInts["vooo"]("a,m,i,j");

    this->T2_RHF("a,b,i,j") = T2_old("a,b,i,j") + this->T2_RHF("a,b,i,j") * this->Dabij_RHF("a,b,i,j");
  }

  template <typename MatsT>
  size_t RCCSD<MatsT>::estimate_mem_peak() const {
    // will need further checking as the TA objects are dynamically allocated and freed at runtime
    TAManager &TAmanager = TAManager::get();

    size_t nDIIS = this->ccSettings_.useDIIS ? this->ccSettings_.nDIIS : 0;
    size_t count = 0;
    //                                             1         2       3       4         5         6      7
    count += 7 * TAmanager.elem_per_TA("oo");   // muX_oo,   muY_oo, muZ_oo, coreH_oo, fock_oo,  Fmi_,  TMPmj
    count += 6 * TAmanager.elem_per_TA("ov");   // muX_ov,   muY_ov, muZ_ov, coreH_oo, fock_ov,  Fme_
    count += 7 * TAmanager.elem_per_TA("vo");   // muX_vo,   muY_vo, muZ_vo, coreH_oo, fock_vo,  Dai,   T1
    count += 7 * TAmanager.elem_per_TA("vv");   // muX_vv,   muY_vv, muZ_vv, coreH_oo, fock_vv,  Fae_,  TMPbe
    count += 2 * TAmanager.elem_per_TA("oooo"); // ERI_oooo, Wmnij_
    count += 1 * TAmanager.elem_per_TA("vooo"); // ERI_vooo
    count += 2 * TAmanager.elem_per_TA("ovvo"); // Wmbej_,   tmp
    count += 3 * TAmanager.elem_per_TA("voov"); // ERI_voov, Wmbej1, Wmbej2
    count += 1 * TAmanager.elem_per_TA("vovo"); // ERI_vovo
    count += 7 * TAmanager.elem_per_TA("vvoo"); // ERI_vvoo, Dabij,  T2,     tau_,     tilde_tau_, tmp,   Pabij
    count += 1 * TAmanager.elem_per_TA("vvvo"); // ERI_vvvo,
    count += 2 * TAmanager.elem_per_TA("vvvv"); // ERI_vvvv, W_abef

    if (nDIIS) {
      count += (nDIIS + 1) * 2 * TAmanager.elem_per_TA("vo");   // T1 DIIS copy?
      count += (nDIIS + 1) * 2 * TAmanager.elem_per_TA("vvoo"); // T2 DIIS copy?
    }

    return count * sizeof(MatsT);
  }

}; // namespace ChronusQ