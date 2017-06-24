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
#ifdef CQ_HAS_TA

#include <tiledarray.h>
#include <particleintegrals/twopints/incore4indextpi.hpp>
#include <particleintegrals/twopints/gtodirecttpi.hpp>
#include <libcint.hpp>
#include <coupledcluster/TAManager.hpp>

namespace ChronusQ {

  template <typename IntsT>
  class TAERI {

    const TwoPInts<IntsT> &aoTPI_;
    TA::TiledRange1 aorange_;

  public:
    TAERI(const TwoPInts<IntsT> &aoTPI, size_t blksize): aoTPI_(aoTPI) {
      aorange_ = generateAOrange(blksize);
    }

    TA::TiledRange1 getAOrange() const { return aorange_; }
    
    static TA::TiledRange1 LinRange(size_t N, size_t blksize) {
      size_t blocks = N % blksize == 0 ? N / blksize : N / blksize + 1;
      std::vector<size_t> blk;
      blk.reserve(blocks + 1);
      
      for (auto i = 0 ; i < blocks; i++)
        blk.push_back(blksize * i);
      blk.push_back(N);
      
      return TA::TiledRange1(blk.begin(), blk.end());
    }

    TA::TiledRange1 generateAOrange(size_t blksize) const {
      try {
        const DirectTPI<IntsT>& directTPI = dynamic_cast<const DirectTPI<IntsT>&>(aoTPI_);
        const BasisSet &basis = directTPI.basisSet().groupGeneralContractionBasis();

        std::vector<size_t> blk(1, 0);

        size_t count = 0, curBlk = 0;
        for (const libint2::Shell &sh : basis.shells) {
          size_t shSize = sh.size();

          if (curBlk + shSize > blksize) {
            if (count > 0)
              blk.push_back(count);
            curBlk = shSize;
          } else {
            curBlk += shSize;
          }
          count += shSize;
        }

        assert(count == basis.nBasis);

        blk.push_back(basis.nBasis);

        return TA::TiledRange1(blk.begin(), blk.end());
      } catch (const std::bad_cast& e) {}

      return LinRange(aoTPI_.nBasis(), blksize);
    }

    template <typename MatsT>
    TA::TArray<MatsT> generateAOERI() const {

      TAManager &TAmanager = TAManager::get();

      // Incore TPI case
      try {
        const InCore4indexTPI<IntsT> &aoIncoreTPI = dynamic_cast<const InCore4indexTPI<IntsT>&>(aoTPI_);

        TA::TArray<MatsT> aoTPIta = TAmanager.template malloc_fresh<MatsT>("aaaa");
        aoTPIta.init_elements([&aoIncoreTPI](const typename TA::TArray<MatsT>::index &i){
          return aoIncoreTPI(i[0], i[1], i[2], i[3]);
        });

        return aoTPIta;
      } catch (const std::bad_cast&) {}

      // Direct TPI case
      try {
        const DirectTPI<IntsT> &directTPI = dynamic_cast<const DirectTPI<IntsT>&>(aoTPI_);
        const BasisSet &basis = directTPI.basisSet().groupGeneralContractionBasis();
        
        if (basis.forceCart)
          CErr("Libcint + cartesian GTO NYI.");
        
        std::set<std::array<double, 3>> shellCenters;
        std::vector<size_t> mapBf2Sh(basis.nBasis, 0);
        for (size_t i = 0, bf = 0; i < basis.shells.size(); i++) {
          std::fill_n(&mapBf2Sh[bf], basis.shells[i].size(), i);
          bf += basis.shells[i].size();
          shellCenters.insert(basis.shells[i].O);
        }

        std::vector<Atom> atoms;
        atoms.reserve(shellCenters.size());
        for (const std::array<double, 3> &center : shellCenters) {
          atoms.emplace_back("HE-4", center);
        }

        Molecule mol(0,1,atoms);
        size_t nAtoms = mol.nAtoms;
        size_t nShells = basis.nShell;

        // ATM_SLOTS = 6; BAS_SLOTS = 8;
        int *atm = CQMemManager::get().malloc<int>(nAtoms * ATM_SLOTS);
        int *bas = CQMemManager::get().malloc<int>(nShells * BAS_SLOTS);
        double *env = CQMemManager::get().malloc<double>(basis.getLibcintEnvLength(mol));

        basis.setLibcintEnv(mol, atm, bas, env);

        TA::TArray<MatsT> aoTPIta = TAmanager.template malloc_fresh<MatsT>("aaaa");
        aoTPIta.init_tiles([&basis, &mapBf2Sh, atm, bas, env, nAtoms, nShells](const TA::Range &range){
          TA::Tensor<MatsT> tile(range, 0.0);
          const auto& lobound = range.lobound();
          const auto& upbound = range.upbound();
          int shls[4];
          std::array<size_t, 4> shellSizes;

          // TODO: implement
          std::size_t x[] = {0,0,0,0};
          for(x[0] = lobound[0]; x[0] < upbound[0]; x[0] += shellSizes[0]) {
            shls[0] = mapBf2Sh[x[0]];
            shellSizes[0] = basis.shells[shls[0]].size();
            for(x[1] = lobound[1]; x[1] < upbound[1]; x[1] += shellSizes[1]) {
              shls[1] = mapBf2Sh[x[1]];
              shellSizes[1] = basis.shells[shls[1]].size();
              for(x[2] = lobound[2]; x[2] < upbound[2]; x[2] += shellSizes[2]) {
                shls[2] = mapBf2Sh[x[2]];
                shellSizes[2] = basis.shells[shls[2]].size();
                for(x[3] = lobound[3]; x[3] < upbound[3]; x[3] += shellSizes[3]) {
                  shls[3] = mapBf2Sh[x[3]];
                  shellSizes[3] = basis.shells[shls[3]].size();

                  std::vector<double> buff(shellSizes[0]*shellSizes[1]*shellSizes[2]*shellSizes[3], 0.0);
                  if(cint2e_sph(buff.data(), shls, atm, nAtoms, bas, nShells, env, nullptr)==0) {
                    continue;
                  }

                  for (size_t s = 0, pqrs = 0; s < shellSizes[3]; s++)
                    for (size_t r = 0; r < shellSizes[2]; r++)
                      for (size_t q = 0; q < shellSizes[1]; q++)
                        for (size_t p = 0; p < shellSizes[0]; p++, pqrs++) {
                          std::size_t xx[] = {x[0]+p,x[1]+q,x[2]+r,x[3]+s};
                          tile[xx] = buff[pqrs];
                        }

                }
              }
            }
          }

          return tile;
        });

        TA::get_default_world().gop.fence();
        CQMemManager::get().free(env, bas, atm);
        return aoTPIta;
      } catch (const std::bad_cast&) {}

      CErr("Two-particle integral type not supported in coupled cluster calculation.");

    }

  }; // class TAERI

}; // namespace ChronusQ
#endif
