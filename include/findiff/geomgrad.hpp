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

#include <findiff.hpp>
#include <singleslater.hpp>
#include <grid/integrator.hpp>


namespace ChronusQ {

  class NumGradient {
 
    CQInputFile& input_;
    std::shared_ptr<SingleSlaterBase>& ref_;
    std::shared_ptr<SingleSlaterBase> curr_;
    EMPerturbation emPert;
    std::shared_ptr<IntegralsBase> aoints;
    Molecule mol;
    std::shared_ptr<BasisSet> basis;
    std::shared_ptr<BasisSet> dfbasis;
    std::shared_ptr<BasisSet> prot_basis;
    
    size_t nAtoms;

    template <typename MatsT, typename IntsT>
    bool goodType() {
      return dynamic_cast<SingleSlater<MatsT,IntsT>*>(ref_.get()) != nullptr;
    };

    // Movement methods (see src/findiff/geomgrad.cxx for details)
    template <typename MatsT, typename IntsT>
    void moveMol(size_t iAtm, size_t iXYZ, double diff);

    // Properties to get (bound to valFunc for the finite differencer)
    //   (see src/findiff/geomgrad.cxx for details)
    std::vector<double> getEnergies();
    std::vector<double> getGridVals(bool doGrad, size_t iA, size_t iXYZ);

    template<typename IntsT>
    std::vector<IntsT> getOneEInts();

    template<typename IntsT>
    std::vector<IntsT> getERI();

    public:

    NumGradient() = delete;

    NumGradient(CQInputFile& in, std::shared_ptr<SingleSlaterBase>& ss, std::shared_ptr<BasisSet> b) :
      input_(in), ref_(ss), basis(b) {
      
      if ( goodType<dcomplex,dcomplex>() )
        nAtoms = std::dynamic_pointer_cast<SingleSlater<dcomplex,dcomplex>>(ref_)->molecule().nAtoms;
      else if ( goodType<dcomplex,double>() )
        nAtoms = std::dynamic_pointer_cast<SingleSlater<dcomplex,double>>(ref_)->molecule().nAtoms;
      else
        nAtoms = std::dynamic_pointer_cast<SingleSlater<double,double>>(ref_)->molecule().nAtoms;

    };

    NumGradient(const NumGradient&) = default;
    NumGradient(NumGradient&&) = default;

    void doGrad(size_t acc = 2) {

      FiniteDifferencer<double,double> diff;

      diff.setDerOrder(1);
      diff.setAccOrder(acc);

      diff.setStepSize(1e-2);

      
      diff.setValFunc(std::bind(&NumGradient::getEnergies, this));

      bool cInt = goodType<dcomplex,dcomplex>();
      bool cMat = cInt or goodType<dcomplex,double>();

      for ( auto iAtm = 0; iAtm < nAtoms; iAtm++ ) {
        for ( auto iXYZ = 0; iXYZ < 3; iXYZ++ ) {

          if(cInt)
            diff.addChangeFunc(
              std::bind(&NumGradient::moveMol<dcomplex,dcomplex>, this,
                        iAtm, iXYZ, std::placeholders::_1));

          else if(cMat)
            diff.addChangeFunc(
              std::bind(&NumGradient::moveMol<dcomplex,double>, this,
                        iAtm, iXYZ, std::placeholders::_1));

          else
            diff.addChangeFunc(
              std::bind(&NumGradient::moveMol<double,double>, this,
                        iAtm, iXYZ, std::placeholders::_1));

        }
      }

      diff.doDifference();

      auto res = diff.getResults();

      printResults("Nuclear Gradient (Eh / Bohr)", res, 0, 1e-12);
      printResults("OneE Gradient (Eh / Bohr)", res, 1);
      printResults("TwoE Gradient (Eh / Bohr)", res, 2);
      printResults("NucRep Gradient (Eh / Bohr)", res, 3);
      if (res.size() > 4)
        printResults("XC Gradient (Eh / Bohr)", res, 4);

    };

    void gridGrad(size_t acc = 2) {
      
      FiniteDifferencer<double,double> diff;

      diff.setDerOrder(1);
      diff.setAccOrder(acc);

      diff.setStepSize(1e-4);

      
      diff.setValFunc(std::bind(&NumGradient::getGridVals, this, false, 0, 0));

      bool cInt = goodType<dcomplex,dcomplex>();
      bool cMat = cInt or goodType<dcomplex,double>();

      for ( auto iAtm = 0; iAtm < nAtoms; iAtm++ ) {
        for ( auto iXYZ = 0; iXYZ < 3; iXYZ++ ) {

          if(cInt)
            diff.addChangeFunc(
              std::bind(&NumGradient::moveMol<dcomplex,dcomplex>, this,
                        iAtm, iXYZ, std::placeholders::_1));

          else if(cMat)
            diff.addChangeFunc(
              std::bind(&NumGradient::moveMol<dcomplex,double>, this,
                        iAtm, iXYZ, std::placeholders::_1));

          else
            diff.addChangeFunc(
              std::bind(&NumGradient::moveMol<double,double>, this,
                        iAtm, iXYZ, std::placeholders::_1));

        }
      }

      diff.doDifference();

      auto res = diff.getResults();

      moveMol<double,double>(0, 0, 0.);

      for ( auto iAtom = 0; iAtom < nAtoms; iAtom++ ) {
        for ( auto iXYZ = 0; iXYZ < 3; iXYZ++ ) {
          std::cout << "=== Weight Gradients [" << iAtom << "," << iXYZ << "]\n";

          auto anaRes = getGridVals(true, iAtom, iXYZ);

          bool allgood = true;

          for ( auto i = 0; i < anaRes.size(); i++ ) {
            if ( std::max( res[iAtom*3 + iXYZ][i], anaRes[i] ) < 1e-12 )
              continue;
            if ( std::abs( res[iAtom*3 + iXYZ][i] - anaRes[i] ) / anaRes[i] > 1e-8 ) {
              std::cout << "  " << res[iAtom*3 + iXYZ][i];
              std::cout << "  " << anaRes[i];
              std::cout << "  " << (res[iAtom*3 + iXYZ][i] - anaRes[i])/anaRes[i];
              std::cout << "\n";
              allgood = false;
            }
          }

          if (allgood)
            std::cout << "All good\n";
          std::cout << std::endl;
        }
      }

    };

    template <typename IntsT>
    void intGrad(size_t acc = 2) {

      FiniteDifferencer<double,IntsT> diff;

      diff.setDerOrder(1);
      diff.setAccOrder(acc);

      diff.setStepSize(1e-3);

      diff.setValFunc(std::bind(&NumGradient::getOneEInts<IntsT>, this));

      bool cInt = goodType<dcomplex,dcomplex>();
      bool cMat = cInt or goodType<dcomplex,double>();

      for ( auto iAtm = 0; iAtm < nAtoms; iAtm++ ) {
        for ( auto iXYZ = 0; iXYZ < 3; iXYZ++ ) {

          if(cInt)
            diff.addChangeFunc(
              std::bind(&NumGradient::moveMol<dcomplex,dcomplex>, this,
                        iAtm, iXYZ, std::placeholders::_1));

          else if(cMat)
            diff.addChangeFunc(
              std::bind(&NumGradient::moveMol<dcomplex,double>, this,
                        iAtm, iXYZ, std::placeholders::_1));

          else
            diff.addChangeFunc(
              std::bind(&NumGradient::moveMol<double,double>, this,
                        iAtm, iXYZ, std::placeholders::_1));

        }
      }

      diff.doDifference();

      auto res = diff.getResults();

      std::vector<std::string> names = {"Numerical Overlap Gradient",
                                        "Numerical Kinetic Gradient",
                                        "Numerical Potential Gradient"};

      auto NB = basis->nBasis;

      std::cout << " ===  Integral Numerical Gradients  === " << std::endl;
      for ( auto iOp = 0; iOp < 3; iOp++ )
      for ( auto iAtom = 0; iAtom < nAtoms; iAtom++ ) {

        for ( auto iXYZ = 0; iXYZ < 3; iXYZ++ ) {

          std::cout << "Atom: " << std::setw(8) << std::left << iAtom;
          std::cout << "Cart: " << std::setw(8) << std::left << iXYZ;
          std::cout << std::endl;

          prettyPrintSmart(std::cout, names[iOp],
                           res[iAtom*3 + iXYZ].data() + iOp*NB*NB,
                           NB, NB, NB);
        }
      }
      std::cout << " ====================================== " << std::endl;

    };

    template <typename IntsT>
    void eriGrad(size_t acc = 2) {

      FiniteDifferencer<double,IntsT> diff;

      diff.setDerOrder(1);
      diff.setAccOrder(acc);

      diff.setStepSize(1e-3);

      diff.setValFunc(std::bind(&NumGradient::getERI<IntsT>, this));

      bool cInt = goodType<dcomplex,dcomplex>();
      bool cMat = cInt or goodType<dcomplex,double>();

      for ( auto iAtm = 0; iAtm < nAtoms; iAtm++ ) {
        for ( auto iXYZ = 0; iXYZ < 3; iXYZ++ ) {

          if(cInt)
            diff.addChangeFunc(
              std::bind(&NumGradient::moveMol<dcomplex,dcomplex>, this,
                        iAtm, iXYZ, std::placeholders::_1));

          else if(cMat)
            diff.addChangeFunc(
              std::bind(&NumGradient::moveMol<dcomplex,double>, this,
                        iAtm, iXYZ, std::placeholders::_1));

          else
            diff.addChangeFunc(
              std::bind(&NumGradient::moveMol<double,double>, this,
                        iAtm, iXYZ, std::placeholders::_1));

        }
      }

      diff.doDifference();

      auto res = diff.getResults();

      auto NB = basis->nBasis;
      auto MB = prot_basis ? prot_basis->nBasis : basis->nBasis;

      std::cout << " ===  Numerical ERI Gradient  === " << std::endl;
      std::cout << std::setprecision(10);
      for ( auto iAtom = 0; iAtom < nAtoms; iAtom++ ) {

        for ( auto iXYZ = 0; iXYZ < 3; iXYZ++ ) {

          std::cout << "Atom: " << std::setw(8) << std::left << iAtom;
          std::cout << "Cart: " << std::setw(8) << std::left << iXYZ;
          std::cout << '\n' << std::endl;

          for ( auto i = 0, ijkl = 0; i < NB; i++ )
          for ( auto j = 0; j < NB; j++ )
          for ( auto k = 0; k < MB; k++ )
          for ( auto l = 0; l < MB; l++, ijkl++ ) {
            std::cout << "    (" << i << "," << j << "," << k << "," << l << ")";
            std::cout << "    ";
            if (std::abs(res[iAtom*3 +iXYZ][i + j*NB + k*NB*NB + l*NB*NB*MB]) < 1e-14)
              std::cout << 0.;
            else 
              std::cout << res[iAtom*3 +iXYZ][i + j*NB + k*NB*NB + l*NB*NB*MB];
            std::cout << std::endl;
          }
        }
      }
      std::cout << " ====================================== " << std::endl;

    };

    void printResults(std::string name, std::vector<std::vector<double>> res, size_t resIndex, double minPrint=1e-8) {
      std::cout << "\n" << name << "\n\n";

      std::cout << std::setw(18) << std::left << "  Atom";
      std::cout << std::setw(18) << std::right << "X";
      std::cout << std::setw(18) << std::right << "Y";
      std::cout << std::setw(18) << std::right << "Z" << std::endl;
      for ( auto iAtom = 0; iAtom < nAtoms; iAtom++ ) {
        std::cout << "  " << std::setw(16) << std::left << iAtom;

        std::cout << std::setprecision(8) << std::scientific << std::right;
        for ( auto iXYZ = 0; iXYZ < 3; iXYZ++ )
          std::cout << std::setw(18) << (std::abs(res[iAtom*3 + iXYZ][resIndex]) > minPrint ? res[iAtom*3+iXYZ][resIndex] : 0.);

        std::cout << std::endl;
      }

    }

  };

} // namespace ChronusQ

