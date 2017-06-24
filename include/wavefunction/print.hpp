/* 
 *  This file is part of the Chronus Quantum (ChronusQ) software package
 *  
 *  Copyright (C) 2014-2022 Li Research Group (University of Washington)
 *  
 *  This program is free software; you ca redistribute it and/or modify
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

#include <wavefunction.hpp>

#include <util/matout.hpp>
#include <util/math.hpp>
#include <physcon.hpp>

namespace ChronusQ {

  template <typename MatsT, typename IntsT>
  void WaveFunction<MatsT,IntsT>::printMO(std::ostream &out) {

    size_t NB = nAlphaOrbital() * this->nC;
   
    prettyPrintSmart(std::cout,"MO 1",mo[0].pointer(),NB,NB,NB);
    if( mo.size() > 1 )
      prettyPrintSmart(std::cout,"MO 2",mo[1].pointer(),NB,NB,NB);

  }; // WaveFunction<T,MatsT,IntsT>::printMO

  template <typename MatsT, typename IntsT>
  void WaveFunction<MatsT,IntsT>::printEPS(std::ostream &out) {

    size_t NB = nAlphaOrbital() * this->nC;
   
    prettyPrintSmart(std::cout,"EPS 1",eps1,NB,NB,NB);
    if( eps2 != nullptr )
      prettyPrintSmart(std::cout,"EPS 2",eps2,NB,NB,NB);

  }; // WaveFunction<T,MatsT,IntsT>::printEPS

  /**
   * \brief Function to print detailed MO coefficients
   *
   */
  template <typename T, typename ValManipOp>
  void prettyMOPrint(std::ostream &out, size_t NB, size_t NOrb, double *EPS,
    T* MO, size_t LDM, Molecule &mol, BasisSet &basis,
    ValManipOp op ){

    constexpr size_t maxLPrint = 6 + 1;

    std::array< std::vector<std::string>, maxLPrint > cartLabel, sphLabel;
    std::array< std::string, maxLPrint > angLabel =
    { "S", "P", "D", "F", "G", "H", "I" };

    std::array< std::string, 3 > axes = { "X", "Y", "Z" };

    // S Functions
    cartLabel[0] = { angLabel[0] };

    for(auto i = 0; i < 3; i++) {

      // P Functions
      cartLabel[1].emplace_back( angLabel[1] + axes[i] );

    for(auto j = i; j < 3; j++) {

      // D Functions
      cartLabel[2].emplace_back( angLabel[2] + axes[i] + axes[j] );

    for(auto k = j; k < 3; k++) {

      // F Functions
      cartLabel[3].emplace_back( angLabel[3] + axes[i] + axes[j] + axes[k] );

    for(auto l = k; l < 3; l++) {

      // G Functions
      cartLabel[4].emplace_back( angLabel[4] + axes[i] + axes[j] + axes[k]
          + axes[l] );

    for(auto m = l; m < 3; m++) {

      // H Functions
      cartLabel[5].emplace_back( angLabel[5] + axes[i] + axes[j] + axes[k]
          + axes[l] + axes[m] );

    for(auto n = m; n < 3; n++) {

      // I Functions
      cartLabel[6].emplace_back( angLabel[6] + axes[i] + axes[j] + axes[k]
          + axes[l] + axes[m] + axes[n] );

    }}}}}}

    // Spherical Labels
    auto pm_string = [](int x) -> std::string {

      if( x > 0 )      return "+" + std::to_string(x) ;
      else if (x < 0 ) return "" + std::to_string(x) ;
      else             return " " + std::to_string(x) ;

    };

    sphLabel[0] = cartLabel[0];
    sphLabel[1] = cartLabel[1];
    for(auto i = 2; i < maxLPrint; i++) {

      for(auto ml = -i; ml <= i; ml++)
        sphLabel[i].emplace_back( angLabel[i] + pm_string(ml) );

    }

    // Printing format setting
    size_t list = 4; // #mo per line
    size_t printWidth = 14;
    size_t eValOff = 24;

    out << std::endl << bannerTop << std::endl;

    out << std::scientific << std::left << std::setprecision(5);
    for(size_t p = 0; p < NOrb; p += list) {

      out << std::left << std::endl;
      int end = list;

      if( (p + list) >= NOrb ) end = NOrb - p;

      out << std::left << std::setw(eValOff) << " ";
      out << std::right;
      for(size_t k = p; k < p + end; k++)
        out << std::setw(printWidth) << k + 1;
      out << std::endl;

      out << std::left << std::setw(eValOff) << " EigV --";
      out << std::right;
      for(size_t k = p; k < p + end; k++)
        out << std::setw(printWidth) << EPS[k];
      out << std::endl << std::endl;

      // Function to get the atom info
      auto getAtmSymb = [&](size_t iAtm) -> std::string {

        std::map<std::string,Atom>::const_iterator it =
        std::find_if(atomicReference.begin(),atomicReference.end(),
          [&](const std::pair<std::string,Atom> &st){
            return (st.second.atomicNumber == mol.atoms[iAtm].atomicNumber) and
                   (st.second.massNumber == mol.atoms[iAtm].massNumber);}
           );

        return (it == atomicReference.end() ? "X" : it->first);

      };

      size_t iAtm = 0;
      std::string atmSymb = getAtmSymb(iAtm);
      size_t iShellAtm = 0;
      size_t nShell = basis.shells.size();

      for(size_t iShell = 0; iShell < nShell ; iShell++, iShellAtm++) {

        size_t bfst = basis.mapSh2Bf[iShell];
        size_t sz   = basis.shells[iShell].size();
        size_t L    = basis.shells[iShell].contr[0].l;

        size_t curCen = basis.mapSh2Cen[iShell];
        bool newAtm = false;
        if( curCen != iAtm ) {
          newAtm = true;
          iAtm++;
          atmSymb = getAtmSymb(iAtm);
          iShellAtm = 0;
        }

        bool lastAtm  = iAtm == (mol.atoms.size() - 1);
        bool firstAtm = iShell == 0;

        for(size_t mu = bfst, ml = 0 ; ml < sz; mu++, ml++) {

          // BF number
          out << " " << std::setw(5) << std::left <<  mu + 1; // 6

          // atom info
          if( (newAtm or firstAtm) ) {
            out << std::setw(3) << iAtm;
            out << std::setw(7) << atmSymb;
          } else out << std::setw(10) << " "; // 16

          // #shell and spherical label
          out << std::setw(2) << iShellAtm ;
          out << std::setw(4) << sphLabel[L][ml]; //22

          out << std::setw(2) << " "; // 24

          // MO coefficients
          for(auto q = p; q < p + end; q++) {

            double VAL = op(MO[mu + q*LDM]);
            out << std::right << std::setw(printWidth);

            if(std::abs(VAL) > PRINT_SMALL)    out << VAL;
            else if(std::isnan(std::abs(VAL))) out << "NAN";
            else if(std::isinf(std::abs(VAL))) out << "INF";
            else                               out << 0.;

          }

          out << std::endl;

          bool nextAtmNew = lastAtm ? false : (mu+1) >= basis.mapCen2BfSt[iAtm+1];

          if( nextAtmNew ) out<<std::endl;
        }
      }
      out << std::endl;
    }
    out << bannerEnd << std::endl;

  }; // prettyMOPrint

  template <typename T>
  void prettyMOPrint(std::ostream &out, size_t NB, size_t NOrb, double *EPS,
    T* MO, size_t LDM, Molecule &mol, BasisSet &basis){

    prettyMOPrint(out,NB,NOrb,EPS,MO,LDM,mol,basis,
        [](T x){ return x; });

  }; // prettyMOPrint

  /**
   * \brief Print out the components of different angular momentum of MOs
   *
   */
  template <typename T, typename IntsT>
  void analyzeMOPrint(std::ostream &out, size_t NB, size_t NOrb, IntsT* S, T* MO,
        size_t LDM, Molecule &mol, BasisSet &basis,
        bool groupAtm = false, T* MO2 = nullptr) {

    out << "\nMO components projected to basis functions of different angular momentum"; 

    constexpr size_t maxLPrint = 6 + 1;
    T* SCR  = CQMemManager::get().malloc<T>(NB*NOrb);
    blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
               NB,NOrb,NB,T(1.),S,NB,MO,LDM,T(0.),SCR,NB);
    // If alpha and beta parts are summed together
    T *SCR2;
    if (MO2) {
      SCR2 = CQMemManager::get().malloc<T>(NB*NOrb);
      blas::gemm(blas::Layout::ColMajor,blas::Op::NoTrans,blas::Op::NoTrans,
               NB,NOrb,NB,T(1.),S,NB,MO2,LDM,T(0.),SCR2,NB);
    }

    std::array< std::string, maxLPrint > angLabel =
                { "S", "P", "D", "F", "G", "H", "I" };

    // Printing format setting
    size_t list = 7;
    size_t printWidth = 10;
    size_t eValOff = 10;

    out << std::endl << bannerTop << std::endl;

    auto getAtmSymb = [&](size_t iAtm) -> std::string {

      std::map<std::string,Atom>::const_iterator it =
      std::find_if(atomicReference.begin(),atomicReference.end(),
        [&](const std::pair<std::string,Atom> &st){
          return (st.second.atomicNumber == mol.atoms[iAtm].atomicNumber) and
                 (st.second.massNumber == mol.atoms[iAtm].massNumber);}
         );

      return (it == atomicReference.end() ? "X" : it->first);

    };


    for(size_t p = 0; p < NOrb; p += list) {
      out << std::left << std::endl;
      int end = list;
      if( (p + list) >= NOrb ) end = NOrb - p;
      std::vector<std::vector<T>> angWeight(end,
                std::vector<T>(maxLPrint, T(0.0)));
      out << std::left << std::setw(eValOff) << " ";
      out << std::right;
      for(size_t k = p; k < p + end; k++)
        out << std::setw(printWidth) << k + 1;
      out << std::endl;

      size_t iAtm = 0;
      std::string atmSymb = getAtmSymb(iAtm);

      size_t nShell = basis.shells.size();
      size_t maxL = 0;
      for(size_t iShell = 0; iShell < nShell ; iShell++) {
        size_t bfst = basis.mapSh2Bf[iShell];
        size_t sz   = basis.shells[iShell].size();
        size_t L    = basis.shells[iShell].contr[0].l;
        if (L > maxL) maxL = L;
        size_t curCen = basis.mapSh2Cen[iShell];
        bool newAtm = false;
        if( curCen != iAtm ) {
          iAtm++;
          if( !groupAtm or (groupAtm and (atmSymb != getAtmSymb(iAtm))) ) {
            newAtm = true;
            atmSymb = getAtmSymb(iAtm);
            for(size_t i = 0; i < end; i++)
              std::fill(angWeight[i].begin(), angWeight[i].end(), 0.0);
          }
        }
        bool lastAtm  = iAtm == (mol.atoms.size() - 1);
        bool firstAtm = iShell == 0;
        bool nextAtmNew;

        for(size_t mu = bfst; mu < bfst + sz; mu++) {
          for(auto q = p; q < p + end; q++) {
            angWeight[q-p][L] += SmartConj(MO[mu + q*LDM]) * SCR[mu + q*NB];
            if (MO2)
              angWeight[q-p][L] += SmartConj(MO2[mu + q*LDM]) * SCR2[mu + q*NB];

          }
          nextAtmNew = lastAtm ? false : (mu+1) >= basis.mapCen2BfSt[iAtm+1];
          if (groupAtm and (atmSymb == getAtmSymb(iAtm+1))) nextAtmNew = false;
        }

        if( (newAtm or firstAtm) ) {
            out << std::setw(6) << " ";
            out << std::setw(3) << std::left << iAtm;
            out << std::setw(7) << std::left << atmSymb;
            out << std::endl;
        }
        if( nextAtmNew or (iShell + 1 == nShell) ) {
          for (size_t i = 0; i <= maxL; i++) {
            out << " " << std::left << std::setw(eValOff-1) << angLabel[i];
            out << std::right;
            for (auto q = p; q < p+end; q++) {
              double VAL = std::real(angWeight[q-p][i]);
              std::cout << std::fixed << std::right<< std::setprecision(4);
              std::cout.fill(' ');

              if(std::abs(VAL) >= 1e-4)    out << std::setw(printWidth) << VAL;
              else if(std::isnan(std::abs(VAL))) out << std::setw(printWidth) << "NAN";
              else if(std::isinf(std::abs(VAL))) out << std::setw(printWidth) << "INF";
              else                               out << std::setw(printWidth) << "-";
            }
            out << std::endl;
          }
        }
      }
      out << std::endl;
    }

    CQMemManager::get().free(SCR);
    out << bannerEnd << std::endl;

  } // analyzeMOPrint


  template <typename MatsT, typename IntsT>
  void WaveFunction<MatsT,IntsT>::printMOInfo(std::ostream &out, size_t printLevel) {

    size_t NB = this->nAlphaOrbital();
    size_t NOrb = NB * this->nC;

    bool MOcoeffs = printLevel % 2 == 1;
    bool analyzeMO = printLevel / 2 > 0;

    if (MOcoeffs) {

    printLevel -= 1;

    // Pretty MO print
    std::function<double(MatsT)> printOp1 = [](MatsT x) { return std::real(x); };
    std::function<double(MatsT)> printOp2 = printOp1;

    if( std::is_same<MatsT,dcomplex>::value ) {
      printOp1 = [](MatsT x) { return std::abs(x); };
      printOp2 = [](MatsT x) { return std::arg(x); };
    }

    if( this->nC >= 2 )
    out << " *** NOTICE: Alpha and Beta Coefficients refer to the SAME "
      << "Canonical MOs ***\n";

    out << "\n\nCanonical Molecular Orbital Coefficients (Alpha)";
    if( std::is_same<MatsT,dcomplex>::value )
      out << " Magnitude";
    if( this->nC == 4 ) out << " for Large component";

    prettyMOPrint(out,NB,NOrb,this->eps1,this->mo[0].pointer(),NOrb,
        molecule(),basisSet(),printOp1);

    if( std::is_same<MatsT,dcomplex>::value ) {
      out << "\n\nCanonical Molecular Orbital Coefficients (Alpha) Phase";
      if( this->nC == 4 ) out << " for Large component"; 
      prettyMOPrint(out,NB,NOrb,this->eps1,this->mo[0].pointer(),NOrb,
          molecule(),basisSet(),printOp2);
    }

    if( this->nC >= 2 or not this->iCS ) {

      out << "\n\nCanonical Molecular Orbital Coefficients (Beta)";
      if( std::is_same<MatsT,dcomplex>::value )
        out << " Magnitude";
      if( this->nC == 4 ) out << " for Large component";

      if( this->nC == 1 )
        prettyMOPrint(out,NB,NOrb,this->eps2,this->mo[1].pointer(),NOrb,
            molecule(),basisSet(),printOp1);
      else
        prettyMOPrint(out,NB,NOrb,this->eps1,this->mo[0].pointer() + (this->nC/2)*NB,
            NOrb,molecule(),basisSet(),printOp1);

      if( std::is_same<MatsT,dcomplex>::value ) {
        out << "\n\nCanonical Molecular Orbital Coefficients (Beta) Phase";
        if( this->nC == 4 ) out << " for Large component";

        if( this->nC == 1 )
          prettyMOPrint(out,NB,NOrb,this->eps2,this->mo[1].pointer(),NOrb,
              molecule(),basisSet(),printOp2);
        else
          prettyMOPrint(out,NB,NOrb,this->eps1,this->mo[0].pointer() + (this->nC/2)*NB,
              NOrb,molecule(),basisSet(),printOp2);
      }
    }

    if( this->nC == 4 ) {

      out << "\n\nCanonical Molecular Orbital Coefficients (Alpha)";
      if( std::is_same<MatsT,dcomplex>::value )
        out << " Magnitude";
      out << " for Small component";
      prettyMOPrint(out,NB,NOrb,this->eps1,this->mo[0].pointer() + NB,NOrb,
        molecule(),basisSet(),printOp1);

      if( std::is_same<MatsT,dcomplex>::value ) {
        out << "\n\nCanonical Molecular Orbital Coefficients (Alpha) Phase for Small component";
        prettyMOPrint(out,NB,NOrb,this->eps1,this->mo[0].pointer() + NB,NOrb,
        molecule(),basisSet(),printOp2);
      }

      out << "\n\nCanonical Molecular Orbital Coefficients (Beta)";
      if( std::is_same<MatsT,dcomplex>::value )
        out << " Magnitude";
      out << " for Small component";
      prettyMOPrint(out,NB,NOrb,this->eps1,this->mo[0].pointer() + 3*NB,NOrb,
        molecule(),basisSet(),printOp1);

      if( std::is_same<MatsT,dcomplex>::value ) {
        out << "\n\nCanonical Molecular Orbital Coefficients (Beta) Phase for Small component";
        prettyMOPrint(out,NB,NOrb,this->eps1,this->mo[0].pointer() + 3*NB,NOrb,
        molecule(),basisSet(),printOp2);
      }

    }
    out << "\n" << BannerEnd << "\n\n";

    }

    if ( analyzeMO ) {

    // Analysis angular momentum components of MOs
    bool groupAtm = (printLevel/2) % 2 == 0;
    bool groupAB = printLevel/2 >= 3;

    MatsT * MO2 = nullptr;
    if( groupAB ) {
      out << "\n\nCanonical Molecular Orbital based Mulliken Population Analysis (Alpha + Beta)";
      if (this->nC == 1 and not this->iCS) MO2 = this->mo[1].pointer();
      else MO2 =  this->mo[0].pointer() + (this->nC/2)*NB;
    }
    else if( this->nC >= 2 )
      out << "\n *** NOTICE: Alpha and Beta Analysis refer to the SAME "
        << "Canonical MOs ***\n";

    out << "\n\nCanonical Molecular Orbital based Mulliken Population Analysis (Alpha)";

    if( this->nC == 4 ) out << " for Large component";
    analyzeMOPrint(out, NB, NOrb, aoints_->overlap->pointer(), this->mo[0].pointer(),
            NOrb, molecule(), basisSet(), groupAtm, MO2);

    if( not groupAB and (this->nC >= 2 or not this->iCS) ) {
      out << "\n\nCanonical Molecular Orbital based Mulliken Population Analysis (Beta)";
      if( this->nC == 4 ) out << " for Large component";
      if( this->nC == 1 )
        analyzeMOPrint(out, NB, NOrb, aoints_->overlap->pointer(),
              this->mo[1].pointer(), NOrb, molecule(), basisSet(), groupAtm);
      else
        analyzeMOPrint(out, NB, NOrb, aoints_->overlap->pointer(),
              this->mo[0].pointer() + (this->nC/2)*NB, NOrb, molecule(), basisSet(),
              groupAtm);

    }

    if( this->nC == 4 ) {
      IntsT* ssOverlap = CQMemManager::get().malloc<IntsT>(NB*NB);
      SetMat('N',NB,NB,1./(2*SpeedOfLight*SpeedOfLight),this->aoints_->kinetic->pointer(),
                NB,ssOverlap,NB);
      out << "\n\nCanonical Molecular Orbital based Mulliken Population Analysis (Alpha) for Small component";
      analyzeMOPrint(out, NB, NOrb, ssOverlap, this->mo[0].pointer() + NB,
            NOrb, molecule(), basisSet(), groupAtm, MO2);
      if( not groupAB ) {
        out << "\n\nCanonical Molecular Orbital based Mulliken Population Analysis (Beta) for Small component";
        analyzeMOPrint(out, NB, NOrb, ssOverlap, this->mo[0].pointer() + 3*NB,
               NOrb, molecule(), basisSet(), groupAtm);
      }
      CQMemManager::get().free(ssOverlap);

    }


    out << "\n" << BannerEnd << "\n\n";

    }

  }; // WaveFunction::printMOInfo



}; // namespace ChronusQ

