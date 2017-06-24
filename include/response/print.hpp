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

#include <response/tbase.hpp>

#include <physcon.hpp>

namespace ChronusQ {

  template <typename T>
  void ResponseTBase<T>::printTMoments(std::ostream &out) {


    out << "\n\n\n* RESIDUE TRANSITION MOMENTS\n\n\n";

    size_t nRoots = resSettings.nRoots;
 
    double xSmall = 1e-6;



    auto printRec = [&](T x) {
      out << std::setprecision(5) << std::scientific << std::setw(15)
          << std::right << ((std::abs(x) > xSmall) ? x : 0. );
    };

    auto printDipole = [&]( T* dipole ) {

      out << "                                { X, Y, Z }\n\n";

      for( auto iSt = 0; iSt < nRoots; iSt++ ) {

        out << "    " << "n = " << std::setw(5) << std::left << iSt + 1;
        printRec(dipole[3*iSt + 0]);
        printRec(dipole[3*iSt + 1]);
        printRec(dipole[3*iSt + 2]);

        out << "\n";

      }

    };


    auto printQuadrupole = [&]( T* quadrupole ) {

      out << "                              { XX, XY, XZ }\n";
      out << "                              { YY, YZ, ZZ }\n\n";

      for( auto iSt = 0; iSt < nRoots; iSt++ ) {

        out << "    " << "n = " << std::setw(5) << std::left << iSt + 1;

        printRec(quadrupole[6*iSt + 0]);
        printRec(quadrupole[6*iSt + 1]);
        printRec(quadrupole[6*iSt + 2]);

        out << "\n" << "    " << "    " << std::setw(5) << " ";

        printRec(quadrupole[6*iSt + 3]);
        printRec(quadrupole[6*iSt + 4]);
        printRec(quadrupole[6*iSt + 5]);

        out << "\n\n";

      }

    };

    auto printOctupole = [&]( T* octupole ) {

      out << "                             { XXX, XXY, XXZ }\n";
      out << "                             { XYY, XYZ, XZZ }\n";
      out << "                             { YYY, YYZ, YZZ }\n";
      out << "                             { ZZZ           }\n\n";

      for( auto iSt = 0; iSt < nRoots; iSt++ ) {

        out << "    " << "n = " << std::setw(5) << std::left << iSt + 1;

        printRec(octupole[10*iSt + 0]);
        printRec(octupole[10*iSt + 1]);
        printRec(octupole[10*iSt + 2]);

        out << "\n" << "    " << "    " << std::setw(5) << " ";

        printRec(octupole[10*iSt + 3]);
        printRec(octupole[10*iSt + 4]);
        printRec(octupole[10*iSt + 5]);

        out << "\n" << "    " << "    " << std::setw(5) << " ";

        printRec(octupole[10*iSt + 6]);
        printRec(octupole[10*iSt + 7]);
        printRec(octupole[10*iSt + 8]);

        out << "\n" << "    " << "    " << std::setw(5) << " ";

        printRec(octupole[10*iSt + 9]);

        out << "\n\n";

      }

    };





    if( resResults.tLenElecDipole_ge ) {

      out << "  Transition Dipole Moments (Length) : < 0 | r | n > (AU)\n\n\n";

      printDipole(resResults.tLenElecDipole_ge);

      out << "\n\n\n\n";

    }


    if( resResults.tVelElecDipole_ge ) {

      out << "  Transition Dipole Moments (Velocity) : i < 0 | p | n > (AU)\n\n\n";

      printDipole(resResults.tVelElecDipole_ge);

      out << "\n\n\n\n";

    }

    if( resResults.tLenElecQuadrupole_ge ) {

      out << "  Transition Quadrupole Moments (Length) : < 0 | r_i r_j | n > (AU)\n\n\n";

      printQuadrupole(resResults.tLenElecQuadrupole_ge);

      out << "\n\n\n\n";

    }


    if( resResults.tVelElecQuadrupole_ge ) {

      out << "  Transition Quadrupole Moments (Velocity) : i < 0 | p_i r_j + r_i p_j | n > (AU)\n\n\n";

      printQuadrupole(resResults.tVelElecQuadrupole_ge);

      out << "\n\n\n\n";

    }


    if( resResults.tLenElecOctupole_ge ) {

      out << "  Transition Octupole Moments (Length) : < 0 | r_i r_j r_k | n > (AU)\n\n\n";


      printOctupole(resResults.tLenElecOctupole_ge);

      out << "\n\n\n\n";

    }



    if( resResults.tVelElecOctupole_ge ) {

      out << "  Transition Octupole Moments (Velocity) : i < 0 | p_i r_j r_k + r_i p_j r_k + r_i r_j p_k | n > (AU)\n\n\n";

      printOctupole(resResults.tVelElecOctupole_ge);

      out << "\n\n\n\n";

    }


    if( resResults.tMagDipole_ge ) {
      out << "  Transition Magnetic Moments (Length) : \n";
      out << "    i < 0 | r x p | n > (AU)\n\n\n";

      printDipole(resResults.tMagDipole_ge);

      out << "\n\n\n\n";
    }


  }


  template <typename T>
  void ResponseTBase<T>::printResObservables(std::ostream &out) {

    size_t nRoots = resSettings.nRoots;
 
    double xSmall = 1e-6;



    auto printRec = [&](T x) {
      out << std::setprecision(8) << std::scientific << std::setw(15)
          << std::right << ((std::abs(x) > xSmall) ? x : 0. );
    };


    if( resObs.rotatory_len_RM ) {

      out << "\n\n\nROTATORY STRENGTH (LENGTH) \n\n";
      out << "  * R(n) = 0.5 * Im[ < 0 | r_k | n > < n | (r x p)_k | 0 > ]\n";

      out << "\n\n\n";

      out << "   " << "    " << std::setw(5) << " " << "R(n) (10^40 erg-esu-cm / Gauss)\n\n";
      for(auto iO = 0; iO < nRoots; iO++) {

        out << "   " << "n = " << std::setw(5) << std::left << iO + 1;
        printRec(Rotatory_CGS_Length * resObs.rotatory_len_RM[iO]);
        out << "\n";

      };
      out << "\n";

    }

  };


  template <typename T>
  template <typename U>
  void ResponseTBase<T>::printRF(
    FDResponseResults<T,U> &results, std::ostream &out) {

    out << "\n\n\n* RESPONSE FUNCTIONS (POLARIZABILITIES)\n\n\n";

    size_t nOmega = fdrSettings.bFreq.size();


    // Printing Lambdas

    double xSmall = 1e-6;

    auto printRecRe = [&](U x) {
      out << std::setprecision(5) << std::scientific << std::setw(15)
          << std::right << ((std::abs(x) > xSmall) ? std::real(x) : 0. );
    };

    auto printRecIm = [&](U x) {
      out << std::setprecision(5) << std::scientific << std::setw(15)
          << std::right << ((std::abs(x) > xSmall) ? std::imag(x) : 0. );
    };



    // Helper strings
    std::string ededHelper = 
     "                         { << X; X >>, << X; Y >>, << X; Z >> }\n"
     "                         { << Y; X >>, << Y; Y >>, << Y; Z >> }\n"
     "                         { << Z; X >>, << Z; Y >>, << Z; Z >> }\n";

    std::string eqedHelper =
     "                        { << XX; X >>, << XX; Y >>, << XX; Z >> }\n"
     "                        { << XY; X >>, << XY; Y >>, << XY; Z >> }\n"
     "                        { << XZ; X >>, << XZ; Y >>, << XZ; Z >> }\n"
     "                        { << YY; X >>, << YY; Y >>, << YY; Z >> }\n"
     "                        { << YZ; X >>, << YZ; Y >>, << YZ; Z >> }\n"
     "                        { << ZZ; X >>, << ZZ; Y >>, << ZZ; Z >> }\n";


    // Print out a D-D Polarizability
    auto printDD = [&]( std::function<void(U)> printRec, U *dd_polar ) -> void {

      out << ededHelper << std::endl;
      for(auto iOmega = 0; iOmega < nOmega; iOmega++) {

        double omega = fdrSettings.bFreq[iOmega];
        U*     dd  = dd_polar + iOmega*9;

        out << "    " << "W(AU) = " << std::setw(8) << std::setprecision(4)
            << std::fixed << std::left << omega;

        printRec(dd[0]); printRec(dd[3]); printRec(dd[6]);

        out << "\n" << "    " << "        " << std::setw(8) << " ";
        printRec(dd[1]); printRec(dd[4]); printRec(dd[7]);

        out << "\n" << "    " << "        " << std::setw(8) << " ";
        printRec(dd[2]); printRec(dd[5]); printRec(dd[8]);

        out << "\n\n";

      }
      out << "\n\n\n\n";
    };



    // Print out a Q-D Polarizability
    auto printQD = [&]( std::function<void(U)> printRec, U* qd_polar ) -> void {
      out << eqedHelper << std::endl;
      for(auto iOmega = 0; iOmega < nOmega; iOmega++) {

        double omega = fdrSettings.bFreq[iOmega];
        U*     qd  = qd_polar + iOmega*3*6;

        out << "    " << "W(AU) = " << std::setw(8) << std::setprecision(4)
            << std::fixed << std::left << omega;

        printRec(qd[0]); printRec(qd[6]); printRec(qd[12]);

        out << "\n" << "    " << "        " << std::setw(8) << " ";
        printRec(qd[1]); printRec(qd[7]); printRec(qd[13]);

        out << "\n" << "    " << "        " << std::setw(8) << " ";
        printRec(qd[2]); printRec(qd[8]); printRec(qd[14]);

        out << "\n" << "    " << "        " << std::setw(8) << " ";
        printRec(qd[3]); printRec(qd[9]); printRec(qd[15]);

        out << "\n" << "    " << "        " << std::setw(8) << " ";
        printRec(qd[4]); printRec(qd[10]); printRec(qd[16]);

        out << "\n" << "    " << "        " << std::setw(8) << " ";
        printRec(qd[5]); printRec(qd[11]); printRec(qd[17]);

        out << "\n\n";

      }
      out << "\n\n\n\n";
    };





    if(results.ed_ed_Polar) {

      out << "  Electric Dipole - Electric Dipole (Length) : ";
      out << "- Re [<< r_i; r_j >>] (AU)\n\n\n";


      printDD(printRecRe, results.ed_ed_Polar);

      if( std::is_same<U,dcomplex>::value ) {
        out << "  Electric Dipole - Electric Dipole (Length) : ";
        out << "- Im[<< r_i; r_j >>] (AU)\n\n\n";

        printDD(printRecIm, results.ed_ed_Polar);
      }

    }

    if(results.eq_ed_Polar) {

      out << "  Electric Quadrupole - Electric Dipole (Length) : ";
      out << "- Re[<< r_i r_j; r_k >>] (AU)\n\n\n";

      printQD(printRecRe, results.eq_ed_Polar);

      if( std::is_same<U,dcomplex>::value ) {
        out << "  Electric Quadrupole - Electric Dipole (Length) : ";
        out << "- Im[<< r_i r_j; r_k >>] (AU)\n\n\n";

        printQD(printRecIm, results.eq_ed_Polar);
      }

    }

    if(results.md_ed_Polar) {

      out << "  Magnetic Dipole - Electric Dipole (Length) : ";
      out << "Im[<< (r x p)_i; r_j >>] (AU)\n\n\n";


      printDD(printRecRe, results.md_ed_Polar);

      if( std::is_same<U,dcomplex>::value ) {
        out << "  Magnetic Dipole - Electric Dipole (Length) : ";
        out << "- Re[<< (r x p)_i; r_j >>] (AU)\n\n\n";


        printDD(printRecIm, results.md_ed_Polar);
      }

    }


    if(results.md_md_Polar) {

      out << "  Magnetic Dipole - Magnetic Dipole : ";
      out << "Im[<< (r x p)_i; (r x p)_j >>] (AU)\n\n\n";


      printDD(printRecRe, results.md_md_Polar);

      if( std::is_same<U,dcomplex>::value ) {
        out << "  Magnetic Dipole - Magnetic Dipole : ";
        out << "- Re[<< (r x p)_i; (r x p)_j >>] (AU)\n\n\n";


        printDD(printRecIm, results.md_md_Polar);
      }

    }
    

  };

  template <typename T>
  void ResponseTBase<T>::printFDObservables(std::ostream &out) {


    // FIXME: This is because only the OPA is calculated now
    if( not fdObs.opaCross_eda ) return;

    size_t nOmega = fdrSettings.bFreq.size();

    double xSmall = 1e-6;

    auto printRec = [&](double x) {
      out << std::setprecision(10) << std::scientific << std::setw(25)
          << std::right << ((std::abs(x) > xSmall) ? x : 0. );
    };

    out << "\n\n\n* OBSERVABLES \n\n\n";

    out << "  ONE-PHOTON ABSORPTION CROSS-SECTION (EDA) \n\n";
    out << "    * SIGMA(W) = 4 * PI * W / C * IM[ALPHA(-W,W)]\n"; 
    out << "    * ALPHA(-W,W) = TR[ << r; r >>(W) ]\n" ;

    out << "\n\n\n";

    out << "    " << std::right << std::setw(25) << "W"
                                << std::setw(25) << "SIGMA(W) (AU)\n";

    for(auto iOmega = 0; iOmega < nOmega; iOmega++) {

      double omega = fdrSettings.bFreq[iOmega];

      out << "    ";
      printRec(omega); printRec(-fdObs.opaCross_eda[iOmega]);
      out << "\n";

    }

    out << "\n\n\n\n";
  }

}; // namespace ChronusQ

